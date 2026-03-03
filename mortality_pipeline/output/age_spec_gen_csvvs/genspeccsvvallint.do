/**************************************************************************
* CSVV maker for Spec 2 age-specific interacted regressions — per product
* 36 gammas: ages [0–4,5–64,65+] × polys 1..4 × {1, LR, loginc}
**************************************************************************/

version 17.0
clear all
set more off

/************ CONFIG ************/
global REPO   "/user/ab5405/summeraliaclimate"
local PANELS  "$REPO/code/regressions/output_panels"
local STERDIR "$REPO/code/regressions/output/age_spec_interacted"
local CSVVDIR "$REPO/code/regressions/output/age_spec_interacted_csvv"
cap mkdir "`CSVVDIR'"

local products "ERA5_025 GMFD MERRA2 JRA_3Q"

/* ---------- Helper: find numeric codes for bins 0–4, 5–64, 65+ ---------- */
program define _age_bin_codes, rclass
    syntax , PRODUCT(name) PANELS(string)

    local try1 "`panels'/panel_prepped_for_regressions_`product'_aligned.dta"
    local try2 "`panels'/panel_prepped_for_regressions_`product'.dta"

    preserve
        capture confirm file "`try1'"
        if !_rc {
            use "`try1'", clear
        }
        else {
            capture confirm file "`try2'"
            if !_rc use "`try2'", clear
            else {
                di as err "!! _age_bin_codes: panel not found for `product'"
                restore
                exit 498
            }
        }

        /* Ensure there is a numeric age var WITH value labels, using tempvars only */
        tempvar a s t
        capture confirm variable agegroup
        if _rc {
            di as err "!! _age_bin_codes: variable agegroup not found."
            restore
            exit 498
        }

        capture confirm numeric variable agegroup
        if _rc {
            /* agegroup is string -> encode to temp numeric `a' */
            encode agegroup, gen(`a')
        }
        else {
            /* agegroup is numeric -> use it as-is */
            gen double `a' = agegroup
        }

        /* Make sure `a' has a value label; if not, rebuild via decode/encode */
        local lbl : value label `a'
        if "`lbl'" == "" {
            decode `a', gen(`s')
            encode `s', gen(`t')
            drop `a'
            gen double `a' = `t'
            local lbl : value label `a'
        }
        if "`lbl'" == "" {
            di as err "!! agegroup has no value label; cannot map bins."
            restore
            exit 498
        }

        /* Find codes whose labels are exactly 0-4, 5-64, 65+ (after dash/space normalization) */
        levelsof `a', local(allcodes)
        local code_0_4 ""
        local code_5_64 ""
        local code_65p ""

        foreach c of local allcodes {
            local L : label (`lbl') `c'
            /* normalize dashes/spaces */
            local L = strtrim("`L'")
            local L = subinstr("`L'", uchar(8211), "-", .)
            local L = subinstr("`L'", uchar(8212), "-", .)
            local L = subinstr("`L'", " ", "", .)
            if "`L'" == "0-4"   local code_0_4  `c'
            if "`L'" == "5-64"  local code_5_64 `c'
            if inlist("`L'","65+","65-+") local code_65p `c'
        }

        if "`code_0_4'"=="" | "`code_5_64'"=="" | "`code_65p'"=="" {
            di as err "!! Could not find all bins (labels must be 0-4, 5-64, 65+)."
            restore
            exit 498
        }

        return local bin_codes "`code_0_4' `code_5_64' `code_65p'"
    restore
end

/* ---------- Helper: pick first existing matrix column by name ---------- */
program define _choose_colname, rclass
    syntax , MAT(name) CANDIDATES(string)
    tempname M
    matrix `M' = `mat'
    local pick ""
    foreach nm of local candidates {
        if (colnumb(`M', "`nm'") > 0) {
            local pick "`nm'"
            continue, break
        }
    }
    return local name "`pick'"
end

/* ======================= MAIN LOOP ======================= */
foreach P of local products {

    local prod_lbl "`P'"
    local prod_var : subinstr local P "-" "_", all

    di as txt "=============================================================="
    di as txt "CSVV build for product: `prod_lbl'  (stem: `prod_var')"
    di as txt "=============================================================="

    /* 0) Age-bin codes in order: 0–4, 5–64, 65+ */
    quietly _age_bin_codes , product("`prod_lbl'") panels("`PANELS'")
    if _rc {
        di as err "!! Skipping `prod_lbl' (age-bin discovery failed)."
        continue
    }
    local AGE3 `r(bin_codes)'

    /* 1) Load ster */
    local s1 "`STERDIR'/agespec_interaction_response_`prod_lbl'.ster"
    local s2 "`STERDIR'/agespec_interaction_response_`prod_var'.ster"
    capture estimates use "`s1'"
    if _rc {
        capture estimates use "`s2'"
        if _rc {
            di as err "!! Could not open ster: `s1' or `s2'"
            continue
        }
        else local STERFILE "agespec_interaction_response_`prod_var'.ster"
    }
    else     local STERFILE "agespec_interaction_response_`prod_lbl'.ster"

    /* 2) Pull matrices */
    matrix B = e(b)
    matrix V = e(V)

    /* 3) Build ordered names (36): for each age in AGE3, k=1..4, slots {1,LR,log} */
    local coefnames
    local bad 0
    foreach a of local AGE3 {
        forvalues k = 1/4 {
            /* base */
            local c1 "`a'.agegroup#c.tavg_poly_`k'_`prod_var'"
            local c2 "c.tavg_poly_`k'_`prod_var'#`a'.agegroup"
            quietly _choose_colname , mat(V) candidates("`c1' `c2'")
            local slot1 `r(name)'

            /* ×LR */
            local d1 "`a'.agegroup#c.tavg_poly_`k'_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg"
            local d2 "c.tavg_poly_`k'_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#`a'.agegroup"
            quietly _choose_colname , mat(V) candidates("`d1' `d2'")
            local slot2 `r(name)'

            /* ×loginc */
            local e1 "`a'.agegroup#c.tavg_poly_`k'_`prod_var'#c.loggdppc_adm1_avg"
            local e2 "c.tavg_poly_`k'_`prod_var'#c.loggdppc_adm1_avg#`a'.agegroup"
            quietly _choose_colname , mat(V) candidates("`e1' `e2'")
            local slot3 `r(name)'

            if "`slot1'"=="" | "`slot2'"=="" | "`slot3'"=="" {
                di as err "!! Missing coef for agecode=`a' poly=`k' (base/LR/loginc)."
                local bad 1
            }
            else {
                local coefnames `"`coefnames' "`slot1'" "`slot2'" "`slot3'""'
            }
        }
    }
    if `bad' {
        di as err "!! Aborting `prod_lbl' — not all 36 coefs found."
        continue
    }
    local ncoef : word count `coefnames'
    if (`ncoef' != 36) {
        di as err "!! Expected 36 coefs; found `ncoef'. Aborting `prod_lbl'."
        continue
    }

    /* 4) Write CSVV */
    local CSVVFILE "agespec_interaction_response_spec2_`prod_var'.csvv"
    cap file close csvv
    file open csvv using "`CSVVDIR'/`CSVVFILE'", write replace

    file write csvv "---" _n
    file write csvv "oneline: Age-specific interacted (Spec 2) response – `prod_lbl'" _n
    file write csvv "version: MORTALITY-AGE-SPEC-INTERACTED-`prod_var'" _n
    file write csvv "dependencies: `STERFILE'" _n
    file write csvv "description: 4th-order poly OLS; age×ADM2 & age×Country×Year FE; temp×income & temp×LR-temp interactions" _n
    file write csvv "csvv-version: auto-`c(current_date)'" _n

    file write csvv "variables:" _n
    file write csvv "  tavg_poly_1: T1_sum [C]" _n
    file write csvv "  tavg_poly_2: T2_sum [C^2]" _n
    file write csvv "  tavg_poly_3: T3_sum [C^3]" _n
    file write csvv "  tavg_poly_4: T4_sum [C^4]" _n
    file write csvv "  lr_tavg: long-run mean annual temperature [C]" _n
    file write csvv "  loggdppc_adm1_avg: log GDP per capita" _n
    file write csvv "  deathrate_w99: mortality rate winsorized at 99th percentile" _n
    file write csvv "..." _n

    file write csvv "observations" _n
    file write csvv "`=e(N)'" _n

    file write csvv "prednames" _n
    file write csvv "tavg_poly_1, tavg_poly_2, tavg_poly_3, tavg_poly_4" _n

    file write csvv "covarnames" _n
    file write csvv "1, lr_tavg, loggdppc_adm1_avg" _n

    /* GAMMA */
    file write csvv "gamma" _n
    local i 0
    foreach nm of local coefnames {
        local ++i
        file write csvv "`=_b[`nm']'"
        if mod(`i',3)!=0 file write csvv ", "
        else             file write csvv "" _n
    }

    /* VCV */
    file write csvv "gammavcv" _n
    forvalues r = 1/36 {
        local rn : word `r' of `coefnames'
        local ir = colnumb(V, "`rn'")
        forvalues c = 1/36 {
            local cn : word `c' of `coefnames'
            local ic = colnumb(V, "`cn'")
            if `c'==1 file write csvv "`=el(V,`ir',`ic')'"
            else      file write csvv ", `=el(V,`ir',`ic')'"
        }
        file write csvv "" _n
    }

    /* Residual variance */
    file write csvv "residvcv" _n
    file write csvv "`=e(rmse)^2'" _n

    file close csvv
    di as res "Wrote `CSVVFILE' to `CSVVDIR'"
}

di as res "All CSVVs generated."
