/**************************************************************************
* CSVVs — Spec 2 (age-interacted, ALLREA base + product deltas) — SIMPLE
* Input : code/regressions/output/age_spec_ALLREA/agespec_prod_interact_spec2_ALLREA.ster
* Output: code/regressions/output/age_spec_interacted_csvv/agespec_interaction_response_spec2_<PRODUCT>.csvv
* Shape : gamma(36) in order [age 1..3 × poly 1..4 × covars (1, lr, loginc)],
*         gammavcv(36×36) in the same flat order
**************************************************************************/

version 17.0
set more off

/**************** CONFIG (EDIT THESE 3 LINES IF NEEDED) ****************/
global REPO "/user/ab5405/summeraliaclimate"
local AGEVAR      "agegroup"        // or "agegroup_"
local AGE_LEVELS  "1 2 3"           // or "0 1 2" if encoded that way
local PROD_MAP    "GMFD=1 ERA5_025=2 MERRA2=3 JRA_3Q=4"  // your i.product map
local PRODUCTS    "ERA5_025 GMFD MERRA2 JRA_3Q"          // which CSVVs to write
/***********************************************************************/

local STERDIR  "$REPO/code/regressions/output/age_spec_ALLREA"
local STERFILE "agespec_prod_interact_spec2_ALLREA.ster"

local OUTDIR   "$REPO/code/regressions/output/age_spec_interacted_csvv"
cap mkdir "`OUTDIR'"

/**************** LOAD .STER ONCE ****************/
capture estimates use "`STERDIR'/`STERFILE'"
if _rc {
    di as err "!! Could not open ster: `STERDIR'/`STERFILE'"
    exit 198
}
matrix B = e(b)
matrix V = e(V)

/**************** HELPER ****************/
cap program drop _pickfirst
program define _pickfirst, rclass
    // Return first existing colname among up to 8 candidates
    syntax , [C1(string) C2(string) C3(string) C4(string) C5(string) C6(string) C7(string) C8(string)]
    local rv ""
    forvalues j=1/8 {
        local nm : word `j' of "`c1' `c2' `c3' `c4' `c5' `c6' `c7' `c8'"
        if "`nm'"=="" continue
        if (colnumb(V,"`nm'")>0) {
            local rv "`nm'"
            continue, break
        }
    }
    return local name "`rv'"
end

/**************** MAIN LOOP (HARD-CODED MAP) ****************/
tokenize "`PROD_MAP'"
while "`1'"!="" {
    local pair "`1'"
    macro shift

    local prod_lbl = substr("`pair'", 1, strpos("`pair'", "=")-1)
    local lvl_str  = substr("`pair'", strpos("`pair'", "=")+1, .)
    local pcode = real("`lvl_str'")
    if missing(`pcode') {
        di as err "!! Invalid i.product level for `prod_lbl': `lvl_str'"
        exit 198
    }

    // only build requested products
    local do_it 0
    foreach P of local PRODUCTS {
        if "`P'"=="`prod_lbl'" local do_it 1
    }
    if !`do_it' continue

    di as txt "=============================================================="
    di as txt "CSVV for product: `prod_lbl'  (i.product = `pcode')"
    di as txt "Using AGEVAR = `AGEVAR' with levels: `AGE_LEVELS'"
    di as txt "=============================================================="

    /******** Resolve 36 slots: ages × polys × [1, lr, loginc] ********/
    local coef_base_names   // 36 entries (must exist)
    local coef_delta_names  // 36 entries ("__ZERO__" if missing)

    foreach a of local AGE_LEVELS {
        forvalues k = 1/4 {
            /* ----- BASE ×1 (no product) ----- */
            local base1 ""
            local c1 "`a'.`AGEVAR'#c.tavg_poly_`k'"
            local c2 "c.tavg_poly_`k'#`a'.`AGEVAR'"
            _pickfirst, c1("`c1'") c2("`c2'")
            if "`r(name)'"!="" local base1 "`r(name)'"
            if "`base1'"=="" {
                di as err "!! Missing BASE ×1 (age=`a', poly=`k') — check AGEVAR/LEVELS"
                exit 459
            }

            /* ----- BASE ×LR (no product) ----- */
            local base2 ""
            local c1 "`a'.`AGEVAR'#c.tavg_poly_`k'#c.lr_tavg_adm1_avg"
            local c2 "c.tavg_poly_`k'#c.lr_tavg_adm1_avg#`a'.`AGEVAR'"
            _pickfirst, c1("`c1'") c2("`c2'")
            if "`r(name)'"!="" local base2 "`r(name)'"
            if "`base2'"=="" {
                di as err "!! Missing BASE ×LR (age=`a', poly=`k')"
                exit 459
            }

            /* ----- BASE ×loginc (no product) ----- */
            local base3 ""
            local c1 "`a'.`AGEVAR'#c.tavg_poly_`k'#c.loggdppc_adm1_avg"
            local c2 "c.tavg_poly_`k'#c.loggdppc_adm1_avg#`a'.`AGEVAR'"
            _pickfirst, c1("`c1'") c2("`c2'")
            if "`r(name)'"!="" local base3 "`r(name)'"
            if "`base3'"=="" {
                di as err "!! Missing BASE ×loginc (age=`a', poly=`k')"
                exit 459
            }

            /* ----- DELTA ×1 (with #`pcode'.product) ----- */
            local dlt1 "__ZERO__"
            local d1 "`a'.`AGEVAR'#c.tavg_poly_`k'#`pcode'.product"
            local d2 "c.tavg_poly_`k'#`pcode'.product#`a'.`AGEVAR'"
            local d3 "`a'.`AGEVAR'#`pcode'.product#c.tavg_poly_`k'"
            local d4 "`pcode'.product#`a'.`AGEVAR'#c.tavg_poly_`k'"
            _pickfirst, c1("`d1'") c2("`d2'") c3("`d3'") c4("`d4'")
            if "`r(name)'"!="" local dlt1 "`r(name)'"

            /* ----- DELTA ×LR ----- */
            local dlt2 "__ZERO__"
            local e1 "`a'.`AGEVAR'#c.tavg_poly_`k'#c.lr_tavg_adm1_avg#`pcode'.product"
            local e2 "c.tavg_poly_`k'#c.lr_tavg_adm1_avg#`pcode'.product#`a'.`AGEVAR'"
            local e3 "`a'.`AGEVAR'#`pcode'.product#c.tavg_poly_`k'#c.lr_tavg_adm1_avg"
            local e4 "`pcode'.product#`a'.`AGEVAR'#c.tavg_poly_`k'#c.lr_tavg_adm1_avg"
            _pickfirst, c1("`e1'") c2("`e2'") c3("`e3'") c4("`e4'")
            if "`r(name)'"!="" local dlt2 "`r(name)'"

            /* ----- DELTA ×loginc ----- */
            local dlt3 "__ZERO__"
            local f1 "`a'.`AGEVAR'#c.tavg_poly_`k'#c.loggdppc_adm1_avg#`pcode'.product"
            local f2 "c.tavg_poly_`k'#c.loggdppc_adm1_avg#`pcode'.product#`a'.`AGEVAR'"
            local f3 "`a'.`AGEVAR'#`pcode'.product#c.tavg_poly_`k'#c.loggdppc_adm1_avg"
            local f4 "`pcode'.product#`a'.`AGEVAR'#c.tavg_poly_`k'#c.loggdppc_adm1_avg"
            _pickfirst, c1("`f1'") c2("`f2'") c3("`f3'") c4("`f4'")
            if "`r(name)'"!="" local dlt3 "`r(name)'"

            // accumulate in Spec-2 order [1, lr, loginc]
            local coef_base_names  `coef_base_names'  `base1'  `base2'  `base3'
            local coef_delta_names `coef_delta_names' `dlt1'   `dlt2'   `dlt3'
        }
    }

    /* Sanity */
    local n1 : word count `coef_base_names'
    local n2 : word count `coef_delta_names'
    if (`n1'!=36 | `n2'!=36) {
        di as err "!! Expected 36 slots; got base=`n1', delta=`n2' for `prod_lbl'"
        exit 459
    }
    foreach nm of local coef_base_names {
        if (colnumb(V,"`nm'")==0 | colnumb(B,"`nm'")==0) {
            di as err "!! Missing base column in V or B: `nm'"
            exit 459
        }
    }

    /**************** WRITE CSVV (structure unchanged) ****************/
    local CSVVFILE "agespec_interaction_response_spec2_`=subinstr("`prod_lbl'","-","_",.)'.csvv"
    cap file close csvv
    file open csvv using "`OUTDIR'/`CSVVFILE'", write replace

    file write csvv "---" _n
    file write csvv "oneline: Age-specific interacted (Spec 2) response – `prod_lbl'" _n
    file write csvv "version: MORTALITY-AGE-SPEC-INTERACTED-`=subinstr("`prod_lbl'","-","_",.)'" _n
    file write csvv "dependencies: `STERFILE'" _n
    file write csvv "description: 4th-order poly OLS; age×ADM2 & age×Country×Year FE; base AGE×Temp(+LR,+loginc) + PRODUCT deltas (ALLREA pooled)" _n
    file write csvv "csvv-version: aliabonanno-`c(current_date)'" _n

    file write csvv "variables:" _n
    file write csvv "  tavg_poly_1: T1_sum [C]" _n
    file write csvv "  tavg_poly_2: T2_sum [C^2]" _n
    file write csvv "  tavg_poly_3: T3_sum [C^3]" _n
    file write csvv "  tavg_poly_4: T4_sum [C^4]" _n
    file write csvv "  lr_tavg: long-run mean annual temperature [C]" _n
    file write csvv "  deathrate_w99: mortality rate winsorized at 99th percentile" _n
    file write csvv "  loggdppc_adm1_avg: log GDP per capita" _n
    file write csvv "..." _n

    file write csvv "observations" _n
    file write csvv " `e(N)'" _n

    file write csvv "prednames" _n
    file write csvv "tavg_poly_1, tavg_poly_2, tavg_poly_3, tavg_poly_4" _n

    file write csvv "covarnames" _n
    file write csvv "1, lr_tavg, loggdppc_adm1_avg" _n

    // --- gamma: base + (delta if present), 36 numbers ---
    file write csvv "gamma" _n
    forvalues i = 1/36 {
        local bn : word `i' of `coef_base_names'
        local dn : word `i' of `coef_delta_names'
        scalar gval = _b[`bn']
        if ("`dn'"!="__ZERO__") scalar gval = gval + _b[`dn']
        if mod(`i',3)!=1 file write csvv ", "
        file write csvv %21.10g (gval)
        if mod(`i',3)==0 file write csvv _n
    }

    // --- gammavcv: Var(base+delta) = Vbb + Vdd + Vbd + Vdb ---
    file write csvv "gammavcv" _n
    forvalues r = 1/36 {
        local bn_r : word `r' of `coef_base_names'
        local dn_r : word `r' of `coef_delta_names'
        local irb = colnumb(V, "`bn_r'")
        local ird = ("`dn_r'"=="__ZERO__" ? . : colnumb(V, "`dn_r'"))

        forvalues c = 1/36 {
            local bn_c : word `c' of `coef_base_names'
            local dn_c : word `c' of `coef_delta_names'
            local icb = colnumb(V, "`bn_c'")
            local icd = ("`dn_c'"=="__ZERO__" ? . : colnumb(V, "`dn_c'"))

            scalar vrc = el(V, `irb', `icb')                 // Vbb
            if (`ird'<.) scalar vrc = vrc + el(V, `ird', `icb')  // Vdb
            if (`icd'<.) scalar vrc = vrc + el(V, `irb', `icd')  // Vbd
            if (`ird'<. & `icd'<.) scalar vrc = vrc + el(V, `ird', `icd') // Vdd

            if (`c'==1) file write csvv %21.10g (vrc)
            else        file write csvv ", " %21.10g (vrc)
        }
        file write csvv _n
    }

    file write csvv "residvcv" _n
    file write csvv %21.10g (e(rmse)^2) _n

    file close csvv
    di as res "✅ Wrote agespec_interaction_response_spec2_`=subinstr("`prod_lbl'","-","_",.)'.csvv"
}

di as res "All CSVVs generated."
