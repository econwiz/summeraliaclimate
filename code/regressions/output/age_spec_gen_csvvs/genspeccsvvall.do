/**************************************************************************
* Make CSVVs from Spec 2 age-specific interacted .ster files — all products
* Products: ERA5_025, GMFD, MERRA2, JRA_3Q
**************************************************************************/

version 17.0
set more off

/************  CONFIG  ************/
global REPO "/user/ab5405/summeraliaclimate"
local STERDIR  "$REPO/code/regressions/output/age_spec_interacted"
local CSVVDIR  "$REPO/code/regressions/output/age_spec_interacted_csvv"

cap mkdir "`CSVVDIR'"

/* Products may be written with hyphens elsewhere; we normalize to underscores
   for variable stems, but we try both names for file loading. */
local products "ERA5_025 GMFD MERRA2 JRA_3Q"

/* Helper program: returns a scalar col index for a named coef in V */
program define _colix, rclass
    syntax , MNAME(name) CNAME(name)
    tempname V
    matrix `V' = `mname'
    local ix = colnumb(`V', "`cname'")
    return scalar ix = `ix'
end

foreach P of local products {

    /* Map to readable product tokens */
    local prod_lbl "`P'"
    local prod_var : subinstr local P "-" "_", all  // variable stem

    di as txt "=============================================================="
    di as txt "CSVV build for product: `prod_lbl'  (var stem: `prod_var')"
    di as txt "=============================================================="

    /************ 1) Load ster ************/
    local ster_try1 "`STERDIR'/agespec_interaction_response_`prod_lbl'_noint.ster"
    local ster_try2 "`STERDIR'/agespec_interaction_response_`prod_var'_noint.ster"

    capture estimates use "`ster_try1'"
    if _rc {
        capture estimates use "`ster_try2'"
        if _rc {
            di as err "!! Could not open ster: tried:"
            di as err "   `ster_try1'"
            di as err "   `ster_try2'"
            continue
        }
        else {
            local STERFILE : di "agespec_interaction_response_`prod_var'_noint.ster"
        }
    }
    else {
        local STERFILE : di "agespec_interaction_response_`prod_lbl'_noint.ster"
    }

    /* Grab e(b) / e(V) so we can read by name */
    matrix B = e(b)
    matrix V = e(V)

    /************ 2) Build ordered coef name list (36 entries) ************
       Order per CSVV:
         covarnames: 1, lr_tavg, loggdppc_adm1_avg
         loop ages 1..3, polys k=1..4, for each write [1, lr, loginc]
       Factor-var column names in e(b)/e(V) typically look like:
         "1.agegroup#c.tavg_poly_1_<prod>"
         "1.agegroup#c.tavg_poly_1_<prod>#c.lr_tavg_<prod>_adm1_avg"
         "1.agegroup#c.tavg_poly_1_<prod>#c.loggdppc_adm1_avg"
       To be robust, we’ll try both "c.tavg...#1.agegroup" and "1.agegroup#c.tavg..."
    ***********************************************************************/
    local coefnames

    forvalues a = 1/3 {
        forvalues k = 1/4 {
            /* Base (×1) */
            local nm1a "1.agegroup#c.tavg_poly_`k'_`prod_var'"
            local nm1b "c.tavg_poly_`k'_`prod_var'#1.agegroup"
            /* × LR */
            local nmla "1.agegroup#c.tavg_poly_`k'_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg"
            local nmlb "c.tavg_poly_`k'_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#1.agegroup"
            /* × log income */
            local ncia "1.agegroup#c.tavg_poly_`k'_`prod_var'#c.loggdppc_adm1_avg"
            local ncib "c.tavg_poly_`k'_`prod_var'#c.loggdppc_adm1_avg#1.agegroup"

            /* Replace '1.' with current age 'a.' in both alternatives */
            foreach TOKEN in nm1a nm1b nmla nmlb ncia ncib {
                local `TOKEN' = subinstr("``TOKEN''","1.agegroup","`a'.agegroup",.)
            }

            /* Resolve to existing column names by checking V colnames */
            local trylist " `nm1a' `nm1b' `nmla' `nmlb' `ncia' `ncib' "
            /* For each of the three slots (1, lr, loginc) pick the first that exists */
            /* Slot 1: base */
            local slot1 ""
            foreach cand in `nm1a' `nm1b' {
                local ix = colnumb(V, "`cand'")
                if (`ix'>0) {
                    local slot1 "`cand'"
                    continue, break
                }
            }
            if "`slot1'"=="" {
                di as err "!! Missing coef for age `a', poly `k' (base)."
                continue, break
            }

            /* Slot 2: ×LR */
            local slot2 ""
            foreach cand in `nmla' `nmlb' {
                local ix = colnumb(V, "`cand'")
                if (`ix'>0) {
                    local slot2 "`cand'"
                    continue, break
                }
            }
            if "`slot2'"=="" {
                di as err "!! Missing coef for age `a', poly `k' (×LR)."
                continue, break
            }

            /* Slot 3: ×loginc */
            local slot3 ""
            foreach cand in `ncia' `ncib' {
                local ix = colnumb(V, "`cand'")
                if (`ix'>0) {
                    local slot3 "`cand'"
                    continue, break
                }
            }
            if "`slot3'"=="" {
                di as err "!! Missing coef for age `a', poly `k' (×loginc)."
                continue, break
            }

            /* Append in desired order */
            local coefnames `"`coefnames' "`slot1'" "`slot2'" "`slot3'""'
        }
    }

    /* Sanity: we expect 36 names */
    local ncoef : word count `coefnames'
    if (`ncoef'!=36) {
        di as err "!! Expected 36 coefficients; found `ncoef'. Stopping for `prod_lbl'."
        continue
    }

    /************ 3) Open CSVV and write header/meta ************/
    local CSVVFILE "agespec_interaction_response_spec2_`prod_var'_noint.csvv"

    cap file close csvv
    file open csvv using "`CSVVDIR'/`CSVVFILE'", write replace

    file write csvv "---" _n
    file write csvv "oneline: Age-specific interacted (Spec 2) response – `prod_lbl'" _n
    file write csvv "version: MORTALITY-AGE-SPEC-INTERACTED-`prod_var'" _n
    file write csvv "dependencies: `STERFILE'" _n
    file write csvv "description: 4th-order poly OLS; age×ADM2 & age×Country×Year FE; temp×income & temp×LR-temp interactions" _n
    file write csvv "csvv-version: aliabonanno-`c(current_date)'" _n

    /************ 4) Variables (generic names as used downstream) ************/
    file write csvv "variables:" _n
    file write csvv "  tavg_poly_1: T1_sum [C]" _n
    file write csvv "  tavg_poly_2: T2_sum [C^2]" _n
    file write csvv "  tavg_poly_3: T3_sum [C^3]" _n
    file write csvv "  tavg_poly_4: T4_sum [C^4]" _n
    file write csvv "  lr_tavg: long-run mean annual temperature [C]" _n
    file write csvv "  deathrate_w99: mortality rate winsorized at 99th percentile" _n
    file write csvv "  loggdppc_adm1_avg: log GDP per capita" _n
    file write csvv "..." _n

    /************ 5) Observations & predictor names ************/
    file write csvv "observations" _n
    file write csvv " `e(N)'" _n

    file write csvv "prednames" _n
    file write csvv "tavg_poly_1, tavg_poly_2, tavg_poly_3, tavg_poly_4" _n

    file write csvv "covarnames" _n
    file write csvv "1, lr_tavg, loggdppc_adm1_avg" _n

    /************ 6) GAMMA (36 numbers, 12 lines x 3 per line) ************/
    file write csvv "gamma" _n
    local i 0
    foreach nm of local coefnames {
        local ++i
        /* write coefficient value */
        file write csvv "`= _b[`nm']'"
        /* comma separate triplets, newline every 3 */
        if mod(`i',3)!=0 {
            file write csvv ", "
        }
        else {
            file write csvv "" _n
        }
    }

    /************ 7) VCV in the same order (36x36) ************/
    file write csvv "gammavcv" _n
    forvalues r = 1/36 {
        /* name for row r */
        local rn : word `r' of `coefnames'
        /* column index for rn */
        local ir = colnumb(V, "`rn'")
        forvalues c = 1/36 {
            local cn : word `c' of `coefnames'
            local ic = colnumb(V, "`cn'")
            /* element */
            local vcf = el(V, `ir', `ic')
            if `c'==1 {
                file write csvv "`vcf'"
            }
            else {
                file write csvv ", `vcf'"
            }
        }
        file write csvv "" _n
    }

    /************ 8) Residual variance ************/
    file write csvv "residvcv" _n
    file write csvv "`= e(rmse)^2'" _n

    file close csvv
    di as res "Wrote `CSVVFILE' to `CSVVDIR'"

}  // end products loop

di as res "All CSVVs generated."
