/**************************************************************************
* CSVV GENERATOR — Age-Combined ALLREA pooled ster → per-product CSVVs
* 1deg version
*
* Product encoding: 1=ERA5_025, 2=GMFD, 3=JRA_3Q, 4=MERRA2
* Base category (product=1b) = ERA5_025
*
* 4 gammas per product: polys 1..4 (age-combined, no covariate interactions)
**************************************************************************/

clear all
set more off

global REPO "/user/ab5405/summeraliaclimate"

local STERDIR "$REPO/code/mortality_pipeline/output/1deg/age_combined_ALLREA"
local CSVVDIR "$REPO/code/mortality_pipeline/output/1deg/age_combined_csvv_ALLREA"
cap mkdir "`CSVVDIR'"

local prod_codes "1 2 3 4"
local prod_names "ERA5_025 GMFD JRA_3Q MERRA2"

forvalues SPEC = 1/5 {

    local STERFILE "pooled_response_spec`SPEC'_public_ALLREA.ster"

    capture estimates use "`STERDIR'/`STERFILE'"
    if _rc {
        di as yellow "Note: couldn't find `STERDIR'/`STERFILE' — skipping Spec `SPEC'."
        continue
    }

    matrix B = e(b)
    matrix V = e(V)

    di ""
    di "========================================================"
    di "Spec `SPEC' — generating per-product CSVVs"
    di "========================================================"

    forvalues pi = 1/4 {

        local prod_name : word `pi' of `prod_names'
        local prod_code : word `pi' of `prod_codes'

        local ncoef 4
        matrix GAMMA = J(1, `ncoef', 0)
        matrix GVCV  = J(`ncoef', `ncoef', 0)

        forvalues k = 1/4 {

            local base_nm "1b.product#c.tavg_poly_`k'"
            local ib = colnumb(B, "`base_nm'")
            if `ib' == 0 di as err "!! base coef not found: `base_nm'"

            if `pi' == 1 {
                local has_delta 0
                local id 0
                matrix GAMMA[1,`k'] = el(B,1,`ib')
            }
            else {
                local has_delta 1
                local delta_nm "`prod_code'.product#c.tavg_poly_`k'"
                local id = colnumb(B, "`delta_nm'")
                if `id' == 0 di as err "!! delta coef not found: `delta_nm'"
                matrix GAMMA[1,`k'] = el(B,1,`ib') + el(B,1,`id')
            }

            local base_idx_`k'  = `ib'
            local delta_idx_`k' = `id'
            local has_delta_`k' = `has_delta'
        }

        forvalues r = 1/4 {
            forvalues c = 1/4 {
                local ib_r = `base_idx_`r''
                local ib_c = `base_idx_`c''
                local id_r = `delta_idx_`r''
                local id_c = `delta_idx_`c''
                local hd_r = `has_delta_`r''
                local hd_c = `has_delta_`c''

                local vcv_val = el(V,`ib_r',`ib_c')
                if `hd_r' & `hd_c'  local vcv_val = `vcv_val' + el(V,`id_r',`id_c') + el(V,`ib_r',`id_c') + el(V,`id_r',`ib_c')
                else if `hd_r'      local vcv_val = `vcv_val' + el(V,`id_r',`ib_c')
                else if `hd_c'      local vcv_val = `vcv_val' + el(V,`ib_r',`id_c')
                matrix GVCV[`r',`c'] = `vcv_val'
            }
        }

        local CSVVFILE "pooled_response_spec`SPEC'_`prod_name'_prodint.csvv"
        cap file close csvv
        file open csvv using "`CSVVDIR'/`CSVVFILE'", write replace

        file write csvv "---" _n
        file write csvv "oneline: Age-combined ALLREA 1deg (Spec `SPEC') response – `prod_name'" _n
        file write csvv "version: MORTALITY-AGE-COMBINED-ALLREA-1DEG-SPEC`SPEC'-`prod_name'" _n
        file write csvv "dependencies: `STERFILE'" _n
        file write csvv "description: 4th-order poly OLS pooled ALLREA 1deg age-combined; per-product coefficients reconstructed from base+delta" _n
        file write csvv "csvv-version: auto-`c(current_date)'" _n

        file write csvv "variables:" _n
        file write csvv "  tavg_poly_1: T1_sum [C]" _n
        file write csvv "  tavg_poly_2: T2_sum [C^2]" _n
        file write csvv "  tavg_poly_3: T3_sum [C^3]" _n
        file write csvv "  tavg_poly_4: T4_sum [C^4]" _n
        file write csvv "  deathrate_w99: mortality rate winsorized at 99th percentile" _n
        file write csvv "..." _n

        file write csvv "observations" _n
        file write csvv "`=e(N)'" _n

        file write csvv "prednames" _n
        file write csvv "tavg_poly_1, tavg_poly_2, tavg_poly_3, tavg_poly_4" _n

        file write csvv "covarnames" _n
        file write csvv "1" _n

        file write csvv "gamma" _n
        forvalues k = 1/4 {
            if `k'==1 file write csvv "`=el(GAMMA,1,`k')'"
            else      file write csvv ", `=el(GAMMA,1,`k')'"
        }
        file write csvv "" _n

        file write csvv "gammavcv" _n
        forvalues r = 1/4 {
            forvalues c = 1/4 {
                if `c'==1 file write csvv "`=el(GVCV,`r',`c')'"
                else      file write csvv ", `=el(GVCV,`r',`c')'"
            }
            file write csvv "" _n
        }

        file write csvv "residvcv" _n
        file write csvv "`=e(rmse)^2'" _n

        file close csvv
        di as res "  ✅ Wrote `CSVVFILE'"
    }
}

di ""
di as res "All 1deg age-combined per-product CSVVs written to `CSVVDIR'"
