/**************************************************************************
* CSVV GENERATOR — ALLREA pooled ster → per-product CSVVs
* 1deg version
*
* Product encoding: 1=ERA5_025, 2=GMFD, 3=JRA_3Q, 4=MERRA2
* Base category (product=1b) = ERA5_025
*
* 36 gammas per product: 3 ages × 4 polys × 3 covariates {1, LR, loginc}
* Age codes: 1b=0-4, 2=5-64, 3=65+, 4=total (skip total in CSVV)
**************************************************************************/

clear all
set more off

global REPO "/user/ab5405/summeraliaclimate"

local STERFILE "$REPO/code/mortality_pipeline/output/1deg/age_spec_ALLREA/agespec_prod_interact_spec2_ALLREA.ster"
local CSVVDIR  "$REPO/code/mortality_pipeline/output/1deg/age_spec_interacted_csvv"
cap mkdir "`CSVVDIR'"

local prod_codes "1 2 3 4"
local prod_names "ERA5_025 GMFD JRA_3Q MERRA2"
local age_codes  "1 2 3"

estimates use "`STERFILE'"
matrix B = e(b)
matrix V = e(V)

di as res "Loaded: `STERFILE'"
di "N = " e(N)

forvalues pi = 1/4 {

    local prod_name : word `pi' of `prod_names'
    local prod_code : word `pi' of `prod_codes'

    di ""
    di "============================================================"
    di "Building CSVV for product: `prod_name' (code=`prod_code')"
    di "============================================================"

    local ncoef 36
    matrix GAMMA = J(1, `ncoef', 0)
    matrix GVCV  = J(`ncoef', `ncoef', 0)

    local idx 0

    foreach a of local age_codes {
        forvalues k = 1/4 {
            foreach slot in base lr loginc {

                local ++idx

                if "`slot'" == "base" {
                    if `a' == 1 local base_nm "`a'b.agegroup#c.tavg_poly_`k'"
                    else        local base_nm "`a'.agegroup#c.tavg_poly_`k'"
                }
                else if "`slot'" == "lr" {
                    if `a' == 1 local base_nm "`a'b.agegroup#c.tavg_poly_`k'#c.lr_tavg_adm1_avg"
                    else        local base_nm "`a'.agegroup#c.tavg_poly_`k'#c.lr_tavg_adm1_avg"
                }
                else {
                    if `a' == 1 local base_nm "`a'b.agegroup#c.tavg_poly_`k'#c.loggdppc_adm1_avg"
                    else        local base_nm "`a'.agegroup#c.tavg_poly_`k'#c.loggdppc_adm1_avg"
                }

                if `pi' == 1 {
                    local has_delta 0
                    local id 0
                }
                else {
                    local has_delta 1
                    if "`slot'" == "base" {
                        if `a' == 1 local delta_nm "`a'b.agegroup#`prod_code'.product#c.tavg_poly_`k'"
                        else        local delta_nm "`a'.agegroup#`prod_code'.product#c.tavg_poly_`k'"
                    }
                    else if "`slot'" == "lr" {
                        if `a' == 1 local delta_nm "`a'b.agegroup#`prod_code'.product#c.tavg_poly_`k'#c.lr_tavg_adm1_avg"
                        else        local delta_nm "`a'.agegroup#`prod_code'.product#c.tavg_poly_`k'#c.lr_tavg_adm1_avg"
                    }
                    else {
                        if `a' == 1 local delta_nm "`a'b.agegroup#`prod_code'.product#c.tavg_poly_`k'#c.loggdppc_adm1_avg"
                        else        local delta_nm "`a'.agegroup#`prod_code'.product#c.tavg_poly_`k'#c.loggdppc_adm1_avg"
                    }
                }

                local ib = colnumb(B, "`base_nm'")
                if `ib' == 0 di as err "!! base coef not found: `base_nm'"

                if `has_delta' {
                    local id = colnumb(B, "`delta_nm'")
                    if `id' == 0 di as err "!! delta coef not found: `delta_nm'"
                    matrix GAMMA[1,`idx'] = el(B,1,`ib') + el(B,1,`id')
                }
                else {
                    matrix GAMMA[1,`idx'] = el(B,1,`ib')
                }

                local base_idx_`idx'  = `ib'
                local delta_idx_`idx' = `id'
                local has_delta_`idx' = `has_delta'

            }
        }
    }

    forvalues r = 1/`ncoef' {
        forvalues c = 1/`ncoef' {
            local ib_r = `base_idx_`r''
            local ib_c = `base_idx_`c''
            local id_r = `delta_idx_`r''
            local id_c = `delta_idx_`c''
            local hd_r = `has_delta_`r''
            local hd_c = `has_delta_`c''

            local vcv_val = el(V, `ib_r', `ib_c')
            if `hd_r' & `hd_c'  local vcv_val = `vcv_val' + el(V,`id_r',`id_c') + el(V,`ib_r',`id_c') + el(V,`id_r',`ib_c')
            else if `hd_r'      local vcv_val = `vcv_val' + el(V,`id_r',`ib_c')
            else if `hd_c'      local vcv_val = `vcv_val' + el(V,`ib_r',`id_c')
            matrix GVCV[`r',`c'] = `vcv_val'
        }
    }

    local CSVVFILE "agespec_interaction_response_spec2_`prod_name'_prodint.csvv"
    cap file close csvv
    file open csvv using "`CSVVDIR'/`CSVVFILE'", write replace

    file write csvv "---" _n
    file write csvv "oneline: Age-specific interacted (Spec 2) ALLREA 1deg response – `prod_name'" _n
    file write csvv "version: MORTALITY-AGE-SPEC-ALLREA-1DEG-`prod_name'" _n
    file write csvv "dependencies: agespec_prod_interact_spec2_ALLREA.ster" _n
    file write csvv "description: 4th-order poly OLS pooled ALLREA 1deg; age×product interactions; per-product coefficients reconstructed from base+delta" _n
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

    file write csvv "gamma" _n
    forvalues i = 1/`ncoef' {
        file write csvv "`=el(GAMMA,1,`i')'"
        if mod(`i',3) != 0 file write csvv ", "
        else               file write csvv "" _n
    }

    file write csvv "gammavcv" _n
    forvalues r = 1/`ncoef' {
        forvalues c = 1/`ncoef' {
            if `c'==1 file write csvv "`=el(GVCV,`r',`c')'"
            else      file write csvv ", `=el(GVCV,`r',`c')'"
        }
        file write csvv "" _n
    }

    file write csvv "residvcv" _n
    file write csvv "`=e(rmse)^2'" _n

    file close csvv
    di as res "✅ Wrote `CSVVFILE'"
}

di ""
di as res "All 4 per-product CSVVs written to `CSVVDIR'"
