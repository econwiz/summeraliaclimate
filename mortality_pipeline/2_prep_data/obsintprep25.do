/**************************************************************************
* stack_all_products_to_output_panels.do
* - Stacks per-product aligned panels into ONE pooled dataset (ALLREA)
* - Saves to $OUTPUT (your output_panels folder)
**************************************************************************/

clear all
set more off

global CODE    "/user/ab5405/summeraliaclimate/code"

* Input/output directories include resolution
global PANELS  "$CODE/mortality_pipeline/prep_panels"
global OUTPUT  "$CODE/mortality_pipeline/output_panels"
capture mkdir "$OUTPUT"


capture mkdir "$OUTPUT"

tempfile base
local first = 1

foreach product of local products {
    local fpath "$PANELS/panel_prepped_for_regressions_`product'_aligned.dta"

    capture confirm file "`fpath'"
    if _rc {
        di as yellow "Note: file not found for `product' → `fpath' (skipping)"
        continue
    }

    di as txt "Loading `product' → `fpath'"
    use "`fpath'", clear

    * Add product label
    gen str20 proj_base = "`product'"

    * ---- Standardize product-specific names to generic ----
    capture confirm variable tavg_poly_1_`product'
    if !_rc {
        rename tavg_poly_1_`product' tavg_poly_1
        capture rename tavg_poly_2_`product' tavg_poly_2
        capture rename tavg_poly_3_`product' tavg_poly_3
        capture rename tavg_poly_4_`product' tavg_poly_4
    }

    capture confirm variable prcp_poly_1_`product'
    if !_rc {
        rename prcp_poly_1_`product' prcp_poly_1
        capture rename prcp_poly_2_`product' prcp_poly_2
    }

    capture confirm variable lr_tavg_`product'_adm1_avg
    if !_rc {
        rename lr_tavg_`product'_adm1_avg lr_tavg_adm1_avg
    }

    * Keep the vars needed downstream
    keep iso iso2 adm0_code adm0_code2 adm1_id adm1_code adm2_id adm2_code ///
     year agegroup CHN_ts deaths population deathrate_w99 ///
     tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
     prcp_poly_1 prcp_poly_2 loggdppc_adm1_avg lr_tavg_adm1_avg ///
     /* keep FE IDs from prep */ ///
     adm0_agegrp_code adm1_agegrp_code ///
     proj_base

    order iso iso2 adm0_code adm0_code2 adm1_id adm1_code adm2_id adm2_code ///
         year agegroup CHN_ts deaths population deathrate_w99 ///
         tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
         prcp_poly_1 prcp_poly_2 loggdppc_adm1_avg lr_tavg_adm1_avg ///
         proj_base

    if `first' {
        save `base'
        local first = 0
    }
    else {
        append using `base'
        save `base', replace
    }
}

capture confirm file `base'
if _rc {
    di as error "No product files were found in $OUTPUT. Nothing to stack."
    exit 1
}

use `base', clear
sort proj_base iso adm1_id adm2_id year agegroup
count
di as result "Stacked observations: " r(N)

* Save the stacked panel inside your output_panels folder
local out_dta "$OUTPUT/panel_prepped_for_regressions_ALLREA_aligned.dta"
local out_csv "$OUTPUT/panel_prepped_for_regressions_ALLREA_aligned.csv"

save "`out_dta'", replace
export delimited using "`out_csv'", replace

di as text "Wrote pooled ALLREA to:"
di as result "  `out_dta'"
di as result "  `out_csv'"
