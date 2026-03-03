/**************************************************************************
* stack_all_products_to_output_panels_FLEXIBLE.do
*
* FLEXIBLE RESOLUTION VERSION
*
* - Stacks per-product aligned panels into ONE pooled dataset (ALLREA)
* - Saves to $OUTPUT (your output_panels folder)
*
* CHANGE LOG:
* - Added resolution parameter
* - Fixed input path to read from OUTPUT (where prep script saves)
* - Added products list definition
**************************************************************************/

clear all
set more off

* ============================================================================
* CONFIG
* ============================================================================

local resolution "1deg"  /* ← CHANGE THIS: "025deg" or "1deg" */

global CODE    "/user/ab5405/summeraliaclimate/code"

* Input/output directories include resolution
global PANELS  "$CODE/mortality_pipeline/prep_panels/`resolution'"
global OUTPUT  "$CODE/mortality_pipeline/output_panels/`resolution'"
capture mkdir "$OUTPUT"

* Products to stack
local products GMFD ERA5_025 JRA_3Q MERRA2

di ""
di "========================================================================"
di "CONFIGURATION:"
di "  Resolution: `resolution'"
di "  Products: `products'"
di "  Input dir:  $OUTPUT (individual product panels)"
di "  Output: panel_prepped_for_regressions_ALLREA_aligned.dta"
di "========================================================================"
di ""

* ============================================================================
* STACK PRODUCTS
* ============================================================================

tempfile base
local first = 1

foreach product of local products {
    
    * Input comes from prep (where merge_panel_ALL_aligned saves)
    local fpath "$PANELS/global_mortality_panel_public_`product'.dta"
    
    capture confirm file "`fpath'"
    if _rc {
        di as yellow ""
        di as yellow "  File not found for `product':"
        di as yellow "    `fpath'"
        di as yellow "    (skipping)"
        di as yellow ""
        continue
    }
    
    di ""
    di "--------------------------------------------------------------------"
    di "Loading `product'"
    di "  → `fpath'"
    di "--------------------------------------------------------------------"
    
    use "`fpath'", clear
    
    * Add product label
    gen str20 proj_base = "`product'"
    
    * ========================================================================
    * Standardize product-specific variable names to generic names
    * ========================================================================
    
    * Temperature polynomials
    capture confirm variable tavg_poly_1_`product'
    if !_rc {
        rename tavg_poly_1_`product' tavg_poly_1
        capture rename tavg_poly_2_`product' tavg_poly_2
        capture rename tavg_poly_3_`product' tavg_poly_3
        capture rename tavg_poly_4_`product' tavg_poly_4
    }
    
    * Precipitation polynomials
    capture confirm variable prcp_poly_1_`product'
    if !_rc {
        rename prcp_poly_1_`product' prcp_poly_1
        capture rename prcp_poly_2_`product' prcp_poly_2
    }
    
    * Long-run temperature average
    capture confirm variable lr_tavg_`product'_adm1_avg
    if !_rc {
        rename lr_tavg_`product'_adm1_avg lr_tavg_adm1_avg
    }
    
    * ========================================================================
    * Keep only variables needed downstream
    * ========================================================================
    
    keep iso iso2 adm0_code adm0_code2 adm1_id adm1_code adm2_id adm2_code ///
         year agegroup CHN_ts deaths population deathrate_w99 ///
         tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
         prcp_poly_1 prcp_poly_2 loggdppc_adm1_avg lr_tavg_adm1_avg ///
         adm0_agegrp_code adm1_agegrp_code ///
         proj_base
    
    order iso iso2 adm0_code adm0_code2 adm1_id adm1_code adm2_id adm2_code ///
          year agegroup CHN_ts deaths population deathrate_w99 ///
          tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
          prcp_poly_1 prcp_poly_2 loggdppc_adm1_avg lr_tavg_adm1_avg ///
          proj_base
    
    di "  Observations: " _N
    
    * ========================================================================
    * Append to base
    * ========================================================================
    
    if `first' {
        save `base'
        local first = 0
    }
    else {
        append using `base'
        save `base', replace
    }
}

* ============================================================================
* CHECK AND SAVE
* ============================================================================

capture confirm file `base'
if _rc {
    di as error ""
    di as error "❌ ERROR: No product files were found!"
    di as error "   Location: $OUTPUT"
    di as error ""
    di as error "Make sure you ran prep_panel_ALL_aligned_FLEXIBLE.do first!"
    exit 1
}

use `base', clear
sort proj_base iso adm1_id adm2_id year agegroup

* ============================================================================
* SUMMARY
* ============================================================================

di ""
di "========================================================================"
di "STACKING SUMMARY"
di "========================================================================"

count
di "Total observations: " r(N)

* Show breakdown by product
di ""
di "Breakdown by product:"
tab proj_base

* Save stacked panel
local out_dta "$OUTPUT/panel_prepped_for_regressions_ALLREA_aligned.dta"
local out_csv "$OUTPUT/panel_prepped_for_regressions_ALLREA_aligned.csv"

save "`out_dta'", replace
export delimited using "`out_csv'", replace

di ""
di "========================================================================"
di "SUCCESS!"
di "========================================================================"
di ""
di "Wrote pooled ALLREA panel to:"
di "   `out_dta'"
di "   `out_csv'"
di ""
di "Resolution: `resolution'"
di "Location: $OUTPUT"
di ""
di "========================================================================"
di ""
