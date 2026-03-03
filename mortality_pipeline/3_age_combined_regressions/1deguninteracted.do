/**************************************************************************
* prep_regressions_ALLREA_interacted_FLEXIBLE.do
* 
* FLEXIBLE RESOLUTION VERSION
*
* - Estimates temperature–mortality response with FULL product interactions
* - Input: stacked ALLREA panel in $OUTPUT (output_panels/{resolution})
* - Saves Spec1–Spec5 .ster files for ALLREA with product interactions
*
* CHANGE LOG:x
* - Added resolution parameter to control input/output directories
* - Reads from mortality_pipeline/output_panels/{resolution}/
* - Saves to mortality_pipeline/output/{resolution}/
**************************************************************************/

clear all
set more off

* ============================================================================
* CONFIG
* ============================================================================

local resolution "1deg"  /* ← CHANGE THIS: "025deg" or "1deg" */

global CODE    "/user/ab5405/summeraliaclimate/code"

* Input/output directories include resolution
global PANELS  "$CODE/mortality_pipeline/output_panels/`resolution'"
global OUTPUT  "$CODE/mortality_pipeline/output/`resolution'"

* Input panel and output directory for sters
local inpanel  "$PANELS/panel_prepped_for_regressions_ALLREA_aligned.dta"
local ster_out "$OUTPUT/age_combined_ALLREA"
capture mkdir "$OUTPUT"
capture mkdir "`ster_out'"

di ""
di "========================================================================"
di "CONFIGURATION:"
di "  Resolution: `resolution'"
di "  Input panel: `inpanel'"
di "  Output dir:  `ster_out'"
di "========================================================================"
di ""

* ============================================================================
* LOAD DATA
* ============================================================================

* Check if input exists
capture confirm file "`inpanel'"
if _rc {
    di as error ""
    di as error "ERROR: Input panel not found!"
    di as error "  `inpanel'"
    di as error ""
    di as error "You need to create the stacked ALLREA panel first."
    di as error "This combines all products (GMFD, ERA5, JRA_3Q, MERRA2)"
    di as error "into a single stacked dataset."
    exit 601
}

di as txt "Loading stacked ALLREA panel → `inpanel'"
use "`inpanel'", clear

* ============================================================================
* DATA PREP
* ============================================================================

* Ensure agegroup is numeric for factor notation
capture confirm numeric variable agegroup
if _rc {
    encode agegroup, gen(agegroup_num)
    drop agegroup
    rename agegroup_num agegroup
}

* Encode product label (proj_base) to a factor var "product"
capture confirm numeric variable product
if _rc {
    encode proj_base, gen(product)
}

* Weights normalized within product × year (mirrors per-product scheme)
bysort product year: egen tot_pop = total(population)
gen double weight = population / tot_pop

* ============================================================================
* SPECIFICATIONS
* ============================================================================

* All specs below:
* - Allow country-specific precip slopes to vary by product:
*     absorb: adm0_code#i.product##c.prcp_poly_1 and prcp_poly_2
* - Make key FE product-specific so levels don't mix across products.
* - Cluster at adm1_code (kept for continuity).

* ----------------------------------------------------------------------------
* SPEC 1: Baseline FE (country×year) and subnational FE, all BY PRODUCT
* ----------------------------------------------------------------------------

di ""
di "========================================================================"
di "ESTIMATING SPEC 1: Baseline FE by product"
di "========================================================================"
di ""

reghdfe deathrate_w99                                                         ///
    /* Temp quartic × product */                                              ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.product             ///
    [pw = weight],                                                            ///
    absorb(                                                                   ///
        /* Country precip slopes by product */                                ///
        adm0_code#i.product##c.prcp_poly_1                                    ///
        adm0_code#i.product##c.prcp_poly_2                                    ///
        /* Subnational FE by product */                                       ///
        i.adm2_code#i.CHN_ts#i.agegroup#i.product                             ///
        /* Country × year FE by product */                                    ///
        i.adm0_code#i.year#i.product                                          ///
    )                                                                         ///
    cluster(adm1_code)

estimates save "`ster_out'/pooled_response_spec1_public_ALLREA.ster", replace
di " Saved Spec 1"

* ----------------------------------------------------------------------------
* SPEC 2 (PREFERRED): country×year×agegroup FE by product
* ----------------------------------------------------------------------------

di ""
di "========================================================================"
di "ESTIMATING SPEC 2: Country×year×agegroup FE by product (PREFERRED)"
di "========================================================================"
di ""

reghdfe deathrate_w99                                                         ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.product             ///
    [pw = weight],                                                            ///
    absorb(                                                                   ///
        adm0_code#i.product##c.prcp_poly_1                                    ///
        adm0_code#i.product##c.prcp_poly_2                                    ///
        i.adm2_code#i.CHN_ts#i.agegroup#i.product                             ///
        i.adm0_code#i.year#i.agegroup#i.product                               ///
    )                                                                         ///
    cluster(adm1_code) residual(e_hat)

estimates save "`ster_out'/pooled_response_spec2_public_ALLREA.ster", replace
di " Saved Spec 2"

* ----------------------------------------------------------------------------
* SPEC 3: Adds adm1×agegroup linear time slope BY PRODUCT
* ----------------------------------------------------------------------------

di ""
di "========================================================================"
di "ESTIMATING SPEC 3: + ADM1×agegroup time trends"
di "========================================================================"
di ""

reghdfe deathrate_w99                                                         ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.product             ///
    [pw = weight],                                                            ///
    absorb(                                                                   ///
        adm0_code#i.product##c.prcp_poly_1                                    ///
        adm0_code#i.product##c.prcp_poly_2                                    ///
        i.adm2_code#i.CHN_ts#i.agegroup#i.product                             ///
        i.adm0_code#i.year#i.agegroup#i.product                               ///
        /* Adm1×agegroup linear time slope by product */                      ///
        adm1_agegrp_code#i.product##c.year                                    ///
    )                                                                         ///
    cluster(adm1_code)

estimates save "`ster_out'/pooled_response_spec3_public_ALLREA.ster", replace
di " Saved Spec 3"

* ----------------------------------------------------------------------------
* SPEC 4: FGLS-style using product-specific precision weights
* ----------------------------------------------------------------------------

di ""
di "========================================================================"
di "ESTIMATING SPEC 4: FGLS precision weights"
di "========================================================================"
di ""

gen double e2 = e_hat^2
bysort adm1_code product: egen double omega = sd(e_hat) if e_hat < .
gen double precisionweight = weight * 1/(omega^2)

sort iso adm1_id adm2_id year agegroup

reghdfe deathrate_w99                                                         ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.product             ///
    [pw = precisionweight],                                                   ///
    absorb(                                                                   ///
        adm0_code#i.product##c.prcp_poly_1                                    ///
        adm0_code#i.product##c.prcp_poly_2                                    ///
        i.adm2_code#i.CHN_ts#i.agegroup#i.product                             ///
        i.adm0_code#i.year#i.agegroup#i.product                               ///
    )                                                                         ///
    cluster(adm1_code) tol(1e-7)

estimates save "`ster_out'/pooled_response_spec4_public_ALLREA.ster", replace
di " Saved Spec 4"

* ----------------------------------------------------------------------------
* SPEC 5: Full interactions (temp×income, temp×LR-mean, all × product)
* ----------------------------------------------------------------------------

di ""
di "========================================================================"
di "ESTIMATING SPEC 5: Full interactions (temp×income×LR-mean)"
di "========================================================================"
di ""

reghdfe deathrate_w99                                                         ///
    /* 1) Temp quartic × product */                                           ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.product             ///
    /* 2) Temp × income × product */                                          ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.product ///
    /* 3) Temp × LR-mean-temp × product */                                    ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.product  ///
    [pw = weight],                                                            ///
    absorb(                                                                   ///
        /* Country precip slopes by product */                                ///
        adm0_code#i.product##c.prcp_poly_1                                    ///
        adm0_code#i.product##c.prcp_poly_2                                    ///
        /* Subnational FE by product */                                       ///
        i.adm2_code#i.CHN_ts#i.agegroup#i.product                             ///
        /* Country × year FE by product */                                    ///
        i.adm0_code#i.year#i.product                                          ///
    )                                                                         ///
    cluster(adm1_code)

estimates save "`ster_out'/pooled_response_spec5_public_ALLREA.ster", replace
di " Saved Spec 5"

* ============================================================================
* SUMMARY
* ============================================================================

di ""
di "========================================================================"
di "ALL SPECS COMPLETE!"
di "========================================================================"
di ""
di "All ALLREA product-interacted specs (1–5) estimated and saved to:"
di "  `ster_out'"
di ""
di "Resolution: `resolution'"
di ""
di "Files created:"
local specs "spec1 spec2 spec3 spec4 spec5"
foreach spec of local specs {
    local file "`ster_out'/pooled_response_`spec'_public_ALLREA.ster"
    capture confirm file "`file'"
    if !_rc {
        di "   pooled_response_`spec'_public_ALLREA.ster"
    }
    else {
        di "   pooled_response_`spec'_public_ALLREA.ster (MISSING!)"
    }
}
di ""
di "========================================================================"
di ""
