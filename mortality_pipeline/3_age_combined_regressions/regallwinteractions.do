/**************************************************************************
* prep_regressions_ALLREA_interacted.do
* - Estimates temperature–mortality response with FULL product interactions
* - Input: stacked ALLREA panel in $OUTPUT (output_panels)
* - Saves Spec1–Spec5 .ster files for ALLREA with product interactions
**************************************************************************/

clear all
set more off

* CONFIG
global CODE    "/user/ab5405/summeraliaclimate/code"
global PANELS  "$CODE/regressions/output_panels"
global OUTPUT  "$CODE/regressions/output"

* INPUT / OUTPUT
local inpanel  "$PANELS/panel_prepped_for_regressions_ALLREA_aligned.dta"
local ster_out "$OUTPUT/age_combined_ALLREA"
capture mkdir "`ster_out'"

di as txt "Loading stacked ALLREA panel → `inpanel'"
use "`inpanel'", clear

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

* 1) Weights normalized within product × year (mirrors per-product scheme)
bysort product year: egen tot_pop = total(population)
gen double weight = population / tot_pop

* ----------------------------- SPECS ----------------------------------- *
* All specs below:
* - Allow country-specific precip slopes to vary by product:
*     absorb: adm0_code#i.product##c.prcp_poly_1 and prcp_poly_2
* - Make key FE product-specific so levels don’t mix across products.
* - Cluster at adm1_code (kept for continuity).
* ----------------------------------------------------------------------- *

/*** Spec 1: Baseline FE (country×year) and subnational FE (adm2×CHN_ts×agegroup), all BY PRODUCT ***/
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

/*** Spec 2 (preferred): country×year×agegroup FE by product + subnational FE by product ***/
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

/*** Spec 3: adds adm1×agegroup linear time (year) slope BY PRODUCT ***/
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

/*** Spec 4: FGLS-style using product-specific precision weights (built from Spec 2 residuals) ***/
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

/*** Spec 5: “Full interactions” — temp quartic, temp×income, temp×LR-mean, ALL × product ***/
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

di as result "All ALLREA product-interacted specs (1–5) estimated and saved to:"
di as txt    "  `ster_out'"
