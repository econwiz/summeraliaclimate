/**************************************************************************
* Age-specific × Product interacted regressions (ALL products pooled)
* - Spec 2: Temp quartic and interactions vary by AGE, with PRODUCT deltas
* - Includes Temp×(Income, LR) heterogeneity by AGE, with PRODUCT deltas
* - Precip controls as country×age FE-slopes, by PRODUCT
* - All key FE are PRODUCT-specific to prevent level mixing
* - Resolution: 1deg (all paths hardcoded for 1deg)
**************************************************************************/

clear all
set more off

* --- CONFIG -------------------------------------------------------------
global CODE    "/user/ab5405/summeraliaclimate/code"

local resolution "1deg"

global PANELS  "$CODE/mortality_pipeline/output_panels/`resolution'"
global OUTPUT  "$CODE/mortality_pipeline/output/`resolution'"

* Create output directory chain (one level at a time)
cap mkdir "$CODE/mortality_pipeline"
cap mkdir "$CODE/mortality_pipeline/output"
cap mkdir "$CODE/mortality_pipeline/output/`resolution'"

local inpanel "$PANELS/panel_prepped_for_regressions_ALLREA_aligned.dta"
local outdir  "$OUTPUT/age_spec_ALLREA"
cap mkdir "`outdir'"

local altspec 0   // 0 = only Spec 2; 1 = also run 1/3/4/5

* --- Load pooled (stacked) panel ---------------------------------------
di as txt "Loading ALLREA stacked panel → `inpanel'"
use "`inpanel'", clear

* --- Ensure factors exist in the right type ----------------------------
capture confirm numeric variable agegroup
if _rc {
    encode agegroup, gen(agegroup_)
    drop agegroup
    rename agegroup_ agegroup
}

capture confirm numeric variable product
if _rc {
    encode proj_base, gen(product)
}

* Rebuild FE IDs if missing (defensive)
capture confirm variable adm0_agegrp_code
if _rc {
    egen adm0_agegrp_code = group(adm0_code agegroup)
}
capture confirm variable adm1_agegrp_code
if _rc {
    egen adm1_agegrp_code = group(adm1_code agegroup)
}

* Ensure geographic IDs are numeric (required for absorb)
capture confirm numeric variable adm2_code
if _rc {
    encode adm2_code, gen(adm2_code_n)
    drop adm2_code
    rename adm2_code_n adm2_code
}

capture confirm numeric variable adm1_code
if _rc {
    encode adm1_code, gen(adm1_code_n)
    drop adm1_code
    rename adm1_code_n adm1_code
}

capture confirm numeric variable adm0_code
if _rc {
    encode adm0_code, gen(adm0_code_n)
    drop adm0_code
    rename adm0_code_n adm0_code
}

* --- Population weights: normalize within PRODUCT x YEAR ---------------
bysort product year: egen double tot_pop = total(population)
gen double weight = population / tot_pop

* =============================== SPEC 2 (Preferred) ==================== *
reghdfe deathrate_w99                                                          ///
    /* Base AGE x Temp quartic (anchors slopes) */                              ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup              ///
    /* PRODUCT deltas for AGE x Temp */                                         ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup#i.product    ///
    /* Base AGE x Temp x log income */                                          ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup ///
    /* PRODUCT deltas for AGE x Temp x log income */                            ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup#i.product ///
    /* Base AGE x Temp x LR */                                                  ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup  ///
    /* PRODUCT deltas for AGE x Temp x LR */                                    ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup#i.product   ///
    [pw = weight],                                                              ///
    absorb(                                                                     ///
        adm0_agegrp_code#i.product##c.prcp_poly_1                               ///
        adm0_agegrp_code#i.product##c.prcp_poly_2                               ///
        i.adm2_code#i.CHN_ts#i.agegroup#i.product                               ///
        i.adm0_code#i.year#i.agegroup#i.product                                 ///
    ) cluster(adm1_code) residual(e_hat)

estimates save "`outdir'/agespec_prod_interact_spec2_ALLREA_`resolution'.ster", replace
di as res "Saved Spec 2 -> `outdir'/agespec_prod_interact_spec2_ALLREA_`resolution'.ster"

* =========================== Alt specs toggle ========================== *
if (`altspec') {

    * ------------------------------ SPEC 1 ----------------------------- *
    reghdfe deathrate_w99                                                      ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup          ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup#i.product ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup#i.product ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup  ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup#i.product  ///
        [pw = weight],                                                          ///
        absorb(                                                                 ///
            adm0_agegrp_code#i.product##c.prcp_poly_1                           ///
            adm0_agegrp_code#i.product##c.prcp_poly_2                           ///
            i.adm2_code#i.CHN_ts#i.agegroup#i.product                           ///
            i.adm0_code#i.year#i.product                                        ///
        ) cluster(adm1_code)
    estimates save "`outdir'/agespec_prod_interact_spec1_ALLREA_`resolution'.ster", replace

    * ------------------------------ SPEC 3 ----------------------------- *
    reghdfe deathrate_w99                                                      ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup          ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup#i.product ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup#i.product ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup  ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup#i.product  ///
        [pw = weight],                                                          ///
        absorb(                                                                 ///
            adm0_agegrp_code#i.product##c.prcp_poly_1                           ///
            adm0_agegrp_code#i.product##c.prcp_poly_2                           ///
            i.adm2_code#i.CHN_ts#i.agegroup#i.product                           ///
            i.adm0_code#i.year#i.agegroup#i.product                              ///
            adm1_agegrp_code#i.product##c.year                                   ///
        ) cluster(adm1_code)
    estimates save "`outdir'/agespec_prod_interact_spec3_ALLREA_`resolution'.ster", replace

    * ------------------------------ SPEC 4 ----------------------------- *
    bysort adm1_code product: egen double omega = sd(e_hat) if e_hat < .
    gen double precisionweight = 1/(omega^2)

    sort iso adm1_id adm2_id agegroup year
    reghdfe deathrate_w99                                                      ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup          ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup#i.product ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup#i.product ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup  ///
        c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup#i.product  ///
        [pw = precisionweight],                                                 ///
        absorb(                                                                 ///
            adm0_agegrp_code#i.product##c.prcp_poly_1                           ///
            adm0_agegrp_code#i.product##c.prcp_poly_2                           ///
            i.adm2_code#i.CHN_ts#i.agegroup#i.product                           ///
            i.adm0_code#i.year#i.agegroup#i.product                              ///
        ) cluster(adm1_code) tol(1e-7)
    estimates save "`outdir'/agespec_prod_interact_spec4_ALLREA_`resolution'.ster", replace

    * ------------------------------ SPEC 5 ----------------------------- *
    capture confirm variable tavg_poly_1_13m
    if !_rc {
        reghdfe deathrate_w99                                                  ///
            c.(tavg_poly_1_13m tavg_poly_2_13m tavg_poly_3_13m tavg_poly_4_13m)#i.agegroup ///
            c.(tavg_poly_1_13m tavg_poly_2_13m tavg_poly_3_13m tavg_poly_4_13m)#i.agegroup#i.product ///
            c.(tavg_poly_1_13m tavg_poly_2_13m tavg_poly_3_13m tavg_poly_4_13m)#c.loggdppc_adm1_avg#i.agegroup ///
            c.(tavg_poly_1_13m tavg_poly_2_13m tavg_poly_3_13m tavg_poly_4_13m)#c.loggdppc_adm1_avg#i.agegroup#i.product ///
            c.(tavg_poly_1_13m tavg_poly_2_13m tavg_poly_3_13m tavg_poly_4_13m)#c.lr_tavg_adm1_avg#i.agegroup  ///
            c.(tavg_poly_1_13m tavg_poly_2_13m tavg_poly_3_13m tavg_poly_4_13m)#c.lr_tavg_adm1_avg#i.agegroup#i.product  ///
            [pw = weight],                                                      ///
            absorb(                                                             ///
                adm0_agegrp_code#i.product##c.prcp_poly_1                       ///
                adm0_agegrp_code#i.product##c.prcp_poly_2                       ///
                i.adm2_code#i.CHN_ts#i.agegroup#i.product                       ///
                i.adm0_code#i.year#i.agegroup#i.product                          ///
            ) cluster(adm1_code)
        estimates save "`outdir'/agespec_prod_interact_spec5_ALLREA_`resolution'.ster", replace
    }
}

di as txt "Done. Results saved in: `outdir'"
