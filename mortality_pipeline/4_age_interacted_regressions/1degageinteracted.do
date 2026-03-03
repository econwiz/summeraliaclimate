/**************************************************************************
* Age-specific × Product interacted regressions (ALL products pooled)
* - Spec 2: Temp quartic and interactions vary by AGE, with PRODUCT deltas
* - Includes Temp×(Income, LR) heterogeneity by AGE, with PRODUCT deltas
* - Precip controls as country×age FE-slopes, by PRODUCT
* - All key FE are PRODUCT-specific to prevent level mixing
**************************************************************************/

clear all
set more off

* --- CONFIG -------------------------------------------------------------
global CODE    "/user/ab5405/summeraliaclimate/code"

local resolution "1deg"   /* ← CHANGE THIS to match stacking script: "025deg" or "1deg" */

* Input path matches where stack_all_products_to_output_panels_FLEXIBLE.do saves
global PANELS  "$CODE/mortality_pipeline/output_panels/`resolution'"
global OUTPUT  "$CODE/mortality_pipeline/output/1deg"

local inpanel "$PANELS/panel_prepped_for_regressions_ALLREA_aligned.dta"
local outdir  "$OUTPUT/age_spec_ALLREA"
capture mkdir "`outdir'"

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

* --- Population weights: normalize within PRODUCT × YEAR ---------------
bysort product year: egen double tot_pop = total(population)
gen double weight = population / tot_pop

* =============================== SPEC 2 (Preferred) ==================== *
* We explicitly include:
*   Base AGE×Temp quartic                     (no product)
*   + PRODUCT deltas for AGE×Temp             (age×temp×product)
*   Base AGE×Temp×log income                  (no product)
*   + PRODUCT deltas for AGE×Temp×log income  (…×product)
*   Base AGE×Temp×LR                          (no product)
*   + PRODUCT deltas for AGE×Temp×LR          (…×product)
* FE (all PRODUCT-specific):
*   - Country×AGE FE-slopes for precip: adm0_agegrp_code#i.product ## c.prcp_poly_1/2
*   - Subnational FE: i.adm2_code # i.CHN_ts # i.agegroup # i.product
*   - Country×Year×AGE FE: i.adm0_code # i.year # i.agegroup # i.product
reghdfe deathrate_w99                                                          ///
    /* Base AGE×Temp quartic (anchors slopes) */                                ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup              ///
    /* PRODUCT deltas for AGE×Temp */                                           ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.agegroup#i.product    ///
    /* Base AGE×Temp×log income */                                              ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup ///
    /* PRODUCT deltas for AGE×Temp×log income */                                ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.loggdppc_adm1_avg#i.agegroup#i.product ///
    /* Base AGE×Temp×LR */                                                      ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup  ///
    /* PRODUCT deltas for AGE×Temp×LR */                                        ///
    c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#c.lr_tavg_adm1_avg#i.agegroup#i.product   ///
    [pw = weight],                                                              ///
    absorb(                                                                     ///
        /* country×age precip slopes BY PRODUCT */                              ///
        adm0_agegrp_code#i.product##c.prcp_poly_1                               ///
        adm0_agegrp_code#i.product##c.prcp_poly_2                               ///
        /* subnational FE BY PRODUCT (China break) */                           ///
        i.adm2_code#i.CHN_ts#i.agegroup#i.product                               ///
        /* country×year×age FE BY PRODUCT */                                    ///
        i.adm0_code#i.year#i.agegroup#i.product                                 ///
    ) cluster(adm1_code) residual(e_hat)

estimates save "`outdir'/agespec_prod_interact_spec2_ALLREA.ster", replace
di as res "✅ Saved Spec 2 → `outdir'/agespec_prod_interact_spec2_ALLREA.ster"

* =========================== Alt specs toggle ========================== *
if (`altspec') {

    * ------------------------------ SPEC 1 ----------------------------- *
    * Like Spec 2 but with country×year FE (no age on that FE)
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
    estimates save "`outdir'/agespec_prod_interact_spec1_ALLREA.ster", replace

    * ------------------------------ SPEC 3 ----------------------------- *
    * Add ADM1×AGE linear time trend BY PRODUCT
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
    estimates save "`outdir'/agespec_prod_interact_spec3_ALLREA.ster", replace

    * ------------------------------ SPEC 4 ----------------------------- *
    * FGLS-style using precision weights within ADM1 × PRODUCT
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
    estimates save "`outdir'/agespec_prod_interact_spec4_ALLREA.ster", replace

    * ------------------------------ SPEC 5 ----------------------------- *
    * If 13-month exposure polynomials exist in the stacked file
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
        estimates save "`outdir'/agespec_prod_interact_spec5_ALLREA.ster", replace
    }
}

di as txt "Done. Results saved in: `outdir'"
