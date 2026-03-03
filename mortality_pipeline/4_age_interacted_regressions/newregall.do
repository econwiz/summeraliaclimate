version 17.0
clear all
set more off

/**************************************************************************
* Age-specific × Product regressions (per-product runs)
* Spec 2 (preferred): Temp quartic × AGE
*      + Temp×log-income × AGE
*      + Temp×LR-mean × AGE
* FE: (AGE×ADM2) & (AGE×Country×Year)
* Precip controls: (AGE×ADM0×prcp_poly_1/2)
* Weights: population normalized within YEAR (per-product panel)
* Cluster: ADM1
**************************************************************************/

* ---- ABSOLUTE PATHS ----
global ROOT   "/user/ab5405/summeraliaclimate"
global CODE   "$ROOT/code/regressions"
global PANELS "$CODE/output_panels"
global OUTPUT "$CODE/output"

* sanity print
di as txt "PANELS = $PANELS"
di as txt "OUTPUT = $OUTPUT"

* Products (use underscore form for variable stems)
local prods "ERA5_025 GMFD MERRA2 JRA_3Q"

* Toggle to also run 1/3/4/5
local altspec 0

* Where .ster files go
local sterdir "$OUTPUT/age_spec_interacted"
cap mkdir "`sterdir'"

foreach P of local prods {

    di as txt "=============================================================="
    di as txt "Running age-specific interacted regressions for: `P'"
    di as txt "=============================================================="

    * File basenames to try (your folder shows both aligned & plain files)
    local try1 "$PANELS/panel_prepped_for_regressions_`P'_aligned.dta"
    local try2 "$PANELS/panel_prepped_for_regressions_`P'.dta"

    capture confirm file "`try1'"
    if _rc {
        capture confirm file "`try2'"
        if _rc {
            di as err "!! Missing panel for `P'. Tried:"
            di as err "   `try1'"
            di as err "   `try2'"
            continue
        }
        use "`try2'", clear
        di as txt "Loaded: `try2'"
    }
    else {
        use "`try1'", clear
        di as txt "Loaded: `try1'"
    }

    * Ensure agegroup is numeric for factor notation
    capture confirm numeric variable agegroup
    if _rc {
        encode agegroup, gen(agegroup_)
        drop agegroup
        rename agegroup_ agegroup
    }

    * Required variables check (suffix uses `P`, e.g., ERA5_025 / JRA_3Q)
    foreach v in ///
        tavg_poly_1_`P' tavg_poly_2_`P' tavg_poly_3_`P' tavg_poly_4_`P' ///
        prcp_poly_1_`P' prcp_poly_2_`P' ///
        lr_tavg_`P'_adm1_avg loggdppc_adm1_avg deathrate_w99 ///
        adm0_code adm1_code adm2_code CHN_ts year population {
        capture confirm variable `v'
        if _rc {
            di as err "!! Missing variable `v' in panel for `P'. Skipping product."
            continue, break
        }
    }

    * Weights: pop normalized within YEAR (per-product panel)
    bysort year: egen double tot_pop = total(population)
    gen double weight = population / tot_pop

    * ====================== SPEC 2 (Preferred) =========================
    reghdfe deathrate_w99                                                     ///
        /* Temp quartic × AGE */                                              ///
        c.tavg_poly_1_`P'#i.agegroup                                         ///
        c.tavg_poly_2_`P'#i.agegroup                                         ///
        c.tavg_poly_3_`P'#i.agegroup                                         ///
        c.tavg_poly_4_`P'#i.agegroup                                         ///
        /* Temp × log income × AGE */                                         ///
        c.tavg_poly_1_`P'#c.loggdppc_adm1_avg#i.agegroup                     ///
        c.tavg_poly_2_`P'#c.loggdppc_adm1_avg#i.agegroup                     ///
        c.tavg_poly_3_`P'#c.loggdppc_adm1_avg#i.agegroup                     ///
        c.tavg_poly_4_`P'#c.loggdppc_adm1_avg#i.agegroup                     ///
        /* Temp × LR-mean × AGE */                                            ///
        c.tavg_poly_1_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup                  ///
        c.tavg_poly_2_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup                  ///
        c.tavg_poly_3_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup                  ///
        c.tavg_poly_4_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup,                 ///
        absorb(                                                               ///
            /* Precip FE-slopes: AGE×ADM0×prcp */                            ///
            i.adm0_code#i.agegroup##c.prcp_poly_1_`P'                        ///
            i.adm0_code#i.agegroup##c.prcp_poly_2_`P'                        ///
            /* Subnational FE with China break: AGE×ADM2×CHN_ts */           ///
            i.adm2_code#i.CHN_ts#i.agegroup                                   ///
            /* AGE×Country×Year FE */                                         ///
            i.adm0_code#i.year#i.agegroup                                    ///
        ) cluster(adm1_code) residual(e_hat)

    estimates save "`sterdir'/agespec_interaction_response_`P'.ster", replace
    di as res "✅ Saved Spec 2 → `sterdir'/agespec_interaction_response_`P'.ster"

    * ======================= OPTIONAL SPECS ============================
    if (`altspec') {

        * Spec 1
        reghdfe deathrate_w99                                                 ///
            c.tavg_poly_1_`P'#i.agegroup                                     ///
            c.tavg_poly_2_`P'#i.agegroup                                     ///
            c.tavg_poly_3_`P'#i.agegroup                                     ///
            c.tavg_poly_4_`P'#i.agegroup                                     ///
            c.tavg_poly_1_`P'#c.loggdppc_adm1_avg#i.agegroup                 ///
            c.tavg_poly_2_`P'#c.loggdppc_adm1_avg#i.agegroup                 ///
            c.tavg_poly_3_`P'#c.loggdppc_adm1_avg#i.agegroup                 ///
            c.tavg_poly_4_`P'#c.loggdppc_adm1_avg#i.agegroup                 ///
            c.tavg_poly_1_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup              ///
            c.tavg_poly_2_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup              ///
            c.tavg_poly_3_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup              ///
            c.tavg_poly_4_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup              ///
            [pw = weight],                                                    ///
            absorb(                                                           ///
                i.adm0_code#i.agegroup##c.prcp_poly_1_`P'                    ///
                i.adm0_code#i.agegroup##c.prcp_poly_2_`P'                    ///
                i.adm2_code#i.CHN_ts#i.agegroup                               ///
                i.adm0_code#i.year                                           ///
            ) cluster(adm1_code)
        estimates save "`sterdir'/agespec_response_spec1_interacted_`P'.ster", replace

        * Spec 3
        reghdfe deathrate_w99                                                 ///
            c.tavg_poly_1_`P'#i.agegroup                                     ///
            c.tavg_poly_2_`P'#i.agegroup                                     ///
            c.tavg_poly_3_`P'#i.agegroup                                     ///
            c.tavg_poly_4_`P'#i.agegroup                                     ///
            c.tavg_poly_1_`P'#c.loggdppc_adm1_avg#i.agegroup                 ///
            c.tavg_poly_2_`P'#c.loggdppc_adm1_avg#i.agegroup                 ///
            c.tavg_poly_3_`P'#c.loggdppc_adm1_avg#i.agegroup                 ///
            c.tavg_poly_4_`P'#c.loggdppc_adm1_avg#i.agegroup                 ///
            c.tavg_poly_1_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup              ///
            c.tavg_poly_2_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup              ///
            c.tavg_poly_3_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup              ///
            c.tavg_poly_4_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup              ///
            [pw = weight],                                                    ///
            absorb(                                                           ///
                i.adm0_code#i.agegroup##c.prcp_poly_1_`P'                    ///
                i.adm0_code#i.agegroup##c.prcp_poly_2_`P'                    ///
                i.adm2_code#i.CHN_ts#i.agegroup                               ///
                i.adm0_code#i.year#i.agegroup                                 ///
                adm1_agegrp_code##c.year                                      ///
            ) cluster(adm1_code)
        estimates save "`sterdir'/agespec_response_spec3_interacted_`P'.ster", replace

        * Spec 4 (FGLS two-step)
        preserve
            reghdfe deathrate_w99                                             ///
                c.tavg_poly_1_`P'#i.agegroup                                 ///
                c.tavg_poly_2_`P'#i.agegroup                                 ///
                c.tavg_poly_3_`P'#i.agegroup                                 ///
                c.tavg_poly_4_`P'#i.agegroup                                 ///
                c.tavg_poly_1_`P'#c.loggdppc_adm1_avg#i.agegroup             ///
                c.tavg_poly_2_`P'#c.loggdppc_adm1_avg#i.agegroup             ///
                c.tavg_poly_3_`P'#c.loggdppc_adm1_avg#i.agegroup             ///
                c.tavg_poly_4_`P'#c.loggdppc_adm1_avg#i.agegroup             ///
                c.tavg_poly_1_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup          ///
                c.tavg_poly_2_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup          ///
                c.tavg_poly_3_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup          ///
                c.tavg_poly_4_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup          ///
                [pw = weight],                                                ///
                absorb(                                                       ///
                    i.adm0_code#i.agegroup##c.prcp_poly_1_`P'                ///
                    i.adm0_code#i.agegroup##c.prcp_poly_2_`P'                ///
                    i.adm2_code#i.CHN_ts#i.agegroup                           ///
                    i.adm0_code#i.year#i.agegroup                             ///
                ) cluster(adm1_code) residual(e_hat)

            bysort adm1_code: egen double omega = sd(e_hat) if e_hat < .
            gen double precisionweight = 1/(omega^2)

            sort iso adm1_id adm2_id agegroup year

            reghdfe deathrate_w99                                             ///
                c.tavg_poly_1_`P'#i.agegroup                                 ///
                c.tavg_poly_2_`P'#i.agegroup                                 ///
                c.tavg_poly_3_`P'#i.agegroup                                 ///
                c.tavg_poly_4_`P'#i.agegroup                                 ///
                c.tavg_poly_1_`P'#c.loggdppc_adm1_avg#i.agegroup             ///
                c.tavg_poly_2_`P'#c.loggdppc_adm1_avg#i.agegroup             ///
                c.tavg_poly_3_`P'#c.loggdppc_adm1_avg#i.agegroup             ///
                c.tavg_poly_4_`P'#c.loggdppc_adm1_avg#i.agegroup             ///
                c.tavg_poly_1_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup          ///
                c.tavg_poly_2_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup          ///
                c.tavg_poly_3_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup          ///
                c.tavg_poly_4_`P'#c.lr_tavg_`P'_adm1_avg#i.agegroup          ///
                [pw = precisionweight],                                       ///
                absorb(                                                       ///
                    i.adm0_code#i.agegroup##c.prcp_poly_1_`P'                ///
                    i.adm0_code#i.agegroup##c.prcp_poly_2_`P'                ///
                    i.adm2_code#i.CHN_ts#i.agegroup                           ///
                    i.adm0_code#i.year#i.agegroup                             ///
                ) cluster(adm1_code) tol(1e-7) maxiter(10000)
            estimates save "`sterdir'/agespec_response_spec4_interacted_`P'.ster", replace
        restore

        * Spec 5 (13-month), if present
        capture confirm variable tavg_poly_1_`P'_13m
        if !_rc {
            reghdfe deathrate_w99                                             ///
                c.tavg_poly_1_`P'_13m#i.agegroup                              ///
                c.tavg_poly_2_`P'_13m#i.agegroup                              ///
                c.tavg_poly_3_`P'_13m#i.agegroup                              ///
                c.tavg_poly_4_`P'_13m#i.agegroup                              ///
                c.tavg_poly_1_`P'_13m#c.loggdppc_adm1_avg#i.agegroup          ///
                c.tavg_poly_2_`P'_13m#c.loggdppc_adm1_avg#i.agegroup          ///
                c.tavg_poly_3_`P'_13m#c.loggdppc_adm1_avg#i.agegroup          ///
                c.tavg_poly_4_`P'_13m#c.loggdppc_adm1_avg#i.agegroup          ///
                c.tavg_poly_1_`P'_13m#c.lr_tavg_`P'_adm1_avg#i.agegroup       ///
                c.tavg_poly_2_`P'_13m#c.lr_tavg_`P'_adm1_avg#i.agegroup       ///
                c.tavg_poly_3_`P'_13m#c.lr_tavg_`P'_adm1_avg#i.agegroup       ///
                c.tavg_poly_4_`P'_13m#c.lr_tavg_`P'_adm1_avg#i.agegroup       ///
                [pw = weight],                                                ///
                absorb(                                                       ///
                    i.adm0_code#i.agegroup##c.prcp_poly_1_`P'                 ///
                    i.adm0_code#i.agegroup##c.prcp_poly_2_`P'                 ///
                    i.adm2_code#i.CHN_ts#i.agegroup                           ///
                    i.adm0_code#i.year#i.agegroup                             ///
                ) cluster(adm1_code)
            estimates save "`sterdir'/agespec_response_spec5_interacted_`P'.ster", replace
        }
    }

    di as res "✅ Finished `P'"

} // products

di as res "All per-product age-specific interacted regressions completed."
