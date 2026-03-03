/**************************************************************************
* Age-specific interacted regressions across products
* Carleton et al. (2022) Spec 2 by default; altspec=1 runs 1/3/4/5.
* - Adds precip controls in FE: (AGE×ADM0×prcp_poly_1/2)
* - Handles product name hyphen vs underscore in file/var names
**************************************************************************/

* --- CONFIG -------------------------------------------------------------
global CODE    "/user/ab5405/summeraliaclimate/code/regressions"
global OUTPUT  "/user/ab5405/summeraliaclimate/code/regressions/output"

* Use human-readable product labels (filenames may use hyphens),
* but variable stems will use underscores.
local prods "ERA5-025 GMFD MERRA2 JRA-3Q"

local altspec 0   // toggle: 0 = only spec2, 1 = run all altspecs

foreach prod of local prods {

    di as txt "=============================="
    di as txt ">>> Running age-spec regression for `prod'"
    di as txt "=============================="

    * Map to variable stem with underscores (e.g., ERA5-025 -> ERA5_025)
    local prod_var : subinstr local prod "-" "_", all

    * Resolve input panel path (try hyphen name, then underscore name)
    local in_hyp "$CODE/output_panels/panel_prepped_for_regressions_`prod'_aligned.dta"
    local in_und "$CODE/output_panels/panel_prepped_for_regressions_`prod_var'_aligned.dta"

    capture confirm file "`in_hyp'"
    if _rc {
        capture confirm file "`in_und'"
        if _rc {
            di as err "!! Could not find input panel for `prod'. Tried:"
            di as err "   `in_hyp'"
            di as err "   `in_und'"
            continue
        }
        else {
            use "`in_und'", clear
        }
    }
    else {
        use "`in_hyp'", clear
    }

    * Output dir
    local sterdir "$OUTPUT/age_spec_interacted"
    cap mkdir "`sterdir'"

    * --- Ensure agegroup is numeric for factor notation ---
    capture confirm numeric variable agegroup
    if _rc {
        encode agegroup, gen(agegroup_)
        drop agegroup
        rename agegroup_ agegroup
    }

    * --- Sanity checks for required variables (fail fast with clear msg) ---
    foreach v in ///
        tavg_poly_1_`prod_var' tavg_poly_2_`prod_var' tavg_poly_3_`prod_var' tavg_poly_4_`prod_var' ///
        prcp_poly_1_`prod_var' prcp_poly_2_`prod_var' ///
        lr_tavg_`prod_var'_adm1_avg loggdppc_adm1_avg deathrate_w99 ///
        adm0_code adm1_code adm2_code CHN_ts year {
        capture confirm variable `v'
        if _rc {
            di as err "!! Missing variable `v' for `prod' (`prod_var')."
            continue, break
        }
    }

    /********************************************************************
     * SPEC 2 (preferred): 4th-order polynomial × AGE,
     * with Temp×Income×AGE and Temp×LR-Mean×AGE;
     * FE: (AGE×ADM2) & (AGE×Country×Year);
     * Precipitation controls: (AGE×ADM0×prcp_poly_1/2)
     ********************************************************************/
    reghdfe deathrate_w99 ///
        /* Temp quartic × AGE */ ///
        c.tavg_poly_1_`prod_var'#i.agegroup ///
        c.tavg_poly_2_`prod_var'#i.agegroup ///
        c.tavg_poly_3_`prod_var'#i.agegroup ///
        c.tavg_poly_4_`prod_var'#i.agegroup ///
        /* Temp × income × AGE */ ///
        c.tavg_poly_1_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
        c.tavg_poly_2_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
        c.tavg_poly_3_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
        c.tavg_poly_4_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
        /* Temp × LR-mean × AGE */ ///
        c.tavg_poly_1_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
        c.tavg_poly_2_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
        c.tavg_poly_3_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
        c.tavg_poly_4_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
        , absorb( ///
            /* Precip controls as FE interactions (AGE×ADM0×prcp) */ ///
            i.adm0_code#c.prcp_poly_1_`prod_var'#i.agegroup ///
            i.adm0_code#c.prcp_poly_2_`prod_var'#i.agegroup ///
            /* AGE×ADM2 FE with China discontinuity */ ///
            i.adm2_code#i.CHN_ts#i.agegroup ///
            /* AGE×Country×Year FE */ ///
            i.adm0_code#i.year#i.agegroup ///
        ) cluster(adm1_code)

    estimates save "`sterdir'/agespec_interaction_response_`prod_var'_noint.ster", replace
    di as res "Saved Spec2 for `prod' (`prod_var')"

    * -----------------------------------------------------------------
    * OPTIONAL: Alternative specs if requested
    * -----------------------------------------------------------------
    if (`altspec') {

        * Spec 1: AGE×ADM2 & ADM0×Year; precip controls AGE×ADM0
        reghdfe deathrate_w99 ///
            c.tavg_poly_1_`prod_var'#i.agegroup ///
            c.tavg_poly_2_`prod_var'#i.agegroup ///
            c.tavg_poly_3_`prod_var'#i.agegroup ///
            c.tavg_poly_4_`prod_var'#i.agegroup ///
            c.tavg_poly_1_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
            c.tavg_poly_2_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
            c.tavg_poly_3_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
            c.tavg_poly_4_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
            c.tavg_poly_1_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
            c.tavg_poly_2_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
            c.tavg_poly_3_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
            c.tavg_poly_4_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
            , absorb( ///
                i.adm0_code#c.prcp_poly_1_`prod_var'#i.agegroup ///
                i.adm0_code#c.prcp_poly_2_`prod_var'#i.agegroup ///
                i.adm2_code#i.CHN_ts#i.agegroup ///
                i.adm0_code#i.year ///
            ) cluster(adm1_code)
        estimates save "`sterdir'/agespec_response_spec1_interacted_`prod_var'.ster", replace

        * Spec 3: add (ADM1×Age) linear trends
        reghdfe deathrate_w99 ///
            c.tavg_poly_1_`prod_var'#i.agegroup ///
            c.tavg_poly_2_`prod_var'#i.agegroup ///
            c.tavg_poly_3_`prod_var'#i.agegroup ///
            c.tavg_poly_4_`prod_var'#i.agegroup ///
            c.tavg_poly_1_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
            c.tavg_poly_2_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
            c.tavg_poly_3_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
            c.tavg_poly_4_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
            c.tavg_poly_1_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
            c.tavg_poly_2_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
            c.tavg_poly_3_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
            c.tavg_poly_4_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
            , absorb( ///
                i.adm0_code#c.prcp_poly_1_`prod_var'#i.agegroup ///
                i.adm0_code#c.prcp_poly_2_`prod_var'#i.agegroup ///
                i.adm2_code#i.CHN_ts#i.agegroup ///
                i.adm0_code#i.year#i.agegroup ///
                adm1_agegrp_code##c.year ///
            ) cluster(adm1_code)
        estimates save "`sterdir'/agespec_response_spec3_interacted_`prod_var'.ster", replace

        * Spec 4: FGLS (two-step), weights from Spec 2 residuals
        preserve
            reghdfe deathrate_w99 ///
                c.tavg_poly_1_`prod_var'#i.agegroup ///
                c.tavg_poly_2_`prod_var'#i.agegroup ///
                c.tavg_poly_3_`prod_var'#i.agegroup ///
                c.tavg_poly_4_`prod_var'#i.agegroup ///
                c.tavg_poly_1_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_2_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_3_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_4_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_1_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_2_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_3_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_4_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                , absorb( ///
                    i.adm0_code#c.prcp_poly_1_`prod_var'#i.agegroup ///
                    i.adm0_code#c.prcp_poly_2_`prod_var'#i.agegroup ///
                    i.adm2_code#i.CHN_ts#i.agegroup ///
                    i.adm0_code#i.year#i.agegroup ///
                ) cluster(adm1_code) residual(e_hat)

            bysort adm1_code: egen omega = sd(e_hat) if e_hat!=.
            gen precisionweight = 1/(omega*omega)

            sort iso adm1_id adm2_id agegroup year

            reghdfe deathrate_w99 ///
                c.tavg_poly_1_`prod_var'#i.agegroup ///
                c.tavg_poly_2_`prod_var'#i.agegroup ///
                c.tavg_poly_3_`prod_var'#i.agegroup ///
                c.tavg_poly_4_`prod_var'#i.agegroup ///
                c.tavg_poly_1_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_2_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_3_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_4_`prod_var'#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_1_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_2_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_3_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_4_`prod_var'#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                [pw = precisionweight], ///
                absorb( ///
                    i.adm0_code#c.prcp_poly_1_`prod_var'#i.agegroup ///
                    i.adm0_code#c.prcp_poly_2_`prod_var'#i.agegroup ///
                    i.adm2_code#i.CHN_ts#i.agegroup ///
                    i.adm0_code#i.year#i.agegroup ///
                ) cluster(adm1_code) tol(1e-7) maxiter(10000)
            estimates save "`sterdir'/agespec_response_spec4_interacted_`prod_var'.ster", replace
        restore

        * Spec 5: 13-month exposure (if available)
        capture confirm variable tavg_poly_1_`prod_var'_13m
        if !_rc {
            reghdfe deathrate_w99 ///
                c.tavg_poly_1_`prod_var'_13m#i.agegroup ///
                c.tavg_poly_2_`prod_var'_13m#i.agegroup ///
                c.tavg_poly_3_`prod_var'_13m#i.agegroup ///
                c.tavg_poly_4_`prod_var'_13m#i.agegroup ///
                c.tavg_poly_1_`prod_var'_13m#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_2_`prod_var'_13m#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_3_`prod_var'_13m#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_4_`prod_var'_13m#c.loggdppc_adm1_avg#i.agegroup ///
                c.tavg_poly_1_`prod_var'_13m#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_2_`prod_var'_13m#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_3_`prod_var'_13m#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                c.tavg_poly_4_`prod_var'_13m#c.lr_tavg_`prod_var'_adm1_avg#i.agegroup ///
                , absorb( ///
                    i.adm0_code#c.prcp_poly_1_`prod_var'#i.agegroup ///
                    i.adm0_code#c.prcp_poly_2_`prod_var'#i.agegroup ///
                    i.adm2_code#i.CHN_ts#i.agegroup ///
                    i.adm0_code#i.year#i.agegroup ///
                ) cluster(adm1_code)
            estimates save "`sterdir'/agespec_response_spec5_interacted_`prod_var'.ster", replace
        }
    }

    di as res "Finished all requested regressions for `prod' (`prod_var')"

}
