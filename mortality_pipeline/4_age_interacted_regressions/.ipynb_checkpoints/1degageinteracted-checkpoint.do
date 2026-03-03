/**************************************************************************
* age_specific_regressions_FLEXIBLE.do
*
* FLEXIBLE RESOLUTION VERSION
*
* Age-specific interacted regressions across products
* Carleton et al. (2022) Spec 2 by default; altspec=1 runs 1/3/4/5.
* 
* CHANGE LOG:
* - Added resolution parameter to control input/output directories
* - Reads from mortality_pipeline/output_panels/{resolution}/
* - Saves to mortality_pipeline/output/{resolution}/age_spec_interacted/
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

* Products to process
local prods "GMFD ERA5 JRA_3Q MERRA2"

local altspec 0   // toggle: 0 = only spec2, 1 = run all altspecs

di ""
di "========================================================================"
di "CONFIGURATION:"
di "  Resolution: `resolution'"
di "  Products: `prods'"
di "  Input dir: $PANELS"
di "  Output dir: $OUTPUT/age_spec_interacted"
di "  Alt specs: `altspec'"
di "========================================================================"
di ""

* ============================================================================
* PRODUCT LOOP
* ============================================================================

foreach prod of local prods {

    di ""
    di "========================================================================"
    di ">>> PROCESSING: `prod'"
    di "========================================================================"
    di ""

    * Map to variable stem with underscores (e.g., ERA5-025 -> ERA5_025)
    local prod_var : subinstr local prod "-" "_", all

    * ========================================================================
    * LOAD DATA
    * ========================================================================
    
    * Try with hyphen first, then underscore
    local in_hyp "$PANELS/panel_prepped_for_regressions_`prod'_aligned.dta"
    local in_und "$PANELS/panel_prepped_for_regressions_`prod_var'_aligned.dta"

    capture confirm file "`in_hyp'"
    if _rc {
        capture confirm file "`in_und'"
        if _rc {
            di as error ""
            di as error "❌ ERROR: Could not find input panel for `prod'"
            di as error "   Tried:"
            di as error "     `in_hyp'"
            di as error "     `in_und'"
            di as error ""
            continue
        }
        else {
            di "Loading: `in_und'"
            use "`in_und'", clear
        }
    }
    else {
        di "Loading: `in_hyp'"
        use "`in_hyp'", clear
    }

    * Output directory
    local sterdir "$OUTPUT/age_spec_interacted"
    capture mkdir "$OUTPUT"
    capture mkdir "`sterdir'"

    * ========================================================================
    * DATA PREP
    * ========================================================================
    
    * Ensure agegroup is numeric for factor notation
    capture confirm numeric variable agegroup
    if _rc {
        encode agegroup, gen(agegroup_)
        drop agegroup
        rename agegroup_ agegroup
    }

    * ========================================================================
    * SANITY CHECKS
    * ========================================================================
    
    di "Checking required variables..."
    
    local required_vars ///
        tavg_poly_1_`prod_var' tavg_poly_2_`prod_var' ///
        tavg_poly_3_`prod_var' tavg_poly_4_`prod_var' ///
        prcp_poly_1_`prod_var' prcp_poly_2_`prod_var' ///
        lr_tavg_`prod_var'_adm1_avg loggdppc_adm1_avg ///
        deathrate_w99 adm0_code adm1_code adm2_code ///
        CHN_ts year
    
    local all_ok = 1
    foreach v of local required_vars {
        capture confirm variable `v'
        if _rc {
            di as error "  ❌ Missing: `v'"
            local all_ok = 0
        }
    }
    
    if !`all_ok' {
        di as error ""
        di as error "Skipping `prod' due to missing variables"
        di as error ""
        continue
    }
    
    di "  ✅ All required variables present"

    * ========================================================================
    * SPEC 2 (PREFERRED): 4th-order polynomial × AGE
    * ========================================================================
    
    di ""
    di "--------------------------------------------------------------------"
    di "ESTIMATING SPEC 2 (Preferred) for `prod'"
    di "--------------------------------------------------------------------"
    di ""
    
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
    di ""
    di "✅ Saved Spec 2 for `prod'"

    * ========================================================================
    * ALTERNATIVE SPECS (if requested)
    * ========================================================================
    
    if (`altspec') {
        
        di ""
        di "--------------------------------------------------------------------"
        di "ESTIMATING ALT SPECS for `prod'"
        di "--------------------------------------------------------------------"

        * --------------------------------------------------------------------
        * SPEC 1: Baseline FE
        * --------------------------------------------------------------------
        
        di ""
        di "  Spec 1..."
        
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
        di "  ✅ Saved Spec 1"

        * --------------------------------------------------------------------
        * SPEC 3: + ADM1×Age time trends
        * --------------------------------------------------------------------
        
        di ""
        di "  Spec 3..."
        
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
        di "  ✅ Saved Spec 3"

        * --------------------------------------------------------------------
        * SPEC 4: FGLS precision weights
        * --------------------------------------------------------------------
        
        di ""
        di "  Spec 4 (FGLS)..."
        
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
        di "  ✅ Saved Spec 4"

        * --------------------------------------------------------------------
        * SPEC 5: 13-month exposure (if available)
        * --------------------------------------------------------------------
        
        capture confirm variable tavg_poly_1_`prod_var'_13m
        if !_rc {
            di ""
            di "  Spec 5 (13-month)..."
            
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
            di "  ✅ Saved Spec 5"
        }
    }

    di ""
    di "========================================================================"
    di "✅ FINISHED: `prod' (`prod_var')"
    di "========================================================================"
    di ""
}

* ============================================================================
* FINAL SUMMARY
* ============================================================================

di ""
di "========================================================================"
di "ALL PRODUCTS COMPLETE!"
di "========================================================================"
di ""
di "Resolution: `resolution'"
di "Output directory: $OUTPUT/age_spec_interacted"
di ""
di "========================================================================"
di ""
