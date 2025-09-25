/**************************************************************************
* prep_regressions_all.do
* - Estimates the temperature-mortality response function using 
*   subnational data across 40 countries pooled across 3 age groups
* - Loops over observational products: ERA5_025, JRA_3Q, MERRA2, GMFD
* - Saves spec1–spec5 regression .ster files for each product
**************************************************************************/

global CODE   "/user/ab5405/summeraliaclimate/code"
global OUTPUT "$CODE/regressions/output"
local products MERRA2

foreach product of local products {

    di "Running regressions for `product'..."

    * Input/output
    local inpanel "$CODE/regressions/output_panels/panel_prepped_for_regressions_`product'_aligned.dta"
    local ster    "$OUTPUT/age_combined_`product'"
    capture mkdir "`ster'"

    use "`inpanel'", clear

    * --- Ensure agegroup numeric for factor notation ---
    capture confirm numeric variable agegroup
    if _rc {
        encode agegroup, gen(agegroup_num)
        drop agegroup
        rename agegroup_num agegroup
    }

    * 1. set weighting scheme
    bysort year: egen tot_pop = total(population)
    gen weight = population / tot_pop

    *****************************************************************************
    * PART 2. OLS Regressions
    *****************************************************************************

    * Spec 1
    reghdfe deathrate_w99 tavg_poly_1_`product' tavg_poly_2_`product' tavg_poly_3_`product' tavg_poly_4_`product' ///
            [pw = weight], ///
            absorb(adm0_code##c.prcp_poly_1_`product' adm0_code##c.prcp_poly_2_`product' ///
                   i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year) ///
            cluster(adm1_code)
    estimates save "`ster'/pooled_response_spec1_public_`product'.ster", replace

    * Spec 2 (preferred)
    reghdfe deathrate_w99 tavg_poly_1_`product' tavg_poly_2_`product' tavg_poly_3_`product' tavg_poly_4_`product' ///
            [pw = weight], ///
            absorb(adm0_code##c.prcp_poly_1_`product' adm0_code##c.prcp_poly_2_`product' ///
                   i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup) ///
            cluster(adm1_code) residual(e_hat)
    estimates save "`ster'/pooled_response_spec2_public_`product'.ster", replace

    * Spec 3
    reghdfe deathrate_w99 tavg_poly_1_`product' tavg_poly_2_`product' tavg_poly_3_`product' tavg_poly_4_`product' ///
            [pw = weight], ///
            absorb(adm0_code##c.prcp_poly_1_`product' adm0_code##c.prcp_poly_2_`product' ///
                   i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup ///
                   adm1_agegrp_code##c.year) ///
            cluster(adm1_code)
    estimates save "`ster'/pooled_response_spec3_public_`product'.ster", replace

    *****************************************************************************
    * PART 3. FGLS Regressions
    *****************************************************************************

    gen e2 = e_hat * e_hat
    bysort adm1_code: egen omega = sd(e_hat) if e_hat != .
    gen precisionweight = weight * 1/(omega*omega)

    sort iso adm1_id adm2_id year agegroup

    reghdfe deathrate_w99 tavg_poly_1_`product' tavg_poly_2_`product' tavg_poly_3_`product' tavg_poly_4_`product' ///
            [pw = precisionweight], ///
            absorb(adm0_code##c.prcp_poly_1_`product' adm0_code##c.prcp_poly_2_`product' ///
                   i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup) ///
            cluster(adm1_code) tol(1e-7)
    estimates save "`ster'/pooled_response_spec4_public_`product'.ster", replace

    *****************************************************************************
    * PART 4. Full Interactions
    *****************************************************************************

    reghdfe deathrate_w99                                                              ///
        /* 1) Temperature quartic */                                                   ///
        c.tavg_poly_1_`product'   c.tavg_poly_2_`product'  c.tavg_poly_3_`product'   c.tavg_poly_4_`product' ///
                                                                                     ///
        /* 2) Temp × income */                                                         ///
        c.tavg_poly_1_`product'#c.loggdppc_adm1_avg   ///
        c.tavg_poly_2_`product'#c.loggdppc_adm1_avg   ///
        c.tavg_poly_3_`product'#c.loggdppc_adm1_avg   ///
        c.tavg_poly_4_`product'#c.loggdppc_adm1_avg   ///
                                                                                     ///
        /* 3) Temp × LR-mean-temp */                                                   ///
        c.tavg_poly_1_`product'#c.lr_tavg_`product'_adm1_avg   ///
        c.tavg_poly_2_`product'#c.lr_tavg_`product'_adm1_avg   ///
        c.tavg_poly_3_`product'#c.lr_tavg_`product'_adm1_avg   ///
        c.tavg_poly_4_`product'#c.lr_tavg_`product'_adm1_avg   ///
                                                                                     ///
        [pw = weight],                                                                ///
        absorb(                                                                       ///
          adm0_code##c.prcp_poly_1_`product'                                          ///
          adm0_code##c.prcp_poly_2_`product'                                          ///
          i.adm2_code#i.CHN_ts                                                        ///
          i.adm0_code#i.year                                                          ///
        )                                                                             ///
        cluster(adm1_code)

    estimates save "`ster'/pooled_response_spec5_public_`product'.ster", replace

    di "Finished regressions for `product'"

    clear
}

di "All products processed successfully."
