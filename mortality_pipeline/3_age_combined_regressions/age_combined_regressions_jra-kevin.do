/**************************************************************************
* Age-combined JRA-3Q pooled response — KEVIN panel
**************************************************************************/

*****************************************************************************
* PART 1. Initializing
*****************************************************************************
* 0. Set your output folder and define where to save .ster files
global OUTPUT "/user/ab5405/summeraliaclimate/code/regressions/output"
local ster "$OUTPUT/age_combined_kev_JRA_3Q"    // <— changed: separate KEVIN results

* create the results folder if it doesn't exist
capture mkdir "`ster'"

* 1. Load the prepped JRA panel (KEVIN version)
use "/user/ab5405/summeraliaclimate/code/regressions/output_panels/panel_prepped_for_regressions_kev_JRA_3Q.dta", clear
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// <— changed: point to your KEVIN-prepped panel

*****************************************************************************
* PART 2. OLS Regressions
*****************************************************************************

* 1. set weighting schemes
bysort year: egen tot_pop = total(population)
gen weight = population / tot_pop

* 2. run regressions

* specification 1
reghdfe deathrate_w99 tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
        [pw = weight]  ///
        , absorb(adm0_code##c.prcp_poly_1_JRA_3Q adm0_code##c.prcp_poly_2_JRA_3Q ///
                i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year ) ///
        cluster(adm1_code)
estimates save "`ster'/pooled_response_spec1_public_kev_JRA_3Q.ster", replace
//                                                     ^^^ add "kev" tag

* specification 2
reghdfe deathrate_w99 tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
        [pw = weight] ///
        , absorb(adm0_code##c.prcp_poly_1_JRA_3Q      adm0_code##c.prcp_poly_2_JRA_3Q ///
                 i.adm2_code#i.CHN_ts#i.agegroup  i.adm0_code#i.year#i.agegroup ) ///
        cluster(adm1_code) residual(e_hat)
estimates save "`ster'/pooled_response_spec2_public_kev_JRA_3Q.ster", replace

* specification 3
reghdfe deathrate_w99 tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
        [pw = weight] ///
        , absorb(adm0_code##c.prcp_poly_1_JRA_3Q adm0_code##c.prcp_poly_2_JRA_3Q ///
                 i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup adm1_agegrp_code##c.year) ///
        cluster(adm1_code)
estimates save "`ster'/pooled_response_spec3_public_kev_JRA_3Q.ster", replace


*****************************************************************************
* PART 3. FGLS Regressions
*****************************************************************************

* 1. generate weighting matrix
gen e2 = e_hat * e_hat
bysort adm1_code: egen omega = sd(e_hat) if e_hat != .
gen precisionweight = weight * 1/(omega*omega)

sort iso adm1_id adm2_id year agegroup

* 2. run regressions
* specification 4
reghdfe deathrate_w99 tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
        [pw = precisionweight] ///
        , absorb(adm0_code##c.prcp_poly_1_JRA_3Q adm0_code##c.prcp_poly_2_JRA_3Q ///
                i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup) ///
        cluster(adm1_code) tol(1e-7)
estimates save "`ster'/pooled_response_spec4_public_kev_JRA_3Q.ster", replace


*****************************************************************************
* PART 4. All cross
*****************************************************************************

reghdfe deathrate_w99                                                             ///
    /* 1) Temperature quartic */                                                  ///
    c.tavg_poly_1_JRA_3Q   c.tavg_poly_2_JRA_3Q  c.tavg_poly_3_JRA_3Q   c.tavg_poly_4_JRA_3Q    ///
                                                                                  ///
    /* 2) Temp × income */                                                        ///
    c.tavg_poly_1_JRA_3Q#c.loggdppc_adm1_avg   ///
    c.tavg_poly_2_JRA_3Q#c.loggdppc_adm1_avg   ///
    c.tavg_poly_3_JRA_3Q#c.loggdppc_adm1_avg   ///
    c.tavg_poly_4_JRA_3Q#c.loggdppc_adm1_avg   ///
                                                                                  ///
    /* 3) Temp × LR-mean-temp */                                                  ///
    c.tavg_poly_1_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg   ///
    c.tavg_poly_2_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg   ///
    c.tavg_poly_3_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg   ///
    c.tavg_poly_4_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg   ///
                                                                                  ///
    [pw = weight],                                                                 ///
    absorb(                                                                         ///
      /* a) country-specific precip polynomials */                                  ///
      adm0_code##c.prcp_poly_1_JRA_3Q                                           ///
      adm0_code##c.prcp_poly_2_JRA_3Q                                              ///
                                                                                  ///
      /* b) regional seasonality (China break) */                                   ///
      i.adm2_code#i.CHN_ts#i.agegroup                                              ///
                                                                                  ///
      /* c) country × year shocks */                                               ///
      i.adm0_code#i.year                                                           ///
    )                                                                              ///
    cluster(adm1_code)

estimates save "`ster'/pooled_response_spec5_public_kev_JRA_3Q.ster", replace
