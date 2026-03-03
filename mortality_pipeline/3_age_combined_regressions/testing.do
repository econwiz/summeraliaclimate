/**************************************************************************
* Pooled temperature–mortality regressions (JRA-3Q, aligned panel)
**************************************************************************/

* 0) Paths
global OUTPUT "/user/ab5405/summeraliaclimate/code/regressions/output"
local ster     "$OUTPUT/age_combined_JRA_3Q"
capture mkdir "`ster'"

* 1) Load the prepped (ALIGNED) JRA panel
use "/user/ab5405/summeraliaclimate/code/regressions/output_panels/panel_prepped_for_regressions_JRA_3Q_aligned.dta", clear

* 1a) Make sure agegroup is numeric so factor notation works
capture confirm numeric variable agegroup
if _rc {
    encode agegroup, gen(agegroup_)
    drop agegroup
    rename agegroup_ agegroup
}

*****************************************************************************
* PART 2. OLS Regressions
*****************************************************************************

* Weights
bysort year: egen tot_pop = total(population)
gen weight = population / tot_pop

* ---------- Spec 1 ----------
reghdfe deathrate_w99 ///
       tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
       [pw = weight], ///
       absorb( ///
         adm0_code##c.prcp_poly_1_JRA_3Q ///
         adm0_code##c.prcp_poly_2_JRA_3Q ///
         i.adm2_code#i.CHN_ts#i.agegroup ///
         i.adm0_code#i.year ///
       ) ///
       cluster(adm1_code)
estimates save "`ster'/pooled_response_spec1_public_JRA_3Q.ster", replace

* ---------- Spec 2 (preferred) ----------
reghdfe deathrate_w99 ///
       tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
       [pw = weight], ///
       absorb( ///
         adm0_code##c.prcp_poly_1_JRA_3Q ///
         adm0_code##c.prcp_poly_2_JRA_3Q ///
         i.adm2_code#i.CHN_ts#i.agegroup ///
         i.adm0_code#i.year#i.agegroup ///
       ) ///
       cluster(adm1_code) residual(e_hat)
estimates save "`ster'/pooled_response_spec2_public_JRA_3Q.ster", replace

*****************************************************************************
* PART 3. FGLS (using residuals from Spec 2)
*****************************************************************************
gen e2 = e_hat*e_hat
bys adm1_code: egen omega = sd(e_hat) if e_hat < .
gen precisionweight = weight * (1/(omega*omega))

sort iso adm1_id adm2_id year agegroup

* ---------- Spec 4 (FGLS) ----------
reghdfe deathrate_w99 ///
       tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
       [pw = precisionweight], ///
       absorb( ///
         adm0_code##c.prcp_poly_1_JRA_3Q ///
         adm0_code##c.prcp_poly_2_JRA_3Q ///
         i.adm2_code#i.CHN_ts#i.agegroup ///
         i.adm0_code#i.year#i.agegroup ///
       ) ///
       cluster(adm1_code) tol(1e-7)
estimates save "`ster'/pooled_response_spec4_public_JRA_3Q.ster", replace

*****************************************************************************
* PART 4. Fully crossed controls (Spec 5)
*****************************************************************************
reghdfe deathrate_w99 ///
    /* 1) Temperature quartic */                                                  ///
    c.tavg_poly_1_JRA_3Q   c.tavg_poly_2_JRA_3Q  c.tavg_poly_3_JRA_3Q   c.tavg_poly_4_JRA_3Q ///
    /* 2) Temp × income */                                                        ///
    c.tavg_poly_1_JRA_3Q#c.loggdppc_adm1_avg ///
    c.tavg_poly_2_JRA_3Q#c.loggdppc_adm1_avg ///
    c.tavg_poly_3_JRA_3Q#c.loggdppc_adm1_avg ///
    c.tavg_poly_4_JRA_3Q#c.loggdppc_adm1_avg ///
    /* 3) Temp × LR-mean-temp */                                                  ///
    c.tavg_poly_1_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg ///
    c.tavg_poly_2_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg ///
    c.tavg_poly_3_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg ///
    c.tavg_poly_4_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg ///
    [pw = weight], ///
    absorb( ///
      /* a) country-specific precip polynomials */                                  ///
      adm0_code##c.prcp_poly_1_JRA_3Q ///
      adm0_code##c.prcp_poly_2_JRA_3Q ///
      /* b) regional seasonality (China break, age-specific) */                     ///
      i.adm2_code#i.CHN_ts#i.agegroup ///
      /* c) country × year shocks (age-specific) */                                 ///
      i.adm0_code#i.year#i.agegroup ///
    ) ///
    cluster(adm1_code)
estimates save "`ster'/pooled_response_spec5_public_JRA_3Q.ster", replace
