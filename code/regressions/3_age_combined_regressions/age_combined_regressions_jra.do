/*

Purpose: Estimates the temperature-mortality response function using 
subnational data across 40 countries pooled across the 3 age groups
(as discussed in Appendix D.2).

Inputs
------

- `regressions/prep_panels/global_mortality_panel_public_JRA_3Q.dta` - Final mortality panel.

Outputs
-------

- `regressions/output_panels`
	- `pooled_response_spec*_JRA_2Q.ster` where * is model 1 through 5. Ster files
	containing uninteracted, age combined regressuion results under various
	fixed effects, estimation, and data construction assumptions.

Notes
------

Summary of models:
    1. 4th-order polynomial OLS (Age x ADM2) & (Age x ADM2) FE
   *2. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE
    3. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) (Age x ADM1 linear trend)
	4. 4th-order polynomial FGLS (Age x ADM2) & (AGE x Country x Year) FE
	5. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE with 13-month climate exposure
* indicates preferred model.

- All regressions are population weighted with clustered standard errors at the
ADM1 level.

- Models also control for AGE x ADM0 precipitation.

- `i.CHN_ts` indicator variable accounts for a discontinuity in China's
mortality/population time series. See 1_estimation/README.md for additional
details.

*/


*****************************************************************************
* 						PART 1. Initializing		 						*
*****************************************************************************
* 0. Set your output folder and define where to save .ster files
global OUTPUT "/user/ab5405/summeraliaclimate/code/regressions/output"
local ster "$OUTPUT/age_combined_JRA_3Q"

* create the results folder if it doesn't exist
capture mkdir "`ster'"

* 1. Load the prepped JRA panel
use "/user/ab5405/summeraliaclimate/code/regressions/output_panels/panel_prepped_for_regressions_JRA_3Q.dta", clear

*****************************************************************************
* 						PART 2. OLS Regressions                 		    *
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
estimates save "`ster'/pooled_response_spec1_public_JRA_3Q.ster", replace

* specification 2
reghdfe deathrate_w99 tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
		[pw = weight] ///
		, absorb(adm0_code##c.prcp_poly_1_JRA_3Q      adm0_code##c.prcp_poly_2_JRA_3Q ///
				 i.adm2_code#i.CHN_ts#i.agegroup  i.adm0_code#i.year#i.agegroup ) ///
		cluster(adm1_code) residual(e_hat)
estimates save "`ster'/pooled_response_spec2_public_JRA_3Q.ster", replace
* save rediduals for FGLS regressions

* specification 3
reghdfe deathrate_w99 tavg_poly_1_JRA_3Q tavg_poly_2_JRA_3Q tavg_poly_3_JRA_3Q tavg_poly_4_JRA_3Q ///
		[pw = weight] ///
		, absorb(adm0_code##c.prcp_poly_1_JRA_3Q adm0_code##c.prcp_poly_2_JRA_3Q ///
				 i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup adm1_agegrp_code##c.year) ///
		cluster(adm1_code)
estimates save "`ster'/pooled_response_spec3_public_JRA_3Q.ster", replace


*****************************************************************************
* 						PART 3. FGLS Regressions                 		    *
*****************************************************************************

* 1. generate weighting matrix
* get residuals from 1st stage and
* calculate diagonal elements of the weighting matrix
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
estimates save "`ster'/pooled_response_spec4_public_JRA_3Q.ster", replace


*****************************************************************************
* 						PART 4. All cross        		    *
*****************************************************************************

reghdfe deathrate_w99                                                             ///
    /* 1) Temperature quartic */                                                  ///
    c.tavg_poly_1_JRA_3Q   c.tavg_poly_2_JRA_3Q  c.tavg_poly_3_JRA_3Q   c.tavg_poly_4_JRA_3Q    ///
                                                                                  ///
    /* 2) Temp × income */                                                         ///
    c.tavg_poly_1_JRA_3Q#c.loggdppc_adm1_avg   ///
    c.tavg_poly_2_JRA_3Q#c.loggdppc_adm1_avg   ///
    c.tavg_poly_3_JRA_3Q#c.loggdppc_adm1_avg   ///
    c.tavg_poly_4_JRA_3Q#c.loggdppc_adm1_avg   ///
                                                                                  ///
    /* 3) Temp × LR‑mean‑temp */                                                   ///
    c.tavg_poly_1_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg   ///
    c.tavg_poly_2_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg   ///
    c.tavg_poly_3_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg   ///
    c.tavg_poly_4_JRA_3Q#c.lr_tavg_JRA_3Q_adm1_avg   ///
                                                                                  ///
    [pw = weight],                                                                 ///
    absorb(                                                                         ///
      /* a) country‑specific precip polynomials */                                  ///
      adm0_code##c.prcp_poly_1_JRA_3Q                                           ///
      adm0_code##c.prcp_poly_2_JRA_3Q                                              ///
                                                                                  ///
      /* b) regional seasonality (China break) */                                   ///
      i.adm2_code#i.CHN_ts                                                         ///
                                                                                  ///
      /* c) country × year shocks */                                                ///
      i.adm0_code#i.year                                                           ///
    )                                                                              ///
    cluster(adm1_code)

estimates save "`ster'/pooled_response_spec5_public_JRA_3Q.ster", replace