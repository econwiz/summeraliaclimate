/*

Purpose: Estimates age specific temperature-mortality response function with
termperature-income and temperature-LR climate interaction terms. 

The `Agespec_interaction_response` file is the main regression that is discussed
in Carleton et al 2022, and is GMFDried through to the projection system 
and beyond:

*2. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE


Notes: In this public release, the `_public` suffix is added on to signify
that the regressions would be run on the publically available mortality 
sample, which excludes USA and China.

Inputs
------

- `output_panels/panel_prepped_for_regression_car.dta` - Final mortality panel for Carleton data.

Outputs
-------

- `output/age_spec_interacted'
	- `agespec_interaction_response_car.ster` - Ster file containing results from
	an age-stacked regression interacted with ADM1 average income and climate. 


Notes
------

Summary of models:
    1. 4th-order polynomial OLS (Age x ADM2) & (Age x ADM2) FE
   *2. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE
    3. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) (Age x ADM1 linear trend)
	4. 4th-order polynomial FGLS (Age x ADM2) & (AGE x Country x Year) FE
	5. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE with 13-month climate exposure
* indicates preferred model.

- Temperature polynomials are interacted with time-invariant ADM1-level 
covariates of average income and climate. Specifically, income and climate are
averages of 13-year and 30-year bartlett kernels, respectively.

- Models also control for AGE x ADM0 precipitation.

- All regressions are unweighted with clustered standard errors at the
ADM1 level.

- `i.CHN_ts` indicator variable accounts for a discontinuity in China's
mortality/population time series. See 1_estimation/README.md for additional
details.

NOTE: As we only GMFDry forward the perferred model (Spec. 2) through the rest of
the analysis, this code only runs that model by default. To run interacted
regressions for other specifications, modify the `altspec` toggle below from
0 to 1.

*/


*****************************************************************************
* 						PART 1. Initializing		 						*
*****************************************************************************

* 1) set paths
global CODE   "/user/ab5405/summeraliaclimate/code/regressions"
global OUTPUT "/user/ab5405/summeraliaclimate/code/regressions/output"

* 2)Load GMFD mortality panel
use "/user/ab5405/summeraliaclimate/code/regressions/output_panels/panel_prepped_for_regression_car.dta", clear


* 4) ——— Define where you’ll save your Ster files ———
local sterdir = "$OUTPUT/age_spec_interacted"
cap mkdir "`sterdir'"

local altspec 0

*****************************************************************************
* 						PART 2. OLS Regressions                 		    *
*****************************************************************************

* 1. run regressions
* specification 

reghdfe deathrate_w99 c.tavg_poly_1_GMFD#i.agegroup c.tavg_poly_2_GMFD#i.agegroup ///
        c.tavg_poly_3_GMFD#i.agegroup c.tavg_poly_4_GMFD#i.agegroup ///
        c.tavg_poly_1_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.loggdppc_adm1_avg#i.agegroup ///
        c.tavg_poly_3_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.loggdppc_adm1_avg#i.agegroup ///
        c.tavg_poly_1_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
        c.tavg_poly_3_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
                , absorb(i.adm2_code#i.CHN_ts#i.agegroup  i.adm0_code#i.year#i.agegroup ) ///
                cluster(adm1_code)

estimates save "`sterdir'/agespec_interaction_response_car.ster", replace 


if (`altspec') {

	*****************************************************************************
	* 							Alternative FE OLS Regressions          		*
	*****************************************************************************

	* 1. set weighting schemes
	bysort year agegroup: egen tot_pop = total(population)
	gen weight = population / tot_pop

	* 2. run regressions
	* specification 1
	reghdfe deathrate_w99 c.tavg_poly_1_GMFD#i.agegroup c.tavg_poly_2_GMFD#i.agegroup ///
		c.tavg_poly_3_GMFD#i.agegroup c.tavg_poly_4_GMFD#i.agegroup ///
		c.tavg_poly_1_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.loggdppc_adm1_avg#i.agegroup ///
		c.tavg_poly_3_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.loggdppc_adm1_avg#i.agegroup ///	
		c.tavg_poly_1_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
		c.tavg_poly_3_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
			, absorb(i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year ) ///
			cluster(adm1_code)
	estimates save "`sterdir'/altspec/agespec_response_spec1_interacted_car.ster", replace

	* specification 3
	reghdfe deathrate_w99 c.tavg_poly_1_GMFD#i.agegroup c.tavg_poly_2_GMFD#i.agegroup ///
		c.tavg_poly_3_GMFD#i.agegroup c.tavg_poly_4_GMFD#i.agegroup ///
		c.tavg_poly_1_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.loggdppc_adm1_avg#i.agegroup ///
		c.tavg_poly_3_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.loggdppc_adm1_avg#i.agegroup ///	
		c.tavg_poly_1_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
		c.tavg_poly_3_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
			, absorb(i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup adm1_agegrp_code##c.year) ///
			cluster(adm1_code)
	estimates save "`sterdir'/altspec/agespec_response_spec3_interacted_car.ster", replace


	*****************************************************************************
	* 								FGLS Regressions                 		    *
	*****************************************************************************

	preserve

		* 1. generate weighting matrix
		bysort adm2_code agegroup: gen n = _N
		drop if n < 5
		* get residuals from 1st stage
		* specification 2
		reghdfe deathrate_w99 c.tavg_poly_1_GMFD#i.agegroup c.tavg_poly_2_GMFD#i.agegroup ///
			c.tavg_poly_3_GMFD#i.agegroup c.tavg_poly_4_GMFD#i.agegroup ///
			c.tavg_poly_1_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.loggdppc_adm1_avg#i.agegroup ///
			c.tavg_poly_3_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.loggdppc_adm1_avg#i.agegroup ///			
			c.tavg_poly_1_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
			c.tavg_poly_3_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
				, absorb(i.adm2_code#i.CHN_ts#i.agegroup   i.adm0_code#i.year#i.agegroup ) ///
				cluster(adm1_code) residual(e_hat)

		* calculate diagonal elements of the weighting matrix
		bysort adm1_code: egen omega = sd(e_hat) if e_hat != .
		gen precisionweight = 1/(omega*omega) 

		sort iso adm1_id adm2_id agegroup year

		* 2. run regressions
		* specification 4
		reghdfe deathrate_w99 c.tavg_poly_1_GMFD#i.agegroup c.tavg_poly_2_GMFD#i.agegroup ///
			c.tavg_poly_3_GMFD#i.agegroup c.tavg_poly_4_GMFD#i.agegroup ///
			c.tavg_poly_1_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.loggdppc_adm1_avg#i.agegroup ///
			c.tavg_poly_3_GMFD#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.loggdppc_adm1_avg#i.agegroup ///	
			c.tavg_poly_1_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_2_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
			c.tavg_poly_3_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_4_GMFD#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
				[pw = precisionweight] ///
				, absorb(i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup)  ///
				cluster(adm1_code) tol(1e-7) maxiter(10000)

		estimates save "`sterdir'/altspec/agespec_response_spec4_interacted_car.ster", replace
	restore

	*****************************************************************************
	* 								13-month ave T Regressions         		    *
	*****************************************************************************

	* specification 5
	reghdfe deathrate_w99 c.tavg_poly_1_GMFD_13m#i.agegroup c.tavg_poly_2_GMFD_13m#i.agegroup ///
		c.tavg_poly_3_GMFD_13m#i.agegroup c.tavg_poly_4_GMFD_13m#i.agegroup ///
		c.tavg_poly_1_GMFD_13m#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_2_GMFD_13m#c.loggdppc_adm1_avg#i.agegroup ///
		c.tavg_poly_3_GMFD_13m#c.loggdppc_adm1_avg#i.agegroup c.tavg_poly_4_GMFD_13m#c.loggdppc_adm1_avg#i.agegroup ///	
		c.tavg_poly_1_GMFD_13m#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_2_GMFD_13m#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
		c.tavg_poly_3_GMFD_13m#c.lr_tavg_GMFD_adm1_avg#i.agegroup c.tavg_poly_4_GMFD_13m#c.lr_tavg_GMFD_adm1_avg#i.agegroup ///
			, absorb(i.adm2_code#i.CHN_ts#i.agegroup i.adm0_code#i.year#i.agegroup ) ///
			cluster(adm1_code)
	estimates save "`sterdir'/altspec/agespec_response_spec5_interacted_public_car.ster", replace

}
*/