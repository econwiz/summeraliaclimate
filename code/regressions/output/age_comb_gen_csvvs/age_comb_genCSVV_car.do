/*

Purpose: Converts ster file into response function configuration file for the
Climate Impact Lab projection system. Commonly referred to in documentation as a
"CSVV".

The CSVV generated is for the main age-spec interacted model (Spec 2), that is 
the one that is carried through to the projection system/through the rest of 
the paper. 

Inputs
------

- `data/1_estimation/1_ster/age_spec_interacted`
	- `Agespec_interaction_response.ster` - Ster file containing results from
	an age-stacked regression interacted with ADM1 average income and climate. 

Outputs
-------

- `output/1_estimation/2_csvv`
	- `Agespec_interaction_response.csvv` - File containing latex output in xlsx format.

Notes
------

Summary of models:
    1. 4th-order polynomial OLS (Age x ADM2) & (Age x ADM2) FE
   *2. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE
    3. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) (Age x ADM1 linear trend)
	4. 4th-order polynomial FGLS (Age x ADM2) & (AGE x Country x Year) FE
	5. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE with 13-month climate exposure
* indicates preferred model.

NOTE: As we only carry forward the perferred model (Spec. 2) through the rest of
the analysis, this code generates a CSVV file only for that model. To produce
CSVVs for the other specifications, specify the ster file in Part 1 below.

*/


*****************************************************************************
* 						PART 0. Initializing		 						*
*****************************************************************************
global REPO "/user/ab5405/summeraliaclimate"

* note: in Stata locals, you don’t use “=”
local STERDIR   "$REPO/code/regressions/output/age_combined_car"
local CSVVDIR   "$REPO/code/regressions/output/age_combined_car_csvv"

cap mkdir "`CSVVDIR'"

//── 1. Loop over your 5 specs ────────────────────────────────────────────────
forvalues SPEC = 1/5 {
    local STERFILE = "pooled_response_spec`SPEC'_public_car.ster"
    local CSVVFILE = "pooled_response_spec`SPEC'_public_car.csvv"

    *****************************************************************************
    *  PART 1. Import & inspect .ster
    *****************************************************************************
    estimates use "`STERDIR'/`STERFILE'"
    ereturn display

    *****************************************************************************
    *  PART 2. Open CSVV & write metadata
    *****************************************************************************
    file open csvv using "`CSVVDIR'/`CSVVFILE'", write replace
    file write csvv "---" _n
    file write csvv "oneline: Mortality global age‑combined spec`SPEC' `PRODUCT' results" _n
    file write csvv "version: MORTALITY-AGE-COMBINED-SPEC`SPEC'-`PRODUCT'" _n
    file write csvv "dependencies: `STERFILE'" _n
    file write csvv "description: 4th-order polynomial OLS, population‑weighted; ADM2 & ADM0×Year FE" _n
    file write csvv "csvv-version: aliabonanno-2025-07-23" _n

    *****************************************************************************
    *  PART 3. Variable descriptions
    *****************************************************************************
    file write csvv "variables:" _n
    file write csvv "  tavg_poly_1: T1_sum [C]" _n
    file write csvv "  tavg_poly_2: T2_sum [C^2]" _n
    file write csvv "  tavg_poly_3: T3_sum [C^3]" _n
    file write csvv "  tavg_poly_4: T4_sum [C^4]" _n
    file write csvv "  lr_tavg: yearly avg tas [C]" _n
    file write csvv "  deathrate_w99: death rate winsorized at 99th percentile" _n
    file write csvv "  loggdppc_adm0_avg: log GDP per capita" _n
    file write csvv "..." _n

    *****************************************************************************
    *  PART 4. Observations & predictors
    *****************************************************************************
    file write csvv "observations" _n
    file write csvv " `e(N)'" _n

    file write csvv "prednames" _n
    file write csvv "tavg_poly_1, tavg_poly_2, tavg_poly_3, tavg_poly_4" _n

    file write csvv "covarnames" _n
    file write csvv "1, lr_tavg, loggdppc_adm0_avg" _n

    *****************************************************************************
    *  PART 5. Coefficients (all 36 gammas across 3 age groups × 4 polys)
    *****************************************************************************
    file write csvv "gamma" _n
    matrix B = e(b)
    forvalues row = 1/36 {
        local coeff = el(B,1,`row')

        // if this is the very first overall (row=1) OR
        // it's the first in a new block of 12 ((row−1)%12==0),
        // write it without a comma prefix
        if (`row'==1) | (mod(`row'-1,12)==0) {
            file write csvv "`coeff'"
        }
        else {
            file write csvv ", `coeff'"
        }

        // after every 12 entries, emit a newline
        if mod(`row',12)==0 {
            file write csvv "" _n
        }
    }
    // (no extra newline needed if 36 is a multiple of 12)


    *****************************************************************************
    *  PART 6. Variance–covariance matrix (36×36)
    *****************************************************************************
    file write csvv "gammavcv" _n
    matrix V = e(V)
    forvalues r = 1/36 {
        forvalues c = 1/36 {
            local vcv = el(V,`r',`c')
            if `c'==1    file write csvv "`vcv'"
            else         file write csvv ", `vcv'"
        }
        file write csvv "" _n
    }

    *****************************************************************************
    *  PART 7. Residual variance
    *****************************************************************************
    file write csvv "residvcv" _n
    local rmse = e(rmse)
    file write csvv "`= `rmse'^2'" _n

    *****************************************************************************
    *  PART 8. Close
    *****************************************************************************
    file close csvv
    di as txt "✅ Wrote `CSVVFILE' in `CSVVDIR'"
}

