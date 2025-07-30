/*

Purpose: Converts ster file into response function configuration file for the
Climate Impact Lab projection system. Commonly referred to in documentation as a
"CSVV".

The CSVV generated is for the main age-spec interacted model (Spec 2), that is 
the one that is carried through to the projection system/through the rest of 
the paper. 

Inputs
------

- `output/age_spec_interacted'
    - `agespec_interaction_response_ERA5.ster`` - Ster file containing results from
    an age-stacked regression interacted with ADM1 average income and climate. 

Outputs
-------

- `output/age_spec_interacted_csvv`
    - `agespec_interaction_response_spec2_ERA5_025.csvv` - File containing latex output in xlsx format.

Notes
------

Summary of models:
    1. 4th-order polynomial OLS (Age x ADM2) & (Age x ADM2) FE
   *2. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE
    3. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) (Age x ADM1 linear trend)
    4. 4th-order polynomial FGLS (Age x ADM2) & (AGE x Country x Year) FE
    5. 4th-order polynomial OLS (Age x ADM2) & (AGE x Country x Year) FE with 13-month climate exposure
* indicates preferred model.

NOTE: As we only carry forward the preferred model (Spec. 2) through the rest of
the analysis, this code generates a CSVV file only for that model. To produce
CSVVs for the other specifications, specify the ster file in Part 1 below.
*/


/*****************************************************************************
*  PART 0. INITIALIZING                                                    *
*****************************************************************************/
global REPO "/user/ab5405/summeraliaclimate"
local PRODUCT  "ERA5_025"
local STERDIR  "$REPO/code/regressions/output/age_spec_interacted"
local CSVVDIR  "$REPO/code/regressions/output/age_spec_interacted_csvv"

cap mkdir "`CSVVDIR'"

* We only carry forward Spec 2 for ERA5_025:
local STERFILE "agespec_interaction_response_ERA5_025.ster"
local CSVVFILE "agespec_interaction_response_spec2_ERA5_025.csvv"


/*****************************************************************************
*  PART 1. LOAD THE SPEC 2 STER FILE                                        *
*****************************************************************************/
estimates use "`STERDIR'/`STERFILE'"
ereturn display


/*****************************************************************************
*  PART 2. OPEN CSVV & WRITE METADATA                                       *
*****************************************************************************/
cap file close csvv
file open csvv using "`CSVVDIR'/`CSVVFILE'", write replace

file write csvv "---" _n
file write csvv "oneline: Age‑specific interacted (Spec 2) response – ERA5_025" _n
file write csvv "version: MORTALITY-AGE-SPEC-INTERACTED-ERA5_025" _n
file write csvv "dependencies: `STERFILE'" _n
file write csvv "description: 4th-order poly OLS; age×ADM2 & age×Country×Year FE; temp×income & temp×LR‑temp interactions" _n
file write csvv "csvv-version: aliabonanno-2025-07-28" _n


/*****************************************************************************
*  PART 3. VARIABLE DESCRIPTIONS                                            *
*****************************************************************************/
file write csvv "variables:" _n
file write csvv "  tavg_poly_1: T1_sum [C]" _n
file write csvv "  tavg_poly_2: T2_sum [C^2]" _n
file write csvv "  tavg_poly_3: T3_sum [C^3]" _n
file write csvv "  tavg_poly_4: T4_sum [C^4]" _n
file write csvv "  lr_tavg: long‑run mean annual temperature [C]" _n
file write csvv "  deathrate_w99: mortality rate winsorized at 99th percentile" _n
file write csvv "  loggdppc_adm1_avg: log GDP per capita" _n
file write csvv "..." _n


/*****************************************************************************
*  PART 4. OBSERVATIONS & PREDICTORS                                         *
*****************************************************************************/
file write csvv "observations" _n
file write csvv " `e(N)'" _n

file write csvv "prednames" _n
file write csvv "tavg_poly_1, tavg_poly_2, tavg_poly_3, tavg_poly_4" _n

file write csvv "covarnames" _n
file write csvv "1, lr_tavg, loggdppc_adm1_avg" _n


/*****************************************************************************
*  PART 5. COEFFICIENTS (36 GAMMAS: 3 AGES × 4 POLYS × 3 AGE GROUPS)         *
*****************************************************************************/
file write csvv "gamma" _n
matrix B = e(b)
forvalues row = 1/36 {
    local coeff = el(B,1,`row')
    if (`row' == 1) | (mod(`row'-1,12) == 0) {
        file write csvv "`coeff'"
    }
    else {
        file write csvv ", `coeff'"
    }
    if mod(`row',12) == 0 {
        file write csvv "" _n
    }
}


****************************************************************************
*  PART 6. VARIANCE–COVARIANCE MATRIX (36×36)                              *
****************************************************************************/
file write csvv "gammavcv" _n
matrix V = e(V)
forvalues r = 1/36 {
    forvalues c = 1/36 {
        local vcf = el(V,`r',`c')
        if `c' == 1 {
            file write csvv "`vcf'"
        }
        else {
            file write csvv ", `vcf'"
        }
    }
    file write csvv "" _n
}


****************************************************************************
*  PART 7. RESIDUAL VARIANCE                                                *
****************************************************************************/
file write csvv "residvcv" _n
local rmse = e(rmse)
file write csvv "`= `rmse'^2'" _n


/*****************************************************************************
*  PART 8. CLOSE CSVV & FINISH                                             *
*****************************************************************************/
file close csvv
di as txt "✅ Wrote `CSVVFILE' to `CSVVDIR'"
