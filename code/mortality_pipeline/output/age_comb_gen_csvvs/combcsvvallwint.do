/*****************************************************************************
*  CSVV GENERATOR — Pooled (ALLREA) Age-Combined, Product-Interacted (Spec 1–5)
*  Reads:  $REPO/code/regressions/output/age_combined_ALLREA/*.ster
*  Writes: $REPO/code/regressions/output/age_combined_csvv_ALLREA/*.csvv
*****************************************************************************/

clear all
set more off

global REPO "/user/ab5405/summeraliaclimate"

local STERDIR "$REPO/code/regressions/output/age_combined_ALLREA"
local CSVVDIR "$REPO/code/regressions/output/age_combined_csvv_ALLREA"
cap mkdir "`CSVVDIR'"

forvalues SPEC = 1/5 {
    cap file close csvv

    local STERFILE  "pooled_response_spec`SPEC'_public_ALLREA.ster"
    local CSVVFILE  "pooled_response_spec`SPEC'_public_ALLREA.csvv"

    * 1) Load .ster (skip nicely if missing)
    capture estimates use "`STERDIR'/`STERFILE'"
    if _rc {
        di as yellow "Note: couldn't find `STERDIR'/`STERFILE' — skipping Spec `SPEC'."
        continue
    }

    * 2) Metadata preview (optional)
    quietly ereturn list
    ereturn display

    * 3) Open CSVV & write header/meta
    file open csvv using "`CSVVDIR'/`CSVVFILE'", write replace
    file write csvv "---" _n
    file write csvv "oneline: ALLREA pooled (age-combined) spec`SPEC' results with product interactions" _n
    file write csvv "version: MORTALITY-AGE-COMBINED-ALLREA-SPEC`SPEC'" _n
    file write csvv "dependencies: `STERFILE'" _n
    file write csvv "description: 4th-order polynomial OLS; country×age precip FE-slopes by product; subnational FE and country×year FE by product (per spec)." _n
    file write csvv "csvv-version: aliabonanno-2025-10-17" _n

    * 4) Variable descriptions (generic)
    file write csvv "variables:" _n
    file write csvv "  tavg_poly_1: T1_sum [C]" _n
    file write csvv "  tavg_poly_2: T2_sum [C^2]" _n
    file write csvv "  tavg_poly_3: T3_sum [C^3]" _n
    file write csvv "  tavg_poly_4: T4_sum [C^4]" _n
    file write csvv "  deathrate_w99: mortality rate winsorized at 99th percentile" _n
    file write csvv "..." _n

    * 5) Observations & predictor names (base poly names kept)
    file write csvv "observations" _n
    file write csvv " `e(N)'" _n

    file write csvv "prednames" _n
    file write csvv "tavg_poly_1, tavg_poly_2, tavg_poly_3, tavg_poly_4" _n

    * 6) Coefficient names (full, as in e(b) — documents product interactions)
    file write csvv "coefnames" _n
    local cn : colnames e(b)
    local cn_commas : subinstr local cn " " ", ", all
    file write csvv "`cn_commas'" _n

    * 7) Covariate/FE description by spec
    file write csvv "covarnames" _n
    local covar ""
    if `SPEC'==1 local covar "FE(by product): ADM2×CHN_ts×Age; ADM0×Year; Precip FE-slopes: (ADM0×Age)##prcp_{1,2}; Weights: pop(norm within product×year); Cluster: ADM1"
    else if `SPEC'==2 local covar "FE(by product): ADM2×CHN_ts×Age; ADM0×Year×Age; Precip FE-slopes: (ADM0×Age)##prcp_{1,2}; Weights: pop(norm within product×year); Cluster: ADM1"
    else if `SPEC'==3 local covar "FE(by product): ADM2×CHN_ts×Age; ADM0×Year×Age; Trends: (ADM1×Age)×Year; Precip FE-slopes: (ADM0×Age)##prcp_{1,2}; Weights: pop(norm within product×year); Cluster: ADM1"
    else if `SPEC'==4 local covar "Estimator: FGLS on Spec2 FE; Precision weights ∝ 1/ω² within ADM1×product; Precip FE-slopes: (ADM0×Age)##prcp_{1,2}; Cluster: ADM1"
    else if `SPEC'==5 local covar "Exposure: 13-month quartic (if applicable in ALLREA run); FE as in Spec2; Precip FE-slopes: (ADM0×Age)##prcp_{1,2}; Cluster: ADM1"
    file write csvv "`covar'" _n

    file write csvv "# note: CSVV gamma is full parameter vector; interactions (product, age) are in coefnames" _n

    * 8) Coefficients (gamma): write entire e(b) row vector
    file write csvv "gamma" _n
    matrix B = e(b)
    local k = colsof(B)
    forvalues j = 1/`k' {
        local coeff = el(B,1,`j')
        if `j'==1 {
            file write csvv "`coeff'"
        }
        else {
            file write csvv ", `coeff'"
        }
    }
    file write csvv "" _n

    * 9) Variance–covariance (gammavcv): full matrix
    file write csvv "gammavcv" _n
    matrix V = e(V)
    local n = colsof(V)
    forvalues r = 1/`n' {
        forvalues c = 1/`n' {
            local vcv = el(V,`r',`c')
            if `c'==1 file write csvv "`vcv'"
            else      file write csvv ", `vcv'"
        }
        file write csvv "" _n
    }

    * 10) Residual variance
    file write csvv "residvcv" _n
    file write csvv "`= (e(rmse))^2'" _n

    file close csvv
    di as txt "Wrote `CSVVFILE'"
}

di as res "CSVV generation finished for ALLREA (age-combined, product-interacted)"
