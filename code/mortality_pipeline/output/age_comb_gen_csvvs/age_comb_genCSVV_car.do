/*****************************************************************************
*  CSVV GENERATOR FOR AGE‑COMBINED GMFD REGRESSIONS (Specs 1–5)             *
*****************************************************************************/

/* PART 0. INITIALIZING */
global REPO   "/user/ab5405/summeraliaclimate"
local PRODUCT "GMFD"
local STERDIR "$REPO/code/regressions/output/age_combined_car"
local CSVVDIR "$REPO/code/regressions/output/age_combined_csvv"

cap mkdir "`CSVVDIR'"

/* PART 1–8: Loop over specs 1 through 5 */
forvalues SPEC = 1/5 {
    /* ensure no stale file handle */
    cap file close csvv

    /* define file names */
    local STERFILE  = "pooled_response_spec`SPEC'_public_car.ster"
    local CSVVFILE  = "pooled_response_spec`SPEC'_public_car.csvv"

    /* 1. load .ster */
    estimates use "`STERDIR'/`STERFILE'"
    ereturn display

    /* 2. open CSVV & write metadata */
    file open csvv using "`CSVVDIR'/`CSVVFILE'", write replace
    file write csvv "---" _n
    file write csvv "oneline: Age‑combined GMFD spec`SPEC' results" _n
    file write csvv "version: MORTALITY-AGE-COMBINED-GMFD-SPEC`SPEC'" _n
    file write csvv "dependencies: `STERFILE'" _n
    file write csvv "description: 4th-order poly OLS; ADM2×Seasonality & ADM0×Year FE" _n
    file write csvv "csvv-version: aliabonanno-2025-07-28" _n

    /* 3. variable descriptions */
    file write csvv "variables:" _n
    file write csvv "  tavg_poly_1: T1_sum [C]" _n
    file write csvv "  tavg_poly_2: T2_sum [C^2]" _n
    file write csvv "  tavg_poly_3: T3_sum [C^3]" _n
    file write csvv "  tavg_poly_4: T4_sum [C^4]" _n
    file write csvv "  deathrate_w99: mortality rate winnowed at 99th percentile" _n
    file write csvv "..." _n

    /* 4. observations & predictors */
    file write csvv "observations" _n
    file write csvv " `e(N)'" _n
    file write csvv "prednames" _n
    file write csvv "tavg_poly_1, tavg_poly_2, tavg_poly_3, tavg_poly_4" _n
   /* coefnames: comma-separated names exactly as in e(b), including _cons last */
    file write csvv "coefnames" _n
    local cn : colnames e(b)
    local cn_commas : subinstr local cn " " ", ", all
    file write csvv "`cn_commas'" _n
    
    /* covarnames: human-readable description (metadata only) */
    file write csvv "covarnames" _n
    if `SPEC'==1 {
        file write csvv "FE: ADM2×CHN_ts×Age; FE: ADM0×Year; Controls: ADM0×prcp_poly_{1,2} (GMFD); Weights: population; Clusters: ADM1" _n
    }
    else if `SPEC'==2 {
        file write csvv "FE: ADM2×CHN_ts×Age; FE: ADM0×Year×Age; Controls: ADM0×prcp_poly_{1,2} (GMFD); Weights: population; Clusters: ADM1" _n
    }
    else if `SPEC'==3 {
        file write csvv "FE: ADM2×CHN_ts×Age; FE: ADM0×Year×Age; Trend: ADM1×Age linear; Controls: ADM0×prcp_poly_{1,2} (GMFD)" _n
    }
    else if `SPEC'==4 {
        file write csvv "Estimator: FGLS (Spec2 FE); Weights: precision (pop×1/ω²); Controls: ADM0×prcp_poly_{1,2} (GMFD)" _n
    }
    else if `SPEC'==5 {
        file write csvv "Exposure: 13-month; FE: as Spec2; Controls: ADM0×prcp_poly_{1,2} (GMFD)" _n
    }
    else {
        file write csvv "" _n   // leave empty if not applicable
    }
    
    /* (optional) make intercept convention explicit */
    file write csvv "# note: gamma includes the intercept (_cons) as the last entry" _n

    /* 5. coefficients (gamma) */
    file write csvv "gamma" _n
    matrix B = e(b)
    local n = colsof(B)
    forvalues row = 1/`n' {
        local coeff = el(B,1,`row')
        if (`row' == 1) | (mod(`row'-1,12) == 0) {
            file write csvv "`coeff'"
        }
        else {
            file write csvv ", `coeff'"
        }
        if mod(`row',12) == 0 & `row' < `n' {
            file write csvv "" _n
        }
    }
    file write csvv "" _n

    /* 6. variance–covariance (gammavcv) */
    file write csvv "gammavcv" _n
    matrix V = e(V)
    local n = colsof(V)
    forvalues r = 1/`n' {
        forvalues c = 1/`n' {
            local vcv = el(V,`r',`c')
            if `c' == 1 {
                file write csvv "`vcv'"
            }
            else {
                file write csvv ", `vcv'"
            }
        }
        file write csvv "" _n
    }

    /* 7. residual variance */
    file write csvv "residvcv" _n
    local rmse = e(rmse)
    file write csvv "`= `rmse'^2'" _n

    /* 8. close CSVV */
    file close csvv
    di as txt "✅ Wrote `CSVVFILE' to `CSVVDIR' for car spec`SPEC'"
}
