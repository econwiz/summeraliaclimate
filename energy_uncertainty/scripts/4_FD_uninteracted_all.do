************************************************************
* 2_FD_uninteracted_all.do
* Replicate CIL stacked "uninteracted" FD spec for all products
* (GMFD, ERA5, JRA_3Q, MERRA2), with pop weights + FGLS,
* AND export coefficients and vcov to CSV.
*
* Saves per product:
*   $OUTDIR/FD_global_`product'.ster       (first stage)
*   $OUTDIR/FD_FGLS_global_`product'.ster  (FGLS second stage)
*
* Also exports per product:
*   $OUTDIR/FD_FGLS_global_`product'_coeff.csv
*   $OUTDIR/FD_FGLS_global_`product'_vcov_long.csv
*
* NOTE on FGLS weights: the uninteracted spec uses the population-weighted
* within-FE variance (matching the CIL original). This differs from the
* interacted spec (3_FD_interacted_all.do), which uses the simpler sd()^2
* approach. Do not make these consistent without checking with your mentor.
************************************************************

clear all
set more off

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global DATA    "$PROJECT/data"
global REGD    "$DATA/regression"
global OUTDIR  "$DATA/regression/sters"
cap mkdir "$OUTDIR"

local products "GMFD ERA5 JRA_3Q MERRA2"

foreach product of local products {

    di "==============================================="
    di "  Uninteracted FD spec for `product'"
    di "==============================================="

    * Load FD dataset for this product
    use "$REGD/`product'_TINV_clim_regsort.dta", clear

    * Population weights: country share of world population within
    * each year × product × flow cell, matching the CIL original.
    bysort year product flow: egen year_product_flow_total_pop = total(pop)
    gen pop_weight = pop / year_product_flow_total_pop

    * Set panel structure (match original CIL: tset not xtset)
    sort region_i year
    tset region_i year

    * ------------------------------------------------------------------
    * Build regressor locals. All FD variables are UNSUFFIXED — the
    * product-specific suffix was stripped when building the FD dataset
    * in 0_build_FD_dataset_all.do, so no _GMFD/_ERA5 etc. here.
    * ------------------------------------------------------------------

    * Precipitation: indp × indf1 × FD_precipk
    local precip_r ""
    forval pg = 1/2 {
        forval k = 1/2 {
            local precip_r "`precip_r' c.indp`pg'#c.indf1#c.FD_precip`k'"
        }
    }

    * Temperature quartic: indp × indf1 × FD_tempk
    * NOTE: forval k = 1(1)4 is the correct Stata syntax for 1,2,3,4.
    * forval k = 1/4 would be parsed as 0.25 and the loop would not run.
    local temp_r ""
    forval pg = 1/2 {
        forval k = 1(1)4 {
            local temp_r "`temp_r' c.indp`pg'#c.indf1#c.FD_temp`k'"
        }
    }

    * ------------------------------------------------------------------
    * First stage: pop-weighted FD regression (match CIL original).
    * Residuals are saved for FGLS weight construction below.
    * ------------------------------------------------------------------
    reghdfe FD_load_pc `temp_r' `precip_r' [pw = pop_weight], ///
        absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        vce(cluster region_i) residuals(resid)

    estimates save "$OUTDIR/FD_global_`product'", replace

    * ------------------------------------------------------------------
    * FGLS weight construction (match CIL original for uninteracted spec):
    * Compute population-weighted within-region variance of residuals,
    * then set FGLS weight = pop_weight / weighted_residual_variance.
    * Drop singletons: Stata rounding can produce nonzero weights when
    * (resid - weighted_mean) should be exactly zero.
    * ------------------------------------------------------------------
    drop if resid == .

    bysort region_i: gen count = _N

    * Relative population weights within each region FE
    bysort region_i: egen sum_of_weights_in_FE = total(pop_weight)
    gen for_variance_weighting = pop_weight / sum_of_weights_in_FE

    * Weighted mean residual within each region FE
    gen weighted_residual = for_variance_weighting * resid
    bysort region_i: egen weighted_mean_resid_FE_level = mean(weighted_residual)

    * Weighted within-FE variance
    gen square_term_weighted = for_variance_weighting * ///
        (resid - weighted_mean_resid_FE_level) ^ 2
    bysort region_i: egen weighted_residual_variance = total(square_term_weighted)

    * Final FGLS weights: population weight divided by within-FE variance
    gen FGLS_weight = pop_weight / weighted_residual_variance

    * Drop singletons (weight numerically unstable due to Stata rounding)
    drop if count == 1

    * ------------------------------------------------------------------
    * Second stage: FGLS regression with population-variance weights
    * ------------------------------------------------------------------
    reghdfe FD_load_pc `temp_r' `precip_r' [pw = FGLS_weight], ///
        absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        vce(cluster region_i)

    estimates save "$OUTDIR/FD_FGLS_global_`product'", replace

    * ------------------------------------------------------------------
    * Export coefficient vector (parm, beta) to CSV
    * ------------------------------------------------------------------
    matrix b = e(b)
    local k = colsof(b)
    local names : colnames b

    preserve
        clear
        set obs `k'
        gen strL parm = ""
        gen double beta = .

        local j = 1
        foreach nm of local names {
            replace parm = "`nm'" in `j'
            replace beta = b[1,`j'] in `j'
            local ++j
        }

        order parm beta
        export delimited using ///
            "$OUTDIR/FD_FGLS_global_`product'_coeff.csv", replace
    restore

    * ------------------------------------------------------------------
    * Export variance-covariance matrix (long format: parm_i, parm_j, v)
    * ------------------------------------------------------------------
    matrix V = e(V)
    local k = rowsof(V)
    local names : rownames V

    preserve
        clear
        set obs `=`k' * `k''
        gen strL parm_i = ""
        gen strL parm_j = ""
        gen double v = .

        local r = 1
        forvalues ii = 1/`k' {
            forvalues jj = 1/`k' {
                local name_i : word `ii' of `names'
                local name_j : word `jj' of `names'
                replace parm_i = "`name_i'" in `r'
                replace parm_j = "`name_j'" in `r'
                replace v = V[`ii',`jj'] in `r'
                local ++r
            }
        }

        order parm_i parm_j v
        export delimited using ///
            "$OUTDIR/FD_FGLS_global_`product'_vcov_long.csv", replace
    restore

    di as res "Finished `product': sters and CSVs saved to $OUTDIR"

} // end foreach product

di "All products processed."

************************************************************
* End 2_FD_uninteracted_all.do
************************************************************
