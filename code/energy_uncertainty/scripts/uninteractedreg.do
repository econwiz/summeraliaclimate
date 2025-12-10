************************************************************
* 2_FD_uninteracted_all.do
* Replicate CIL stacked "uninteracted" FD spec for all products
* (GMFD, ERA5, JRA_3Q, MERRA2), with pop weights + FGLS,
* AND export coeffs / vcov to CSV.
************************************************************
clear all
set more off

* Paths
global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global DATA    "$PROJECT/data"
global REGD    "$DATA/regression"
global OUTDIR  "$DATA/regression/sters"
cap mkdir "$OUTDIR"

* Products to loop over
local products "GMFD ERA5 JRA_3Q MERRA2"

foreach product of local products {

    di "==============================================="
    di "  Uninteracted FD spec for `product'"
    di "==============================================="

    * suffix for climate variables
    local suf = "`product'"

    * 1) Load FD dataset for this product
    use "$REGD/`product'_TINV_clim_regsort.dta", clear

    * 2) Population weights (same as your original)
    bysort year product flow: egen year_product_flow_total_pop = total(pop)
    gen pop_weight = pop / year_product_flow_total_pop

    * 3) Set panel structure
    sort region_i year
    tset region_i year

    * 4) Precip block: indp × indf1 × FD_precipk  (NO product suffix on FD vars)
    local precip_r ""
    forval pg = 1/2 {
        forval k = 1/2 {
            local precip_r "`precip_r' c.indp`pg'#c.indf1#c.FD_precip`k'"
        }
    }
    
    * 5) Temp block: indp × indf1 × FD_tempk  (NO product suffix on FD vars)
    local temp_r ""
    forval pg = 1/2 {
        forval k = 1/4 {
            local temp_r "`temp_r' c.indp`pg'#c.indf1#c.FD_temp`k'"
        }
    }

    * ------------------------------------------------------
    * First-stage: pop-weighted FD regression
    * ------------------------------------------------------
    reghdfe FD_load_pc `temp_r' `precip_r' [pw = pop_weight], ///
        absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        cluster(region_i) residuals(resid)

    estimates store FD_global_`product'
    estimates save "$OUTDIR/FD_global_`product'", replace

    * ------------------------------------------------------
    * FGLS weights (your original variance weighting logic)
    * ------------------------------------------------------
    drop if resid == .

    bysort region_i: gen count = _N

    * relative pop weights within FE
    bysort region_i: egen sum_of_weights_in_FE = total(pop_weight)
    gen for_variance_weighting = pop_weight / sum_of_weights_in_FE

    * weighted residual mean within FE
    gen weighted_residual = for_variance_weighting * resid
    bysort region_i: egen weighted_mean_resid_FE_level = mean(weighted_residual)

    * weighted within-FE variance
    gen square_term_weighted = for_variance_weighting * ///
        (resid - weighted_mean_resid_FE_level)^2
    bysort region_i: egen weighted_residual_variance = total(square_term_weighted)

    * FGLS weights
    gen FGLS_weight = pop_weight / weighted_residual_variance

    * drop singletons (weight numerically unstable)
    drop if count == 1

    * ------------------------------------------------------
    * Second-stage FGLS regression
    * ------------------------------------------------------
    reghdfe FD_load_pc `temp_r' `precip_r' [pw = FGLS_weight], ///
        absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        cluster(region_i)

    estimates store FD_FGLS_global_`product'
    estimates save "$OUTDIR/FD_FGLS_global_`product'", replace

    * ------------------------------------------------------
    * Export coeffs and vcov for this PRODUCT (long format)
    * ------------------------------------------------------

    * 1) Coefficient vector: parm, beta
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

    * 2) Variance–covariance matrix: parm_i, parm_j, v
    matrix V = e(V)
    local k = rowsof(V)
    local names : rownames V

    preserve
        clear
        * We’ll create k^2 rows: one per (i,j) pair
        set obs `=`k'*`k''

        gen int i = .
        gen int j = .
        gen strL parm_i = ""
        gen strL parm_j = ""
        gen double v = .

        local r = 1
        forvalues ii = 1/`k' {
            forvalues jj = 1/`k' {
                * store indices
                replace i = `ii' in `r'
                replace j = `jj' in `r'

                * look up parameter names
                local name_i : word `ii' of `names'
                local name_j : word `jj' of `names'
                replace parm_i = "`name_i'" in `r'
                replace parm_j = "`name_j'" in `r'

                * vcov entry
                replace v = V[`ii',`jj'] in `r'

                local ++r
            }
        }

        drop i j
        order parm_i parm_j v
        export delimited using ///
            "$OUTDIR/FD_FGLS_global_`product'_vcov_long.csv", replace
    restore

}

************************************************************
* End 2_FD_uninteracted_all.do
************************************************************
