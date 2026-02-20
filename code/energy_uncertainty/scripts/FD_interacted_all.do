/************************************************************
* 3_FD_interacted_all.do
* Replicate CIL stacked "interacted" FD spec for all products
* (GMFD, ERA5, JRA_3Q, MERRA2), matching original CIL regression
* logic (unweighted first stage, sd-based FGLS weights).
*
* Saves per product:
*   $OUTDIR/FD_inter_`model_name'_`product'.ster       (first stage)
*   $OUTDIR/FD_FGLS_inter_`model_name'_`product'.ster  (FGLS second stage)
*
* Also exports per product:
*   $OUTDIR/FD_FGLS_inter_`model_name'_`product'_coeff.csv
*   $OUTDIR/FD_FGLS_inter_`model_name'_`product'_vcov_long.csv
*
* Optional:
*   complete_for_plot = 1 reposts e(b), e(V) to include omitted
*   terms as 0 so predictnl/_b[...] works even if collinear/omitted.
************************************************************/

clear all
set more off

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global DATA    "$PROJECT/data"
global REGD    "$DATA/regression"
global OUTDIR  "$DATA/regression/sters"
cap mkdir "$OUTDIR"

* Model naming: submodel controls whether year and year^2 interactions
* are included. Options: "" (base), "lininter" (+ year), "quadinter" (+ year^2)
local model    "TINV_clim"
local submodel "lininter"

if "`submodel'" != "" local model_name "`model'_`submodel'"
else                   local model_name "`model'"

* 1 = pad e(b)/e(V) with zeros for omitted terms (safe for predictnl/_b[...])
* 0 = save default Stata estimates
local complete_for_plot 0

local products "GMFD ERA5 JRA_3Q MERRA2"

foreach product of local products {

    di "==============================================="
    di "  Interacted FD spec for `product' (`model_name')"
    di "==============================================="

    * Load FD dataset for this product
    use "$REGD/`product'_TINV_clim_regsort.dta", clear

    * Set panel structure (match original CIL: tset not xtset)
    sort region_i year
    tset region_i year

    * ------------------------------------------------------------------
    * Build regressor locals. All FD variables are UNSUFFIXED — the
    * product-specific suffix was stripped when building the FD dataset
    * in 0_build_FD_dataset_all.do, so no _GMFD/_ERA5 etc. here.
    * ------------------------------------------------------------------

    * Long-run income × income group
    local lgdppc_MA15_r ""
    forval pg = 1/2 {
        forval lg = 1/2 {
            local lgdppc_MA15_r ///
                "`lgdppc_MA15_r' c.indp`pg'#c.indf1#c.FD_I`lg'lgdppc_MA15"
        }
    }

    * Large income group dummies (FD of indicator × product group interactions)
    forval pg = 1/2 {
        forval lg = 1/2 {
            gen DumIncG`lg'F1P`pg' = FD_largeind`lg' * indf1 * indp`pg'
        }
    }

    * Precipitation: indp × indf1 × FD_precipk
    local precip_r ""
    forval pg = 1/2 {
        forval k = 1/2 {
            local precip_r "`precip_r' c.indp`pg'#c.indf1#c.FD_precip`k'"
        }
    }

    * Temperature quartic: indp × indf1 × FD_tempk
    local temp_r ""
    forval pg = 1/2 {
        forval k = 1/4 {
            local temp_r "`temp_r' c.indp`pg'#c.indf1#c.FD_temp`k'"
        }
    }

    * Temperature × long-run climate (TINV cross terms):
    * indp × indf1 × FD_hdd20_TINVtemp`k' and FD_cdd20_TINVtemp`k'.
    * The redundant lg loop matches the CIL original exactly.
    local climate_r ""
    forval pg = 1/2 {
        forval lg = 1/2 {
            forval k = 1/2 {
                local climate_r "`climate_r' c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k'"
                local climate_r "`climate_r' c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k'"
            }
        }
    }

    * Temperature × income spline: indp × indf1 × FD_dc1_lgdppc_MA15I`lg'temp`k'
    local income_spline_r ""
    forval pg = 1/2 {
        forval lg = 1/2 {
            forval k = 1/2 {
                local income_spline_r ///
                    "`income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`lg'temp`k'"
            }
        }
    }

    * Temperature × year and temperature × year × income spline
    * (added for lininter and quadinter submodels)
    local year_temp_r ""
    local year_income_spline_r ""

    if ("`submodel'" == "lininter" | "`submodel'" == "quadinter") {

        * temp × year
        forval pg = 1/2 {
            forval k = 1/2 {
                local year_temp_r ///
                    "`year_temp_r' c.indp`pg'#c.indf1#c.FD_yeartemp`k'"
            }
        }

        * temp × year × income spline
        forval pg = 1/2 {
            forval lg = 1/2 {
                forval k = 1/2 {
                    local year_income_spline_r ///
                        "`year_income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`lg'temp`k'"
                }
            }
        }
    }

    * Temperature × year^2 and temperature × year^2 × income spline
    * (added for quadinter submodel only)
    if ("`submodel'" == "quadinter") {

        * temp × year^2
        forval pg = 1/2 {
            forval k = 1/2 {
                local year_temp_r ///
                    "`year_temp_r' c.indp`pg'#c.indf1#c.FD_year2temp`k'"
            }
        }

        * temp × year^2 × income spline
        forval pg = 1/2 {
            forval lg = 1/2 {
                forval k = 1/2 {
                    local year_income_spline_r ///
                        "`year_income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`lg'temp`k'"
                }
            }
        }
    }

    * ------------------------------------------------------------------
    * First stage: unweighted FD regression (match CIL original —
    * no population weights at this stage, only in FGLS second stage).
    * Residuals are saved for FGLS weight construction below.
    * ------------------------------------------------------------------
    capture drop resid
    reghdfe FD_load_pc ///
        `temp_r' `precip_r' `climate_r' ///
        `lgdppc_MA15_r' `income_spline_r' ///
        `year_temp_r' `year_income_spline_r' ///
        DumInc*, ///
        absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        vce(cluster region_i) residuals(resid)

    estimates save "$OUTDIR/FD_inter_`model_name'_`product'", replace

    * ------------------------------------------------------------------
    * FGLS weight construction (match CIL original):
    * omega = sd(resid)^2 within each region; FGLS weight = 1/omega.
    * Drop singletons (weight undefined for regions with only one obs).
    * ------------------------------------------------------------------
    drop if resid == .

    bysort region_i: gen count = _N
    drop if count == 1

    bysort region_i: egen sd_resid = sd(resid)
    gen omega = sd_resid ^ 2
    drop sd_resid
    gen FGLS_weight = 1 / omega

    * ------------------------------------------------------------------
    * Second stage: FGLS regression with 1/omega weights
    * ------------------------------------------------------------------
    reghdfe FD_load_pc ///
        `temp_r' `precip_r' `climate_r' ///
        `lgdppc_MA15_r' `income_spline_r' ///
        `year_temp_r' `year_income_spline_r' ///
        DumInc* [pw = FGLS_weight], ///
        absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        vce(cluster region_i)

    di "e(V) trace: " trace(e(V))
    di "e(V) rows: " rowsof(e(V))
    * ------------------------------------------------------------------
    * Optional: pad e(b) and e(V) with zeros for omitted/collinear terms
    * so that predictnl and _b[...] lookups work for all terms.
    * Only applies to the saved .ster — does not affect the regression.
    * ------------------------------------------------------------------
    if (`complete_for_plot') {

        local want_all ""
        local want_all "`want_all' `temp_r' `precip_r' `climate_r'"
        local want_all "`want_all' `lgdppc_MA15_r' `income_spline_r'"
        local want_all "`want_all' `year_temp_r' `year_income_spline_r'"
        unab dumlist : DumInc*
        local want_all "`want_all' `dumlist'"

        * Deduplicate
        local want ""
        foreach nm of local want_all {
            if strpos(" `want' ", " `nm' ") == 0 local want "`want' `nm'"
        }
        local K_want : word count `want'

        tempname b V bfull Vfull
        matrix `b' = e(b)
        matrix `V' = e(V)

        matrix `bfull' = J(1, `K_want', 0)
        matrix colnames `bfull' = `want'

        matrix `Vfull' = J(`K_want', `K_want', 0)
        matrix rownames `Vfull' = `want'
        matrix colnames `Vfull' = `want'

        local i = 0
        foreach ni of local want {
            local ++i
            local ci = colnumb(`b', "`ni'")
            if (`ci' < .) matrix `bfull'[1,`i'] = `b'[1,`ci']

            local ri = rownumb(`V', "`ni'")
            if (`ri' < .) {
                local j = 0
                foreach nj of local want {
                    local ++j
                    local cj = colnumb(`V', "`nj'")
                    if (`cj' < .) matrix `Vfull'[`i',`j'] = `V'[`ri',`cj']
                }
            }
        }

        ereturn repost b=`bfull' V=`Vfull', rename
    }

    estimates save "$OUTDIR/FD_FGLS_inter_`model_name'_`product'", replace

    * ------------------------------------------------------------------
    * Export coefficient vector (parm, beta) to CSV
    * ------------------------------------------------------------------
    matrix b = e(b)
    local kb = colsof(b)
    local bnames : colnames b

    preserve
        clear
        set obs `kb'
        gen strL parm = ""
        gen double beta = .
        local j = 1
        foreach nm of local bnames {
            replace parm = "`nm'" in `j'
            replace beta = b[1,`j'] in `j'
            local ++j
        }
        order parm beta
        export delimited using ///
            "$OUTDIR/FD_FGLS_inter_`model_name'_`product'_coeff.csv", replace
    restore

    * ------------------------------------------------------------------
    * Export variance-covariance matrix (long format: parm_i, parm_j, v)
    * ------------------------------------------------------------------
    matrix V = e(V)
    local kv = rowsof(V)
    local vnames : rownames V

    preserve
        clear
        set obs `=`kv' * `kv''
        gen strL parm_i = ""
        gen strL parm_j = ""
        gen double v = .
        local r = 1
        forvalues ii = 1/`kv' {
            forvalues jj = 1/`kv' {
                local name_i : word `ii' of `vnames'
                local name_j : word `jj' of `vnames'
                replace parm_i = "`name_i'" in `r'
                replace parm_j = "`name_j'" in `r'
                replace v = V[`ii',`jj'] in `r'
                local ++r
            }
        }
        order parm_i parm_j v
        export delimited using ///
            "$OUTDIR/FD_FGLS_inter_`model_name'_`product'_vcov_long.csv", replace
    restore

    di as res "Finished `product': sters and CSVs saved to $OUTDIR"

} // end foreach product

di "All products processed."

/************************************************************
* End 3_FD_interacted_all.do
************************************************************/
