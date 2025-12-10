************************************************************
* 3_FD_interacted_all.do
* Replicate CIL stacked "interacted" FD spec for all products,
* dropping TINV HDD/CDD terms, and export coeff/vcov.
************************************************************
clear all
set more off

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global DATA    "$PROJECT/data"
global REGD    "$DATA/regression"
global OUTDIR  "$DATA/regression/sters"
cap mkdir "$OUTDIR"

* Choose model / submodel names as you like;
* these only really affect ster FILE naming (not est names).
local model    "TINV_clim"
local submodel "lininter"   // or "quadinter" or "" as needed

if "`submodel'" != "" local model_name = "`model'_`submodel'"
else local model_name = "`model'"

local products "GMFD ERA5 JRA_3Q MERRA2"

foreach product of local products {

    di "==============================================="
    di "  Interacted FD spec for `product' (`model_name')"
    di "==============================================="

    * suffix only used for FILE names, NOT variable names
    local suf = "`product'"

    * Load FD dataset
    use "$REGD/`product'_TINV_clim_regsort.dta", clear

    * time set
    sort region_i year
    tset region_i year

    * -------- long run income x income group (FD_I*lgdppc_MA15) ----------
    local lgdppc_MA15_r ""
    forval pg = 1/2 {
        forval lg = 1/2 {
            local lgdppc_MA15_r = ///
                "`lgdppc_MA15_r' c.indp`pg'#c.indf1#c.FD_I`lg'lgdppc_MA15"
        }
    }

    * -------- large income group dummies ----------
    forval pg = 1/2 {
        forval lg = 1/2 {
            gen DumIncG`lg'F1P`pg' = FD_largeind`lg' * indf1 * indp`pg'
        }
    }

    * -------- precip: indp × indf1 × FD_precipk (NO product suffix) ----------
    local precip_r ""
    forval pg = 1/2 {
        forval k  = 1/2 {
            local precip_r = ///
                "`precip_r' c.indp`pg'#c.indf1#c.FD_precip`k'"
        }
    }

    * -------- temp: indp × indf1 × FD_tempk (NO product suffix) ----------
    local temp_r ""
    forval pg = 1/2 {
        forval k  = 1/2 {   // interacted spec uses k=1/2
            local temp_r = ///
                "`temp_r' c.indp`pg'#c.indf1#c.FD_temp`k'"
        }
    }

    * -------- temp x long-run climate (TINV HDD/CDD) ----------
    * You decided NOT to use these, so leave empty.
    local climate_r ""

    * -------- temp x income spline (FD_dc1_lgdppc_MA15I`lg'temp`k') ----------
    * These exist without product suffix: FD_dc1_lgdppc_MA15I1temp1, ...
    local income_spline_r ""
    forval pg = 1/2 {
        forval lg = 1/2 {
            forval k  = 1/2 {
                local income_spline_r = ///
                    "`income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`lg'temp`k'"
            }
        }
    }

    * -------- temp x year (optional, for lininter / quadinter) ----------
    local year_temp_r ""
    if ("`submodel'" == "lininter" | "`submodel'" == "quadinter") {
        forval pg = 1/2 {
            forval k  = 1/2 {
                * Your FD builder created FD_yeartemp1, FD_yeartemp2 (no suffix)
                local year_temp_r = ///
                    "`year_temp_r' c.indp`pg'#c.indf1#c.FD_yeartemp`k'"
            }
        }
    }

    * -------- temp x year x income spline ----------
    local year_income_spline_r ""
    if ("`submodel'" == "lininter" | "`submodel'" == "quadinter") {
        forval pg = 1/2 {
            forval lg = 1/2 {
                forval k  = 1/2 {
                    local year_income_spline_r = ///
                        "`year_income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`lg'temp`k'"
                }
            }
        }
    }

    * -------- temp x year^2 parts if quadinter ----------
    if ("`submodel'" == "quadinter") {

        * add year^2 × temp: FD_year2temp`k'
        forval pg = 1/2 {
            forval k  = 1/2 {
                local year_temp_r = ///
                    "`year_temp_r' c.indp`pg'#c.indf1#c.FD_year2temp`k'"
            }
        }

        * add year^2 × temp × income spline
        forval pg = 1/2 {
            forval lg = 1/2 {
                forval k  = 1/2 {
                    local year_income_spline_r = ///
                        "`year_income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`lg'temp`k'"
                }
            }
        }
    }

    * ------------------------------------------------------
    * First-stage interacted regression
    * ------------------------------------------------------
    reghdfe FD_load_pc ///
        `temp_r' `precip_r' `climate_r' ///
        `lgdppc_MA15_r' `income_spline_r' ///
        `year_temp_r' `year_income_spline_r' ///
        DumInc*, absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        cluster(region_i) residuals(resid)

    * Use SHORT estimation name so Stata is happy
    estimates store FDint_`product'
    estimates save "$OUTDIR/FD_inter_`model_name'_`product'", replace

    * ------------------------------------------------------
    * FGLS weights (omega-by-region_i)
    * ------------------------------------------------------
    drop if resid == .

    bysort region_i: egen sd_resid = sd(resid)
    gen omega = sd_resid^2
    drop sd_resid

    gen weight = 1/omega
    drop resid

    * ------------------------------------------------------
    * Second-stage FGLS interacted regression
    * ------------------------------------------------------
    reghdfe FD_load_pc ///
        `temp_r' `precip_r' `climate_r' ///
        `lgdppc_MA15_r' `income_spline_r' ///
        `year_temp_r' `year_income_spline_r' ///
        DumInc* [pw = weight], ///
        absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        cluster(region_i)

    estimates store FDintF_`product'
    estimates save "$OUTDIR/FD_FGLS_inter_`model_name'_`product'", replace

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
            "$OUTDIR/FD_FGLS_inter_`model_name'_`product'_coeff.csv", replace
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
                replace i = `ii' in `r'
                replace j = `jj' in `r'

                local name_i : word `ii' of `names'
                local name_j : word `jj' of `names'
                replace parm_i = "`name_i'" in `r'
                replace parm_j = "`name_j'" in `r'

                replace v = V[`ii',`jj'] in `r'

                local ++r
            }
        }

        drop i j
        order parm_i parm_j v
        export delimited using ///
            "$OUTDIR/FD_FGLS_inter_`model_name'_`product'_vcov_long.csv", replace
    restore

}

************************************************************
* End 3_FD_interacted_all.do
************************************************************
