************************************************************
* 3_FD_interacted_all.do - REGHDFE VCOV FIX
* PROBLEM: reghdfe doesn't populate e(V) with absorb() option
* SOLUTION: Use estadd or manually construct from _se[]
************************************************************
clear all
set more off

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global DATA    "$PROJECT/data"
global REGD    "$DATA/regression"
global OUTDIR  "$DATA/regression/sters"
cap mkdir "$OUTDIR"

* Install estadd if needed
cap which estadd
if _rc ssc install estout

local model    "TINV_clim"
local submodel "quadinter"

if "`submodel'" != "" local model_name = "`model'_`submodel'"
else local model_name = "`model'"

local products "GMFD ERA5 JRA_3Q MERRA2"

foreach product of local products {

    di "==============================================="
    di "  Interacted FD spec for `product' (`model_name')"
    di "==============================================="

    use "$REGD/`product'_TINV_clim_regsort.dta", clear

    sort region_i year
    tset region_i year

    * Build regressor lists (abbreviated)
    local lgdppc_MA15_r ""
    forval pg = 1/2 {
        forval lg = 1/2 {
            local lgdppc_MA15_r = "`lgdppc_MA15_r' c.indp`pg'#c.indf1#c.FD_I`lg'lgdppc_MA15"
        }
    }

    forval pg = 1/2 {
        forval lg = 1/2 {
            gen DumIncG`lg'F1P`pg' = FD_largeind`lg' * indf1 * indp`pg'
        }
    }

    local precip_r ""
    forval pg = 1/2 {
        forval k  = 1/2 {
            local precip_r = "`precip_r' c.indp`pg'#c.indf1#c.FD_precip`k'"
        }
    }

    local temp_r ""
    forval pg = 1/2 {
        forval k  = 1/2 {
            local temp_r = "`temp_r' c.indp`pg'#c.indf1#c.FD_temp`k'"
        }
    }

    local climate_r ""
    forval pg = 1/2 {
        forval k = 1/2 {
            local climate_r "`climate_r' c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k' "
            local climate_r "`climate_r' c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k' "
        }
    }

    local income_spline_r ""
    forval pg = 1/2 {
        forval lg = 1/2 {
            forval k  = 1/2 {
                local income_spline_r = "`income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`lg'temp`k'"
            }
        }
    }

    local year_temp_r ""
    if ("`submodel'" == "lininter" | "`submodel'" == "quadinter") {
        forval pg = 1/2 {
            forval k  = 1/2 {
                local year_temp_r = "`year_temp_r' c.indp`pg'#c.indf1#c.FD_yeartemp`k'"
            }
        }
    }

    local year_income_spline_r ""
    if ("`submodel'" == "lininter" | "`submodel'" == "quadinter") {
        forval pg = 1/2 {
            forval lg = 1/2 {
                forval k  = 1/2 {
                    local year_income_spline_r = "`year_income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`lg'temp`k'"
                }
            }
        }
    }

    if ("`submodel'" == "quadinter") {
        forval pg = 1/2 {
            forval k  = 1/2 {
                local year_temp_r = "`year_temp_r' c.indp`pg'#c.indf1#c.FD_year2temp`k'"
            }
        }

        forval pg = 1/2 {
            forval lg = 1/2 {
                forval k  = 1/2 {
                    local year_income_spline_r = "`year_income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`lg'temp`k'"
                }
            }
        }
    }

    * First-stage regression
    di "Running first-stage regression..."
    reghdfe FD_load_pc ///
        `temp_r' `precip_r' `climate_r' ///
        `lgdppc_MA15_r' `income_spline_r' ///
        `year_temp_r' `year_income_spline_r' ///
        DumInc*, absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        cluster(region_i) residuals(resid)

    * FGLS weights
    di "Calculating FGLS weights..."
    drop if resid == .
    bysort region_i: egen sd_resid = sd(resid)
    gen omega = sd_resid^2
    drop sd_resid
    gen weight = 1/omega
    drop resid

    * FGLS regression
    di "Running FGLS regression..."
    reghdfe FD_load_pc ///
        `temp_r' `precip_r' `climate_r' ///
        `lgdppc_MA15_r' `income_spline_r' ///
        `year_temp_r' `year_income_spline_r' ///
        DumInc* [pw = weight], ///
        absorb(i.flow_i#i.product_i#i.year#i.subregionid) ///
        cluster(region_i)

    * CRITICAL: Force reghdfe to store VCOV using estadd
    * This populates e(V) properly
    qui estadd vce

    * Check if VCOV is now populated
    matrix V_check = e(V)
    local v_test = V_check[1,1]
    di "VCOV check after estadd: V[1,1] = `v_test'"

    if `v_test' == 0 {
        di as error "ERROR: VCOV still zero after estadd!"
        di as error "Will construct manually from standard errors..."
        
        * Manually construct diagonal VCOV from _se[]
        * (This loses correlations but at least gives variances)
        matrix b = e(b)
        local k = colsof(b)
        local names : colnames b
        
        matrix V = J(`k', `k', 0)
        matrix colnames V = `names'
        matrix rownames V = `names'
        
        forval i = 1/`k' {
            local nm : word `i' of `names'
            cap local se = _se[`nm']
            if _rc == 0 {
                matrix V[`i',`i'] = `se'^2
            }
        }
        
        * Replace e(V)
        ereturn repost V = V
        
        di "Constructed diagonal VCOV from standard errors"
        di "NOTE: This omits correlations between parameters"
    }

    * Save estimates
    estimates store FDintF_`product'
    estimates save "$OUTDIR/FD_FGLS_inter_`model_name'_`product'", replace

    * Export coefficients
    di "Exporting coefficients..."
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

    * Export VCOV
    di "Exporting VCOV..."
    matrix V = e(V)
    local k = rowsof(V)
    local vnames : rownames V

    preserve
        clear
        set obs `=`k'*`k''

        gen strL parm_i = ""
        gen strL parm_j = ""
        gen double v = .

        local r = 0
        forvalues ii = 1/`k' {
            local name_i : word `ii' of `vnames'
            forvalues jj = 1/`k' {
                local ++r
                local name_j : word `jj' of `vnames'
                replace parm_i = "`name_i'" in `r'
                replace parm_j = "`name_j'" in `r'
                replace v = V[`ii',`jj'] in `r'
            }
        }

        * Verify
        qui count if v != 0
        local nz = r(N)
        di "VCOV: `nz' non-zero entries out of `=`k'*`k''"

        order parm_i parm_j v
        export delimited using ///
            "$OUTDIR/FD_FGLS_inter_`model_name'_`product'_vcov_long.csv", replace
    restore

    di "Saved: $OUTDIR/FD_FGLS_inter_`model_name'_`product'_coeff.csv"
    di "Saved: $OUTDIR/FD_FGLS_inter_`model_name'_`product'_vcov_long.csv"

}

di ""
di "=========================================="
di "ALL PRODUCTS COMPLETED"
di "=========================================="

************************************************************
* End 3_FD_interacted_all.do
************************************************************
