************************************************************
* 3_FD_interacted_all.do
* Replicate CIL stacked "interacted" FD spec for all products,
* dropping TINV HDD/CDD terms.
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

local products "GMFD ERA5 JRA_3Q"

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

    * Again, short name:
    estimates store FDintF_`product'
    estimates save "$OUTDIR/FD_FGLS_inter_`model_name'_`product'", replace
}

************************************************************
* End 3_FD_interacted_all.do
************************************************************
