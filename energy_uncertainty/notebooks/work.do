/******************************************************************
plot_3x3_CIL_exact_fix.do

- Exact CIL plotting logic:
  * kink-at-20 basis
  * use predictnl with _b[...]
  * BEFORE predictnl: repost e(b), e(V) so every needed _b[c....] exists
    (missing or omitted -> set to 0)
******************************************************************/

version 16.0
clear all
set more off
set scheme s1color

* -----------------------------
* PATHS
* -----------------------------
global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global EDR     "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/DATA/regression"
global DATA    "$PROJECT/data"
global STERD   "$DATA/regression/sters"
global FIGDIR  "$PROJECT/figures"
cap mkdir "$FIGDIR"

* If something upstream left us in a preserved state, unwind it.
cap restore


* -----------------------------
* CONTROLS
* -----------------------------
local model_main "TINV_clim"
local submodel   "quadinter"
local model_name "`model_main'_`submodel'"

local fuels    "electricity other_energy"
local products "GMFD ERA5 JRA_3Q MERRA2"
local year 2099

local col_electricity  "dknavy"
local col_other_energy "dkorange"

cap log close _all
log using "$FIGDIR/plot_3x3_CIL_exact_fix_`model_name'.log", text replace

di "=== START `c(current_date)' `c(current_time)' ==="
di "model_name=`model_name'"

* ===============================================================
* Precompute tercile means ONCE (NO preserve/restore)
* ===============================================================
use "$EDR/break_data_TINV_clim.dta", clear

preserve
    keep tpid avgCDD_tpid avgHDD_tpid
    duplicates drop tpid, force
    forval tr=1/3 {
        qui su avgCDD_tpid if tpid==`tr'
        local CDD`tr' = r(mean)
        qui su avgHDD_tpid if tpid==`tr'
        local HDD`tr' = r(mean)
    }
restore

preserve
    keep tgpid avgInc_tgpid
    duplicates drop tgpid, force
    forval lg=1/3 {
        qui su avgInc_tgpid if tgpid==`lg'
        local INC`lg' = r(mean)
    }
restore

* after this, you can safely "clear" before making the T grid
clear


* ===============================================================
* MAIN LOOP
* ===============================================================
foreach product of local products {
    foreach var of local fuels {

        di "-------------------------------------------------------------"
        di "PLOT: product=`product' fuel=`var' model=`model_name' year=`year'"
        di "-------------------------------------------------------------"

        * pg for fuel
        local pg = .
        if ("`var'"=="electricity") local pg = 1
        else if ("`var'"=="other_energy") local pg = 2
        else continue

        * income knot (ibar)
        preserve
            use "$EDR/break_data_TINV_clim.dta", clear
            qui summ maxInc_largegpid_`var' if largegpid_`var' == 1
            local ibar = r(max)
        restore
        di "ibar=`ibar'"

        * load ster
        local sterfile "$STERD/FD_FGLS_inter_`model_name'_`product'.ster"
        di "ster=`sterfile'"
        cap estimates use "`sterfile'"
        if (_rc) {
            di as error "FAILED to load ster rc=`=_rc'"
            continue
        }

        * coefficient names
        matrix b0 = e(b)
        matrix V0 = e(V)
        local bnames0 : colnames b0

        * detect suffix for temp/precip/climate/year terms
        local SUF ""
        local candidates "_`product' _GMFD _ERA5 _JRA_3Q _MERRA2"
        local found 0
        foreach s of local candidates {
            foreach nm of local bnames0 {
                if ("`nm'"=="c.indp`pg'#c.indf1#c.FD_temp1`s'") {
                    local SUF "`s'"
                    local found 1
                }
            }
            if (`found') continue, break
        }
        di "SUF=`SUF'"

        * ===========================================================
        * CIL EXACT STEP: ensure all plotting coefficients EXIST as c....
        * - If only co.... exists -> create c.... = 0
        * - If neither exists -> create c.... = 0
        * ===========================================================
        local need ""

        * For plotting we only need terms for this pg, indf1, k=1/2, ig=1/2.
        foreach k of num 1/2 {
            local need "`need' c.indp`pg'#c.indf1#c.FD_temp`k'`SUF'"
            local need "`need' c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k'`SUF'"
            local need "`need' c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k'`SUF'"

            * time interactions (suffix may apply in THEIR old code)
            local need "`need' c.indp`pg'#c.indf1#c.FD_yeartemp`k'`SUF'"
            local need "`need' c.indp`pg'#c.indf1#c.FD_year2temp`k'`SUF'"

            * income spline interactions (UNSUFFIXED in their code)
            foreach ig of num 1/2 {
                local need "`need' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`ig'temp`k'"
                local need "`need' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`ig'temp`k'"
                local need "`need' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`ig'temp`k'"
            }
        }

        * Expand b and V with missing coefs as zeros (CIL “complete_for_plot” idea)
        matrix b = b0
        matrix V = V0
        local bnames : colnames b

        foreach coef of local need {

            * check if coef exists already
            local has_c 0
            foreach nm of local bnames {
                if ("`nm'"=="`coef'") local has_c 1
            }

            if (!`has_c') {
                * even if there is a co. version, we still create the c. version for plotting
                local oldK = colsof(b)
                matrix b = b , J(1,1,0)
                matrix V = ( V , J(`oldK',1,0) ) \ ( J(1,`oldK',0) , J(1,1,0) )

                * assign col/row names
                local newnames "`bnames' `coef'"
                matrix colnames b = `newnames'
                matrix colnames V = `newnames'
                matrix rownames V = `newnames'
                local bnames "`newnames'"

                di "REPOST: added missing coef as 0 -> `coef'"
            }
        }

        * repost into e()
        ereturn repost b=b V=V

        * ===========================================================
        * Build temperature grid and CIL kink bases
        * ===========================================================
        clear
        set obs 41
        gen double T  = _n - 6           // -5..35
        gen double T1 = T
        gen double T2 = T^2
        gen byte ABOVE20 = (T >= 20)
        gen byte BELOW20 = (T < 20)

        gen double dT1_ABV = (T1 - 20^1)
        gen double dT2_ABV = (T2 - 20^2)
        gen double dT1_BLW = (20^1 - T1)
        gen double dT2_BLW = (20^2 - T2)

        local col_main "`col_`var''"
        local graphicM ""

        * ===========================================================
        * 3x3 panels (rows: income 3->1, cols: climate 3->1)
        * ===========================================================
        forval lg=3(-1)1 {
            forval tr=3(-1)1 {

                local cellid = `lg' + `tr'*100

                local subCDD = `CDD`tr''
                local subHDD = `HDD`tr''
                local subInc = `INC`lg''

                local deltacut = `subInc' - `ibar'
                local ig = 1
                if (`subInc' > `ibar') local ig = 2

                di "CELL `cellid': lg=`lg' tr=`tr' subCDD=`subCDD' subHDD=`subHDD' subInc=`subInc' ig=`ig' deltacut=`deltacut'"

                * CIL exact expression (safe now because we reposted missing coefs)
                local line "0"
                foreach k of num 1/2 {

                    * base temp (same coef, different side basis)
                    local line "`line' + ABOVE20*_b[c.indp`pg'#c.indf1#c.FD_temp`k'`SUF']*dT`k'_ABV"
                    local line "`line' + BELOW20*_b[c.indp`pg'#c.indf1#c.FD_temp`k'`SUF']*dT`k'_BLW"

                    * climate heterogeneity
                    local line "`line' + ABOVE20*_b[c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k'`SUF']*`subCDD'*dT`k'_ABV"
                    local line "`line' + BELOW20*_b[c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k'`SUF']*`subHDD'*dT`k'_BLW"

                    * income spline
                    local line "`line' + ABOVE20*_b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`ig'temp`k']*`deltacut'*dT`k'_ABV"
                    local line "`line' + BELOW20*_b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`ig'temp`k']*`deltacut'*dT`k'_BLW"

                    * year x temp (will be 0 if omitted/missing)
                    local line "`line' + ABOVE20*_b[c.indp`pg'#c.indf1#c.FD_yeartemp`k'`SUF']*`year'*dT`k'_ABV"
                    local line "`line' + BELOW20*_b[c.indp`pg'#c.indf1#c.FD_yeartemp`k'`SUF']*`year'*dT`k'_BLW"

                    * year x income spline
                    local line "`line' + ABOVE20*_b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`ig'temp`k']*`deltacut'*`year'*dT`k'_ABV"
                    local line "`line' + BELOW20*_b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`ig'temp`k']*`deltacut'*`year'*dT`k'_BLW"

                    * year^2 x temp
                    local line "`line' + ABOVE20*_b[c.indp`pg'#c.indf1#c.FD_year2temp`k'`SUF']*(`year'^2)*dT`k'_ABV"
                    local line "`line' + BELOW20*_b[c.indp`pg'#c.indf1#c.FD_year2temp`k'`SUF']*(`year'^2)*dT`k'_BLW"

                    * year^2 x income spline
                    local line "`line' + ABOVE20*_b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`ig'temp`k']*`deltacut'*(`year'^2)*dT`k'_ABV"
                    local line "`line' + BELOW20*_b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`ig'temp`k']*`deltacut'*(`year'^2)*dT`k'_BLW"
                }

                cap predictnl yhat`cellid' = `line', se(se`cellid') ci(lo`cellid' hi`cellid')
                if (_rc) {
                    di as error "predictnl FAILED cell `cellid' rc=`=_rc'"
                    di as error "`line'"
                    continue
                }

                tw ///
                    rarea hi`cellid' lo`cellid' T, col(`col_main'%30) || ///
                    line yhat`cellid' T, lc(`col_main') lw(medthick) || ///
                    , yline(0, lwidth(vthin) lc(black) lp(dash)) ///
                      xline(20, lwidth(vthin) lc(gray) lp(dot)) ///
                      xlabel(-5(10)35, labsize(vsmall)) ///
                      ylabel(, labsize(vsmall) nogrid format(%5.2f)) ///
                      legend(off) xtitle("") ytitle("") ///
                      plotregion(color(white)) graphregion(color(white)) ///
                      nodraw name(M`cellid', replace)

                local graphicM "`graphicM' M`cellid'"
            }
        }

        graph combine `graphicM', imargin(zero) ycommon rows(3) ///
            title("`product': `var'", size(medsmall)) ///
            subtitle("`model_name' | Year=`year' | SUF=`SUF'", size(vsmall)) ///
            note("Rows: High→Mid→Low income | Cols: Hot→Mid→Cold climate", size(vsmall)) ///
            plotregion(color(white)) graphregion(color(white)) ///
            name(comb_`product'_`var', replace)

        local figname "$FIGDIR/fig_3x3_`var'_`model_name'_`product'_CILexact_fix.pdf"
        graph export "`figname'", replace
        di "SAVED: `figname'"

        graph drop _all
    }
}

di "=== DONE ==="
log close
