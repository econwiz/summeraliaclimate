/*
Working version - treat missing/omitted coefficients as zero
*/

set scheme s1color
clear all
set more off

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global EDR "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/DATA/regression"
global DATA    "$PROJECT/data"
global STERD   "$DATA/regression/sters"
global OUTPUT  "$PROJECT/figures"

cap mkdir "$OUTPUT"

local model_name "TINV_clim_quadinter"
local fuels "electricity other_energy"
local products "GMFD ERA5 JRA_3Q MERRA2"
local year 2099

local col_electricity "dknavy"
local col_other_energy "dkorange"

* Precompute tercile means
tempfile tpid_means tgpid_means
preserve
    use "$EDR/break_data_TINV_clim.dta", clear
    collapse (mean) avgCDD_tpid avgHDD_tpid, by(tpid)
    save `tpid_means'
restore

preserve
    use "$EDR/break_data_TINV_clim.dta", clear  
    collapse (mean) avgInc_tgpid, by(tgpid)
    save `tgpid_means'
restore

foreach product of local products {
    foreach var of local fuels {

        di ""
        di "================================================"
        di "`product' - `var'"
        di "================================================"

        * pg assignment
        local pg = cond("`var'"=="electricity", 1, 2)

        * Get income knot
        preserve
            use "$EDR/break_data_TINV_clim.dta", clear
            qui summ maxInc_largegpid_`var' if largegpid_`var' == 1
            local ibar = r(max)
        restore

        * Load estimates
        estimates use "$STERD/FD_FGLS_inter_`model_name'_`product'"
        
        * Get all coefficient names
        matrix b = e(b)
        local allnames : colnames b

        * Create temp grid
        clear
        set obs 41
        gen double temp1 = _n - 6
        gen double temp2 = temp1^2
        gen byte above20 = (temp1 >= 20)
        gen byte below20 = (temp1 < 20)

        local col_main "`col_`var''"
        local graphicM ""

        forval lg=3(-1)1 {
            forval tr=3(-1)1 {

                local cellid = `lg' + `tr'*100

                * Get cell values from precomputed means
                preserve
                    use `tpid_means', clear
                    qui su avgCDD_tpid if tpid==`tr'
                    local subCDD = r(mean)
                    qui su avgHDD_tpid if tpid==`tr'
                    local subHDD = r(mean)
                restore
                
                preserve
                    use `tgpid_means', clear
                    qui su avgInc_tgpid if tgpid==`lg'
                    local subInc = r(mean)
                restore

                local deltacut = `subInc' - `ibar'
                local ig = cond(`subInc' > `ibar', 2, 1)

                * Build expression - only use coefficients that EXIST
                local line "0"
                
                foreach k of num 1/2 {
                    * Base temp (always exists)
                    local line "`line' + _b[c.indp`pg'#c.indf1#c.FD_temp`k']*(temp`k'-20^`k')"
                    
                    * Climate
                    local line "`line' + above20*_b[c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k']*`subCDD'*(temp`k'-20^`k')"
                    local line "`line' + below20*_b[c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k']*`subHDD'*(20^`k'-temp`k')"
                    
                    * Income spline
                    local line "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`ig'temp`k']*`deltacut'*(temp`k'-20^`k')"
                    
                    * Year terms - check if they exist (not omitted)
                    local yt "c.indp`pg'#c.indf1#c.FD_yeartemp`k'"
                    local found 0
                    foreach nm of local allnames {
                        if "`nm'"=="`yt'" local found 1
                    }
                    if `found' {
                        local line "`line' + _b[`yt']*(temp`k'-20^`k')*`year'"
                    }
                    
                    * Year x income
                    local yinc "c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`ig'temp`k'"
                    local found 0
                    foreach nm of local allnames {
                        if "`nm'"=="`yinc'" local found 1
                    }
                    if `found' {
                        local line "`line' + _b[`yinc']*`deltacut'*`year'*(temp`k'-20^`k')"
                    }
                    
                    * Year²
                    local y2t "c.indp`pg'#c.indf1#c.FD_year2temp`k'"
                    local found 0
                    foreach nm of local allnames {
                        if "`nm'"=="`y2t'" local found 1
                    }
                    if `found' {
                        local line "`line' + _b[`y2t']*(temp`k'-20^`k')*`year'^2"
                    }
                    
                    * Year² x income
                    local y2inc "c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`ig'temp`k'"
                    local found 0
                    foreach nm of local allnames {
                        if "`nm'"=="`y2inc'" local found 1
                    }
                    if `found' {
                        local line "`line' + _b[`y2inc']*`deltacut'*`year'^2*(temp`k'-20^`k')"
                    }
                }

                * Predict
                cap predictnl yhat`cellid' = `line', se(se`cellid') ci(lo`cellid' hi`cellid')
                
                if _rc {
                    di as error "  Error in cell `cellid', rc=`=_rc'"
                    continue, break
                }

                * Plot
                tw ///
                    rarea hi`cellid' lo`cellid' temp1, col(`col_main'%30) || ///
                    line yhat`cellid' temp1, lc(`col_main') lw(medthick) || ///
                    , yline(0, lwidth(vthin) lc(black) lp(dash)) ///
                    xline(20, lwidth(vthin) lc(gray) lp(dot)) ///
                    xlabel(-5(10)35, labsize(vsmall)) ///
                    ylabel(, labsize(vsmall) nogrid) ///
                    legend(off) xtitle("") ytitle("") ///
                    plotregion(color(white)) graphregion(color(white)) ///
                    nodraw name(M`cellid', replace)

                local graphicM "`graphicM' M`cellid'"
            }
        }

        * Combine
        graph combine `graphicM', imargin(zero) ycommon rows(3) ///
            title("`product': `var'", size(medsmall)) ///
            subtitle("`model_name' | Year=`year'", size(vsmall)) ///
            note("Rows: High→Mid→Low income | Cols: Hot→Mid→Cold climate", size(vsmall)) ///
            plotregion(color(white)) graphregion(color(white))

        graph export "$OUTPUT/fig_3x3_`var'_`model_name'_`product'.pdf", replace
        di "  ✓ Saved: fig_3x3_`var'_`model_name'_`product'.pdf"

        graph drop _all
    }
}

di ""
di "Done!"
