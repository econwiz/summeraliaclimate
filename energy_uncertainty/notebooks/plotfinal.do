/*
Purpose: 3x3 arrays energy-temperature response
FIXED to handle omitted co.* coefficients
*/

set scheme s1color
clear all
set more off

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global EDR "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/DATA/regression"

global DATA    "$PROJECT/data"
global OUTPUT  "$PROJECT"
global STERD   "$DATA/regression/sters"

local model_main "TINV_clim"
local submodel   "quadinter"
local model_name "`model_main'_`submodel'"

local fuels "electricity other_energy"
local products "GMFD ERA5 JRA_3Q MERRA2"

local year 2099

local col_electricity "dknavy"
local col_other_energy "dkorange"

foreach product of local products {
    foreach var of local fuels {

        di ""
        di "================================================"
        di "Plotting product=`product' fuel=`var' model=`model_name'"
        di "================================================"

        * Create temperature grid
        clear
        set obs 41
        
        gen double temp1 = _n - 6
        gen double temp2 = temp1^2
        gen byte above20 = (temp1 >= 20)
        gen byte below20 = (temp1 < 20)

        * Get income knot
        preserve
            use "$EDR/break_data_TINV_clim.dta", clear
            qui summ maxInc_largegpid_`var' if largegpid_`var' == 1
            local ibar = `r(max)'
        restore

        di "  ibar = `ibar'"

        local col_main "`col_`var''"

        * Create 3x3 plots
        local graphicM ""

        forval lg=3(-1)1 {
            forval tr=3(-1)1 {

                local cellid = `lg' + `tr'*100

                * Get tercile covariates
                preserve
                    use "$EDR/break_data_TINV_clim.dta", clear
                    duplicates drop tpid tgpid, force
                    sort tpid tgpid
                    local tr_index = `tr' * 3
                    local subCDD = avgCDD_tpid[`tr_index']
                    local subHDD = avgHDD_tpid[`tr_index']
                    local subInc = avgInc_tgpid[`lg']
                restore

                di "  Cell (lg=`lg', tr=`tr'): Inc=`subInc', CDD=`subCDD', HDD=`subHDD'"

                * pg assignment - MATCH YOUR DATA
                if "`var'"=="electricity" local pg=1
                else local pg=2

                * Income group
                local deltacut = `subInc' - `ibar'
                if (`subInc' > `ibar') local ig=2
                else local ig=1

                * Load estimates FIRST (so we can check for omitted terms)
                estimates use "$STERD/FD_FGLS_inter_`model_name'_`product'.ster"

                * Build predictnl expression - AVOID omitted co.* terms
                local line ""
                local add ""
                
                foreach k of num 1/2 {
                    * Base temperature effects
                    local line = "`line' `add' _b[c.indp`pg'#c.indf1#c.FD_temp`k'] * (temp`k' - 20^`k')"
                    
                    * Climate heterogeneity
                    local line = "`line' + above20*_b[c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k']*`subCDD' * (temp`k' - 20^`k')"
                    local line = "`line' + below20*_b[c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k']*`subHDD' * (20^`k' - temp`k')"
                    
                    * Income spline
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`ig'temp`k']*`deltacut'*(temp`k' - 20^`k')"

                    * Year × temp (check if omitted first)
                    cap _coef_table
                    matrix b = e(b)
                    local bnames : colnames b
                    
                    * Check if year terms exist and are NOT omitted
                    local has_year`k' = 0
                    foreach nm of local bnames {
                        if "`nm'" == "c.indp`pg'#c.indf1#c.FD_yeartemp`k'" {
                            local has_year`k' = 1
                        }
                    }
                    
                    if `has_year`k'' {
                        local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_yeartemp`k'] * (temp`k' - 20^`k')*`year'"
                    }
                    
                    * Year × income × temp
                    local has_dcyear`k' = 0
                    foreach nm of local bnames {
                        if "`nm'" == "c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`ig'temp`k'" {
                            local has_dcyear`k' = 1
                        }
                    }
                    
                    if `has_dcyear`k'' {
                        local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`ig'temp`k']*`deltacut'*`year'*(temp`k' - 20^`k')"
                    }

                    * Year² × temp
                    local has_year2`k' = 0
                    foreach nm of local bnames {
                        if "`nm'" == "c.indp`pg'#c.indf1#c.FD_year2temp`k'" {
                            local has_year2`k' = 1
                        }
                    }
                    
                    if `has_year2`k'' {
                        local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_year2temp`k'] * (temp`k' - 20^`k')*`year'*`year'"
                    }
                    
                    * Year² × income × temp
                    local has_dc2year`k' = 0
                    foreach nm of local bnames {
                        if "`nm'" == "c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`ig'temp`k'" {
                            local has_dc2year`k' = 1
                        }
                    }
                    
                    if `has_dc2year`k'' {
                        local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`ig'temp`k']*`deltacut'*`year'*`year'*(temp`k' - 20^`k')"
                    }

                    local add " + "
                }

                * Generate predictions
                cap predictnl yhat`cellid' = `line', se(se`cellid') ci(lo`cellid' hi`cellid')
                
                if _rc {
                    di as error "  ERROR: predictnl failed for cell `cellid'"
                    di as error "  Expression: `line'"
                    continue, break
                }

                * Create subplot
                tw ///
                    rarea hi`cellid' lo`cellid' temp1, col(`col_main'%30) || ///
                    line yhat`cellid' temp1, lc(`col_main') lw(medthick) || ///
                    , yline(0, lwidth(vthin) lc(black) lp(dash)) ///
                    xline(20, lwidth(vthin) lc(gray) lp(dot)) ///
                    xlabel(-5(10)35, labsize(vsmall)) ///
                    ylabel(, labsize(vsmall) nogrid format(%5.2f)) ///
                    legend(off) ///
                    xtitle("") ytitle("") ///
                    plotregion(color(white)) graphregion(color(white)) ///
                    nodraw name(Maddgraph`cellid', replace)

                local graphicM "`graphicM' Maddgraph`cellid'"
            }
        }

        * Combine and save
        graph combine `graphicM', imargin(zero) ycommon rows(3) ///
            title("`product': `var' energy-temperature response", size(medsmall)) ///
            subtitle("Model: `model_name' | Year: `year'", size(vsmall)) ///
            note("Rows: High→Mid→Low income | Columns: Hot→Mid→Cold climate", size(vsmall)) ///
            plotregion(color(white)) graphregion(color(white)) ///
            name(comb_`product'_`var', replace)

        cap mkdir "$OUTPUT/figures"
        local figname "$OUTPUT/figures/fig_3x3_`var'_`model_name'_`product'.pdf"
        graph export "`figname'", replace
        
        di "  SUCCESS: Saved `figname'"
        
        graph drop _all
    }
}

di ""
di "================================================"
di "ALL FIGURES COMPLETED"
di "================================================"
