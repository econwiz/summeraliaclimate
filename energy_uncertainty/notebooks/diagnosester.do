/*
Purpose: 3x3 arrays energy-temperature response (income x climate)
Loop over climate products using YOUR .ster outputs.
UNSUFFIXED coefficients (no _GMFD on _b[] names).
*/

set scheme s1color
clear all
set more off

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global EDR "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/DATA/regression" 

global DATA    "$PROJECT/data"
global OUTPUT  "$PROJECT"              // for figures output path if you want
global STERD   "$DATA/regression/sters"

* ---- USER CONTROLS ----
local model_main "TINV_clim"
local submodel   "quadinter"
local model_name "`model_main'_`submodel'"

local fuels "electricity other_energy"
local products "GMFD ERA5 JRA_3Q MERRA2"

local year 2099

* colors (optional)
local col_electricity "dknavy"
local col_other_energy "dkorange"

foreach product of local products {
    foreach var of local fuels {

        di "================================================"
        di "Plotting product=`product' fuel=`var' model=`model_name'"
        di "================================================"

        * -----------------------------
        * Step 1: temp grid
        * -----------------------------
        * use any regsort just to anchor environment; then replace with empty grid
        use "$DATA/regression/`product'_`model_main'_regsort.dta", clear

               * -----------------------------
        * Step 1: temp grid (CLEAN)
        * -----------------------------
                * -----------------------------
        * Step 1: temp grid (NUCLEAR SAFE)
        * -----------------------------
        capture restore
        clear all
        set obs 41   // -5..35 inclusive
        
        gen double temp1 = _n - 6
        gen double temp2 = temp1^2

        gen byte above20 = (temp1 >= 20)
        gen byte below20 = (temp1 < 20)

        gen double FD_load_pc = .

        * -----------------------------
        * Step 2: knot + plotting setup
        * -----------------------------
        preserve
            use "$EDR/break_data_TINV_clim.dta", clear
            summ maxInc_largegpid_`var' if largegpid_`var' == 1
            local ibar = `r(max)'
        restore

        local col_main "`col_`var''"

        * -----------------------------
        * Step 3: 3x3 plots
        * -----------------------------
        local graphicM ""

        forval lg=3(-1)1 {
            forval tr=3(-1)1 {

                local cellid = `lg' + `tr'*100

                * pull tercile covariates like their code
                preserve
                    use "$EDR/break_data_TINV_clim.dta", clear
                    duplicates drop tpid tgpid, force
                    sort tpid tgpid
                    local tr_index = `tr' * 3
                    local subCDD = avgCDD_tpid[`tr_index']
                    local subHDD = avgHDD_tpid[`tr_index']
                    local subInc = avgInc_tgpid[`lg']
                restore

                * pg assignment (matches their code)
                if "`var'"=="electricity" local pg=1
                else local pg=2

                * income group at knot
                local deltacut = `subInc' - `ibar'
                if (`subInc' > `ibar') local ig=2
                else local ig=1

                * build predictnl expression (UNSUFFIXED)
                local line ""
                local add ""
                foreach k of num 1/2 {

                    local line = "`line' `add' _b[c.indp`pg'#c.indf1#c.FD_temp`k'] * (temp`k' - 20^`k')"
                    local line = "`line' + above20*_b[c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k']*`subCDD' * (temp`k' - 20^`k')"
                    local line = "`line' + below20*_b[c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k']*`subHDD' * (20^`k' - temp`k')"
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`ig'temp`k']*`deltacut'*(temp`k' - 20^`k')"

                    * lininter/quadinter pieces
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_yeartemp`k'] * (temp`k' - 20^`k')*`year'"
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`ig'temp`k']*`deltacut'*`year'*(temp`k' - 20^`k')"

                    * quadinter only
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_year2temp`k'] * (temp`k' - 20^`k')*`year'*`year'"
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`ig'temp`k']*`deltacut'*`year'*`year'*(temp`k' - 20^`k')"

                    local add " + "
                }

                * load your product-specific ster
                estimates use "$STERD/FD_FGLS_inter_`model_name'_`product'.ster"

                * trace out curve
                predictnl yhat`cellid' = `line', se(se`cellid') ci(lo`cellid' hi`cellid')

                * subplot
                tw ///
                    rarea hi`cellid' lo`cellid' temp1, col(`col_main'%30) || ///
                    line yhat`cellid' temp1, lc(`col_main') || ///
                    yline(0, lwidth(vthin)) ///
                    xlabel(-5(10)35, labsize(vsmall)) ///
                    ylabel(, labsize(vsmall) nogrid) legend(off) ///
                    plotregion(color(white)) graphregion(color(white)) nodraw ///
                    name(Maddgraph`cellid', replace)

                local graphicM "`graphicM' Maddgraph`cellid'"
            }
        }

        * combine + export
        graph combine `graphicM', imargin(zero) ycomm rows(3) ///
            title("Split Degree Days Poly 2 Interaction Model `var' | `product'", size(small)) ///
            subtitle("`model_name' (`col_main')", size(vsmall)) ///
            plotregion(color(white)) graphregion(color(white)) ///
            name(comb, replace)

        cap mkdir "$OUTPUT/figures"
        graph export "$OUTPUT/figures/fig_3x3_`var'_interacted_`model_name'_`product'.pdf", replace
        graph drop _all
    }
}
