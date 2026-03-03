/******************************************************************
Plot 3x3 energy response (income x climate) using YOUR outputs
- Uses the original plotting logic (predictnl) but repoints paths
- Saves PDF into: energy_uncertainty/figures
******************************************************************/

clear all
set more off
set scheme s1color

* -----------------------------
* USER PATHS: EDIT ONLY THIS
* -----------------------------
global PROJECT "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21"
global DATA    "$PROJECT/DATA"
global OUTPUT  "$PROJECT"                 // so $OUTPUT/figures exists
global FIGDIR  "$PROJECT/figures"

cap mkdir "$FIGDIR"

* -----------------------------
* USER CONTROLS: EDIT THESE
* -----------------------------
* model_main should match the suffix in your DTA names:
*   $DATA/regression/GMFD_`model_main'_regsort.dta
*   $DATA/regression/break_data_`model_main'.dta
local model_main "TINV_clim"    // <<< change if needed
local var        "electricity"            // "electricity" or "other_energy"
local submodel_ov ""                      // leave blank for now
local year = 2099

* plotting colors
local col_electricity "dknavy"
local col_electricity_ov "red"
local col_other_energy "dkorange"
local col_other_energy_ov "black"

local col_main "`col_`var''"
local col_ov "`col_`var'_ov'"

* -----------------------------
* Step 1: Load Data and Clean for Plotting
* -----------------------------
use "$DATA/regression/GMFD_`model_main'_regsort.dta", clear

local obs = 35 + abs(-5) + 1

drop if _n > 0
set obs `obs'

replace temp1_GMFD = _n - 6

foreach k of num 1/2 {
    rename temp`k'_GMFD temp`k'
    replace temp`k' = temp1 ^ `k'
}

gen above20 = (temp1 >= 20)
gen below20 = (temp1 < 20)

* -----------------------------
* Step 2: Find knot (ibar) using YOUR break_data
* -----------------------------
preserve
use "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/DATA/regression/break_data_TINV_clim.dta", clear
summ maxInc_largegpid_`var' if largegpid_`var' == 1
local ibar_main = `r(max)'
restore

local ibar_ov = `ibar_main'

local colorGuide " Model Spec: `model_main' (`col_main') "
local plot_title "`model_main'"

local type_list " _main "

* -----------------------------
* Step 3: Plot 3x3
* -----------------------------
local graphicM=""

forval lg=3(-1)1 {       // income tercile (rich -> poor)
    forval tr=3(-1)1 {   // climate tercile (cold -> hot) per their code order

        local cellid = `lg' + `tr'*100

        * pull cell covariates
        preserve
        use "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/DATA/regression/break_data_TINV_clim.dta", clear
        duplicates drop tpid tgpid, force
        sort tpid tgpid

        local tr_index = `tr' * 3
        local subCDD = avgCDD_tpid[`tr_index']
        local subHDD = avgHDD_tpid[`tr_index']
        local subInc = avgInc_tgpid[`lg']
        restore

        * choose fuel index
        if "`var'"=="electricity" local pg=1
        else if "`var'"=="other_energy" local pg=2

        local SE ""

        foreach type in `type_list' {

            local plot_model "`model_main'"

            * income spline
            local deltacut_subInc = `subInc' - `ibar`type''

            * income group side
            if `subInc' > `ibar`type'' local ig = 2
            else if `subInc' <= `ibar`type'' local ig = 1

            local line ""
            local add ""

            foreach k of num 1/2 {

                local line = " `line' `add' _b[c.indp`pg'#c.indf1#c.FD_temp`k'_GMFD] * (temp`k' - 20^`k')"
                local line = "`line' + above20*_b[c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k'_GMFD]*`subCDD' * (temp`k' - 20^`k')"
                local line = "`line' + below20*_b[c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k'_GMFD]*`subHDD' * (20^`k' - temp`k')"
                local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`ig'temp`k']*`deltacut_subInc'*(temp`k' - 20^`k')"

                if ((strpos("`plot_model'", "lininter") > 0) | (strpos("`plot_model'", "quadinter") > 0)) {
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_yeartemp`k'] * (temp`k' - 20^`k')*`year'"
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`ig'temp`k']*`deltacut_subInc'*`year'*(temp`k' - 20^`k')"
                }

                if (strpos("`plot_model'", "quadinter") > 0) {
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_year2temp`k'] * (temp`k' - 20^`k')*`year'*`year'"
                    local line = "`line' + _b[c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`ig'temp`k']*`deltacut_subInc'*`year'*`year'*(temp`k' - 20^`k')"
                }

                local add " + "
            }

            * --- load YOUR model estimates ---
            * <<< EDIT THIS IF NEEDED >>>
            * This assumes you saved a ster named:
            *   $DATA/regression/sters/FD_FGLS_inter_`plot_model'_quadinter_GMFD.ster
            estimates use "$DATA/regression/sters/FD_FGLS_inter_`plot_model'_quadinter_GMFD.ster"

            predictnl yhat`cellid'`type' = `line', se(se`cellid'`type') ci(lower`cellid'`type' upper`cellid'`type')

            local SE = "`SE' rarea upper`cellid'`type' lower`cellid'`type' temp1, col(`col`type''%30) || line yhat`cellid'`type' temp1, lc(`col`type'') ||"
        }

        tw `SE' , yline(0, lwidth(vthin)) xlabel(-5(10)35, labsize(vsmall)) ///
            ylabel(, labsize(vsmall) nogrid) legend(off) ///
            subtitle("", size(vsmall) color(dkgreen)) ///
            ytitle("", size(vsmall)) xtitle("", size(small)) ///
            plotregion(color(white)) graphregion(color(white)) nodraw ///
            name(Maddgraph`cellid', replace)

        local graphicM="`graphicM' Maddgraph`cellid'"
    }
}

graph combine `graphicM', imargin(zero) ycomm rows(3) ///
    title("Split Degree Days Poly 2 Interaction Model `var'", size(small)) ///
    subtitle("`colorGuide'", size(vsmall)) ///
    plotregion(color(white)) graphregion(color(white)) name(comb, replace)

graph export "$FIGDIR/fig_3x3_`var'_`model_main'.pdf", replace
graph drop _all
di as result "Saved: $FIGDIR/fig_3x3_`var'_`model_main'.pdf"
