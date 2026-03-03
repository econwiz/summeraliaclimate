/*
DEBUG version - prints the expression
*/

set scheme s1color
clear all
set more off

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global EDR "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/DATA/regression"
global STERD   "$PROJECT/data/regression/sters"

local var "electricity"
local product "GMFD"
local model_name "TINV_clim_quadinter"
local year 2099

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

di "ibar = `ibar'"

* Just do ONE cell for testing
local lg = 3
local tr = 3
local cellid = 303

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

di "Cell values: Inc=`subInc', CDD=`subCDD', HDD=`subHDD'"

* pg assignment
if "`var'"=="electricity" local pg=1
else local pg=2

di "Using pg = `pg'"

* Income group
local deltacut = `subInc' - `ibar'
if (`subInc' > `ibar') local ig=2
else local ig=1

di "Income group: ig = `ig', deltacut = `deltacut'"

* Load estimates
estimates use "$STERD/FD_FGLS_inter_`model_name'_`product'"

* Get coefficient names
matrix b = e(b)
local bnames : colnames b

* Build expression with detailed checking
local line ""
local add ""

foreach k of num 1/2 {
    di ""
    di "Building terms for k=`k'..."
    
    * Base temperature
    local term "c.indp`pg'#c.indf1#c.FD_temp`k'"
    di "  Checking: `term'"
    cap scalar test = _b[`term']
    if _rc {
        di "    ERROR: Not found!"
    }
    else {
        di "    ✓ Found, value = " _b[`term']
        local line = "`line' `add' _b[`term'] * (temp`k' - 20^`k')"
    }
    
    * Climate - CDD
    local term "c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k'"
    di "  Checking: `term'"
    cap scalar test = _b[`term']
    if _rc {
        di "    ERROR: Not found!"
    }
    else {
        di "    ✓ Found, value = " _b[`term']
        local line = "`line' + above20*_b[`term']*`subCDD' * (temp`k' - 20^`k')"
    }
    
    * Climate - HDD
    local term "c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k'"
    di "  Checking: `term'"
    cap scalar test = _b[`term']
    if _rc {
        di "    ERROR: Not found!"
    }
    else {
        di "    ✓ Found, value = " _b[`term']
        local line = "`line' + below20*_b[`term']*`subHDD' * (20^`k' - temp`k')"
    }
    
    * Income spline
    local term "c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`ig'temp`k'"
    di "  Checking: `term'"
    cap scalar test = _b[`term']
    if _rc {
        di "    ERROR: Not found!"
    }
    else {
        di "    ✓ Found, value = " _b[`term']
        local line = "`line' + _b[`term']*`deltacut'*(temp`k' - 20^`k')"
    }
    
    local add " + "
}

di ""
di "=========================================="
di "FINAL EXPRESSION:"
di "`line'"
di "=========================================="

* Try to predict
di ""
di "Attempting predictnl..."
cap predictnl yhat`cellid' = `line', se(se`cellid') ci(lo`cellid' hi`cellid')

if _rc {
    di as error "ERROR: predictnl failed with code `=_rc'"
}
else {
    di "SUCCESS!"
    list temp1 yhat`cellid' in 1/5
}
