/*
Purpose: Check what coefficients are in YOUR new .ster files
*/

clear all

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global STERD   "$PROJECT/data/regression/sters"

local model_name "TINV_clim_quadinter"
local product "GMFD"

* Load YOUR new estimates
estimates use "$STERD/FD_FGLS_inter_`model_name'_`product'"

di ""
di "=========================================="
di "COEFFICIENTS IN YOUR NEW .STER FILE"
di "File: $STERD/FD_FGLS_inter_`model_name'_`product'.ster"
di "=========================================="
di ""

* Get coefficient names
matrix b = e(b)
local names : colnames b
local k = colsof(b)

di "Total coefficients: `k'"
di ""

* Show first 30 coefficients
di "FIRST 30 COEFFICIENTS:"
di "-------------------------------------------"
local count = 0
foreach nm of local names {
    local ++count
    di "`count'. `nm'"
    if `count' >= 30 {
        continue, break
    }
}

di ""
di "TEMPERATURE COEFFICIENTS (looking for FD_temp):"
di "-------------------------------------------"
foreach nm of local names {
    if strpos("`nm'", "FD_temp") > 0 & strpos("`nm'", "year") == 0 & strpos("`nm'", "dc1") == 0 {
        di "`nm'"
    }
}

di ""
di "CLIMATE HETEROGENEITY (looking for cdd20_TINV, hdd20_TINV):"
di "-------------------------------------------"
foreach nm of local names {
    if (strpos("`nm'", "cdd20_TINV") > 0 | strpos("`nm'", "hdd20_TINV") > 0) {
        di "`nm'"
    }
}

di ""
di "INCOME SPLINE (looking for dc1_lgdppc_MA15I):"
di "-------------------------------------------"
local count2 = 0
foreach nm of local names {
    if strpos("`nm'", "dc1_lgdppc_MA15I") > 0 & strpos("`nm'", "year") == 0 {
        local ++count2
        di "`nm'"
        if `count2' >= 10 {
            di "... (showing first 10 only)"
            continue, break
        }
    }
}

di ""
di "CHECKING EXPECTED COEFFICIENT NAMES:"
di "-------------------------------------------"

* Check electricity (pg=1) temperature coefficients
local test_names "c.indp1#c.indf1#c.FD_temp1 c.indp1#c.indf1#c.FD_temp2"

foreach tn of local test_names {
    local found = 0
    foreach nm of local names {
        if "`nm'" == "`tn'" {
            local found = 1
        }
    }
    
    if `found' {
        di "✓ FOUND: `tn'"
    }
    else {
        di "✗ MISSING: `tn'"
    }
}

* Check with indp2 (in case electricity uses pg=2)
local test_names2 "c.indp2#c.indf1#c.FD_temp1 c.indp2#c.indf1#c.FD_temp2"

di ""
di "Checking with indp2 instead:"
foreach tn of local test_names2 {
    local found = 0
    foreach nm of local names {
        if "`nm'" == "`tn'" {
            local found = 1
        }
    }
    
    if `found' {
        di "✓ FOUND: `tn'"
    }
    else {
        di "✗ MISSING: `tn'"
    }
}

di ""
di "=========================================="
di "DONE"
di "=========================================="
