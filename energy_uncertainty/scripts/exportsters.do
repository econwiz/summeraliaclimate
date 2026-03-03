************************************************************
* export_original_lininter_coeff.do
* Loads the original CIL lininter .ster and exports the
* coefficient vector to CSV for use in Python plotting.
* Does NOT rerun any regression — read-only.
************************************************************

local orig_ster "/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/OUTPUT/sters/FD_inter_TINV_clim_lininter_reg"
local orig_dir  "/user/ab5405/summeraliaclimate/code/energy_uncertainty/output_original/sters"
cap mkdir "`orig_dir'"

estimates use "`orig_ster'"

* Export coefficient vector
matrix b = e(b)
local kb = colsof(b)
local bnames : colnames b

preserve
    clear
    set obs `kb'
    gen strL parm = ""
    gen double beta = .
    local j = 1
    foreach nm of local bnames {
        replace parm = "`nm'" in `j'
        replace beta = b[1,`j'] in `j'
        local ++j
    }
    order parm beta
    export delimited using ///
        "`orig_dir'/FD_FGLS_inter_TINV_clim_lininter_coeff.csv", replace
restore

di "Done. Coefficients written to `orig_dir'"