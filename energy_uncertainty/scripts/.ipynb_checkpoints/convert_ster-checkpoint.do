*******************************************************
* export_from_ster.do
* Usage:
*   do export_from_ster.do "path/to/model.ster" "path/to/coeff.csv" "path/to/vcov.csv"
*******************************************************

version 15

* --- read arguments ---
args ster_path coeff_out vcov_out

di ">>> Loading ster file: `ster_path'"
estimates use "`ster_path'", clear

* ----------------------------------------------------
* 1) Export coefficient vector e(b)
* ----------------------------------------------------
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
    di ">>> Writing coeffs to: `coeff_out'"
    export delimited using "`coeff_out'", replace
restore

* ----------------------------------------------------
* 2) Export variance–covariance matrix e(V)
* ----------------------------------------------------
matrix V = e(V)
local k = rowsof(V)
local names : rownames V

preserve
    clear
    set obs `k'
    gen strL rowname = ""
    local j = 1
    foreach nm of local names {
        replace rowname = "`nm'" in `j'
        local ++j
    }

    local j = 1
    foreach nm of local names {
        gen double v_`nm' = .
        forvalues i = 1/`k' {
            replace v_`nm' = V[`i',`j'] in `i'
        }
        local ++j
    }

    order rowname
    di ">>> Writing vcov to: `vcov_out'"
    export delimited using "`vcov_out'", replace
restore

di ">>> Done exporting from `ster_path'"

