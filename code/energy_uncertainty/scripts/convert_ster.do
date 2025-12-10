clear all
set more off

* paths
global OLD "/user/ab5405/summeraliaclimate/code/energy_uncertainty/output_original/sters"

* helper program: load ster, dump e(b) to CSV
capture program drop est_to_csv
program define est_to_csv
    syntax using/, SAVE(string)

    estimates clear        // <-- this clears estimation results safely
    estimates use "`using'"   // <-- no "clear" allowed here

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
            replace parm = "`nm'"   in `j'
            replace beta = b[1,`j'] in `j'
            local ++j
        }

        order parm beta
        export delimited using "`save'", replace
    restore
end


* ---- export each original ster ----
est_to_csv using "$OLD/FD_FGLS_global_TINV_clim.ster", ///
    save("$OLD/FD_FGLS_global_TINV_clim_coeff.csv")

est_to_csv using "$OLD/FD_global_TINV_clim.ster", ///
    save("$OLD/FD_global_TINV_clim_coeff.csv")

est_to_csv using "$OLD/FD_FGLS_inter_TINV_clim_quadinter.ster", ///
    save("$OLD/FD_FGLS_inter_TINV_clim_quadinter_coeff.csv")

est_to_csv using "$OLD/FD_inter_TINV_clim_quadinter.ster", ///
    save("$OLD/FD_inter_TINV_clim_quadinter_coeff.csv")
