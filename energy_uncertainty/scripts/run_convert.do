
global PROJ "/user/ab5405/summeraliaclimate/code/energy_uncertainty"

* New GMFD FGLS global:
do export_from_ster.do ///
    "$PROJ/data/regression/sters/FD_FGLS_global_GMFD.ster" ///
    "$PROJ/data/regression/sters/FD_FGLS_global_GMFD_coeff.csv" ///
    "$PROJ/data/regression/sters/FD_FGLS_global_GMFD_vcov.csv"

* Original TINV global (in output_original):
do export_from_ster.do ///
    "$PROJ/output_original/sters/FD_FGLS_global_TINV_clim.ster" ///
    "$PROJ/output_original/sters/FD_FGLS_global_TINV_clim_coeff.csv" ///
    "$PROJ/output_original/sters/FD_FGLS_global_TINV_clim_vcov.csv"

