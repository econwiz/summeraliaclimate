************************************************************
* 0_build_FD_dataset_all.do
* Takes merged energy + climate panels (one per product) and constructs
* all first-differenced regressors needed for the Rode et al. energy
* consumption regressions. For each product:
*   - Generates decade dummies and income spline
*   - First-differences load_pc, income group indicators, precipitation,
*     temperature polynomials, polyAbove/polyBelow, and all interactions
*   - Adds TINV cross terms (pixel-level if available, panel-level otherwise)
*
* Inputs:  $MERGED/<product>_energy_panel.dta
* Outputs: $REGD/<product>_TINV_clim_regsort.dta
*
************************************************************
clear all
set more off

* Root + data dirs
global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty"
global DATA    "$PROJECT/data"
global MERGED  "$DATA/merged_panels"
global REGD    "$DATA/regression"

cap mkdir "$REGD"

* Climate products to process
local products "ERA5 JRA_3Q GMFD MERRA2"

foreach product of local products {

    di "==============================================="
    di "  Building FD dataset for `product'"
    di "==============================================="

    * Suffix for climate variables in this file
    local suf "`product'"

    * ---- load merged energy + climate panel for this product (.dta) ----
    * (you must have created these with your Python merge script)
    use "$MERGED/`product'_energy_panel.dta", clear

    * ---- sanity check: key variables must exist ----
    capture confirm variable lgdppc_MA15
    if _rc {
        di as error "lgdppc_MA15 missing in `product'_energy_panel.dta"
        describe
        exit 111
    }

    capture confirm variable region_i
    if _rc {
        di as error "region_i missing in `product'_energy_panel.dta"
        describe
        exit 111
    }

    capture confirm variable year
    if _rc {
        di as error "year missing in `product'_energy_panel.dta"
        describe
        exit 111
    }

    capture confirm variable load_pc
    if _rc {
        di as error "load_pc missing in `product'_energy_panel.dta"
        describe
        exit 111
    }

    * Drop any previous FD variables / spline / decade dummies
    capture drop FD_*
    capture drop dc1_lgdppc_MA15
    capture drop decade decind*

    * Panel setup
    sort region_i year
    xtset region_i year

    * Quick describe
    describe iso year load_pc lgdppc_MA15

    ***********************
    * 1. Generate decade fixed effect indicators (1=1971-80, 2=1981-90, etc.).
    *    tabulate with gen() creates binary indicators decind1–decind4,
    *    used later to interact with temperature in the FD regressions. 
    ***********************
    qui gen decade = .
    qui replace decade = 1 if inrange(year,1971,1980)
    qui replace decade = 2 if inrange(year,1981,1990)
    qui replace decade = 3 if inrange(year,1991,2000)
    qui replace decade = 4 if inrange(year,2001,2012)
    qui tab decade, gen(decind)

    *********************************
    * 2. Construct the income spline dc1_lgdppc_MA15 = lgdppc_MA15 - max(lgdppc_MA15),
    *    where the max is taken within large income group 1 (lgid==1) separately
    *    for other_energy and electricity. This de-means income at the top of the
    *    lower income group, creating a spline that is zero at the group boundary
    *    and negative below it, consistent with the CIL two-segment income specification.
    *********************************
    gen dc1_lgdppc_MA15 = .

    foreach p in "other_energy" "electricity" {
        qui summ lgdppc_MA15 if largegpid == 1 & product == "`p'"
        replace dc1_lgdppc_MA15 = lgdppc_MA15 - r(max) if product == "`p'"
    }

    ************************************
    * 3. First difference load_pc      *
    ************************************
    gen FD_load_pc = load_pc - L1.load_pc

    ************************************************************
    * 4. FD income group and income × income group             *
    * Income group variables (defined upstream in the merge script):
    *   largegpid  = large income group ID (1=lower, 2=upper)
    *   largeind1/2 = binary indicators for each large income group
    *   ind1–ind10  = binary indicators for income decile membership
    ************************************************************
    forval lg = 1/2 {
        qui gen FD_largeind`lg' = largeind`lg' - L1.largeind`lg'
        qui gen double FD_I`lg'lgdppc_MA15 = ///
            (  lgdppc_MA15 * largeind`lg'   ) - ///
            (L1.lgdppc_MA15 * L1.largeind`lg')
    }

    ************************************
    * 5. First difference precip       *
    ************************************
    forval i = 1/2 {
        qui gen double FD_precip`i' = ///
            precip`i'_`suf' - L1.precip`i'_`suf'
    }

    ******************************************************************
    * 6. FD temp, temp × year, temp × year^2, temp × decade,        *
    *    polyBelow / polyAbove interactions with year               *
    * Note on year variables used in interactions:
    *   year   = calendar year (raw)
    *   cyear  = centered year (year - some baseline, defined upstream in merge script)
    *   pyear  = post-period year indicator or polynomial (defined upstream)
    *   p80yr  = year interacted with a post-1980 indicator; captures differential
    *             trends in the cold-side response after 1980
    ******************************************************************
    forval i = 1/4 {

        * temp & polyBelow
        qui gen double FD_temp`i'       = temp`i'_`suf'      - L1.temp`i'_`suf'
        qui gen double FD_polyBelow`i'  = polyBelow`i'_`suf' - L1.polyBelow`i'_`suf'

        foreach yr in year cyear pyear p80yr {

            * temp × year
            qui gen double FD_`yr'temp`i' = ///
                (`yr' * temp`i'_`suf') - (L1.`yr' * L1.temp`i'_`suf')

            * temp × year^2
            qui gen double FD_`yr'2temp`i' = ///
                (`yr' * `yr' * temp`i'_`suf') - ///
                (L1.`yr' * L1.`yr' * L1.temp`i'_`suf')
        }

        * temp × decade
        forval dg = 1/2 {
            qui gen double FD_D`dg'temp`i' = ///
                (decind`dg' * temp`i'_`suf') - (L1.decind`dg' * L1.temp`i'_`suf')
        }

        * polyBelow × year post-1980 (coldside p80)
        qui gen double FD_p80yr_polyBelow`i' = ///
            (p80yr * polyBelow`i'_`suf') - (L1.p80yr * L1.polyBelow`i'_`suf')

        * polyBelow × year (coldside)
        qui gen double FD_year_polyBelow`i' = ///
            (year * polyBelow`i'_`suf') - (L1.year * L1.polyBelow`i'_`suf')

        * polyAbove × year post-1980 (two-sided p80)
        qui gen double FD_p80yr_polyAbove`i' = ///
            (p80yr * polyAbove`i'_`suf') - (L1.p80yr * L1.polyAbove`i'_`suf')
    }

    *********************************************************
    * 7. FD temp × income decile                            *
    *********************************************************
    forval lg = 1/10 {
        forval i = 1/4 {
            gen double FD_I`lg'temp`i' = ///
                ( temp`i'_`suf' * ind`lg' ) - ///
                ( L1.temp`i'_`suf' * L1.ind`lg' )
        }
    }

    ******************************************************************
    * 8. FD temp × year × income spline                              *
    ******************************************************************
    forval lg = 1/2 {
        forval i = 1/4 {

            foreach yr in year cyear pyear p80yr {

                qui gen double FD_dc1_lgdppc_MA15`yr'I`lg'temp`i' = ///
                    ( dc1_lgdppc_MA15 * temp`i'_`suf' * largeind`lg' * `yr' ) ///
                  - ( L1.dc1_lgdppc_MA15 * L1.temp`i'_`suf' * L1.largeind`lg' * L1.`yr' )
            }

            * polyBelow interactions with income spline
            qui gen double FD_lgdppc_MA15p80yrI`lg'polyBelow`i' = ///
                ( dc1_lgdppc_MA15 * polyBelow`i'_`suf' * largeind`lg' * p80yr ) ///
              - ( L1.dc1_lgdppc_MA15 * L1.polyBelow`i'_`suf' * L1.largeind`lg' * L1.p80yr )

            qui gen double FD_lgdppc_MA15yearI`lg'polyBelow`i' = ///
                ( dc1_lgdppc_MA15 * polyBelow`i'_`suf' * largeind`lg' * year ) ///
              - ( L1.dc1_lgdppc_MA15 * L1.polyBelow`i'_`suf' * L1.largeind`lg' * L1.year )

            * polyAbove interactions with income spline
            qui gen double FD_lgdppc_MA15p80yrI`lg'polyAbove`i' = ///
                ( dc1_lgdppc_MA15 * polyAbove`i'_`suf' * largeind`lg' * p80yr ) ///
              - ( L1.dc1_lgdppc_MA15 * L1.polyAbove`i'_`suf' * L1.largeind`lg' * L1.p80yr )
        }
    }

    *********************************************************************
    * 9. FD temp × year^2 × income spline                               *
    *********************************************************************
    forval lg = 1/2 {
        forval i = 1/4 {
            foreach yr in year cyear pyear p80yr {
                qui gen double FD_dc1_lgdppc_MA15`yr'2I`lg'temp`i' = ///
                    ( dc1_lgdppc_MA15 * temp`i'_`suf' * largeind`lg' * `yr' * `yr' ) ///
                  - ( L1.dc1_lgdppc_MA15 * L1.temp`i'_`suf' * L1.largeind`lg' * L1.`yr' * L1.`yr' )
            }
        }
    }

    *********************************************************
    * 10. FD income spline × temp (and polyBelow)            *
    *********************************************************
    * income spline × temp
    forval lg = 1/2 {
        forval i = 1/4 {
            qui gen double FD_dc1_lgdppc_MA15I`lg'temp`i' = ///
                ( dc1_lgdppc_MA15 * temp`i'_`suf' * largeind`lg' ) - ///
                ( L1.dc1_lgdppc_MA15 * L1.temp`i'_`suf' * L1.largeind`lg' )
        }
    }

    * income spline × polyBelow
    forval lg = 1/2 {
        forval i = 1/4 {
            qui gen double FD_dc1_lgdppc_MA15I`lg'polyBelow`i' = ///
                ( dc1_lgdppc_MA15 * polyBelow`i'_`suf' * largeind`lg' ) - ///
                ( L1.dc1_lgdppc_MA15 * L1.polyBelow`i'_`suf' * L1.largeind`lg' )
        }
    }

    *********************************************************
    * 11. FD TINV cross terms: long-run HDD/CDD × polyAbove/polyBelow.
    *     Two sources are handled with fallback logic:
    *     (a) Panel-level TINV (cdd20_TINV_*, hdd20_TINV_*): constructed as
    *         country-year mean HDD/CDD from the aggregated CSV. Used if present.
    *     (b) Pixel-level cross terms (polyAbove*_x_cdd_*, polyBelow*_x_hdd_*):
    *         constructed in Python before spatial aggregation, so they correctly
    *         capture within-country covariance between temperature exposure and
    *         long-run climate. These are preferred and will OVERWRITE (a) if present.
    *     This fallback ensures the script runs even if the Python pipeline was
    *     run without the pixel-level cross-term step.
    *********************************************************
    forval i = 1/4 {

        * --- Panel-level TINV scalars (always available if Python added them) ---
        capture confirm variable cdd20_TINV_`suf'
        if !_rc {
            qui gen double FD_cdd20_TINVtemp`i' = ///
                ( cdd20_TINV_`suf' * polyAbove`i'_`suf' ) - ///
                ( cdd20_TINV_`suf' * L1.polyAbove`i'_`suf' )
        }

        capture confirm variable hdd20_TINV_`suf'
        if !_rc {
            qui gen double FD_hdd20_TINVtemp`i' = ///
                ( hdd20_TINV_`suf' * polyBelow`i'_`suf' ) - ///
                ( hdd20_TINV_`suf' * L1.polyBelow`i'_`suf' )
        }

        * --- If pixel-level cross terms exist, prefer those (CIL-style) ---
        *     polyAbove`i'_x_cdd_`suf' = ∑ w_z polyAbove_i(T_z) · CDD_z^TINV, etc.
        capture confirm variable polyAbove`i'_x_cdd_`suf'
        if !_rc {
            capture drop FD_cdd20_TINVtemp`i'
            qui gen double FD_cdd20_TINVtemp`i' = ///
                polyAbove`i'_x_cdd_`suf' - L1.polyAbove`i'_x_cdd_`suf'
        }

        capture confirm variable polyBelow`i'_x_hdd_`suf'
        if !_rc {
            capture drop FD_hdd20_TINVtemp`i'
            qui gen double FD_hdd20_TINVtemp`i' = ///
                polyBelow`i'_x_hdd_`suf' - L1.polyBelow`i'_x_hdd_`suf'
        }
    }


    ********************************************
    * 12. Compress and save regression-ready dataset*
    ********************************************
    compress
    save "$REGD/`product'_TINV_clim_regsort.dta", replace
    di "Saved: $REGD/`product'_TINV_clim_regsort.dta"
}

************************************************
* End of 0_build_FD_dataset_all.do
************************************************
