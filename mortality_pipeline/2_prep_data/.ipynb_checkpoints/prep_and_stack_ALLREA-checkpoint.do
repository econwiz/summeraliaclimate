/**************************************************************************
* prep_and_stack_ALLREA.do
* 
* COMBINED: Preps per-product panels AND stacks into ALLREA
*
* Step 1: For each product, reads R-merged panel, applies filters,
*         winsorizes, constructs FE IDs → saves panel_prepped_for_regressions_PRODUCT_aligned.dta
* Step 2: Stacks all products into ALLREA panel
*
**************************************************************************/

clear all
set more off

* ============================================================================
* CONFIG — change resolution here only
* ============================================================================

local resolution "1deg"   /* ← CHANGE: "025deg" or "1deg" */

global CODE "/user/ab5405/summeraliaclimate/code"

* R merge script outputs (input for Step 1)
global ROUT   "$CODE/mortality_pipeline/prep_panels/`resolution'"

* Step 1 outputs / Step 2 inputs
global OUTPUT "$CODE/mortality_pipeline/output_panels/`resolution'"
capture mkdir "$OUTPUT"

local products GMFD ERA5 JRA_3Q MERRA2

di ""
di "========================================================================"
di "CONFIGURATION:"
di "  Resolution:  `resolution'"
di "  R panels in: $ROUT"
di "  Output dir:  $OUTPUT"
di "========================================================================"
di ""

* ============================================================================
* STEP 1 — Prep each product panel
* ============================================================================

foreach product of local products {

    local inpanel  "$ROUT/global_mortality_panel_public_`product'.dta"
    local outpanel "$OUTPUT/panel_prepped_for_regressions_`product'_aligned.dta"

    capture confirm file "`inpanel'"
    if _rc {
        di as yellow "⚠️  Not found: `inpanel' (skipping `product')"
        continue
    }

    di ""
    di "--------------------------------------------------------------------"
    di "STEP 1: Prepping `product'"
    di "--------------------------------------------------------------------"

    use "`inpanel'", clear

    * 0a) Coerce agegroup to string
    capture confirm numeric variable agegroup
    if !_rc {
        capture decode agegroup, gen(age_s)
        if !_rc {
            drop agegroup
            rename age_s agegroup
        }
    }
    replace agegroup = strtrim(agegroup)

    * 0b) Build 'total' age group
    preserve
        keep if inlist(agegroup,"0-4","5-64","65+")
        collapse (sum) deaths population, by(iso adm1_id adm2_id year)
        gen agegroup = "total"
        tempfile total
        save `total'
    restore
    append using `total'

    * 1) Recompute raw deathrate
    capture drop deathrate
    gen deathrate = (deaths / population) * 100000

    * 2) Sample filters
    gen sample = 1
    replace sample = 0 if year > 2010
    replace sample = 0 if mi(tavg_poly_1_`product')
    replace sample = 0 if mi(prcp_poly_1_`product')
    replace sample = 0 if mi(loggdppc_adm1_avg)
    replace sample = 0 if mi(lr_tavg_`product'_adm1_avg)
    keep if sample==1
    drop sample

    * 3) Winsorize deathrate at 99th pctile within iso×agegroup
    bysort iso agegroup: egen deathrate_p99 = pctile(deathrate), p(99)
    gen deathrate_w99 = cond(deathrate > deathrate_p99 & !mi(deathrate), deathrate_p99, deathrate)
    drop deathrate_p99

    * 4) FE IDs
    egen adm0_code        = group(iso)
    egen adm1_code        = group(iso adm1_id)
    replace adm2_id       = adm1_id if iso == "JPN"
    egen adm2_code        = group(iso adm1_id adm2_id)
    egen adm0_agegrp_code = group(iso agegroup)
    egen adm1_agegrp_code = group(iso adm1_id agegroup)

    * 5) adm0_code2 (EU pooled)
    gen adm0_code2 = 41
    replace adm0_code2 = 4  if adm0_code == 4
    replace adm0_code2 = 6  if adm0_code == 6
    replace adm0_code2 = 7  if adm0_code == 7
    replace adm0_code2 = 15 if adm0_code == 15
    replace adm0_code2 = 23 if adm0_code == 23
    replace adm0_code2 = 27 if adm0_code == 27
    replace adm0_code2 = 40 if adm0_code == 40

    gen iso2 = iso
    replace iso2 = "EUR" if adm0_code2==41

    * 6) China time series flag
    gen CHN_ts = 1
    replace CHN_ts = 2 if year >= 2004 & iso=="CHN"

    * 7) Rename product-specific vars to generic names
    capture confirm variable tavg_poly_1_`product'
    if !_rc {
        rename tavg_poly_1_`product' tavg_poly_1
        capture rename tavg_poly_2_`product' tavg_poly_2
        capture rename tavg_poly_3_`product' tavg_poly_3
        capture rename tavg_poly_4_`product' tavg_poly_4
    }
    capture confirm variable prcp_poly_1_`product'
    if !_rc {
        rename prcp_poly_1_`product' prcp_poly_1
        capture rename prcp_poly_2_`product' prcp_poly_2
    }
    capture confirm variable lr_tavg_`product'_adm1_avg
    if !_rc {
        rename lr_tavg_`product'_adm1_avg lr_tavg_adm1_avg
    }

    * 8) Save
    save "`outpanel'", replace
    di as res "✅ Saved: `outpanel'"
    di "   Observations: " _N
}

di ""
di "========================================================================"
di "STEP 1 COMPLETE — all products prepped"
di "========================================================================"

* ============================================================================
* STEP 2 — Stack into ALLREA
* ============================================================================

di ""
di "========================================================================"
di "STEP 2: Stacking into ALLREA"
di "========================================================================"

tempfile base
local first = 1

foreach product of local products {

    local fpath "$OUTPUT/panel_prepped_for_regressions_`product'_aligned.dta"

    capture confirm file "`fpath'"
    if _rc {
        di as yellow "⚠️  Not found: `fpath' (skipping)"
        continue
    }

    di "Loading `product' → `fpath'"
    use "`fpath'", clear

    gen str20 proj_base = "`product'"

    keep iso iso2 adm0_code adm0_code2 adm1_id adm1_code adm2_id adm2_code ///
         year agegroup CHN_ts deaths population deathrate_w99 ///
         tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
         prcp_poly_1 prcp_poly_2 loggdppc_adm1_avg lr_tavg_adm1_avg ///
         adm0_agegrp_code adm1_agegrp_code ///
         proj_base

    order iso iso2 adm0_code adm0_code2 adm1_id adm1_code adm2_id adm2_code ///
          year agegroup CHN_ts deaths population deathrate_w99 ///
          tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
          prcp_poly_1 prcp_poly_2 loggdppc_adm1_avg lr_tavg_adm1_avg ///
          proj_base

    di "  Observations: " _N

    if `first' {
        save `base'
        local first = 0
    }
    else {
        append using `base'
        save `base', replace
    }
}

capture confirm file `base'
if _rc {
    di as error "❌ ERROR: No prepped product panels found. Check Step 1 output."
    exit 1
}

use `base', clear
sort proj_base iso adm1_id adm2_id year agegroup

di ""
di "Breakdown by product:"
tab proj_base

local out_dta "$OUTPUT/panel_prepped_for_regressions_ALLREA_aligned.dta"
local out_csv "$OUTPUT/panel_prepped_for_regressions_ALLREA_aligned.csv"

save "`out_dta'", replace
export delimited using "`out_csv'", replace

di ""
di "========================================================================"
di "SUCCESS — ALLREA panel written to:"
di "  ✅ `out_dta'"
di "  ✅ `out_csv'"
di "========================================================================"
