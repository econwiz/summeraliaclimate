/**************************************************************************
* prep_panel_ALL_aligned.do
* - Builds 'total' agegroup from bins (important!)
* - Applies sample filters
* - Winsorizes deathrate within iso×agegroup (AFTER filtering, correct order)
* - Constructs FE IDs
* - Loops over observational products: ERA5_025, JRA_3Q, MERRA2, GMFD
**************************************************************************/

* CONFIG
global CODE    "/user/ab5405/summeraliaclimate/code"
global OUTPUT  "$CODE/regressions/output_panels"
local products MERRA2

foreach product of local products {

    * Define input/output paths
    local inpanel  "$CODE/regressions/prep_panels/global_mortality_panel_public_`product'.dta"
    local outpanel "$OUTPUT/panel_prepped_for_regressions_`product'_aligned.dta"

    di "Processing `product' ..."
    use "`inpanel'", clear

    * 0a) Ensure agegroup is string
    capture confirm numeric variable agegroup
    if !_rc {
        capture decode agegroup, gen(age_s)
        if !_rc {
            drop agegroup
            rename age_s agegroup
        }
    }
    replace agegroup = strtrim(agegroup)

    * 0b) Build 'total' agegroup from bins, then append
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

    * 2) Apply sample filters (product-specific vars)
    gen sample = 1
    replace sample = 0 if year > 2010
    replace sample = 0 if mi(tavg_poly_1_`product')
    * If precip missing is allowed, comment out the next line
    replace sample = 0 if mi(prcp_poly_1_`product')
    replace sample = 0 if mi(loggdppc_adm1_avg)
    replace sample = 0 if mi(lr_tavg_`product'_adm1_avg)

    keep if sample==1
    drop sample

    * 3) Winsorize AFTER drops (iso × agegroup, p99)
    bysort iso agegroup: egen deathrate_p99 = pctile(deathrate), p(99)
    gen deathrate_w99 = cond(deathrate > deathrate_p99 & !mi(deathrate), deathrate_p99, deathrate)
    drop deathrate_p99

    * 4) FE IDs
    egen adm0_code          = group(iso)
    egen adm1_code          = group(iso adm1_id)
    replace adm2_id         = adm1_id if iso == "JPN"
    egen adm2_code          = group(iso adm1_id adm2_id)

    egen adm0_agegrp_code   = group(iso agegroup)
    egen adm1_agegrp_code   = group(iso adm1_id agegroup)

    * 5) EU combined code (but keep iso intact)
    gen adm0_code2 = 41
    replace adm0_code2 = 4  if adm0_code == 4   // BRA
    replace adm0_code2 = 6  if adm0_code == 6   // CHL
    replace adm0_code2 = 7  if adm0_code == 7   // CHN
    replace adm0_code2 = 15 if adm0_code == 15  // FRA
    replace adm0_code2 = 23 if adm0_code == 23  // JPN
    replace adm0_code2 = 27 if adm0_code == 27  // MEX
    replace adm0_code2 = 40 if adm0_code == 40  // USA

    gen iso2 = iso
    replace iso2 = "EUR" if adm0_code2==41

    * 6) China split
    gen CHN_ts = 1
    replace CHN_ts = 2 if year >= 2004 & iso=="CHN"

    * 7) Save aligned panel
    save "`outpanel'", replace
    di "Saved aligned panel for `product' → `outpanel'"

    clear
}

di "All products processed successfully."
