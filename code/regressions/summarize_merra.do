*************************************************************************
*           PART 3. Summarize & Export All Specs for MERRA2 (w/ constants)
*************************************************************************

* 1. Point to your output folder and ster directory
global OUTPUT "/user/ab5405/summeraliaclimate/output"
local sterdir "$OUTPUT/age_combined_MERRA2"

* 2. Clear any existing estimates
estimates clear

* 3. Load each spec-file and store it
foreach i of numlist 1/4 {
    estimates use "`sterdir'/pooled_response_spec`i'_public_MERRA2.ster"
    estimates store spec`i'
}

* 4. Quick console check of all four
estimates dir
esttab spec1 spec2 spec3 spec4, ///
    se label ///
    star(* 0.1 ** 0.05 *** 0.01) ///
    keep(_cons tavg_poly_*) ///
    order(_cons tavg_poly_*) ///
    title("MERRA2: Spec 1–4 (with intercept)") ///
    nonotes compress

* 5. Export to CSV (including _cons)
esttab spec1 spec2 spec3 spec4 using "$OUTPUT/age_combined_MERRA2/merra2_results_table.csv", ///
    se label ///
    star(* 0.1 ** 0.05 *** 0.01) ///
    keep(_cons tavg_poly_*) ///
    order(_cons tavg_poly_*) ///
    replace

* 6. Export to LaTeX (including _cons)
esttab spec1 spec2 spec3 spec4 using "$OUTPUT/age_combined_MERRA2/merra2_results_table.tex", ///
    se label ///
    star(* 0.1 ** 0.05 *** 0.01) ///
    title("MERRA2 Mortality Regression Results (Specs 1–4)") ///
    keep(_cons tavg_poly_*) ///
    order(_cons tavg_poly_*) ///
    replace

di "✅ MERRA2 regression tables (with constants) saved in $OUTPUT/age_combined_MERRA2/"

