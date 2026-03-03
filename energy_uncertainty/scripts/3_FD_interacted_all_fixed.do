/************************************************************ * 3_FD_interacted_all.do * Replicate CIL stacked "interacted" FD spec for all products * using UNSUFFIXED regressors (no _GMFD/_ERA5 in variable names). * * Saves: * $OUTDIR/FD_inter_`model_name'_`product'.ster * $OUTDIR/FD_FGLS_inter_`model_name'_`product'.ster * * Optional: * complete_for_plot = 1 reposts e(b), e(V) to include omitted terms as 0 * so predictnl/_b[...] works even if terms were collinear/omitted. ************************************************************/ 
clear all 
set more off 

global PROJECT "/user/ab5405/summeraliaclimate/code/energy_uncertainty" 
global DATA "$PROJECT/data" 
global REGD "$DATA/regression" 
global OUTDIR "$DATA/regression/sters" 
cap mkdir "$OUTDIR" 

* naming only 
local model "TINV_clim" 
local submodel "quadinter" // "" | "lininter" | "quadinter" 

if "`submodel'" != "" local model_name "`model'_`submodel'" 
else local model_name "`model'" 

local products "GMFD ERA5 JRA_3Q MERRA2" 

* 1 = make saved .ster plotting-safe (include omitted terms as 0) 
* 0 = save Stata default e(b), e(V) 
local complete_for_plot 0 

foreach product of local products { 
	di "===============================================" 
	di " Interacted FD spec for `product' (`model_name')" 
	di "===============================================" 

	* Load FD dataset 
	use "$REGD/`product'_TINV_clim_regsort.dta", clear 

	* time set 
	sort region_i year 
	tset region_i year 

	* -------- long run income x income group ---------- 
	local lgdppc_MA15_r "" 
	forval pg = 1/2 { 
		forval lg = 1/2 { 
			local lgdppc_MA15_r /// 
				"`lgdppc_MA15_r' c.indp`pg'#c.indf1#c.FD_I`lg'lgdppc_MA15" 
		} 
	} 

	* -------- large income group dummies ---------- 
	forval pg = 1/2 { 
		forval lg = 1/2 { 
			gen DumIncG`lg'F1P`pg' = FD_largeind`lg' * indf1 * indp`pg' 
		} 
	} 

	* -------- precip: indp × indf1 × FD_precipk (UNSUFFIXED) ---------- 
	local precip_r "" 
	forval pg = 1/2 { 
		forval k = 1/2 { 
			local precip_r /// 
				"`precip_r' c.indp`pg'#c.indf1#c.FD_precip`k'" 
		} 
	} 

	* -------- temp: indp × indf1 × FD_tempk (UNSUFFIXED) ---------- 
	local temp_r "" 
	forval pg = 1/2 { 
		forval k = 1/2 { 
			local temp_r /// 
				"`temp_r' c.indp`pg'#c.indf1#c.FD_temp`k'" 
		} 
	} 

	* -------- temp x long-run climate (MATCH CIL construction: redundant lg loop) ---------- 
	local climate_r "" 
	forval pg = 1/2 { 
		forval lg = 1/2 { 
			forval k = 1/2 { 
				local climate_r "`climate_r' c.indp`pg'#c.indf1#c.FD_hdd20_TINVtemp`k'" 
				local climate_r "`climate_r' c.indp`pg'#c.indf1#c.FD_cdd20_TINVtemp`k'" 
			} 
		} 
	} 

	* -------- temp x income spline (UNSUFFIXED) ---------- 
	local income_spline_r "" 
	forval pg = 1/2 { 
		forval lg = 1/2 { 
			forval k = 1/2 { 
				local income_spline_r /// 
					"`income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15I`lg'temp`k'" 
			} 
		} 
	} 

	* -------- temp x year (for lininter/quadinter) ---------- 
	local year_temp_r "" 
	local year_income_spline_r "" 

	if ("`submodel'" == "lininter" | "`submodel'" == "quadinter") { 
		forval pg = 1/2 { 
			forval k = 1/2 { 
				local year_temp_r /// 
					"`year_temp_r' c.indp`pg'#c.indf1#c.FD_yeartemp`k'" 
			} 
		} 

		forval pg = 1/2 { 
			forval lg = 1/2 { 
				forval k = 1/2 { 
					local year_income_spline_r /// 
						"`year_income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15yearI`lg'temp`k'" 
				} 
			} 
		} 
	} 

	* -------- temp x year^2 (for quadinter) ---------- 
	if ("`submodel'" == "quadinter") { 
		forval pg = 1/2 { 
			forval k = 1/2 { 
				local year_temp_r /// 
					"`year_temp_r' c.indp`pg'#c.indf1#c.FD_year2temp`k'" 
			} 
		} 

		forval pg = 1/2 { 
			forval lg = 1/2 { 
				forval k = 1/2 { 
					local year_income_spline_r /// 
						"`year_income_spline_r' c.indp`pg'#c.indf1#c.FD_dc1_lgdppc_MA15year2I`lg'temp`k'" 
				} 
			} 
		} 
	} 

	* ------------------------------------------------------ 
	* Pop weights (MATCH CIL) 
	* ------------------------------------------------------ 
	bysort year product flow: egen year_product_flow_total_pop = total(pop) 
	gen double pop_weight = pop / year_product_flow_total_pop 

	* ------------------------------------------------------ 
	* First-stage (MATCH CIL: pop weights) + residuals for FGLS 
	* ------------------------------------------------------ 
	capture drop resid 
	reghdfe FD_load_pc /// 
		`temp_r' `precip_r' `climate_r' /// 
		`lgdppc_MA15_r' `income_spline_r' /// 
		`year_temp_r' `year_income_spline_r' /// 
		DumInc* [pw=pop_weight], /// 
		absorb(i.flow_i#i.product_i#i.year#i.subregionid) /// 
		vce(cluster region_i) residuals(resid) 

	estimates save "$OUTDIR/FD_inter_`model_name'_`product'", replace 

	* ------------------------------------------------------ 
	* FGLS weights (MATCH CIL exactly) 
	* ------------------------------------------------------ 
	drop if missing(resid) 
	bysort region_i: gen count = _N 
	bysort region_i: egen sum_of_weights_in_FE = total(pop_weight) 
	gen double for_variance_weighting = pop_weight / sum_of_weights_in_FE 
	gen double weighted_residual = for_variance_weighting * resid 
	bysort region_i: egen weighted_mean_resid_FE_level = mean(weighted_residual) 
	gen double square_term_weighted = for_variance_weighting * (resid - weighted_mean_resid_FE_level)^2 
	bysort region_i: egen weighted_residual_variance = total(square_term_weighted) 
	gen double FGLS_weight = pop_weight / weighted_residual_variance 
	drop if count == 1 // singleton FE like CIL 

	* ------------------------------------------------------ 
	* Second-stage (MATCH CIL: FGLS weights) 
	* ------------------------------------------------------ 
	reghdfe FD_load_pc /// 
		`temp_r' `precip_r' `climate_r' /// 
		`lgdppc_MA15_r' `income_spline_r' /// 
		`year_temp_r' `year_income_spline_r' /// 
		DumInc* [pw=FGLS_weight], /// 
		absorb(i.flow_i#i.product_i#i.year#i.subregionid) /// 
		vce(cluster region_i) 

	* ------------------------------------------------------ 
	* Optional: make .ster safe for predictnl/_b[...] (omitted -> 0) 
	* ------------------------------------------------------ 
	if (`complete_for_plot') { 
		local want_all "" 
		local want_all "`want_all' `temp_r' `precip_r' `climate_r'" 
		local want_all "`want_all' `lgdppc_MA15_r' `income_spline_r'" 
		local want_all "`want_all' `year_temp_r' `year_income_spline_r'" 
		unab dumlist : DumInc* 
		local want_all "`want_all' `dumlist'" 

		local want "" 
		foreach nm of local want_all { 
			if strpos(" `want' ", " `nm' ") == 0 local want "`want' `nm'" 
		} 

		local K_want : word count `want' 
		tempname b V bfull Vfull 
		matrix `b' = e(b) 
		matrix `V' = e(V) 

		matrix `bfull' = J(1, `K_want', 0) 
		matrix colnames `bfull' = `want' 
		matrix `Vfull' = J(`K_want', `K_want', 0) 
		matrix rownames `Vfull' = `want' 
		matrix colnames `Vfull' = `want' 

		local i = 0 
		foreach ni of local want { 
			local ++i 
			local ci = colnumb(`b', "`ni'") 
			if (`ci' < .) matrix `bfull'[1,`i'] = `b'[1,`ci'] 

			local ri = rownumb(`V', "`ni'") 
			if (`ri' < .) { 
				local j = 0 
				foreach nj of local want { 
					local ++j 
					local cj = colnumb(`V', "`nj'") 
					if (`cj' < .) matrix `Vfull'[`i',`j'] = `V'[`ri',`cj'] 
				} 
			} 
		} 

		ereturn repost b=`bfull' V=`Vfull', rename 
	} 

	estimates save "$OUTDIR/FD_FGLS_inter_`model_name'_`product'", replace 

	* ------------------------------------------------------ 
	* Export coeff vector (parm, beta) 
	* ------------------------------------------------------ 
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
		"$OUTDIR/FD_FGLS_inter_`model_name'_`product'_coeff.csv", replace 
	restore 

	* ------------------------------------------------------ 
	* Export vcov matrix (long format: parm_i, parm_j, v) 
	* ------------------------------------------------------ 
	matrix V = e(V) 
	local kv = rowsof(V) 
	local vnames : rownames V 

	preserve 
	clear 
	set obs `=`kv'*`kv'' 
	gen strL parm_i = "" 
	gen strL parm_j = "" 
	gen double v = . 

	local r = 1 
	forvalues ii = 1/`kv' { 
		forvalues jj = 1/`kv' { 
			local name_i : word `ii' of `vnames' 
			local name_j : word `jj' of `vnames' 
			replace parm_i = "`name_i'" in `r' 
			replace parm_j = "`name_j'" in `r' 
			replace v = V[`ii',`jj'] in `r' 
			local ++r 
		} 
	} 

	order parm_i parm_j v 
	export delimited using /// 
		"$OUTDIR/FD_FGLS_inter_`model_name'_`product'_vcov_long.csv", replace 
	restore 

	di "" 
	di "Exported coefficient and vcov CSVs for `product'" 
	di "" 

} 
* end foreach product 

************************************************************ 
* End 3_FD_interacted_all.do 
************************************************************/
