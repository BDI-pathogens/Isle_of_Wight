// Code to perform Synthetic control analysis for "Epidemiological changes on the Isle of
// Wight after the launch of the NHS Test and Trace programme: a preliminary analysis", 
// Kendall et al., The Lancet Digital Health, 2020

// Date: 11/07/2020
// Author: Luke Milsom

clear all

{/* Initialize */ 
	set more off
	set rmsg on

	// use [dis "`c(username)'"] to find your username 
	if c(username) == "lukem" {
		*global folder = "C:\Users\lukem\Dropbox\App survey\Evaluation\IoW" // laptop
		global folder = "E:\Dropbox\App survey\Evaluation\IoW" // pc
	}
	
	// set folder path
	cd "${folder}"
	
	// globals
	global fig = "ylab(, angle(0) nogrid) graphregion(color(white))"

}

{/* Upload and clean latest cases data */

// sort out pop estimates data (from ONS: https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland )
	clear
	import excel "uk_pop_estimates_2020.xlsx", sheet("Sheet1") firstrow
	rename Geography1 type
	rename Allages pop
	rename Code area_code
	rename Name area_name
	sa "uk_pop_estimates_2020.dta", replace


// import data from GIT
import delimited "https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-cases-uk.csv", varnames(1) clear

	rename area							area_name
	rename areacode 					area_code

// only consider English counties
keep if country == "England" 
drop country 

// merge on populations
merge m:1 area_code using "uk_pop_estimates_2020.dta"
drop if _merge == 2
drop if _merge == 1 // Buckinghamshire 
drop _merge type

// set date
gen date_new = date(date, "YMD")
format date_new %td
drop date 
rename date_new date

// gen numerical ID
drop area_code
encode area_name , gen(area_code)

// transform to cases per 10,000
destring totalcases, replace
replace totalcases = ( totalcases / pop )*10000

sa uk_cases_updated.dta , replace

// fill in 0's 
	// create a blank dataset to merge on
	preserve 
		keep date
		duplicates drop
		gen n = _n
		local N = _N
		sa "temp1.dta" , replace
	restore 
		
	keep area_name
	duplicates drop
		
	expand `N'
	
	bysort area_name : gen n = _n
	
	merge m:1 n using "temp1.dta"
	assert _merge == 3
	drop _merge 
	
	merge 1:1 date area_name using "uk_cases_updated.dta"
	assert _merge != 2
	replace totalcases = 0 if missing(totalcases) 
	drop n
	drop _merge 
	
	sort date
	gen days_since_first_case = date - date[1]
	
	// backfill
	bysort area_name (day) : replace totalcases = totalcases[_n-1] if totalcases < totalcases[_n-1] & day !=0 
	
	// replace missing time invarient data
	bysort area_name (pop) : replace pop = pop[1]
	bysort area_name (area_code) : replace area_code = area_code[1]
	
	// make change in cases variable 
	xtset area_code days_since_first_case
	gen diff_cases = d.totalcases
	
	// for IoW analysis 
	// labelling
	label variable totalcases "Cumulative cases per 10,000"
	label variable days_since "Days since first case"
	label variable diff_cases "Daily new cases per 10,000"
	
	// trim
	drop if days_since < 30 // data is patchy for the first month
	
	// isle of wight identifier 
	gen IoW = ( area_code == 58 )
	
	// drop city of london
	drop if area_code == 25
	
	sa cases_isle_of_wight.dta , replace

}

{/* Add demographic etc. data */

// data from here (ONS, 2011 census) : https://www.nomisweb.co.uk/

// age structure 
	import excel "nomis_age_structure_census.xls", sheet("Data") cellrange(A9:AI564) firstrow clear
	drop in 1

	rename localau area_name 
	rename E pc_0to4
	rename G pc_5to7
	rename I pc_8to9
	rename K pc_10to14
	rename M pc_15
	rename O pc_16to17
	rename Q pc_18to19
	rename S pc_20to24
	rename U pc_25to29
	rename W pc_30to44
	rename Y pc_45to59
	rename AA pc_60to64
	rename AC pc_65to74
	rename AE pc_75to84
	rename AG pc_85to89
	rename AI pc_90plus
	
	keep area_name pc_*
	
	// numeric
	drop if _n > 174	
	destring(pc_*) , replace
	
	gen pc_0to19 = pc_0to4 + pc_5to7 + pc_8to9 + pc_10to14 + pc_15 + pc_16to17 + pc_18to19
	gen pc_20to44 = pc_20to24 + pc_25to29 + pc_30to44 
	gen pc_45to64 = pc_45to59 + pc_60to64
	gen pc_75plus = pc_75to84 + pc_85to89 + pc_90plus
	
	keep pc_0to19 pc_20to44 pc_45to64 pc_65to74 pc_75plus area_name
	
	// make names conformative 
	replace area_name = "Bournemouth, Christchurch and Poole" if area_name == "Bournemouth"
	replace area_name = "Bournemouth, Christchurch and Poole" if area_name == "Poole"
	replace area_name = "Cornwall and Isles of Scilly" if area_name == "Cornwall"
	replace area_name = "Cornwall and Isles of Scilly" if area_name == "Isles of Scilly"
	
	collapse (mean) pc_* , by(area_name)
	
	sa age_structure.dta , replace
	
// ethnicities
	import excel "nomis_ethnicity_census.xls", sheet("Data") cellrange(A9:M186) firstrow clear
	drop in 1
	
	rename E pc_white
	rename G pc_multiple
	rename I pc_asian
	rename K pc_black_african_carribbean
	rename M pc_other
	
	rename localau area_name
	
	keep pc_* area_name 
	
	drop if _n > 174
	
	destring(pc_*) , replace
	
	// make names conformative 
	replace area_name = "Bournemouth, Christchurch and Poole" if area_name == "Bournemouth"
	replace area_name = "Bournemouth, Christchurch and Poole" if area_name == "Poole"
	replace area_name = "Cornwall and Isles of Scilly" if area_name == "Cornwall"
	replace area_name = "Cornwall and Isles of Scilly" if area_name == "Isles of Scilly"
	
	collapse (mean) pc_* , by(area_name)
	
	sa ethnicities.dta , replace
	
// index of multiple deprivation
	// source: Ministry of housing, communities and local government 
	// https://www.gov.uk/government/statistics/english-indices-of-deprivation-2019
	
	import delimited "index_of_mult_deprivation.csv", clear 
	
	rename uppertierlocalauthoritydistrictn area_name
	rename imdrankofaveragerank imd_rank
	rename imdproportionoflsoasinmostdepriv prop_lsoa_10pc_most_deprived
	
	keep area_name imd_rank prop_lsoa_10pc_most_deprived
	
	replace area_name = "Cornwall and Isles of Scilly" if area_name == "Cornwall"
	replace area_name = "Cornwall and Isles of Scilly" if area_name == "Isles of Scilly"
	
	collapse (mean) imd_rank prop , by(area_name)

	sa mult_dep_index.dta ,replace
	
// population density (same source as overall population)
	import delimited "area_and_pop_density.csv", rowrange(7) clear 
	
	rename v2 area_name 
	rename v6 pop_density
	
	keep area_name pop_density
	
	replace pop_density = subinstr(pop_density, "," , "" , .)
	destring pop_density , replace
	
	drop if missing(area_name)
	
	replace area_name = "Cornwall and Isles of Scilly" if area_name == "Cornwall"
	replace area_name = "Cornwall and Isles of Scilly" if area_name == "Isles of Scilly"
	
	collapse (mean) pop_density , by(area_name)

	sa pop_density.dta , replace
	
	
// merge onto main dataset	
	use cases_isle_of_wight.dta , clear 
	
	merge m:1 area_name using age_structure.dta
	drop if _merge == 1 // city of lnd
	drop if _merge == 2 // welsh
	drop _merge 
	
	merge m:1 area_name using ethnicities.dta
	drop if _merge == 1
	drop if _merge == 2
	drop _merge 
	
	merge m:1 area_name using mult_dep_index.dta
	drop if _merge == 1
	drop if _merge == 2
	drop _merge 
	
	merge m:1 area_name using pop_density.dta 
	drop if _merge == 1
	drop if _merge == 2
	drop _merge 	
	
	sa cases_isle_of_wight_demo.dta , replace	
}

{/* Get R backdated approach data from outside models */

	import delimited "time_series_150_regions_backcalculated.csv" , clear
	
	drop v1 iow
	
	gen date_new = date(date, "YMD")
	format date_new %td
	drop dates 
	rename date_new date
	
	sum date , d
	local min_date = r(min)
	
	rename area area_name
	
	sa R_data_michelle_backdated.dta, replace
	
	// make average R series (not inc. IoW)
	drop if area == "Isle of Wight"
	
	collapse (mean) r , by(date)
	
	sa av_rest_of_eng_R_by_date_backdated.dta , replace
	
	// merge on
	use cases_isle_of_wight_demo.dta , clear
	
	*drop if date < `min_date'
	
	// now merge 
	merge 1:1 date area_name using R_data_michelle_backdated.dta 
	drop if _merge != 3 // Buckinghamshire, City of Lnd, 11,12,13 June from R data
	drop _merge 
	
	// start on the 1st of April
	drop if date <= 22001

	sa cases_isle_of_wight_r_backdated.dta , replace 
	
}

// set matching and validation variables for each model	
	#delimit ;

		local R_train "r(58(1)62) r(63(1)69) r(70(1)77)" ;
		local R_validation "r(78(1)82) r(83(1)89) r(90(1)96)" ;
		local age "pc_0to19 pc_20to44 pc_45to64 pc_65to74 pc_75plus" ;
		local ethnicity "pc_white pc_multiple pc_asian pc_black_african_carribbean pc_other" ;
		local poverty_density "pc_other imd_rank prop_lsoa_10pc_most_deprived pop_density" ;

		global matching_vars_training  = `" 
			"`R_train'" 
			"`R_train' `age'" 
			"`R_train' `ethnicity'"
			"`R_train' `poverty_density'"
			"`R_train' `age' `ethnicity'"
			"`R_train' `age' `poverty_density'"
			"`R_train' `age' `ethnicity' `poverty_density'"
			"`R_train' `ethnicity' `poverty_density'"
			"' ;
			
		global matching_vars_validation  = `" 
			"`R_validation'" 
			"`R_validation' `age'" 
			"`R_validation' `ethnicity'"
			"`R_validation' `poverty_density'"
			"`R_validation' `age' `ethnicity'"
			"`R_validation' `age' `poverty_density'"
			"`R_validation' `age' `ethnicity' `poverty_density'"
			"`R_validation' `ethnicity' `poverty_density'"
			"' ;
		
	#delimit cr

{/* Construct synthetic controls */
/* Notes
	- 96 days after first case = 5th May
	- area number 55 = IoW
*/

	use cases_isle_of_wight_r_backdated , clear
	
	// need numeric identifier for SC program
	egen area_number = group(area_name)
	
	// make crosswalk from area_number to area_name
	preserve 
		keep area_number area_name
		duplicates drop
		sa cross_walk_area_num_to_area_name.dta , replace	
	restore 
	
	// xtset for SC
	xtset area_number day
		

	// run SC for each model
	local num_models : word count $matching_vars_training
	forv i = 1(1)`num_models' {	
	
		local training_vars : word `i' of $matching_vars_training
		local validation_vars : word `i' of $matching_vars_validation
	
		levelsof area_number , local(levels)
		local n = 0
		
		foreach area of local levels {
		
			local n = `n' + 1
			
			di "**** `i' ; `n' ****"
			
			qui {					
					// find V
					preserve 
						keep if day < 96 // only consider the pre-TTI period
						synth r `training_vars' , trunit(`area') trperiod(78) mspeperiod(78(1)95)
						mat V = e(V_matrix)
					restore 
				
					// get V in the right form for input into SC program 
					preserve 
						clear
						svmat V
						gen n = _n 
						sum n
						local n_max = r(max)
						local n_min = r(min)
						
						local v_`area' = " "
						forv j = `n_min'(1)`n_max' {
							sum V`j' if n == `j'
							local m = r(mean)
							local v_`area' = " `v_`area'' `m' "
						}
					restore 
			
				synth r `validation_vars' , trunit(`area') trperiod(96) customV(`v_`area'')
					
				mat Y_treat_`area'_`i' = e(Y_treated)
				mat Y_synthetic_`area'_`i' = e(Y_synthetic)
				mat gap_`area'_`i' = Y_treat_`area'_`i' - Y_synthetic_`area'_`i'
				
			}
		}
	}

}

{/* Save and export constructed datasets */

	local num_models : word count $matching_vars_training
	forv i = 1(1)`num_models' {	

		levelsof area_number , local(levels)
		
		preserve
		clear
		foreach area of local levels {
		
			// read in data
			svmat Y_treat_`area'_`i'
			rename Y_treat_`area'_`i' actual_data`area'
			
			svmat Y_synthetic_`area'_`i'
			rename Y_synthetic_`area'_`i' synthetic_cntrl`area'
			
			svmat gap_`area'_`i'
			rename gap_`area'_`i' difference`area'
			
		}
		
		// reshape
		gen n = _n
		reshape long actual_data synthetic_cntrl difference , i(n) j(area_number)
		
		// create area names and dates
		merge m:1 area_number using cross_walk_area_num_to_area_name.dta
		assert _merge == 3
		drop _merge 
		
		gen date = n + 22001
		format date %td
		
		// saving and export to .csv
		save sc_dataset_scenario_`i'.dta , replace
		export delimited using sc_dataset_scenario_`i'.csv , replace
	
		restore
	}

}

{/* Create figure 6 left hand panels */	
	
	// run SC for each model
	local num_models : word count $matching_vars_training
	
	forv i = 1(1)`num_models' {
	preserve
	
		// load data from previous step
		clear
		svmat Y_treat_55_`i'
		svmat Y_synthetic_55_`i'
		
		// create date
		gen n = _n
		gen date = n + 22001
		
		// merge on average R
		merge 1:1 date using av_rest_of_eng_R_by_date_backdated.dta 
		drop if _merge != 3
		drop _merge 
		rename r r_av
		
		// cut off the last 10 days due to incomplete data
		sum date , d
		local max = r(max)
		drop if date > `max' - 10
		sum date ,d 
		local max = r(max)
			
		tw 	( line Y_treat_55 date , lc(red*0.8) ) (line Y_syn date , lp(dash) lc(gs8) ) ///
			( line r_av date , lc(gs4) lp(dot) ) ///
			, $fig  xline(96 , lp(dash) lc(black))  ///
			ytitle("Estimated R") ///
			xlab( #6 , format("%tdmd")) ///
			xline(22040 , lc(black) lp(dash) lw(medthick) ) ///
			text(1.3 22040.2 "Start of TT pilot" , place(e) size(small)) ///
			legend(order(1 2 3) label(2 "Isle of Wight synthetic control") ///
			label(1 "Isle of Wight") label(3 "Rest of England average") ///
			pos(7) ring(0) c(1) region(lc(white)) ) xtitle("") 
		graph display , xsize(6.5)	scale(1)	
		
		graph export "iow_sc_xv_backcalc_basic_`i'.png" , width(2000) replace

	restore
	}
}

{/* Create figure 6 right hand side panels */
	
	// run SC for each model
	local num_models : word count $matching_vars_training
	
	forv i = 1(1)`num_models' {
		
		// create confidence bands 
		preserve
			levelsof area_number , local(levels)
			
			foreach level of local levels {
				svmat gap_`level'_`i'
				rename gap_`level'_`i'1 gap`level'
				cap gen n = _n
			}
		
			reshape long gap , i(n) j(area_num)
			
			bysort n (gap) : gen ci_max_95 = gap[_N - 5]
			bysort n (gap) : gen ci_min_95 = gap[5]
			
			keep n ci_max_95 ci_min_95
			duplicates drop
			
			sa iow_ci_`i'.dta , replace			
		restore
	
		// create lines
		preserve
			levelsof area_number , local(levels)

			clear
			local graph = ""
			foreach level of local levels {
				svmat gap_`level'_`i'
				cap gen n = _n
				local lc = "lc(gs4)"
				local graph = "`graph' ( line gap_`level'_`i'1 date , `lc' lw(vvthin) ) "
			}
			
			// gen empirical mean
			egen mean = rowmean(gap_*)
		
			// find dates
			gen date = n + 22001
			
			// cut off the last 10 days 
			sum date , d
			local temp = r(max)
			drop if date > `temp' - 10

			// merge on confidence bands
			merge m:1 n using iow_ci_`i'.dta
			
			// trick for graph
			rename gap_55_`i'1 true 
			gen gap_55_`i'1 = true 
		
			tw `graph'  ///
				(rarea ci_min_95 ci_max_95 date , lw(none) col(gs8%50) ) ( line true date , lc(red*0.8) ) ///
				, $fig yline(0, lstyle(forground)) xline(99 , lp(dash) lc(black))  ///
					ytitle("Difference in estimated R") ///
					xlab( #6 , format("%tdmd")) ///
					xline(22040 , lc(black) lp(dash) lw(medthick) ) ///
					text(0.9 22040.2 "Start of TT pilot" , place(e) size(small)) ///
					legend(order(148 146 147) label(146 "All other UTLAs") label(148 "Isle of Wight") ///
					label(147 "Non rejection region") ///
					pos(7) ring(0) c(1) region(lc(white)) ) xtitle("") yscale(range(-1 1)) ylab( -1(0.5)1 ) 

			graph display , xsize(6.5)	scale(1)
			
			graph export "iow_sc_xv_backcalc_nonrejection_`i'.png" , width(2000) replace
			
		restore
	
	}
		
}

{/* Create figure S1 - Google mobility data */


	import delimited "Global_Mobility_Report.csv" , clear

	// first just look at IoW on it's own
	keep if sub_region_1 == "Isle of Wight"
	keep date *_percent_*
	
	gen date_new = date(date, "YMD")
	format date_new %td
	drop date 
	rename date_new date

	// may 5th = 22040
	
	tw	(line retail_and date) /// 	
		(line grocery_and_pharmacy_percent_cha date) ///
		(line parks_percent_change_from_baseli date) ///
		(line transit_stations_percent_change_ date) ///
		(line workplaces_percent_change_from_b date) ///
		(line residential_percent_change_from_ date) , ///
		xline(22040 , lp(dash) lw(medthick) lc(black) ) $fig ytitle("% change compared to baseline") ///
		legend( label(1 "Retail") label(2 "Grocery") label(3 "Parks") label(4 "Transit") ///
			label(5 "Workplaces") label(6 "Residential") region(lc(white)) r(1) )  ///
		text(100 22039 "Start of TT pilot" , place(w)) ///
		xtitle("") xlab(#4) title("Google Mobility data for the Isle of Wight" " ")
	
	graph display , xsize(7) scale(0.8)
	
	graph export "iow_google_mobility_data.png" , width(2000) replace
	
	sa mobility_iow.dta , replace

	// plot UK overall
	import delimited "Global_Mobility_Report.csv" , clear
	keep if country_region_code == "GB"
	keep if missing(sub_region_1)

	keep date *_percent_*
	
	gen date_new = date(date, "YMD")
	format date_new %td
	drop date 
	rename date_new date

	// may 5th = 22040
	
		tw 	(line retail_and date) /// 	
			(line grocery_and_pharmacy_percent_cha date) ///
			(line parks_percent_change_from_baseli date) ///
			(line transit_stations_percent_change_ date) ///
			(line workplaces_percent_change_from_b date) ///
			(line residential_percent_change_from_ date) , ///
			xline(22040 , lp(dash) lw(medthick) lc(black) ) $fig ytitle("% change compared to baseline") ///
			legend( label(1 "Retail") label(2 "Grocery") label(3 "Parks") label(4 "Transit") ///
				label(5 "Workplaces") label(6 "Residential") region(lc(white)) r(1) )  ///
			text(100 22039 "Start of TT pilot" , place(w)) ///
			xtitle("") xlab(#4) title("Google Mobility data for the UK overall" " ")
		
		graph display , xsize(7) scale(0.8)
		
		graph export "iow_google_mobility_data_uk_overall.png" , width(2000) replace

	// make comparison figures 
	rename retail iow_retail
	rename grocery iow_grocery
	rename parks iow_parks
	rename transit iow_transit
	rename workplaces iow_work
	rename residential iow_residential
	
	merge 1:1 date using mobility_iow.dta 
	
	rename retail uk_retail
	rename grocery uk_grocery
	rename parks uk_parks
	rename transit uk_transit
	rename workplaces uk_work
	rename residential uk_residential
	
	tw	(line iow_retail date , lc(gs6) ) /// 	
		(line iow_grocery date , lc(navy) ) ///
		(line iow_parks date , lc(maroon) ) ///
		(line iow_transit date , lc(forest_green) ) ///
		(line iow_work date , lc(dkorange) ) ///
		(line iow_residential date , lc(teal) ) ///
		(line uk_retail date , lc(gs6) lp(dash) ) ///
		(line uk_grocery date , lc(navy) lp(dash) ) ///
		(line uk_parks date , lc(maroon) lp(dash) ) ///
		(line uk_transit date , lc(forest_green) lp(dash) ) ///
		(line uk_work date , lc(dkorange) lp(dash) ) ///
		(line uk_residential date , lc(teal) lp(dash)) , ///
			xline(22040 , lp(dash) lw(medthick) lc(black) ) $fig ytitle("% change compared to baseline") ///
			legend( order(1 2 3 4 5 6) label(1 "Retail") label(2 "Grocery") label(3 "Parks") label(4 "Transit") ///
				label(5 "Workplaces") label(6 "Residential") region(lc(white)) r(1) )  ///
			text(100 22039 "Start of TT pilot" , place(w)) ///
			text( 150 21960 "Solid line = IoW" , place(e) ) text( 140 21960 "Dashed line = UK" , place(e)) /// 
			xtitle("") xlab(#4) title("Google Mobility data for the UK overall" " ")
		
		graph display , xsize(7) scale(0.8)
		
	graph export "iow_google_mobility_data_uk_and_iow.png" , width(2000) replace


}

