
*--------------------------------------------------------*
*  Prepare WHO Excess Mortality dataset                  *
*--------------------------------------------------------*

cd "/Users/zixuandong/Desktop/HP434 summative d"

import excel "WHO_COVID_Excess_Deaths_EstimatesByCountry 21.43.46.xlsx", firstrow clear
gen date = ym(year, month)
format date %tm
keep Country date excessmean
rename excessmean ExcessMortality

duplicates report Country date
duplicates list Country date
collapse (mean) ExcessMortality, by(Country date)

save "WHO_monthly_cleaned.dta", replace

*--------------------------------------------------------*
* Prepare OxCGRT Dataset with Instruments                *
*--------------------------------------------------------*

use "OxCGRT_compact_national_v1.dta", clear

* Convert to usable time format
gen date_daily = date(string(Date, "%8.0f"), "YMD")
format date_daily %td
gen date = mofd(date_daily)
format date %tm

* Collapse to country-month
collapse (mean) C6M_Stay_at_home_requirements ///
         StringencyIndex PopulationVaccinated ///
         H1_Public_information_campaigns H2_Testing_policy ///
         C1M_School_closing C2M_Workplace_closing ///
         H3_Contact_tracing H6M_Facial_Coverings ///
         C5M_Close_public_transport C7M_Restrictions_on_internal_mov ///
         E1_Income_support E3_Fiscal_measures E2_Debt_contract_relief ///
         H4_Emergency_investment_in_healt H5_Investment_in_vaccines ///
         H8M_Protection_of_elderly_people V2E_Education ///
         V3_Vaccine_Financial_Support__su ConfirmedCases ///
		 ConfirmedDeaths GovernmentResponseIndex_Average ///
		 ContainmentHealthIndex_Average EconomicSupportIndex, by(CountryName date)

* Standardized lockdown measure
gen Lockdown_Level = round(C6M_Stay_at_home_requirements)
replace Lockdown_Level = 0 if Lockdown_Level < 0
replace Lockdown_Level = 3 if Lockdown_Level > 3

* Clean country name
rename CountryName Country

*--------------------------------------------------------*
* Create Instruments — Early Adopter & Volatility.       *
*--------------------------------------------------------*

* Tag early adopters: lockdown >= 2 before Mar 2020
gen strong_lockdown = C6M_Stay_at_home_requirements >= 2 if !missing(C6M_Stay_at_home_requirements)
bysort Country (date): egen first_lockdown = min(cond(strong_lockdown, date, .))
gen early_adopter = first_lockdown <= tm(2020m3)
label define early 0 "Late" 1 "Early"
label values early_adopter early

encode Country, gen(Country_ID)
xtset Country_ID date

gen lockdown_change = C6M_Stay_at_home_requirements != L.C6M_Stay_at_home_requirements if !missing(C6M_Stay_at_home_requirements, L.C6M)
bysort Country_ID: egen lockdown_volatility = total(lockdown_change)

* Save with instruments
save "OxCGRT_monthly.dta", replace

*--------------------------------------------------------*
* Merge WHO + OxCGRT + IVs                               *
*--------------------------------------------------------*

use "WHO_monthly_cleaned.dta", clear
merge 1:1 Country date using "OxCGRT_monthly.dta"
keep if _merge == 3
drop _merge

save "merged_dataset.dta", replace

*--------------------------------------------------------*
* Filter Countries with Sufficient Observations          *
*--------------------------------------------------------*
use "merged_dataset.dta", clear
summarize ExcessMortality Lockdown_Level ConfirmedCases PopulationVaccinated ///
      GovernmentResponseIndex_Average ContainmentHealthIndex_Average ///
	  StringencyIndex E1_Income_support EconomicSupportIndex ///
	  H8M_Protection_of_elderly_people

* Drop if missing key vars (including IVs now)
drop if missing(ExcessMortality, Lockdown_Level, StringencyIndex, PopulationVaccinated, early_adopter, lockdown_volatility, ConfirmedCases, ConfirmedDeaths, C2M_Workplace_closing, H6M_Facial_Coverings)

* Keep only countries with 24+ observations
bysort Country: gen n_months = _N
drop if n_months < 24
drop n_months

save "final_analysis_data_with_IV.dta", replace

*--------------------------------------------------------*
* Data Preparation                                       *
*--------------------------------------------------------*
use "final_analysis_data_with_IV.dta", clear
xtset Country_ID date, monthly

gen L1_lockdown = L1.Lockdown_Level
gen L2_lockdown = L2.Lockdown_Level
gen L3_lockdown = L3.Lockdown_Level

*--------------------------------------------------------*
* Fixed & Random Effects + Hausman                       *
*--------------------------------------------------------*
xtreg ExcessMortality L2_lockdown ConfirmedCases PopulationVaccinated ///
      GovernmentResponseIndex_Average ContainmentHealthIndex_Average ///
	  StringencyIndex E1_Income_support EconomicSupportIndex ///
	  H8M_Protection_of_elderly_people, fe
estimates store fe

xtreg ExcessMortality L2_lockdown ConfirmedCases PopulationVaccinated ///
      GovernmentResponseIndex_Average ContainmentHealthIndex_Average ///
	   StringencyIndex E1_Income_support EconomicSupportIndex ///
	   H8M_Protection_of_elderly_people, re
estimates store re

hausman fe re

*--------------------------------------------------------*
* IV Regression: 2SLS with Lagged Lockdowns              *
*--------------------------------------------------------*
* Loop over lags 1, 2, 3
local lags 1 2 3
foreach lag of local lags {
    di as txt "—— 2SLS with `lag'-month lag ——"
    
    ivregress 2sls ExcessMortality ///
        (L`lag'_lockdown = early_adopter lockdown_volatility) ///
        PopulationVaccinated ConfirmedCases H4_Emergency_investment_in_healt ///
        StringencyIndex_Average E1_Income_support, first

    * Capture first-stage F-stat and p-value
    quietly estat firststage
    di as res "  Lag `lag': first-stage partial-F = " %6.2f r(F) ///
         ", p-value = " %5.3f r(p)

    * Capture over-identification p-value
    quietly estat overid
    di as res "  Lag `lag': Sargan χ² p = " %5.3f r(p_sargan) _n
}
*--------------------------------------------------------*
* IV Regression: Robust Diagnostics                      *
*--------------------------------------------------------*
foreach lag in 1 2 3 {
    di as txt "—— IV Regression with `lag'-month Lockdown Lag ——"
    
    ivregress 2sls ExcessMortality ///
        (L`lag'_lockdown = early_adopter lockdown_volatility) ///
        PopulationVaccinated ConfirmedCases ///
        H4_Emergency_investment_in_healt ///
        StringencyIndex_Average E1_Income_support, ///
        first vce(robust)

    estat firststage, all
    estat overid
}
*--------------------------------------------------------*
* Instrument Relevance and Endogeneity Tests             *
*--------------------------------------------------------*
reg L2_lockdown early_adopter lockdown_volatility ///
    PopulationVaccinated ConfirmedCases ///
    H4_Emergency_investment_in_healt ///
    StringencyIndex E1_Income_support
test early_adopter lockdown_volatility

ivregress 2sls ExcessMortality ///
    (L2_lockdown = early_adopter lockdown_volatility) ///
    PopulationVaccinated ConfirmedCases ///
    H4_Emergency_investment_in_healt ///
    StringencyIndex E1_Income_support

estat firststage, all 
estat overid
estat endogenous
estimates store iv_2month_lag

*-----------------------------------------------------------*
* Advanced Analysis of Lockdown Effects on Excess Mortality *
*-----------------------------------------------------------*

use "final_analysis_data_with_IV.dta", clear
xtset Country_ID date


* Generate fine-grained lags for better temporal resolution
forvalues i = 1/6 {
    gen L`i'_lockdown = L`i'.Lockdown_Level
}

* Create average lockdown intensity over different time periods
gen lockdown_avg_1to3 = (L1_lockdown + L2_lockdown + L3_lockdown)/3
gen lockdown_avg_4to6 = (L4_lockdown + L5_lockdown + L6_lockdown)/3

* Estimate immediate vs delayed effects
xtreg ExcessMortality lockdown_avg_1to3 lockdown_avg_4to6 ConfirmedCases PopulationVaccinated ///
	   H4_Emergency_investment_in_healt ///
	  StringencyIndex E1_Income_support, fe
	  
* More sophisticated lag selection using information criteria
quietly forvalues i = 1/6 {
    xtreg ExcessMortality L`i'_lockdown ConfirmedCases PopulationVaccinated H4_Emergency_investment_in_healt StringencyIndex E1_Income_support, fe
    estimates store lag`i'
}
estimates stats lag1 lag2 lag3 lag4 lag5 lag6

* Add robust standard errors to optimal 6-month lag model
xtreg ExcessMortality L6_lockdown ConfirmedCases PopulationVaccinated ///
      H4_Emergency_investment_in_healt ///
	  StringencyIndex E1_Income_support, fe vce(robust)

* Test cumulative effects with distributed lag model
xtreg ExcessMortality L1_lockdown L2_lockdown L3_lockdown L4_lockdown L5_lockdown L6_lockdown ConfirmedCases PopulationVaccinated H4_Emergency_investment_in_healt ///
	  StringencyIndex E1_Income_support, fe vce(robust)

* Calculate and test the cumulative effect
lincom L1_lockdown + L2_lockdown + L3_lockdown + L4_lockdown + L5_lockdown + L6_lockdown

*--------------------------------------------------------*
* Non-linear&Interaction Effects  Effects of Lockdown.   *
*--------------------------------------------------------*

* Quadratic term
gen L6_lockdown_sq = L6_lockdown^2
xtreg ExcessMortality L6_lockdown L6_lockdown_sq ConfirmedCases PopulationVaccinated, fe

* Threshold effect
gen high_stringency = L6_lockdown >= 2 if !missing(L6_lockdown)
xtreg ExcessMortality high_stringency ConfirmedCases PopulationVaccinated, fe


* Calculate optimal lockdown level from quadratic model
xtreg ExcessMortality L6_lockdown L6_lockdown_sq ConfirmedCases PopulationVaccinated, fe
* Calculate the turning point (optimal lockdown level)
nlcom -(_b[L6_lockdown]/(2*_b[L6_lockdown_sq]))

* Generate predicted values for a range of lockdown levels to visualize the relationship
predict xtreg_residuals, e
gen pred_exmort = ExcessMortality - xtreg_residuals
gen lockdown_effect = _b[L6_lockdown]*L6_lockdown + _b[L6_lockdown_sq]*L6_lockdown_sq

* Graph the relationship between lockdown and predicted excess mortality
* First create a cleaned dataset with expected values at each lockdown level
preserve
clear
set obs 4
gen lockdown = _n-1
gen lockdown_sq = lockdown^2
gen effect = 7051.533*lockdown - 3503.236*lockdown_sq
twoway (connected effect lockdown, sort), ytitle("Predicted Effect on Excess Mortality") xtitle("Lockdown Stringency Level") title("Non-linear Effect of Lockdown on Excess Mortality") xscale(range(0 3)) xlabel(0(1)3) name(quadratic_effect, replace)
restore

*--------------------------------------------------------*
* Heterogeneous Effects by Vaccination, Case Burden      *
*--------------------------------------------------------*
* Test if lockdown effectiveness depends on vaccination rates
gen L6_lockdown_X_vacc = L6_lockdown * PopulationVaccinated
xtreg ExcessMortality L6_lockdown L6_lockdown_sq PopulationVaccinated L6_lockdown_X_vacc ConfirmedCases, fe

* Test if effectiveness varies with case burden
gen L6_lockdown_X_cases = L6_lockdown * ConfirmedCases
xtreg ExcessMortality L6_lockdown L6_lockdown_sq ConfirmedCases L6_lockdown_X_cases PopulationVaccinated, fe

* Test if countries with early adoption vs. late adoption show different patterns
* First merge back the early_adopter variable if needed
xtreg ExcessMortality L6_lockdown L6_lockdown_sq ConfirmedCases PopulationVaccinated if early_adopter==1, fe
estimates store early
xtreg ExcessMortality L6_lockdown L6_lockdown_sq ConfirmedCases PopulationVaccinated if early_adopter==0, fe
estimates store late
estimates table early late, b(%9.3f) p(%9.3f) stats(N r2_w)

*--------------------------------------------------------*
* Heterogeneity by Policy Support & Health Investment    *
*--------------------------------------------------------*
* Standardize continuous variables to make interaction coefficients more interpretable
egen std_income_support = std(E1_Income_support)
egen std_health_invest = std(H4_Emergency_investment_in_healt)
egen std_lockdown = std(L6_lockdown)

* Create interaction terms
gen lockdown_X_income = std_lockdown * std_income_support
gen lockdown_X_health = std_lockdown * std_health_invest

* Run models with interactions
xtreg ExcessMortality std_lockdown lockdown_X_income std_income_support ConfirmedCases PopulationVaccinated, fe
estimates store m_income

xtreg ExcessMortality std_lockdown lockdown_X_health std_health_invest ConfirmedCases PopulationVaccinated, fe
estimates store m_health

* Categorical analysis based on above/below median policy support
sum E1_Income_support, detail
gen high_income_support = E1_Income_support > r(p50) if !missing(E1_Income_support)
sum H4_Emergency_investment_in_healt, detail
gen high_health_invest = H4_Emergency_investment_in_healt > r(p50) if !missing(H4_Emergency_investment_in_healt)

* Run quadratic models by income support groups
xtreg ExcessMortality L6_lockdown L6_lockdown_sq ConfirmedCases PopulationVaccinated if high_income_support==1, fe
estimates store high_inc
nlcom -(_b[L6_lockdown]/(2*_b[L6_lockdown_sq]))

xtreg ExcessMortality L6_lockdown L6_lockdown_sq ConfirmedCases PopulationVaccinated if high_income_support==0, fe
estimates store low_inc
nlcom -(_b[L6_lockdown]/(2*_b[L6_lockdown_sq]))

* Compare coefficients across models
estimates table high_inc low_inc, b(%9.1f) p(%9.3f) stats(N r2_w)


sum L6_lockdown, meanonly
gen L6_lockdown_c = L6_lockdown - r(mean)
gen L6_lockdown_c_sq = L6_lockdown_c^2

reg ExcessMortality L6_lockdown_c L6_lockdown_c_sq PopulationVaccinated ConfirmedCases i.Country_ID
vif

*-------------------------------*
* visualisation                 *
*-------------------------------

* Create prediction for range of values
preserve
clear
set obs 100
gen Lockdown_Level = _n/33
gen Lockdown_Level_sq = Lockdown_Level^2
gen pred = 7051.533*Lockdown_Level - 3503.236*Lockdown_Level_sq
twoway line pred Lockdown_Level, ///
    xtitle("Lockdown Level") ytitle("Predicted Excess Mortality") ///
    title("Non-linear Effect of Lockdown on Mortality")
restore


