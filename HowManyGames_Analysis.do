clear all
//adopath + "F:\Uni work\Thesis\Data"
set more off
import delimited "STATAdata.csv", delimiter(comma) case(preserve) 
sort uid

replace a = a+d1
replace x = x+d1
gen w = d1
replace c = c+d2
replace y = y+d2 
gen z = d2


local ifRestriction =  "Session ~= 0" // set to Session ~= 0 to include everything

xtset uid Period

qui gen A = 1-ActionAB
qui gen C = 1-ActionCD
qui gen invPeriod = 1/Period

qui gen ActComb = A*C + A*(1-C)*2 + (1-A)*C*3 + (1-A)*(1-C)*4

//bysort ugroupid: tab A C, exact cell nofreq

label define LABELActComb 1 "AC" 2 "AD" 3 "BC" 4 "BD"
label values ActComb LABELActComb

label define LABELTreatment 1 "[1] Baseline" 2 "[2] Add 50 to Game 1" 3 "[3] Add 50 to Game 2"
label values Treatment LABELTreatment

label variable invPeriod "Period$^{-1}$"

gen invPeriodXT2 = 0
replace invPeriodXT2 = invPeriod if Treatment==2
label variable invPeriodXT2 "Period$^{-1} \times $ (Treatment==2)"

gen invPeriodXT3 = 0
replace invPeriodXT3 = invPeriod if Treatment==3
label variable invPeriodXT3 "Period$^{-1} \times $ (Treatment==3)"


save STATAdata, replace



local spec1 = "i.Treatment if  Period ==1 "
	local cluster1 = " "
	local PD1   = "N"
	local PR1   = "1 only"
	local CL1   = "none"
	local invP1  = "N"
	local TxinvP1  = "N"
local spec2 = "i.Treatment invPeriod if  Period >=11 "
	local cluster2 = "vce(cluster ugroupid)"
	local PD2   = "N"
	local PR2   = "11-20"
	local CL2   = "Group"
	local invP2  = "N"
	local TxinvP2  = "N"
local spec3 = "i.Treatment if Period==Period"
	local cluster3 = "vce(cluster ugroupid)"
	local PD3   = "N"
	local PR3   = "none"
	local CL3   = "Group"
	local invP3  = "N"
	local TxinvP3  = "N"
local spec4 = "i.Treatment invPeriod if Period==Period "
	local cluster4 = "vce(cluster ugroupid)"
	local PD4   = "N"
	local PR4   = "none"
	local CL4   = "Group"
	local invP4  = "Y"
	local TxinvP4  = "N"
local spec5 = "i.Treatment invPeriod* if Period==Period "
	local cluster5 = "vce(cluster ugroupid)"
	local PD5   = "N"
	local PR5   = "none"
	local CL5   = "Group"
	local invP5  = "Y"
	local TxinvP5  = "Y"
local spec6 = "i.Treatment i.Period if Period==Period "
	local cluster6 = "vce(cluster ugroupid)"
	local PD6   = "Y"
	local PR6   = "none"
	local CL6   = "Group"
	local invP6  = "N"
	local TxinvP6  = "N"

local numspecs = 6

forvalues ss =  2/`numspecs' {
	
	
	//disp "eststo T23Broad`ss':  xi: mlogit ActComb `spec`ss'', `cluster`ss'' base(1)"
	mlogit ActComb `spec`ss'' & Treatment~=1, `cluster`ss'' base(1)
	eststo T23Broad`ss'
			estadd local PD "`PD`ss''"
			estadd local PR "`PR`ss''"
			estadd local CL "`CL`ss''"
			qui test [AD]3.Treatment=0
			qui test [BC]3.Treatment=0, accum
			 test [BD]3.Treatment=0, accum
			//local test23 = `r(chi2)'
			estadd local test23 = round(`r(chi2)',0.1)
			estadd local test23p = round(`r(p)',0.0001)
			estadd local invP = "`invP`ss''"
			estadd local TxinvP = "`TxinvP`ss''"
			estpost margins , dydx(i.Treatment)
			eststo margins_T23Broad`ss'
	
	mlogit ActComb `spec`ss'', `cluster`ss'' base(1)
	eststo Broad`ss'   
		estadd local PD "`PD`ss''"
		estadd local PR "`PR`ss''"
		estadd local CL "`CL`ss''"
		// Test for no difference between Treatment 2 and Treatment 3
		// i.e. H0: Broad Bracketing
			qui test [AD]3.Treatment=[AD]2.Treatment
			qui test [BC]3.Treatment=[BC]2.Treatment, accum
			 test [BD]3.Treatment=[BD]2.Treatment, accum
			//local test23 = `r(chi2)'
			estadd local test23 = round(`r(chi2)',0.1)
			estadd local test23p = round(`r(p)',0.0001)
	
	 
	
	logit A `spec`ss'', `cluster`ss''
	eststo ANarrow`ss'
		estadd local PD "`PD`ss''"
		estadd local PR "`PR`ss''"
		estadd local CL "`CL`ss''"
	
	logit C `spec`ss'', `cluster`ss''
	eststo CNarrow`ss'
		estadd local PD "`PD`ss''"
		estadd local PR "`PR`ss''"
		estadd local CL "`CL`ss''"
		
	logit A `spec`ss'' & Treatment~=3, `cluster`ss''
	eststo T12ANarrow`ss'
		estadd local PD "`PD`ss''"
		estadd local PR "`PR`ss''"
		estadd local CL "`CL`ss''"
		
	logit C `spec`ss'' & Treatment~=3, `cluster`ss''
	eststo T12CNarrow`ss'
		estadd local PD "`PD`ss''"
		estadd local PR "`PR`ss''"
		estadd local CL "`CL`ss''"	
	
	logit A `spec`ss'' & Treatment~=2, `cluster`ss''
	eststo T13ANarrow`ss'
		estadd local PD "`PD`ss''"
		estadd local PR "`PR`ss''"
		estadd local CL "`CL`ss''"
		
	logit C `spec`ss'' & Treatment~=2, `cluster`ss''
	eststo T13CNarrow`ss'
		estadd local PD "`PD`ss''"
		estadd local PR "`PR`ss''"
		estadd local CL "`CL`ss''"	
	
	
	 biprobit A C `spec`ss'', `cluster`ss''
	 eststo ACBiProbit`ss'
		estadd local PD "`PD`ss''"
		estadd local PR "`PR`ss''"
		estadd local CL "`CL`ss''"
		qui test [A]3.Treatment = 0
		qui test [C]2.Treatment = 0, accum
		estadd local test23 = round(`r(chi2)',0.1)
		estadd local test23p = round(`r(p)',0.0001)
		estadd local invP = "`invP`ss''"
		estadd local TxinvP = "`TxinvP`ss''"
		estpost margins, dydx(i.Treatment) predict(pmarg1) predict(pmarg2)
		eststo margins_ACBiProbit`ss'
		
		
}

//label variable _ITreXinvPe_2 "(Treatment == 1) $ \times $ Period$^{-1}$"

local reglists "ANarrow CNarrow Broad T23Broad T12ANarrow T12CNarrow T13ANarrow T13CNarrow ACBiProbit"
local extensions "txt rtf tex"

 foreach ff of local extensions {
 foreach rr of local reglists {
 disp "`rr'"
 // remove omitted variables from relevant regressions 
	 if "`rr'" =="Broad" {
	 local omitlist = " "
	 local testlist1 = "test23 H0 : T2 = T3, $\chi^2(3)$"
	 local testlist2 = "test23p $\quad p$-value"
	 }
	 else if "`rr'" =="T23Broad" {
	 local omitlist = "" //" _ITreatment_3 invPeriodXT3 AC: "
	 local testlist1 = "test23 H0 : T2 = T3, $\chi^2(3)$"
	 local testlist2 = "test23p $\quad p$-value"
	 }
	 else if "`rr'" == "T12ANarrow" {
	 local omitlist = "" // "  _ITreatment_3 invPeriodXT3 "
	 local testlist1 = " "
	 local testlist2 = " "
	 }
	 else if "`rr'" == "T12CNarrow" {
	 local omitlist = "" //" _ITreatment_3 invPeriodXT3 "
	 local testlist1 = " "
	 local testlist2 = " "
	 }
	 else if "`rr'" == "T13ANarrow" {
	 local omitlist = "" // " _ITreatment_2 invPeriodXT2 "
	 local testlist1 = " "
	 local testlist2 = " "
	 }
	 else if "`rr'" == "T13CNarrow" {
	 local omitlist = "" // " _ITreatment_2 invPeriodXT2 "
	 local testlist1 = " "
	 local testlist2 = " "
	 }
	 else if "`rr'" =="ACBiProbit" {
	 local omitlist = " "
	 local testlist1 = "test23 H0 : Narrow, $\chi^2(2)$"
	 local testlist2 = "test23p $\quad p$-value"
	 }
	 else {
	 local omitlist = " "
	 local testlist1 = " "
	 local testlist2 = " "
	 }
	 
	 
	 
	esttab `rr'* using outputs/`rr'Levels.`ff', replace se label nomtitles noomitted nobaselevels drop(*.Period `omitlist') scalars("CL Cluster level" "N_clust Number of clusters" "PR Period restriction" "PD Period dummies" "`testlist1'" "`testlist2'") coeflabels(_cons "Constant") nogaps
	
	if "`rr'" =="T23Broad" {
		esttab margins_`rr'* using outputs/`rr'Margins.`ff', labcol2("$=0$" "$=0$" "$=0$" "$=0$"  , title("H0")) coeflabels(1._predict "AC" 2._predict "AD" 3._predict "BC" 4._predict "BD") replace se noomitted nomtitles scalars("CL Cluster level" "N_clust Number of clusters" "invP Period$^{-1}$" "TxinvP Period$^{-1}\times$ Treatment" "PR Period restriction" "PD Period dummies" "`testlist1'" "`testlist2'") nogaps
	}
	if "`rr'" =="ACBiProbit" {
		esttab margins_`rr'* using outputs/`rr'Margins.`ff', labcol2(" " "$=0$" "$=0$" " "  , title("H0")) label coeflabels(1._predict "A" 2._predict "C") replace se noomitted nomtitles scalars("CL Cluster level" "N_clust Number of clusters" "invP Period$^{-1}$" "TxinvP Period$^{-1}\times$ Treatment" "PR Period restriction" "PD Period dummies" "`testlist1'" "`testlist2'") nogaps
	}
 }
 
}

esttab margins_ACBiProbit* , labcol2(" " "$=0$" "$=0$" " "  , title(" H0 ")) label coeflabels(1._predict "A" 2._predict "C") replace se noomitted nomtitles scalars("CL Cluster level" "N_clust Number of clusters" "invP Period$^{-1}$" "TxinvP Period$^{-1}\times$ Treatment" "PR Period restriction" "PD Period dummies" "`testlist1'" "`testlist2'") nogaps
