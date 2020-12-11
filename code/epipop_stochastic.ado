program define epipop_stochastic

    version 16.0
	
	epipop nothing
	
	local agegroups "0,10,20,30,40,50,60,130"
	syntax , agevar(string) ///
	         sexvar(string) malecode(integer) femalecode(integer) ///
			 r0(real) [theta(real 0.0) c1(real 0.0) c2(real 0.0) c3(real 1.0)] ///
			 [repeat(integer 100)] [report(string)]

	if (`repeat'<0) {
	  display as error "Number of repetitions must be positive."
	  error 101
	}
	
	fixages , agevar(`agevar') agegroups("`agegroups'") ///
	      sexvar(`sexvar') malecode(`malecode') femalecode(`femalecode') 

	local N=r(N)
	matrix c = r(c)
	matrix F=r(F)

    simulate , n(`N') r0(`r0') theta(`theta') ///
           c1(`c1') c2(`c2') c3(`c3') cname("c") fname("F") ///
		   repeat(`repeat') report(`"`report'"')
end


program stochpois, rclass

    version 16.0
	
	local verbose=0
	
    syntax , n(integer) r0(real) theta(real) ///
	         c1(real) c2(real) c3(real) cname(string) ///
	         frame(string) bedframec(string) bedframeh(string) iteration(integer)
			 
	local steps_day=24
			 
    // Todo: add check that the sum of the frequencies must be 100

    capture assert inrange(`c1', 0.00, 100.00)
	if _rc {
	    display as error "Error! Parameter C1 was expected in percent [0%..100%]"
		error 101
	}
    capture assert inrange(`c2', 0.00, 100.00)
	if _rc {
	    display as error "Error! Parameter C2 was expected in percent [0%..100%]"
		error 101
	}	
	
	// Adjust R0 to reflect the effects C1 and C2:
    local r0=`r0'*(1-`c1'/100)*(1-`c2'/100)
	tempname ccc
    mata `ccc'=st_matrix("`cname'")	
	quietly mata trash=stochPois(`n', `ccc', `r0', `theta', `c3')
	
	local vars "t S Z I R RT H C D"
	matrix colnames RESULTS = `vars'
	quietly svmat RESULTS, names(col)
	
	// now need to apply the variable labels for all these variables.
	foreach v in `vars' {
	    mata st_local("vl", epi_siz_compname("`v'"))
		label variable `v' `"`vl'"'
	}
	
		
	summarize H, meanonly
	local beddaysH=ceil(r(sum)/`steps_day')
	summarize C, meanonly
	local beddaysC=ceil(r(sum)/`steps_day')
		
	if ((`beddaysH'==0) & (`beddaysC'==0)) {
	  return local pass=1
	  exit
	}
	else {
	  return local pass=0
	}
	if `verbose' {
		display "{text}Time to first death: {result:`=string(`tFirstDeath')'}" , 
		display "Time to first hospitalization: {result:`=string(`tFirstHosp')'}"
		display "Total hospitalized: {result:`=string(`TotHosp')'}"
		display "Total ICU: {result:`=string(`TotICU')'}" 
		display "Total deaths: {result:`=string(`TotDeaths')'}"

		display "Number of hospital bed-days: {result:`beddaysH'}"
		display "Number of ICU bed-days: {result:`beddaysC'}"
	}
	local N=S[1]+Z[1]+I[1]+R[1]+RT[1]+H[1]+C[1]+D[1]
	assert (`N'==`n')
	generate double infected=`N'-S
	local Tinfected=infected[_N]
	generate pctinf=infected/`Tinfected'
	summarize t if pctinf>0.25, meanonly
	local p25=floor(r(min))
	summarize t if pctinf>0.50, meanonly
	local p50=floor(r(min))
	summarize t if pctinf>0.75, meanonly
	local p75=floor(r(min))
	
	summarize H, meanonly
	local hmax=r(max)
	summarize t if H>=`hmax', meanonly
	local thmax=r(min)
	
	summarize C, meanonly
	local cmax=r(max)
	summarize t if C>=`cmax', meanonly
	local tcmax=r(min)
	if `verbose' {
		display "{text}Time to occurrence of"
		display "{text}25% of infections: {result:`p25'} days"
		display "50% of infections: {result:`p50'} days"
		display "75% of infections: {result:`p75'} days"
	}
	// Calculate how many beds cover 0,1,2,...100% of bed-days for patients
	// in hospitals and ICUs

	foreach ltr in c h {
	    local compartment=strupper("`ltr'")
		preserve
		    generate int day = floor(t)
			collapse (max) C = `compartment' , by(day)  // number of beds in a day is the max during that day across all the hours
			summarize C, meanonly
			local CBEDDAYS = r(sum)
			local maxbeds=r(max)
			
			tempname BEDS
			
			quietly frame create `BEDS'
			frame `BEDS' {
				generate beds=.
				generate covered=.
			}
			
			forval beds=0/`maxbeds' {
				generate ttt=cond(C>`beds',C-`beds',0)
				quietly summarize ttt
				local bedsx=r(sum)
				local coveredpcnt=100*(1-`bedsx'/`CBEDDAYS')
				frame post `BEDS' (`beds') (`coveredpcnt')
				drop ttt
			}
			forval p=0/100 {
				frame post `BEDS' (.) (`p') // placeholder for interpolation
			}
			frame `BEDS' {
			  ipolate beds covered, generate(ibeds)
			  quietly replace ibeds=ceil(ibeds)
			  quietly keep if missing(beds)
			  drop beds
			  rename ibeds beds

			  assert _N==101
			  forval p=1/101 {
				  frame post `bedframe`ltr'' (`iteration') (`=covered[`p']') (`=beds[`p']')
			  }
			}
			
		restore
	}
	
	quietly {
		pwf
		local cframe `"`r(currentframe)'"'
		frame change `frame'
		set obs `=_N+1'
		replace p25=`p25' in L
		replace p50=`p50' in L
		replace p75=`p75' in L
		replace bedhosp=`beddaysH' in L
		replace bedicu=`beddaysC' in L
		replace totinfected=`Tinfected' in L
		replace percpopinf=`Tinfected'/`N'*100.00 in L
		replace tothosp=`TotHosp' in L
		replace toticu=`TotICU' in L		
		replace totdeath=`TotDeaths' in L
		replace t1death=`tFirstDeath' in L
		replace t1hosp=`tFirstHosp' in L
		replace maxhosp=`hmax' in L
		replace tmaxhosp=`thmax' in L
		replace maxicu=`cmax' in L
		replace tmaxicu=`tcmax' in L
		frame change `cframe'
	}
end

program create_empty
	
	version 16.0
	
	syntax , frame(string) bedframeh(string) bedframec(string)
	
	quietly {
		frame create `frame'
		quietly pwf
		local cframe `"`r(currentframe)'"'
		frame change `frame'
			generate double p25=.
			generate double p50=.
			generate double p75=.
			generate long bedhosp=.
			generate long bedicu=.
			generate long totinfected=.
			generate double percpopinf=.
			generate long tothosp=.
			generate long toticu=.
			generate long totdeath=.
			generate int t1death=.
			generate int t1hosp=.
			generate long maxhosp=.
			generate int tmaxhosp=.
			generate long maxicu=.
			generate int tmaxicu=.
			format * %10.0f
			
			label variable p25         "Time to 25% of infections, days"
			label variable p50         "Time to 50% of infections, days"
			label variable p75         "Time to 75% of infections, days"
			label variable bedhosp     "Number of hospital bed-days, bed-days"
			label variable bedicu      "Number of ICU bed-days, bed-days"
			label variable totinfected "Total infected, persons"
			label variable percpopinf  "Percent of population infected"
			label variable tothosp     "Total hospitalized, persons"
			label variable totdeath    "Total dead, persons"
			label variable toticu      "Total ICU, persons"
			label variable t1death     "Time to 1st death, days"
			label variable t1hosp      "Time to 1st hospitalization, days"
			label variable tmaxhosp    "Time to max hospitalization, days"
			label variable tmaxicu     "Time to max ICUs, days"
			label variable maxhosp     "Maximum hospitalization, persons"
			label variable maxicu      "Maximum ICUs, persons"
			
			quietly frame create `bedframeh'
			frame change `bedframeh'
			generate int iteration=.
			generate int percent=.
			generate int beds=.
			label variable iteration "Simulation iteration"
			label variable percent   "Percent of bed-days covered"
			label variable beds      "Beds for the patients in hospitals"

			quietly frame create `bedframec'
			frame change `bedframec'
			generate int iteration=.
			generate int percent=.
			generate int beds=.
			label variable iteration "Simulation iteration"
			label variable percent   "Percent of bed-days covered"			
			label variable beds      "Beds for the patients in ICUs"
			
		frame change `cframe'
	}
	
end		 
		 
program define frappend
	version 16.0
	syntax , to(string)
	
	quietly count
	if (r(N)==0) {
	    display in yellow "Nothing to append. Exiting."
		exit
	}
	
	capture confirm frame `to'
	if (_rc) {
	  display in green "Frame `to' created."
	  quietly frame create `to'
	}

	capture {
	    tempfile tmp
	    save `"`tmp'"'
		frame `to' : append using `"`tmp'"'
	}
end		 
		 
program define simulate

    version 16.0
	
	syntax , n(integer) cname(string) fname(string) ///
	         r0(real) theta(real) ///
	         c1(real) c2(real) c3(real) ///
			 [repeat(integer 100)] [report(string)] 
			 
	if (`repeat'<=0) {
	  display as error "Number of repetitions must be positive."
	  error 101
	}
			 
	tempname bedframec bedframeh frame dataframe
	
    create_empty, frame(`frame') bedframec(`bedframec') bedframeh(`bedframeh') 
	
	tempfile populog
	set logtype text
	quietly log using `"`populog'"', replace nomsg name(populog)
	  tempname poppercent
	  matrix `poppercent'=100*`fname'
	  matrix list `poppercent', noheader format(%8.4f)
	log close populog
/*
    // not reference lethality
	
	tempfile reflethality
	set logtype text
	quietly log using `"`reflethality'"', replace nomsg name(reflethality)
	  matrix list `cname', noheader format(%12.8f)
	log close reflethality
*/
	
	tempfile simparams
	set logtype text
	local c=69
	quietly log using `"`simparams'"', replace nomsg name(simparams)
	    display as text " "
    	display as text "Population size, N" _col(`c') as result `"`=strtrim(string(`n',"%10.0fc"))'"'
		display as text "Basic reproduction number, R0" _col(`c') as result `"`=strtrim(string(`r0',"%8.5f"))'"'
		display as text "Efficiency in contact trace and removal [0,1], theta" _col(`c') as result %8.6f `theta'
	    display as text "Reduction in contact rate by social distance, proportion [0,1], C1" _col(`c') as result %8.6f `c1'
		display as text "Reduction in contact rate by use of mask, proportion [0,1], C2" _col(`c') as result %8.6f `c2'
		display as text "Increase severity for population conditions, C3" _col(`c') as result %8.6f `c3'
		
		display as text "Number of repetitions" _col(`c') as result `"`=strtrim(string(`repeat',"%10.0fc"))'"'
	log close simparams
	
	display "{text}Running stochastic simulation of SIZ-model ({result:`repeat'} repetitions)"
	forval i=1/`repeat' {
		local flag=1
		display "{text} Iteration: {result:`i'}"
		while (`flag') {
			clear
			stochpois, n(`n') r0(`r0') theta(`theta') ///
					   c1(`c1') c2(`c2') c3(`c3') ///
					   cname("`cname'") frame(`frame') ///
					   bedframec(`bedframec') bedframeh(`bedframeh') iteration(`i')
			local flag=r(pass)
		}
		if (`"`dataframe'"'!="") {
		  generate int iteration=`i'
		  quietly frappend, to(`dataframe') // todo : this can be made general purpose
		}
	}
	display "{text}Simulation complete.{break}"
	
    dographs, frame(`frame') cname("`cname'") ///
	          dataframe("`dataframe'") n(`n') ///
			  bedframec("`bedframec'") bedframeh("`bedframeh'") ///
			  populog(`"`populog'"') reflethality(`reflethality') ///
			  simparams(`"`simparams'"') report(`"`report'"')
end

program define newassign, rclass
    version 16.0
	syntax , at(string) cmatrix(string)
	// Option AT() characterizes the desired division
	
	// The following constants describe the dimensions of age in CMATRIX:
	// Row of 0..STEP-1
	// Row of STEP..2*STEP-1
	// .....
	// Last row X..TOP, where X is where the rows exhaust the steps.
	local TOP=120
	local STEP=10
	tempname T
	
	matrix `T'=J(`TOP'+1,1,.)
	assert colsof(`cmatrix')==7
	forval i=1/`=colsof(`cmatrix')' {
	    local low=(`i'-1)*`STEP'
		local high=`i'*`STEP'-1
		if (`i'==colsof(`cmatrix')) local high=`TOP'
		local w = `high' - `low' + 1
		forval j=`low'/`high' {
			matrix `T'[`j'+1,1]=`cmatrix'[1,`i'] / `w'
		}
	}

	local lst = subinstr(`"`at'"',","," ",.)
	local lst = subinstr(`"`lst'"',"  "," ",.)
	
	local nw=`:word count `lst''
	
	tempname Z ZZ
	local rnames ""
	mata `Z'=J(`nw'-1,1,.)
	forvalue i=1/`=`nw'-1' {
	    local low= `:word `i' of `lst''
		local high= `:word `=`i'+1' of `lst''-1
		if (`i'==`nw'-1) local high=`TOP'
		local rnames `"`rnames' "`low'-`high'""'
		mata `Z'[`i',1]=sum(st_matrix("`T'")[`=`low'+1'..`=`high'+1',1])
	}
	mata st_matrix("`ZZ'",`Z')
	matrix rownames `ZZ' = `rnames'
	return matrix Z=`ZZ'
end


program dographs
    version 16.0
	syntax , n(integer)  ///
	[frame(string)] cname(string) dataframe(string) ///
	[bedframec(string) bedframeh(string)] ///
	[populog(string) reflethality(string) simparams(string) report(string)]

	foreach b in c h {
	
	local bf=`"`bedframe`b''"'
    if (`"`bf'"'!="") {
	    frame change `bf'
		local varlbl `"`:variable label beds'"'
		
		tempvar high2 high mean sd low low2
		egen double `mean'=mean(beds), by(percent)
		egen double `sd'=sd(beds), by(percent)
		generate double `high'=`mean'+`sd'
		generate double `high2'=`mean'+2*`sd'
		generate double `low'=`mean'-`sd'
		generate double `low2'=`mean'-2*`sd'

		preserve
		contract `high2' `high' `mean' `low' `low2' percent
		
		sort percent
		local gname `"Beds_`=strupper("`b'")'"'
		capture graph drop `gname'
		twoway (area `high2' `high' `mean' `low' `low2' percent,  ///
				  fc(gs9 gs3 gs3 gs9 white) lc(white%0 white%0 white%0 white%0 white%0) ///
				  title(`"`varlbl'"') ytitle("Beds required") scale(0.75) ///
				  legend(cols(2) order( 2 1 ) label(6 "`varlbl'") ///
				  label(2 "68%-Confidence area (mean ± 1 S.D.)") ///
				  label(1 "95%-Confidence area (mean ± 2 S.D.)") size(small)))  ///
				  (line `mean' percent, lc(red) name(`"`gname'"')) 
				  
		restore
		drop `high2' `high' `mean' `low' `low2'
	}
	}

	local varlabels `""Time""'
	local varlabels `"`varlabels' "(S) Susceptible""'
	local varlabels `"`varlabels' "(Z) Asympt. inf. that will not require hospitalization""'
	local varlabels `"`varlabels' "(I) Sympt. inf. that will need hospitalization eventually""'
	local varlabels `"`varlabels' "(H) Hospitalized not in intensive care""'
	local varlabels `"`varlabels' "(R) Recovered""'
	local varlabels `"`varlabels' "(RT) Removed temporarily""'
	local varlabels `"`varlabels' "(C) Intensive care""'	
	local varlabels `"`varlabels' "(D) Dead""' 
	
	local mcolnames "t S Z I H R RT C D"
	
	tempfile resultslog

	log using `"`resultslog'"', text replace name(outlog) nomsg // ##############
			
	if (`"`dataframe'"'!="") {
		frame change `dataframe'
		
		generate IZ=I+Z
		preserve
		  tempname mean
		  egen `mean' = mean(IZ), by(t)
		  contract `mean' t
		  sort t, stable
		  quietly summarize `mean'
		  quietly summarize t if `mean'>=r(max)
		  local peak = r(min)
		  //noisily display " {break}{text}Peak of infections (I+Z) at t=" as result %6.2f `peak'
		restore
			
		preserve
		  quietly {
			  generate DAY=floor(t)
			  sort DAY iteration, stable
			  by DAY iteration: generate ismax=(_n==_N) 
			  by DAY iteration: keep if ismax==1  // leave one obs per day
			  replace t=DAY
			  drop DAY
			  summarize t, meanonly
			  local mx=r(max)		  
			  tempvar wt
			  quietly by iteration (t), sort: generate long `wt' = cond(_n == _N, `mx'+1-t, 1)
			  quietly expand `wt'
			  drop `wt'
			  quietly by iteration t, sort: replace t = t + _n - 1
			  
			  collapse S Z I H R RT C D IZ ///
					  (sd) SD_S=S (sd)SD_Z=Z (sd) SD_I=I (sd) SD_H=H (sd) SD_R=R ///
					  (sd) SD_RT=RT (sd) SD_C=C (sd) SD_D=D (sd) SD_IZ=IZ, by(t)
			  
			  summarize t, meanonly
			  local days=round(r(max))
		  }
		  epimodels_util ditable t, days(`days') datefmt("%dCY-N-D") ///
		  mcolnames("`mcolnames'") ivar("IZ") varlabels(`"`varlabels'"') ///
		  ylabel("Population") modeltitle("SIZ-Model") digits(0) stdev comma("c")
		restore
		
		foreach z in `mcolnames' {
		    if (`"`z'"'=="t") continue
		    local varlbl `"`: variable label `z''"'
			preserve
			    quietly replace t=round(t*24) // special knowledge			
			    summarize t, meanonly
				local mx=r(max)
				tempvar wt
				quietly by iteration (t), sort: generate long `wt' = cond(_n == _N, `mx'+1-t, 1)
				quietly expand `wt'
				drop `wt'
				quietly by iteration t, sort: replace t = t + _n - 1
				quietly replace t=t/24 // special knowledge
				
			    tempname mean sd low high low2 high2
				quietly {
					egen `mean'=mean(`z'), by(t)
					label variable `mean' `"[`: variable label `z'']"'
					egen `sd'=sd(`z'), by(t)
					generate double `low2' =max( 0, `mean'-`sd'*2)
					generate double `low'  =max( 0, `mean'-`sd')
					generate double `high' =min(`n',`mean'+`sd')
					generate double `high2'=min(`n',`mean'+`sd'*2)
					contract `mean' `low' `high' `low2' `high2' t
				}
				sort t, stable
				generate float peak=`peak'
				
				capture graph drop "Group_`z'"

				twoway (area `high2' `high' `mean' `low' `low2' t,  ///
				  fc(gs9 gs3 gs3 gs9 white) ///
				  lc(white%0 white%0 white%0 white%0 white%0) ///
				  title(`"`varlbl'"') scale(0.75) ///
				  legend( ///
				      cols(2) order( 2 1 ) ///
					  label(6 "`varlbl'") ///
				      label(2 "68%-Confidence area (mean ± 1 S.D.)") ///
				      label(1 "95%-Confidence area (mean ± 2 S.D.)") ///
					  size(small)))  ///
				  (line `mean' t, lc(red) name("Group_`z'")) ///
				  (dropline `high2' `low2' peak, mcolor(none none) ///
				     lpattern("--" "--") lcolor(green green) ///
					 lwidth(vthin vthin) xline(`peak', lpattern("--") ///
					 lcolor(green) lwidth(vthin))) 
			restore
		}
		

	}
	
	
	
	// --------------------------------------------
	
	
	if (`"`frame'"'!="") frame change `frame'
	
	
	quietly summarize totdeath
	local death_mean=r(mean)
	local death_sd=r(sd)
	quietly summarize tothosp
	local hosp_mean=r(mean)
	local hosp_sd=r(sd)
	quietly summarize toticu
	local icu_mean=r(mean)
	local icu_sd=r(sd)	
	
	tempname c2
	
	mata st_matrix("`c2'", st_matrix("`cname'")/sum(st_matrix("`cname'")))
	// matrix list `c2'
	
	newassign , at(0, 5, 12, 18, 60, 120) cmatrix("`c2'")
	matrix `c2'=r(Z)'
	
	tempname EX
	matrix `EX' = J(`=colsof(`c2')',3,.)

	forval i=1/`=colsof(`c2')' {
	  matrix `EX'[`i',1] = round(`death_mean'*`c2'[1,`i'])
	  matrix `EX'[`i',2] = round(`hosp_mean'*`c2'[1,`i'])
	  matrix `EX'[`i',3] = round(`icu_mean'*`c2'[1,`i'])
	}
	
	matrix colnames `EX' = "Deaths" "Hosp" "ICU"
	//matrix rownames `EX' = "0-9" "10-19" "20-29" "30-39" "40-49" "50-59" "60-.."
	matrix rownames `EX' = `: colnames `c2''
	tempname EXT
	mata st_matrix("`EXT'", colsum(st_matrix("`EX'"))) // this is totals
	matrix rownames `EXT' = "Total"
	matrix `EX'=`EX' \ `EXT'
	
	display as text ""
	display as text ""
	display as text "BURDEN OF EPIDEMICS BY AGE GROUP (PERSONS)"
	display as text "------------------------------------------"
	matrix list `EX', noheader
	display as text ""

	
	local fmt "%20.0fc"
	
	display as text "    TIME TO OCCURRENCE OF:"
	display as text "------------------------------------"
    foreach p in 25 50 75 {
	  quietly summarize p`p'
	  display `"{text}`p'% of infections: {result:`=string(`r(mean)',"`fmt'")' (`=string(`r(sd)',"`fmt'")')} days"'
	}
	display as text "------------------------------------{break}"
	
	
	display as text "    OTHER STATISTICS"
	display as text "----------------------------------------------------------"
	quietly summarize totinfected
	display `"{text}Total number of infections: {col 38}{result:`=string(`r(mean)',"`fmt'")' (`=string(`r(sd)',"`fmt'")')}"'
	quietly summarize percpopinf
	display `"{text}Percent of population infected: {col 38}{result:`=string(`r(mean)',"%20.2f")' (`=string(`r(sd)',"%20.2f")')}"'
	
	// max num hosp
	quietly summarize maxhosp 
	local r1=r(mean)
	local r2=r(sd)
	quietly summarize tmaxhosp
	local r3=r(mean)
	local r4=r(sd)
	display `"{text}Max number of hospitalized in a day: {col 38}{result:`=string(`r1',"`fmt'")' (`=string(`r2',"`fmt'")')} by day {result:`=string(`r3',"`fmt'")'}"'
	
	// max num ICU
	quietly summarize maxicu
	local r1=r(mean)
	local r2=r(sd)
	quietly summarize tmaxicu
	local r3=r(mean)
	local r4=r(sd)
	display `"{text}Max number of ICUs in a day: {col 38}{result:`=string(`r1',"`fmt'")' (`=string(`r2',"`fmt'")')} by day {result:`=string(`r3',"`fmt'")'}"'

	
	quietly summarize t1hosp
	display `"{text}Time to first hospitalization: {col 38}{result:`=string(`r(mean)',"`fmt'")' (`=string(`r(sd)',"`fmt'")')}"'
	quietly summarize t1death
	display `"{text}Time to first death: {col 38}{result:`=string(`r(mean)',"`fmt'")' (`=string(`r(sd)',"`fmt'")')}"'
	
	display "Number of bed-days required:"
	quietly summarize bedhosp
	display `"{text}      - hospitalized: {col 38}{result:`=string(`r(mean)',"`fmt'")' (`=string(`r(sd)',"`fmt'")')}"'
	quietly summarize bedicu
	display `"{text}      - ICUs: {col 38}{result:`=string(`r(mean)',"`fmt'")' (`=string(`r(sd)',"`fmt'")')}"'
	log close outlog
	
	
    quietly {
	    capture graph drop totinfected
		kdensity totinfected, name(totinfected)
		capture graph drop t1death
		kdensity t1death, name(t1death)
		capture graph drop toticu
		kdensity toticu, name(toticu)
	}
	
	popreport /*t S Z I H R RT C D*/, ///
		 modelname("siz") modelparams("`modelparams'") ///
		 appendixgraphs(`appendixfiles') ///
		 results(`"`resultslog'"') reflethality(`reflethality') ///
		 popstruct(`"`populog'"') simparams(`"`simparams'"') ///
		 save(`report')
	
    clear	
end



program define popreport
    version 16.0
	
	local ffontface "Consolas"
	local fontface "Helvetica"
	local fontcolor "steelblue"
	local tfont `""`fontface'",28,`fontcolor'"'
	local sfont `""`fontface'",12,`fontcolor'"'
	local imagewidth=1920

	syntax [varlist(default=none)], modelname(string) [modelparams(string)] ///
	                [modelgraph(string)] [appendixgraphs(string)] ///
					[results(string)] [popstruct(string)] ///
					[simparams(string)] [reflethality(string)] ///
					save(string)

	if (strlen(`"`modelparams'"')>40) local br `"`=char(10)'"'
	
	putpdf begin
	putpdf paragraph
	putpdf text ("EPIPOP Population Report"), bold font(`tfont')
	putpdf paragraph
	putpdf text ("This report was generated on `c(current_date)'"), font(`sfont')
	local s `"`modelname' model."'
	if (`"`modelparams'"'!="") local s `"`modelname' model with parameters:`br'`modelparams'."'
	putpdf paragraph
	putpdf text (`"`s'"'), bold
	
	epipop pdfreport putstaticimage "epimodels_sch_`modelname'.png" // model scheme
	epipop pdfreport putstaticimage "epimodels_eq_`modelname'.png", width(4.0) // model equations

	if (`"`modelgraph'"'!="") {
		putpdf paragraph
		putpdf image `"`modelgraph'"'
	}

	putpdf paragraph
	putpdf text ("Generated with"), font(`sfont')
	mata epipop_about()
	putpdf paragraph
	putpdf text ("EPIPOP version `epipop_version' from `compile_date' built for Stata v`compile_version'"), linebreak(1)
	putpdf text ("For more information, visit EPIPOP' homepage: http://www.radyakin.org/stata/epipop/")
	mata epimodels_about()
	putpdf paragraph
	putpdf text ("EPIMODELS version `epimodels_version' from `compile_date' built for Stata v`compile_version'"), linebreak(1)
	putpdf text ("For more information, visit EPIMODELS' homepage: http://www.radyakin.org/stata/epimodels/")
	
	putpdf paragraph
	putpdf text ("Details:"), font(`sfont')
	putpdf paragraph
	putpdf text ("Carlos Hernandez-Suarez, Paolo Verme, Sergiy Radyakin, Efrén Murillo-Zamora. 2020"), linebreak(1)
	putpdf text ("COVID-19 Outbreaks in Refugee Camps. A Simulation study"), linebreak(1)
	putpdf text ("https://www.medrxiv.org/content/10.1101/2020.10.02.20204818v1.full.pdf")
	
	if (`"`popstruct'"'!="") | (`"`simparams'"'!="") | (`"`reflethality'"'!=""){
		putpdf pagebreak		
		epipop pdfreport putoptionalparagraph `"`popstruct'"', ///
		  title("Population Structure, %") font(`tfont')
		epipop pdfreport putoptionalparagraph `"`reflethality'"', ///
		  title("Reference Lethality Probabilities") font(`tfont')
		epipop pdfreport putoptionalparagraph `"`simparams'"', ///
		  title("Simulation parameters") font(`tfont')
	}

	putpdf pagebreak
	epipop pdfreport putoptionalparagraph `"`results'"', ///
		  title("Simulation results") font(`tfont')

	if (`"`varlist'"'!="") {
		putpdf table t=data(`varlist'), varnames
		putpdf table t(.,.), halign(center) valign(center) nformat(%12.0fc) font("`ffontface'",8)
		putpdf table t(1,.), bgcolor("`fontcolor'")
		putpdf paragraph
	}

	foreach v in `varlist' {
		putpdf text ("`v': `:variable label `v''`=char(10)'")
	}

	// optional appendix	
	local files ""
	foreach f in Group_D Group_I Group_RT Group_Z Group_C Group_H Group_R Group_S Beds_H Beds_C totinfected toticu t1death {
	  tempfile tmp
	  while 1 {
		  capture confirm file `"`tmp'.png"'
		  if _rc continue, break
		  tempfile tmp
	  }
	  local tname `"`tmp'.png"'
	  quietly graph export "`tname'", name(`f') as(png) width(`imagewidth')
	  local files `"`files' "`tname'""'
	}
	epipop pdfreport putallgraphfiles, graphs(`files') landscape

	putpdf save `"`save'"' , replace
	foreach f in `files' {
	  erase `f'
	}
end

program define fixages, rclass

    version 16.0
	
	syntax , agevar(string) agegroups(string) ///
	         sexvar(string) malecode(integer) femalecode(integer)

	preserve
	
	epimodels_util popmatrix , agevar(`agevar') sexvar(`sexvar') ///
				   malecode(`malecode') femalecode(`femalecode') ///
				   at(`agegroups')

	tempname F M T C
	matrix `F'= r(F)
	local N=r(N)
	matrix `M'=`F'*`N'

	tempname ATFRAME
	frame create `ATFRAME'
	
	frame `ATFRAME' {
		epipop_siz_at, clear
		egen group=cut(age), at(`agegroups')
		quietly collapse (sum) male (sum) female, by(group)
		format male female %20.10f
		
		matrix `C'=J(1,`=_N',.)
		forval i=1/`=_N' {
		  matrix `C'[1,`i']=male[`i']*`M'[`i',1]/`N' + female[`i']*`M'[`i',2]/`N' 
		}
	}

	return matrix c=`C'
	return scalar N=`N'
	return matrix F=`F'

end

// END OF FILE