program define epipop_deterministic
    version 16.0
	
	epipop nothing
	
    syntax , ///
    	popsize(integer)  /* Population size N                                    */ ///
		popstruct(string) /* Relative frequencies F                               */ ///
		r0(real)          /* Basic reproduction number                            */ ///
		theta(real)       /* Efficiency in contact trace and removal [0,1]        */ ///
	    [c1(real 0.00)    /* %-reduction in contact rate by social distance [0,1] */ ///
		c2(real 0.00)     /* %-reduction in contact rate by use of mask [0,1]     */ ///
		c3(real 1.00)]    /* Increase severity for population conditions (default=1 (no increase)) */ ///
		tmax(integer)     /* Max time of the simulation (in days)                 */ ///
		agpop(string)     /* Population aggregation parameters (row numbers)      */ ///
		[refdata(string)] /* Reference data for mortality adjustment and rest of parameters */ ///
		[report(string)]  /* Optional PDF report name (to be created)             */

	local imgwidth=1920   // pixels

	confirm matrix `popstruct'
	assert inrange(`theta',0,1)
	assert inrange(`c1',   0,1)
	assert inrange(`c2',   0,1)

	if missing(`"`refdata'"') local refdata="IMSS_Mexico"

	tempname AR
	epi_refdata_`refdata'
	// these parameters are peeked by the epi_simRC() from the mata code
	// todo: perhaps load directly to mata??
	scalar mu1=r(mu1)
	scalar mu2=r(mu2)
	scalar mu3=r(mu3)
	scalar mu4=r(mu4)
	scalar q=r(q)
	scalar w=r(w)
	scalar alpha=r(alpha)
	matrix `AR'=r(M)
	local modelparams `"`r(modelparams)'"'

	tempfile reflethality
	set logtype text
	quietly log using `"`reflethality'"', replace nomsg name(reflethality)
	  matrix list `AR', noheader format(%12.8f)
	log close reflethality

	tempfile populog
	set logtype text
	quietly log using `"`populog'"', replace nomsg name(populog)
	  tempname poppercent
	  matrix `poppercent'=100*`popstruct'
	  matrix list `poppercent', noheader format(%8.4f)
	log close populog

	tempfile simparams
	set logtype text
	local c=69
	quietly log using `"`simparams'"', replace nomsg name(simparams)
	    display as text " "
    	display as text "Population size, N" _col(`c') as result `"`=strtrim(string(`popsize',"%10.0fc"))'"'
		display as text "Basic reproduction number, R0" _col(`c') as result `"`=strtrim(string(`r0',"%8.5f"))'"'
		display as text "Efficiency in contact trace and removal [0,1], theta" _col(`c') as result %8.6f `theta'
	    display as text "Reduction in contact rate by social distance, proportion [0,1], C1" _col(`c') as result %8.6f `c1'
		display as text "Reduction in contact rate by use of mask, proportion [0,1], C2" _col(`c') as result %8.6f `c2'
		display as text "Increase severity for population conditions, C3" _col(`c') as result %8.6f `c3'
	log close simparams

    tempfile resultslog
    set logtype text
	quietly log using `"`resultslog'"', replace nomsg name(resultslog)

	mata epi_sim_siz("`popstruct'",`popsize',`r0',`theta',`c1',`c2',`c3',`tmax',"`agpop'")

	// ------------ BEDS ----------------------

	summarize H, meanonly
	local hosp_bed_days=string(r(sum), "%10.0fc")
	summarize C, meanonly
	local icu_bed_days=string(r(sum), "%10.0fc")

	local beds_files=""
	foreach v in H C {
		summarize `v', meanonly
		local total=r(sum)
		local maxV=r(max)

		tempname ftmp
		frame create `ftmp' int(beds) double(days)

		tempvar ttt
		forval b=0/`maxV' {
			generate `ttt'=min(`v',`b')
			summarize `ttt', meanonly
			local under=r(sum)
			frame post `ftmp' (`b') (`=100*`under'/`total'')
			drop `ttt'
		}
		frame change `ftmp'
		local title = ""
		if (`"`v'"'=="H") local title="Hospital"
		if (`"`v'"'=="C") local title="ICU"

		capture graph drop bed_days_`v'
		graph twoway line days beds, lc(maroon) ///
		       ytitle("% of bed-days fulfilled")  ///
			   xtitle("number of beds") ///
			   title("`title'") ///
		       graphregion(fc(white)) ///
			   name("bed_days_`v'")
		tempfile tmp`v'
		quietly graph export `"`tmp`v''.png"', as(png) width(`imgwidth') replace
		local beds_files=`"`beds_files' "`tmp`v''.png""'

		frame change default // todo: must return to the previous frame which is not necessarily 'default'
		frame drop `ftmp'
	}
	// ------------------------------- END OF BEDS ----------------

	matrix rownames B = `burdenrownames'
	matrix colnames B = `burdencolnames'

    display "{text} BURDEN OF EPIDEMICS BY AGE GROUP"
	display "{text}            (PERSONS)"
	display "{text}-----------------------------------" _continue
	matrix list B, noheader
	display "{text}-----------------------------------"

	display "{text}"
	display "    TIME TO OCCURRENCE OF "
	display "------------------------------------------------------------"
	display " 25% of infections {col 45}{result:`=T[1,1]'} days"
	display " 50% of infections {col 45}{result:`=T[1,2]'} days"
	display " 75% of infections {col 45}{result:`=T[1,3]'} days"
	display "------------------------------------------------------------"

	display "{text}"
	display "    OTHER STATISTICS"
	display "------------------------------------------------------------"
	display `" Total number of infections: {col 45}"' ///
	           `"{result:`=string(tni,"%10.0fc")'} persons"'
	display `" Percent of population infected: {col 45}"' ///
	           `"{result:`=string(100*tni/`popsize',"%6.2f")+"%" '}"'
	display `" Max. number of hospitalized in a day: {col 45}"' ///
	           `"{result:`=string(mhosp,"%10.0fc")'} at day {result:`=tmaxH'}"'
	display `" Max. number of ICUs in a day: {col 45}"' ///
	           `"{result:`=string(micu,"%10.0fc")'} at day {result:`=tmaxC'}"'
	display `" Time to first hospitalization: {col 45}"' ///
	           `"{result:`=thosp1'} days"'
	display `" Time to first death: {col 45}{result:`=tdeath1'} days"'

	display `" Number of bed-days required:"'
	display `" - hospitalized: {col 45}{result:`hosp_bed_days'} bed-days"'
	display `" - ICUs: {col 45}{result:`icu_bed_days'} bed-days"'
	display "------------------------------------------------------------"

	log close resultslog
	confirm file `"`resultslog'"'

	local st `""S Z I H" "R RT C D""'
	local page=1
	local appendixfiles ""
	foreach stages in `st' {
		foreach v in `stages' {
		    capture graph drop `v'
			graph twoway line `v' t, ///
			    name("`v'") graphregion(fc(white)) lc(maroon) ///
				ytitle("Persons") title(`"`:variable label `v''"') ///
				ylabel(,format(%9.0fc))
		}
		capture graph drop p`page'
		graph combine `stages', cols(2) xsize(11.0) ysize(8.0) ///
								scale(0.66) graphregion(fc(white)) ///
								name("p`page'")
		tempfile f
		mata st_local("f", pathrmsuffix(`"`f'"')+".png")
		quietly graph export `"`f'"', name("p`page'") as(png) width(`imgwidth') replace
		graph drop p`page'
		local appendixfiles `"`appendixfiles' `f'"'
		local page=`page'+1
	}

	confirm file `"`resultslog'"'
	local appendixfiles `"`appendixfiles' `beds_files'"'

	if (`"`report'"'!="") {
		popreport /*t S Z I H R RT C D*/, ///
				 modelname("siz") modelparams("`modelparams'") ///
				 appendixgraphs(`appendixfiles') ///
				 results(`"`resultslog'"') reflethality(`reflethality') ///
				 popstruct(`"`populog'"') simparams(`"`simparams'"') ///
				 save(`report')
	}

	foreach appfile in `appendixfiles' {
	    capture erase `"`appfile'"'
	}
end

program define popreport
    version 16.0
	
	local tfont `""Helvetica",28,steelblue"'
	local sfont `""Helvetica",12,steelblue"'

	syntax [varlist(default=none)], modelname(string) modelparams(string) ///
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
	epipop pdfreport putstaticimage "epimodels_eq_`modelname'.png" // model equations

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
		putpdf table t(.,.), halign(center) valign(center) nformat(%12.0fc) font("Consolas",8)
		putpdf table t(1,.), bgcolor("aliceblue")
		putpdf paragraph
	}

	foreach v in `varlist' {
		putpdf text ("`v': `:variable label `v''`=char(10)'")
	}

	// optional appendix	
	epipop pdfreport putallgraphfiles, graphs(`appendixgraphs') landscape

	putpdf save `"`save'"' , replace
end

// END OF FILE
