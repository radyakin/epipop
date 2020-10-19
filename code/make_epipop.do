
clear all
local vers "16.0"
version `vers'
local compdate "`c(current_date)'"

cd "..\code"

mata

real matrix epipop_about() {
    st_local("compile_date", "`compdate'")
	st_local("compile_version", "`vers'")
	st_local("epimodels_version", "2.1.1")
}


void function epi_sim_siz(string fmatname, real N, real R0, real theta, ///
                    real C1, real C2, real C3, real tmax, string agpop) {

    // Clears the data
    // Relies on fixed names as produced by SIZ
					
    // --------------------------------STABLE PARAMETERS------------------------
   
	AR=st_matrix(st_local("AR")) 
	// Results aggregation
	AR2=AR[1,.] \ colsum(AR[2..6,.]) \ AR[7,.] // reduce to ages 0-9, 10-59, 60+
	// AR2=AR // for no aggregation
	
	st_local("burdenrownames",`""0-9" "10-59" "60-129" "Total""')
	st_local("burdencolnames",`""Deaths" "Hosp." "ICU" "')

	mu1   = st_numscalar("mu1")
	mu2   = st_numscalar("mu2")
	mu3   = st_numscalar("mu3")
	mu4   = st_numscalar("mu4")
	q     = st_numscalar("q")
	w     = st_numscalar("w")
	alpha = st_numscalar("alpha")
   
    // -------------------------------------------------------------------------

	F=st_matrix(fmatname)
	TH=AR:*F
	pmod=sum(TH) // .00042468 == 0.00042468 (Matlab) = OK
	p=pmod/(q*w) // .002535403 == 0.002535402985075 (Matlab) = OK
	LAMBDA=R0/(p/mu2 + (1-p)/mu1)
	LAMBDA=LAMBDA*(1-C1)*(1-C2) // .2731049624 == 0.273104962438366 (Matlab) = OK
	p=p*C3
	
	stata(sprintf(`"epimodels simulate siz , "' ///
      + `"lambda(%g) mu1(%g) mu2(%g) mu3(%g) mu4(%g) pi(%g) psi(%g) omega(%g) theta(%g) alpha(%g) "' ///
      + `"pops(%g) popi(%g) steps(1) days(%g) clear nograph "' ///
      + `"xsize(16) ysize(8) legend(size(2) pos(2) cols(1)) graphregion(fc(white)) "', ///
	  LAMBDA, mu1, mu2, mu3, mu4, p, q, w, theta, alpha, N-1, 1, tmax ///
	  ))
	  
	dead=st_data(st_nobs(), "D") // 33.50166 vs 33.3344 (Matlab) ~ OK
	AR2=AR2:/sum(AR2)
    Burden1=round(dead*AR2)
    Burden1=rowsum(Burden1\colsum(Burden1))
	Burden3=Burden1/w
	Burden2=Burden1/(q*w)
	Burden=round((Burden1,Burden2,Burden3)) // (vs Matlab) ~ OK
	
	Total_num_of_infections=round(N-st_data(st_nobs(),"S")) // 78484 vs 78496 (Matlab) ~ OK
	
	pos=epi_obsmax("H")
	Max_num_of_hospitalized_in_a_day=st_data(pos,"H")
	tmaxH=st_data(pos,"t") // 68 vs 67 (Matlab)
	
	pos=epi_obsmax("C")
	Max_num_of_ICU_in_a_day=st_data(pos,"C")
	tmaxC=st_data(pos,"t") // 74 vs 72 (Matlab)
	
	stata("summarize t if H>=1, meanonly")
	Time_to_first_hosp=st_numscalar("r(min)") // ...... (Matlab)
	
	stata("summarize t if D>=1, meanonly")
	Time_to_first_death=st_numscalar("r(min)") // 53 vs 52 (Matlab)

	infe=st_tempname()
	stata(sprintf("generate double %s = (%g-S)/%g",infe,N,Total_num_of_infections))
	
	stata(sprintf("summarize t if %s>0.25, meanonly",infe))
	t25=st_numscalar("r(min)")
	
	stata(sprintf("summarize t if %s>0.50, meanonly",infe))
	t50=st_numscalar("r(min)")
	
	stata(sprintf("summarize t if %s>0.75, meanonly",infe))
	t75=st_numscalar("r(min)")
	
	st_matrix("B",Burden)
	st_matrix("T",(t25,t50,t75))
	st_numscalar("tni",Total_num_of_infections)
	st_numscalar("mhosp",Max_num_of_hospitalized_in_a_day)
	st_numscalar("tmaxH",tmaxH)
	st_numscalar("micu",Max_num_of_ICU_in_a_day)
	st_numscalar("tmaxC",tmaxC)
	st_numscalar("thosp1",Time_to_first_hosp)
	st_numscalar("tdeath1",Time_to_first_death)
}

mata mlib create lepipop, replace
mata mlib add lepipop epipop_about()
mata mlib add lepipop epi_sim_siz()
mata mlib index

end

// END OF FILE
