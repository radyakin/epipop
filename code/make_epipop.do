
clear all
local vers "16.0"
version `vers'
local compdate "`c(current_date)'"

cd "..\code"

mata

real matrix epipop_about() {
    st_local("compile_date", "`compdate'")
	st_local("compile_version", "`vers'")
	st_local("epipop_version", "1.1.1")
}

real function versioncompare(s1,s2) {
     /* Compare two versions with arbitrary number of subnumbers.
	    For example:
		  5.34 vs 3.75.02
		  or
		  7.135.01 vs 7.135.04
		
		Different number of subcomponents is permitted. Zeros are implied,
		where the value is shorter, e.g. comparing 5 vs 4.89 is same as 
		comparing 5.00 vs 4.89.
		
		Each component is assumed a number and compared numerically, 
		e.g.: 5.11>5.2 but 5.11<5.20
		
		An empty version is treated as version 0.
		  
	    Returns: 
	      -1 when v1 is smaller than v2
           0 when v1 is same as v2
           1 when v1 is larger than v2
		   
		Minimal error handling, e.g. for the following situations:
		   - version starts or ends with a dot;
		   - version contains any non-numeric characters except a dot;
		   - two dots are one after another.
     */
	 
	 delim="."   // delimiter used to separate the subcomponents of the versions.
	 
     v1=tokens(s1,delim)
	 v2=tokens(s2,delim)
	 n=min((cols(v1),cols(v2)))
	 for(i=1;i<=n;i++) {
	   if (v1[i]==delim & v2[i]==delim) continue
	   
	   n1=strtoreal(v1[i])
	   n2=strtoreal(v2[i])
	   
	   if (n1==. | n2==.) exit(error(109))
	   
	   if (n1<n2) {
		 return(-1)
	   }
	   if (n1>n2) {
		 return(1)
	   }
	 }
	 
	 if (cols(v1)>cols(v2)) {
		 return(1)
	 }
	 
	 if (cols(v1)<cols(v2)) {
		 return(-1)
	 }

	 return(0)
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

real rowvector multinomial(real N, real rowvector P) {
      /* The probability inputs should be normalized. As an implementation 
	     detail, the value of the last entry is ignored and assumed to take 
		 up any leftover probability mass.

		 ValueError: pvals < 0, pvals > 1 or pvals contains NaNs
		 ValueError: sum(pvals[:-1]) > 1.0
		 (0,0  ,0,0) is treated as (0,0,0,1)
		 (0,0.2,0,0) is treated as (0,0.2,0,0.8)
	  */

      if (N!=floor(N)) {
	    printf("{error}Error! First parameter N must be integer, received: %g{break}",N)
		exit(error(3001))
	  }
	  s=sum(P)
	  if ((s<0.0000) | (s>1.0000)) {
	    printf("{error}Error! Probabilities must sum to 1.0: but they sum to %10.6g{break}",s)
		exit(error(3001))
	  }
	  if (s==0 & N>0) {
	    N, s
		exit(error(3001))
	  }
	  P[1,cols(P)] = P[1,cols(P)] + (1.0000-s)
	  	  
	  result=J(1, cols(P), 0)

	  for(i=1;i<=N;i++) {
	    Z=rdiscrete(1,1,P)
		result[1,Z]=result[1,Z]+1
	  }
	  return(result)  
}

real matrix GetP(p, pRt, pH, q, w) {
	P1=0, (1-p),	0,	0,	0,	0,	0,	p,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0 \ ///   /* 1 */
       0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0     \ ///   /* 2 */
       0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0     \ ///   /* 3 */
       0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0             /* 4 */
	P2=0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0     \ ///   /* 5 */
       0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0     \ ///   /* 6 */
       0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0     \ ///   /* 7 */
       0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0             /* 8 */
	P3=0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	pRt,	0,	pH,	0,	0,	0,	0,	0,	0 \ ///   /* 9 */
       0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0     \ ///   /* 10 */
       0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0     \ ///   /* 11 */ 
       0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0             /* 12 */
	P4=0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0     \ ///   /* 13 */
       0,	0,	0,	0,	0,	0,	0,	0,	0,	(1-q),	0,	0,	0,	0,	q,	0,	0,	0,	0 \ ///   /* 14 */
       0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0     \ ///   /* 15 */
       0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0             /* 16 */
	P5=0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0     \ ///   /* 17 */
       0,	0,	0,	0,	0,	0,	0,	0,	0,	(1-w),	0,	0,	0,	0,	0,	0,	0,	0,	w \ ///   /* 18 */
       0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0             /* 19 */
    
	return(P1\P2\P3\P4\P5)
}

real matrix stochPois(real N, real matrix c, real R0, real theta, real C3) {

	    // R0 is pre-adjusted to capture the effects of C1 and C2
		
		// CHECK PARAMETERS
		
		// N may not be negative
		if (N<0 | missing(N)) {
		    printf("{error}Error! Parameter N may not be negative or missing.")
			exit(error(3001))
		}

		// R0 may not be negative
		if (R0<0.00 | missing(R0)) {
		    printf("{error}Error! Parameter R0 may not be negative or missing.")
			exit(error(3001))
		}
		
		// Theta must be percentages
		if (! ((theta>=0) & (theta<=100))) {
		    printf("{error}Error! Parameter theta was expected in percent [0%..100%].")
			exit(error(3001))
		}
		theta=theta/100.00
		
		// C3 must me 1 or larger and no greater than 5
		if (!((C3>=1.00) & (C3<=5.00)))  {
		    printf("{error}Error! Parameter C3 must be in the range [0..5].")
			exit(error(3001))
		}
		
		if (min(c)<0) {
		    printf("{error}Error! Negative coefficients are not allowed in lethality matrix.")
			exit(error(3001))
		}
		
		// SET PARAMETERS AS IN THE PICKLE FILE
		q=0.25
		w=0.67
		alpha=0
		
		laI=0.4
		muI=0.4
		laZ=0.5454644629902362
		muZ=0.5454644629902362
		laRT=0.3333333333333333
		muRT=0.5
		laH=0.5264432441538478 
		muH=1.0169836265636123
		laC=0.4672897196261682
		muC=0.735348187366718

		/*   NOT USED??
		AdjustedThetas=
		1.58011500e-05, 3.74938000e-06 \ ///
		7.81733000e-06, 7.31053000e-06 \ ///
		8.95337900e-05, 7.83425400e-05 \ ///
		5.44984700e-04, 2.40660680e-04 \ ///
		1.71881400e-03, 9.95424160e-04 \ ///
		4.39487883e-03, 2.25199399e-03 \ ///
		6.67454122e-03, 4.57952314e-03
        */
		
		/* AT=st_matrix("AT") */ // NOT USED???

// ----------------------------------------------------------------------------

	pmod=quadsum(c)   // This is pqw
    // #####  parameter increase due to conditions of camp ##########
	// pmod: Stata: .0003934154 = Python: 0.00039341542228
    p = pmod / (q*w)
	p=p*C3 // Stata: .0070462464 = Python: 0.0070462463691940294
	
	LA=R0 / (p*(1/laI+1/muI) + (1-p)*(3/laZ+3/muZ))
	LA // Stata: .2737844808 ~ Python: 0.27377951887682483
	pH = ((1 - theta)/(1 + alpha + theta + theta^2))
    pRt = 1-pH

	P=GetP(p, pRt, pH, q, w)
	
	Rates = LA*(1 - p)/N , laZ*(1+theta), laZ*(1+theta), laZ*(1+theta), ///
	        muZ*(1+theta), muZ*(1+theta), muZ*(1+theta), ///
			laI*(1+alpha)*(1+theta),muI*(1+alpha)*(1+theta), 0, laRT, ///
			muRT, laH, muH, laC, laC, muC, muC, 0
			
	Y= J(1,19,0)
    Y[1]=N-1
    Y[2]=1
 
    deltaT = 1/24 // time increment of 1 hour
	row=0         // step counter
    act=0         // actual time  
    active = 1    // any positive value to enter the loop
    TotHosp = 0
    TotICU = 0    
	
	tFirstDeath=.
	tFirstHosp=.
	RESULTS=.

	while (active > 0) {
        SI=Y[8]+Y[9]
        SZ=Y[2]+Y[3]+Y[4]+Y[5]+Y[6]+Y[7]
        
        Rates[1]=LA*(SI+SZ)/N   // ???? (1-p) ???

        pr=1:-exp(-Rates*deltaT)        
        TH=Y:*pr // TH is w in pseudocode description
        rV = rpoisson(1,1,TH) // rowvector-19
		rV=mm_cond(rV:!=.,rV,0)
		quita = colmin(Y \ rV) // // quita is U in pseudocode description (rowvector-19)
		
		x=J(19,19,0)
		
		for(i=1;i<=rows(P);i++) {
		    if (sum(P[i,.]:!=0)==1) {
			    x[i,.]=quita[i]*P[i,.]
			}
			else {
			    // https://www.stata.com/statalist/archive/2008-02/msg01145.html
				x[i,.]=multinomial(quita[i],P[i,.])
			}
		}

		pon = colsum(x)
		Y=Y-quita+pon
				
        TotHosp = TotHosp + pon[13]
        TotICU = TotICU + pon[15]
                
        act = act + deltaT

        S = Y[1]
        Z = sum(Y[2..7])
        I = sum(Y[8..9])
        R = Y[10]
        RT= sum(Y[11..12])
        H = sum(Y[13..14])
        C = sum(Y[15..18])
        D = Y[19]
		
        if (H>0 & tFirstHosp==.)  tFirstHosp=round(act)
		if (D>0 & tFirstDeath==.) tFirstDeath=round(act)
			
        row=row+1
        active = Z+I+RT+H+C
		if (RESULTS==.) RESULTS=(act, S,Z,I,R,RT,H,C,D)
		else RESULTS=RESULTS \ (act, S,Z,I,R,RT,H,C,D)
	}

	st_local("tFirstDeath", strofreal(tFirstDeath))
	st_local("tFirstHosp", strofreal(tFirstHosp))
	st_local("TotHosp", strofreal(TotHosp))
	st_local("TotICU", strofreal(TotICU))
	st_local("TotDeaths", strofreal(D))
	
	st_matrix("RESULTS",RESULTS)
	return(RESULTS)
}


mata mlib create lepipop, replace
mata mlib add lepipop epipop_about()
mata mlib add lepipop versioncompare()
mata mlib add lepipop epi_sim_siz()

mata mlib add lepipop multinomial()
mata mlib add lepipop GetP()
mata mlib add lepipop stochPois()
mata mlib index

end


mata assert(versioncompare("2.34","5.11")==-1)
mata assert(versioncompare("5.34.11","5.11.72")==1)
mata assert(versioncompare("19.02.133","19.02.133")==0)
mata assert(versioncompare("19.02.133","19.02.133.99")==-1)
mata assert(versioncompare("02.1.1","2.1.1")==0)
capture noisily mata versioncompare("a.b.c","2.6.7.8.9")
assert _rc==109

mata assert(versioncompare("","1.0")==-1)

display "All OK"

// END OF FILE
