program define epipop
	version 16.0
	`0'
end

program define simulate
	version 16.0
	gettoken model 0 : 0
	local cmd = "epipop_"+strlower(`"`model'"')
	`cmd' `0'
end

program define dialog
	version 16.0
	syntax , repeat(integer) [*]
	
	if (`repeat'==1) {
	  // deterministic
	  epipop_dbdeterministic, `options'
	  
	}
	else {
	  // stochastic
	  epipop_dbstochastic, repeat(`repeat') `options'	  
	}
end

program define epipop_dbdeterministic

	version 16.0
	
	syntax , agevar(string) sexvar(string) ///
	         malecode(integer) femalecode(integer) [*]
	
	epimodels_util popmatrix , ///
	  agevar(dem_age) sexvar(dem_sex) ///
	  malecode(`malecode') femalecode(`femalecode') ///
	  at(0,10,20,30,40,50,60,130) /* for the moment don't change this */
	  
	tempname F
	matrix `F'=r(F)
	local popsize=r(N)

	epipop simulate deterministic, popsize(`popsize') popstruct("`F'") ///
	           agpop(1 2 7) /* for the moment don't change this */ ///
			   `options'
end

program define epipop_dbstochastic

	version 16.0

	syntax , [*] [refdata(string)] [tmax(integer 1)]
	
	// tmax is ignored in stochastic model
	// refdata at the moment doesn't have any alternative.
	
	epipop simulate stochastic, `options'
end	

// end of file
