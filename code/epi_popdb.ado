program define epi_popdb

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

	epi_pop , popsize(`popsize') popstruct("`F'") ///
	           agpop(1 2 7) /* for the moment don't change this */ ///
			   `options'
end

// END OF FILE