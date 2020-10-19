
	clear all
	version 16.0
	/* This is a data version using microdata */
	
	local zaatarifile="C:\data\zaatari.dta" // adjust this to where the data is.
	local outfolder="C:\temp\"              // adjust this if necessary.
	
	frame pwf
	local cf `"`r(currentframe)'"'
	
	tempname tframe
	frame create `tframe'
	frame change `tframe'
	
    use `"`zaatarifile'"'
	
	epimodels_util popmatrix , agevar(dem_age) sexvar(dem_sex) ///
	                           malecode(2) femalecode(1) ///
		                       at(0,10,20,30,40,50,60,130)

	matrix popF=r(F)
	local popsize=_N
	clear
	frame change `cf'
	frame drop `tframe'	

    epi_pop , popsize(`popsize') popstruct("popF") ///
	           r0(3.0) theta(0.00) c3(3) ///
	           agpop(1 2 7) ///
               tmax(200) report("`outfolder'\zaatari_report3d.pdf") 
    clear
			   
    epi_pop , popsize(`popsize') popstruct("popF") ///
	           r0(3.0) theta(0.00) c3(1) ///
	           agpop(1 2 7) ///
               tmax(200) report("`outfolder'\zaatari_report1d.pdf") 
			   
// END OF FILE			   
