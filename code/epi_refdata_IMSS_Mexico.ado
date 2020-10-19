
/* This is a reference file that is not supposed to be edited by the end-user.
If you need to modify any parameters of the model, create a separate reference 
file with epi_refdata_NEWNAME.ado and place the parameter values there. */

program define epi_refdata_IMSS_Mexico, rclass

	version 16.0

	// these parameters are peeked by the epi_simRC() from the mata code
	local mu1   "1/11"
	local mu2   "1/5"
	local mu3   "1/3"
	local mu4   "1/7"
	local q     "0.25"
	local w     "0.67"
	local alpha "0.00"
	
	mata st_local("modelparams", epi_greek(`"mu1 =`mu1'; mu2 =`mu2'; mu3 =`mu3'; mu4 =`mu4'; q =`q'; w =`w'; alpha =`alpha'"'))
	
	tempname M
	matrix `M'= ///
		0.00003246,	0.00000747 \ ///
		0.00001623,	0.00001495 \ ///
		0.00015150,	0.00014200 \ ///
		0.00081702,	0.00039611 \ ///
		0.00199114,	0.00126306 \ ///
		0.00352238,	0.00198054 \ ///
		0.00585981,	0.00461876
	matrix rownames `M' = "0-9" "10-19" "20-29" "30-39" "40-49" "50-59" "60-129"
	matrix colnames `M' = "Males" "Females"
	
	return matrix M=`M'
	
	return scalar mu1   = `mu1'
	return scalar mu2   = `mu2'
	return scalar mu3   = `mu3'
	return scalar mu4   = `mu4'
	return scalar q     = `q'
	return scalar w     = `w'
	return scalar alpha = `alpha'
	
	return local description = "Data from Instituto Mexicano de Seguro Social (Mexican Institute for Social Insurance)"
	return local modelparams=`"`modelparams'"'
end

// END OF FILE
