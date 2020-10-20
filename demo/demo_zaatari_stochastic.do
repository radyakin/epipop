clear all
set seed 1 //33 // 12345678
version 16.0

adopath ++ "..\..\epimodels\code\"
adopath ++ "..\code\"
mata mata mlib index

use "C:\data\zaatari.dta"

epipop simulate stochastic, agevar(dem_age) sexvar(dem_sex) ///
               malecode(2) femalecode(1) ///
		       r0(3.0) theta(0.0) ///
			   c1(0.0) c2(0.0) c3(3.0) ///
			   repeat(10)

// EOF