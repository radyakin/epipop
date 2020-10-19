
// This demo uses example data located online in the specified distribution repository

clear all
import delimited "https://raw.githubusercontent.com/reachjor/reachjor.github.io/master/pop_count/data/pop_count.csv", varnames(1) clear

local outfolder="C:\temp\"              // adjust this if necessary.

generate byte sex=cond("Male"==gender,1,cond("Female"==gender,2,.))
generate byte age=real(substr(age_brakets,1,2))
tabulate age sex, miss

epimodels_util popmatrix , agevar(age) sexvar(sex) ///
	                       malecode(1) femalecode(2) ///
		                   at(0,10,20,30,40,50,60,130)
matrix F=r(F)
local popsize=r(N)
display "Population:" string(`popsize',"%8.0f")

epi_pop , popsize(`popsize') popstruct("F") ///
	           r0(3.0) theta(0.00) c3(3) ///
	           agpop(1 2 7) ///
               tmax(200) report("`outfolder'\zaatari_report_open.pdf")

/* Data should be available for at least 2019 according to this:
   https://reliefweb.int/sites/reliefweb.int/files/resources/71531.pdf
   but not clear whether can get access to this data. */
			   
// END OF FILE