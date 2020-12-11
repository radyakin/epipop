program define epipop
	version 16.0
	check_epimodels
	`0'
end

program define nothing
    // this code does nothing
	if (0) display "nothing"
end

program define check_epimodels
    // check version of epimodels is as required

	local minversion="2.1.1"
	capture mata epimodels_about()
	mata st_local("r",strofreal(versioncompare("0`epimodels_version'","`minversion'")))
	
	if (`r'==-1) {
	    display as error "Incompatible version of -epimodels-"
		display as error "Please install or update the -epimodels- package before using -epipop- package."
		error 6
	}
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

program define pdfreport
	version 16.0
	`0'
end

program define putoptionalparagraph
        version 16.0
		syntax [anything], title(string) font(string)
		
		if ((`"`anything'"'!="") & (`"`anything'"'!=`"`""'"')) {
			putpdf paragraph , halign(center)
			if (`"`title'"'!="") putpdf text (`"`title'"'), font(`font')
			puttextfile `anything'
		}
end

program define puttextfile
    version 16.0
	syntax [anything]

	if (`"`anything'"'!="") {
	    putpdf paragraph
	    file open fh using `anything', read text
		file read fh line
		while r(eof)==0 {
		    putpdf text (`"`line'`=char(10)' "'), font("Consolas", 9)
			file read fh line
		}
		file close fh
	}
end

program define putstaticimage

	syntax [anything], [width(real 5) linebreak(real 2)]

	if (`"`anything'"'=="") exit

	capture findfile `anything'
	if !_rc {
		putpdf paragraph, halign(center)
		putpdf text ("          ")
		putpdf image "`=r(fn)'" , width(`width') linebreak(`linebreak')
    }
end

program define putallgraphfiles
	version 16.0
	syntax , [graphs(string) landscape]
	
	if (`"`graphs'"'!="") {
	    putpdf sectionbreak, `landscape'
		local first=1
	    foreach f in `graphs' {
		    capture confirm file `"`f'"'
			if !_rc {
				if (!`first') {
					putpdf pagebreak
					local first=0
				}
				putpdf paragraph
				putpdf image `"`f'"'
			}
		}
	}
end

program define putallgraphs

	// this doesn't seem to work in Stata
	// Reported here: 
	//  https://www.statalist.org/forums/forum/general-stata-discussion/general/1584316
	// No responses so far.
	
	version 16.0
		
    syntax , graphs(string) [landscape]
	if (`"`graphs'"'!="") {
	    putpdf sectionbreak, `landscape'
		local first=1
	    foreach f in `graphs' {
		    tempfile tmp
			graph export `"`tmp'.png"', name(`f') replace as(png) width(1920)
			
			if (!`first') {
				putpdf pagebreak
				local first=0
			}
			putpdf paragraph
			putpdf image `"`tmp'.png"'

			capture erase `"`tmp'.png"'
		}
	}
end


// end of file
