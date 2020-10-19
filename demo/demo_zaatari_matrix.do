    clear all
    version 16.0
    /* This is a matrix version not using microdata */
	
	local outfolder="C:\temp\"              // adjust this if necessary.

    matrix popF = ///
    /*               MALES  FEMALES        */   ///
    /*  0- 9 */     0.1934,  0.1823         \   ///
    /* 10-19 */     0.1190,  0.1116         \   ///
    /* 20-29 */     0.0714,  0.0794         \   ///
    /* 30-39 */     0.0574,  0.0611         \   ///
    /* 40-49 */     0.0337,  0.0348         \   ///
    /* 50-59 */     0.0127,  0.0171         \   ///
    /* 60-++ */     0.0110,  0.0150

    matrix rownames popF = "0-9" "10-19" "20-29" "30-39" "40-49" "50-59" "60-129"
    matrix colnames popF = "Males" "Females"

    local popsize = 83464

    epi_pop , popsize(`popsize') popstruct("popF") ///
               r0(3.0) theta(0.00) c3(3) ///
               agpop(1 2 7) ///
               tmax(200) report("`outfolder'\zaatari_report3m.pdf") 

    clear

    epi_pop , popsize(`popsize') popstruct("popF") ///
               r0(3.0) theta(0.00) c3(1) ///
               agpop(1 2 7) ///
               tmax(200) report("`outfolder'\zaatari_report1m.pdf") 

// END OF FILE
