{smcl}
{* *! version 2.0.0  02oct2020}{...}
{cmd:help epipop}
{hline}

{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{hi: epi_sir} {hline 2} Simulation of COVID-19 epidemics in a closed population.}
{p_end}
{p2colreset}{...}



{title:Syntax}

{p 8 12 2}
{cmd: epi_pop ,}
{it:options}


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}

{syntab :Model conditions}
{synopt :{opt popsize(#)}} Population size{p_end}
{synopt :{opt popstruct(string)}}  Relative frequencies matrix{p_end}
{synopt :{opt r0(#)}}        Basic reproduction number{p_end}
{synopt :{opt theta(#)}} Efficiency in contact trace and removal [0,1]{p_end}
{synopt :{opt c1(#)}}    %-reduction in contact rate by social distance [0,1]{p_end}
{synopt :{opt c2(#)}}    %-reduction in contact rate by use of mask [0,1]{p_end}
{synopt :{opt c3(#)}}    Increase severity for population conditions (default=1.0 meaning no increase){p_end}
{synopt :{opt refdata(string)}} Reference data for mortality adjustment and rest of parameters {p_end}

{syntab :Reporting options}
{synopt :{opt tmax(#)}}  Max time of the simulation (in days){p_end}
{synopt :{opt agpop(string)}} Population aggregation parameters (row numbers){p_end}
{synopt :{opt report(string)}} Optional PDF report name (to be created){p_end}


{synoptline}
{p2colreset}{...}

{title:Description}

{pstd}
{cmd: epi_pop} 

{pstd} The command is used to simulate the outbreak of the epidemics in a 
closed population (of size {opt popsize}) taking into account the population 
structure (represented by matrix {opt popstruct}) using the SIZ model. 
Reference data ({opt refdata}) on morthality from Mexico and estimates of
other parameters from the current literature are used to calibrate the model.
{p_end}

{pstd} Taking the model of spread of disease as given the command focuses on
the conditions of the population, giving the policy analyst control over the
policy-relevant response parameters, such as measures to increase social 
distancing {opt c1}, use of masks {opt c2}, efficiency of testing {opt theta}, 
etc.{p_end}

{pstd} Detailed information about the assumptions, structure, and calibration
of the model is available in {it: C. M. Hernandes-Suarez et al. (2020)}, see the reference 
below. {p_end}

{pstd} {cmd: epipop} requires {cmd: epimodels}, see {browse "http://www.radyakin.org/stata/epimodels/":epimodels homepage}.
{p_end}


{title:Dialog window}

{pstd}{cmd: epipop} comes with an interactive dialog window, which simplifies 
the use of the model by novice users. Notably the dialog window is to be used
in conjunction with the microdata file (person level file), from which the 
relevant proportions by sex and age groups will be calculated. Hence, the dialog
window does not ask for the population size and population sex-age proportions, 
but instead asks for the age and sex variables and numerical codes corresponding
to males and females. The data file to be used is assumed to be loaded in 
Stata's memory.{p_end}

{pstd}Only the specified age and sex variables are used from the microdata file,
no names, household IDs, or other attributes are needed.{p_end}

{pstd}To open the dialog window type the following command in Stata's command 
line: {p_end}

{phang2}{cmd:db epi_pop}{p_end}

{title:Reference data and parameters}

{pstd}{cmd: epipop} uses data derived from the Mexican IMSS data on mortality and
parameters from the literature for various coefficients of the model. They are
set up in the file {it:epi_refdata_IMSS_Mexico.ado}. It is not recommended for the end
user to edit this file. In case other values of parameters are deemed more 
relevant for the study, the user shall:{p_end}
- create a copy of this file;
- rename it appropriately, to hint at the source, for example, {it:epi_refdata_Italy2021.ado};
- rename the sole program defined by this file to match the name of the file;
- replace any of the values {it:mu1, mu2, mu3, mu4, q, w, alpha} or matrix {it:M} in this file to custom values;
- edit the description line to indicate the source of the reference data;
- register the new set of reference parameters to be selectable in the dialog, by editing the section: 

LIST refdatalist
BEGIN
  IMSS_Mexico
END

to include the name of so defined program (only mention the part after 'epi_refdata_').




{title:Population structure matrix}

{pstd}{cmd: epipop} uses a matrix describing the population structure. This 
matrix is expected in the following form:{p_end}

    ------------------------
          |  MALES  FEMALES
    ------------------------
    00-09 |   F11     F12
    10-19 |   F21     F22
    20-29 |   F31     F32
    30-39 |   F41     F42
    40-49 |   F51     F52
    50-59 |   F61     F62    
    60-++ |   F71     F72 
    ------------------------

{pstd}Values Fij are proportions of the corresponding age-sex group in the 
population meaning the {it:total sum} of Fij is 1.0 (not per column or per row).
{p_end}

{pstd}Certain age groups can be combined into larger groups in the reporting,
if desired. This is done by mentioning the indices of each group (from the 
above) to form a new output group, with (1,2,3,4,5,6,7) implied by default.
For example, (1,2,7) will report on three age groups: 0-9, 10-59, and 60+ years
old.{p_end}



{title:Examples}

    {hline}
{pstd}Simulation{p_end}

{phang2} {cmd:. epipop simulate deterministic} , popsize(100000) popstruct(F) /// {p_end}
	            r0(3.0) theta(0.00) agpop(1 2 7) ///
                    tmax(200) report("C:\out\camp_report.pdf")  

{pstd}Perform SIZ model simulation for a population of 100000 individuals, with 
population structure given by the matrix F, with no testing, basic reproduction
number 3.0 over 200 days, saving the report to the file camp_report.pdf and
grouping age groups 2..6 together in the report into one age group.{p_end}

{phang2} {cmd:. epipop simulate deterministic} , popsize(100000) popstruct(F) /// {p_end}
	            r0(3.0) c3(3.0) theta(0.00) agpop(1 2 7) ///
                    tmax(200) report("C:\out\camp_report.pdf")  

{pstd}Same as above, but with population conditions 3 times worse than for 
reference population.{p_end}

{title:References}

{pstd}CARLOS M Hernandez-Suarez, Paolo Verme, Sergiy Radyakin, Efren Murillo-Zamora (2020).{p_end}
{pstd} 
COVID-19 Outbreaks in Refugee Camps. A simulation study. {browse "https://www.medrxiv.org/content/10.1101/2020.10.02.20204818v1":online}{p_end}



{title:Authors}

{phang}
{it:Carlos M Hernandez-Suarez}, Universidad de Colima
{p_end}

{phang}
{it:Sergiy Radyakin}, The World Bank
{p_end}

{phang}
{it:Paolo Verme}, The World Bank
{p_end}

{title:Also see}

{psee}
Online: {browse "http://www.radyakin.org/stata/epipop/": epipop homepage}
{p_end}


