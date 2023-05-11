# tradeoff_pop_dynamics

This repository hosts the R scripts to reproduce results and figures from the article *Context-dependent tradeoffs and population dynamic consequences across life histories*

By: Louis Bliard, Maria Paniw, Dylan Z. Childs, Arpat Ozgul 

Please contact louis.bliard@evobio.eu if you have any questions, comments, or concerns.


## GENERAL INFORMATION

1. Title: Simulation code from "Context-dependent tradeoffs and population dynamic consequences across life histories".

2. Author Information:
	
        A.  Name: Louis Bliard
		Institution: University of Zurich
		Address: Winterthurerstrasse 190, 8057 Zurich, Switzerland
		Email: louis.bliard@uzh.ch / louis.bliard@evobio.eu
	
        B.  Name: Maria Paniw
		Institution: Estación Biológica de Doñana
	
        C.  Name: Dylan Childs
		Institution: University of Sheffield
    
        D.  Name: Arpat Ozgul
		Institution: University of Zurich
    

3. Information about funding sources that supported the work: Swiss National Science Foundation Grant "31003A_182286".


## SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC0

2. Links to publications that cite or use the data: ADD DOI ARTICLE

3. Links to other publicly accessible locations of the data: ADD DOI ZENODO

4. Was data derived from another source? No

5. Recommended citation for this dataset: Bliard Louis, Paniw Maria, Childs Dylan, Ozgul Arpat. (2023). Simulation code from "Context-dependent tradeoffs and population dynamic consequences across life histories". [Data set]. Zenodo. ADD DOI


## DATA & FILE OVERVIEW

1. `main simulations` folder: 
This folder contains the R scripts to reproduce the main simulations and figures present in the main text.
- `sa_f.r`
Can be used for the simulation with a tradeoff between adult survival and fecundity.
- `sj_f.r`
Can be used for the simulation with a tradeoff between juvenile survival and fecundity.
- `sa_sj - intergenerational.r`
Can be used for the simulation with an intergenerational tradeoff between parental adult survival and offspring juvenile survival.
- `f_sj.R – intergenerational.r`
Can be used for the simulation with an intergenerational tradeoff between parental fecundity and offspring juvenile survival.

2. `appendix B` folder: 
This folder contains the R scripts to reproduce the simulations and figures present in the appendix B.
- `sa_f_alternative.r`
Can be used for the simulation with an alternative method to change the covariation in Appendix 1, for a tradeoff between adult survival and fecundity.
- `sj_f_alternative.r`
Can be used for the simulation with an alternative method to change the covariation in Appendix 1, for a tradeoff between juvenile survival and fecundity.

3. `appendix D` folder: 
This folder contains the R script to reproduce the simulations and figures present in the appendix D.
- `appendix_density.r`
Can be used for the simulation with decreased density dependence.

4. `appendix E` folder: 
This folder contains the R script to reproduce the simulations and figures present in the appendix E.
- `a4_sa_f.r`
Can be used for the simulation of Appendix 4 with a tradeoff between adult survival and fecundity.
- `a4_sj_f.r`
Can be used for the simulation of Appendix 4 with a tradeoff between juvenile survival and fecundity.
- `a4_sa_sj.r`
Can be used for the simulation of Appendix 4 with an intergenerational tradeoff between parental adult survival and offspring juvenile survival.
- `a4_f_sj.r`
Can be used for the simulation of Appendix 4 with an intergenerational tradeoff between parental fecundity and offspring juvenile survival.

5. `methods animation` folder: 
- `method_main.gif`
Animation displaying how individuals 'move' in the demographic-rate space as the environment changes, following the methods presented in the main text for intraindividual tradeoffs.
- `method_appendix.gif`
Animation displaying how individuals 'move' in the demographic-rate space as the environment changes, following the methods presented in the appendix for intraindividual tradeoffs.
