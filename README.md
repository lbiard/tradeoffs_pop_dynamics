# tradeoff_pop_dynamics

This repository hosts the R scripts to reproduce results and figures from the article *Population dynamic consequences of context-dependent tradeoffs across life histories*

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

2. Links to publications that cite or use the data:

3. Links to other publicly accessible locations of the data:

4. Was data derived from another source? No

5. Recommended citation for this dataset: Bliard Louis, Paniw Maria, Childs Dylan, Ozgul Arpat. (2024). Simulation code from "Population dynamic consequences of context-dependent tradeoffs across life histories". [Data set]. Zenodo.


## DATA & FILE OVERVIEW

1. `main simulations` folder: 
This folder contains the R scripts to reproduce the main simulations and Figures 3, 4, and 5 present in the main text, as well as Figure S3.
- `sa_f.R`
Can be used for the simulation with a tradeoff between adult survival and fecundity.
- `sj_f.R`
Can be used for the simulation with a tradeoff between juvenile survival and fecundity.
- `sa_sj - intergenerational.R`
Can be used for the simulation with an intergenerational tradeoff between parental adult survival and offspring juvenile survival.
- `f_sj.R – intergenerational.R`
Can be used for the simulation with an intergenerational tradeoff between parental fecundity and offspring juvenile survival.

2. `Appendix B` folder: 
This folder contains the R scripts to reproduce the simulations present in the supplementary materials, section B, and Figure S2.
- `sa_f_alternative.r`
Can be used for the simulation with an alternative method to change the covariation in the supplementary materials section B, for a tradeoff between adult survival and fecundity (Figure S2).
- `sj_f_alternative.r`
Can be used for the simulation with an alternative method to change the covariation in the supplementary materials section B, for a tradeoff between juvenile survival and fecundity (Figure S2).

3. `Appendix D` folder: 
This folder contains the R scripts to reproduce the simulations present in the supplementary materials, section D, and Figure S4 and S5.
- `appendix_density_sa_f.R`
Can be used for the simulation with decreased density dependence and a tradeoff between adult survival and fecundity (Figure S4).
- `appendix_density_sj_f.R`
Can be used for the simulation with decreased density dependence and a tradeoff between juvenile survival and fecundity (Figure S4).
- `density_sa_f_decrease.R`
Can be used for the simulation with decreased density dependence (Figure S5).
- `density_sa_f_increase.R`
Can be used for the simulation with increased density dependence (Figure S5).

4. `Appendix E` folder: 
This folder contains the R script to reproduce the simulations and figures present in the supplementary materials, section E, and Figure 6.
- `a5_sa_f.r`
Can be used for the simulation of the supplementary materials, section E, with a tradeoff between adult survival and fecundity.
- `a5_sj_f.r`
Can be used for the simulation of the supplementary materials, section E, with a tradeoff between juvenile survival and fecundity.
- `a5_sa_sj.r`
Can be used for the simulation of the supplementary materials, section E, with an intergenerational tradeoff between parental adult survival and offspring juvenile survival.
- `a5_f_sj.r`
Can be used for the simulation of the supplementary materials, section E, with an intergenerational tradeoff between parental fecundity and offspring juvenile survival.

4. `Appendix F` folder:

- `a6_sa_sj.R`
Can be used for the simulation and regression analysis of the supplementary materials, section F, with an intergenerational tradeoff between parental adult survival and offspring juvenile survival (Table S1).
- `a6_f_sj.R`
Can be used for the simulation and regression analysis of the supplementary materials, section F, with an intergenerational tradeoff between parental fecundity and offspring juvenile survival (Table S2).


6. `methods animation` folder: 
- `method_main.gif`
Animation displaying how individuals 'move' in the demographic-rate space as the environment changes, following the methods presented in the main text for intraindividual tradeoffs.
- `method_appendix.gif`
Animation displaying how individuals 'move' in the demographic-rate space as the environment changes, following the methods presented in the appendix for intraindividual tradeoffs.
