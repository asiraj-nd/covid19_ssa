# COVID-19 SSA, Early estimates of COVID-19 infections

This repository contains codes used in the preprint:

[Siraj AS, Worku A, Berhane K, Berhane Y, Siraj DS. 2020] (https://www.medrxiv.org/content/10.1101/2020.04.07.20053421v2) - also in review.

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/). 

====================

## code folder

The scripts in this folder run model simulations, generate results on a daily basis, summarize results into manageable tables, and postprocess results to perform further analyses. The simulation model is an implementation of a process based model structured in to Susceptible, Exposed, Infectious and Removed (SEIR) states.
The main model uses the Partially Observed Markov Process (POMP), with stochastic process model [King et al. 2016], and most parameters obtained from recent studies. 

* 0_numfunctions.R
* 0_snippet_setup.R
* 1_main_model_setup_simulation_output
* 2_generate_cumulative_infections_by_day_t
* 3_postprocess_cumulative_mortality_by_day_t


Each of the scripts in the folder should run in the order indicated starting with #1, #2 and #3. The two initial scripts header files, and are loaded automatically from within the three other.


## data folder

Data included here relate to the age pyramid obtained from Ethiopia [CSA, 2012], and age specific case fatality rates for China [The Novel Coronavirus Pneumonia Emergency Response Epidemiology Team, 2020].


## outputs folder

This is where simulation and other processed data are stored. 


## plots folder

This folder contains figures output based on visualization scripts.  


##References
1. Siraj AS, Worku A, Berhane K, Berhane Y, Siraj DS. Early estimates of COVID-19 infections in small, medium and large population clusters. MedRxiv, doi: https://doi.org/10.1101/2020.04.07.20053421.
2. King AA, Nguyen D, Ionides EL. Statistical Inference for Partially Observed Markov Processes via the R Package pomp. J Stat Soft 2016; 69. DOI:10.18637/jss.v069.i12. 
3. Central Statistics Authority, CSA. The 2007 Population and Housing Census of Ethiopia. Addis Ababa, Ethiopia, 2012.CSA
4. The Novel Coronavirus Pneumonia Emergency Response Epidemiology Team. The Epidemiological Characteristics of an Outbreak of 2019 Novel Coronavirus Diseases (COVID-19) — China, 2020[J]. 2020 http://weekly.chinacdc.cn/en/article/id/e53946e2-c6c4-41e9-9a9b-fea8db1a8f51.
# covid19_ssa
# covid19_ssa
