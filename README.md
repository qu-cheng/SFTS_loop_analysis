# SFTS_loop_analysis

This repository contains data and code for reproducing the analyses in Identifying important transmission loops to the natural persistence of severe fever with thrombocytopenia syndrome virus

## Data

- all_tick_dat.rds: daily counts of questing H. longicornis larvae, nymphs and adults (13,876 rows, 6 columns: clean.var, t, Nymph, Larva, Adult, Lat) for 38 location-years. clean.var identifies the location-year, t is the day index within the year and Lat is latitude.
- base_para.rds:the baseline parameter set used to build the NGM. Values are the column means of 1,000 Latin hypercube draws of the temperature-independent parameters.
- 1000_sets_tick_abundance.rds: 1,000 Latin hypercube draws of total tick abundance (three columns: low, medium, high), used for the three abundance scenarios.


## Code

This folder contains code for the analyses and figures:

- 01_packages_data_functions.R: R script to load packages and data, draw the 1,000 parameter sets, define functions.
- 02_estimate_R0_elasticty.R: R script to evaluate the 27 scenarios in parallel(3 larval : nymphal : adult ratios x 3 co-feeding parameters kappa x 3 abundance levels). For each of 1,000 draws and each of the 38 location-years it returns R0, the loop table with the elasticity of every loop, and the full 7 x 7 elasticity matrix eij.
- 03_process_results.R: R script to flatten the 27 per-scenario files into analysis-ready tables.
- 04_main_figs.R: R script to reproduce main-text Figures.
  
