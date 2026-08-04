# MDSresearch_Final_Report_Ahmed_R-

## Simulating Microfossil Samples 

A Monte Carlo simulation, written in R, that tests whether the total error
formula for microfossil concentration estimation (the linear counting method)
remains accurate when realistic marker-tablet variability is introduced.

`Report_Ahmed_R.Rmd` is the main report
`references.bib` Bibliography.
 `data/` contains CSV files computed at load times 
`numbers_increased.R`  The linear counting implementation; varies the stopping threshold and the target-to-marker ratio separately (the individual sweeps). |
`varying_both_ratio_threhsold.R`  The two-dimensional parameter grid (both parameters varied together). 
 `simulation_stability.R` Convergence run (RMSE against the number of slides). 
 `finalizing_code.R` Slide generation and Poisson verification. 
 `All_plots` contain all the plots used in the report 
