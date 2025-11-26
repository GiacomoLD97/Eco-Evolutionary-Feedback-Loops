# ReadMe Eco-Evolutionary Feedback Loops

## Overview

### Summary
This code aims to explore the potential for feedbacks between intra and interspecific intransitivity and their impacts on community diversity. The simulations are adapted from: 
[Maynard, D. S., C. A. Serván, J. A. Capitán, and S. Allesina. 2019. Phenotypic variability promotes diversity and stability in competitive communities. Ecology Letters 22:1776–1786.](https://onlinelibrary.wiley.com/doi/10.1111/ele.13356)

### Paper Citation
Awaiting Publication

### Authors
Giacomo Delgado 
Daniel Maynard
Jukka Jokela
Thomas Crowther 

### Author Contributions
DM wrote the code and created the figures
GD helped conceptualize and edit the code
JJ and TC provided feedback on the code

## Layout

### File Architecture

- Code
  - Simulated Experiments
    - 0_functions_replicator.R
    - 1_simulation_code.R
    - 2_combine_files.R
    - 3_analyze_plot_results.R
  - Path Analysis
    - 4_path_analysis.R
- Data
  - Simulation_results
  - Simulation_results.csv (>1GB file, not saved in GitHub)
- R Project
- LICENSE
- README

### R script descriptions

#### 0_functions_replicator.R
This script provides 13 functions which are necessary to run the simulations. It is not necessary to run this script, the functions are called using `source()` in 1_simulation_code.R

#### 1_simulation_code.R
This script uses a modification of the replicator-mutator equation to simulate millions of communities across a wide range of model parameters (specifically, rho, p and intransitivity).
The output of this script is a series of csv files saved to Data -> Simulation_results

#### 2_combine_files.R
This script combines the csv files saved as an output of 1_simulation_code.R and combines them, resaving them as Data --> Simulation_results.csv. The resulting dataframe is large (>1GB)

#### 3_analyze_plot_results.R
This script spreads the results of the simulations to calculate relative impact across two extremes of intransitivity, it then creates Figure 4 from the manuscript. 

#### 4_path_analysis.R
This script spreads the results of the simulations to calculate relative impact across two extremes of intransitivity, it then uses this data to explore the relationships between variables 
using path analysis/SEM using the `lavaan` package. 

### Data descriptions

#### Simulation_results.csv
A dataframe of 6,784,640 observations and 36 variables including. Each row corresponds to a single simulation run, with columns containing either simulation unique identifiers (e.g., id, sim, iter, etc.)
or model parameters for that simulation run (e.g., p, rho, m, etc.) and final simulation results (e.g., nm_final, final_cv, final_mean, etc.)






