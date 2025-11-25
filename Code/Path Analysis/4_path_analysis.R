
##### This files plots the relative benefit between two levels of intraspecific intransitivity

### 1. Load required packages ####

library(tidySEM)
library(lavaan)
library(tidyverse)
library(semPlot)
library(here)


### 2. Set Relative Directories to execute script ####

here()

here::i_am("Code/Path Analysis/4_path_analysis.R")

# read in the data
dt <- read_csv(
  here("Data", "Simulation_results.csv") 
)

### 3. Prep the data and take a subset for later modelling ####

# spread the data and calculate the difference in final intransitivity between the two scenarios

relben <- dt %>% 
	select(id, sim, iter, p, m, n_init, rho, model, n_final, init_bab_sum, final_bab_sum) %>% 
    mutate(intran = ifelse(model == "intran", 1, 0)) %>% 
    rename(intran_within = intran) %>% 
    rename(intran_between = final_bab_sum)

set.seed(123)

relben_thous <- relben[sample(nrow(relben), 100000), ]   

relben_mill <- relben[sample(nrow(relben), 1000000), ]   

### 4. Path Analysis ####

#Define model structure
m1 <-   '
  # regressions
  intran_between ~ intran_within + n_init + m + p + rho
  n_final ~ intran_between + n_init
'

#Fit an SEM with a million rows of data
fit <- sem(m1, data = relben_mill)

#Path diagram
graph_sem(model = fit)

#Model fit and variable coefficients
summary(fit, fit.measures = TRUE, rsquare = TRUE, stand = TRUE)

#Refined path diagram
semPaths( 
  object = fit, 
  what = "path", 
  whatLabels = "par",
  style = "ram",
  layout = "tree", 
  rotation = 2, 
  sizeMan = 7, 
  sizeLat = 7, 
  color = "lightgray", 
  edge.label.cex = 1.2, 
  label.cex = 1.3
  )


