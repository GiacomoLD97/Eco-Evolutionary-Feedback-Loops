
##### This files plots the relative benefit between two levels of intraspecific intransitivity

### 1. Load required packages ####

# Install packages before loading if needed

library(viridis)
library(tidyverse)
library(here)

### 2. Set Relative Directories to execute script ####

here()

here::i_am("Code/Simulated Experiments/3_analyze_plot_results.R")

# read in the data
dt <- read_csv(
  here("Data", "Simulation_results.csv") 
  )


### 3. Plot the results ####

# spread the data and calculate the difference in final intransitivity between the two scenarios
relben <- dt %>% 
	select(id, sim, iter, p, m, n_init, rho, model, init_bab_sum, final_bab_sum) %>% 
	spread(model, final_bab_sum) %>% 
	mutate(intran =(intran - tran))

# plot
(g1 <- ggplot(relben, aes(x = p, y = rho, z = intran)) +
		stat_summary_2d(bins =  c(23,23),color = NA, fun = mean) +
		labs(title = "", x = "Phenotypic memory (p)", y = expression(paste("Phenotypic similarity (",tau,")")), fill = "Relative\nBenefit\n") + 
		theme_classic()+
		facet_grid(n_init~m)+
		scale_fill_viridis(option = "A"))
