################################################################################################
##### This files plots the relative benefit between two levels of intraspecific intransitivity

rm(list=ls())
gc()

library(viridis)
library(tidyverse)

setwd("~/Git/ecoevo_intran/")

# read in the data
dt <- read_csv("data/simulation_results.csv") 

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
