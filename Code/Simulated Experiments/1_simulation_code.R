
##### This files runs the replicator model simulations at two levels of intraspecific intransitivity


### 1. Load required packages ####

# Install packages before loading if needed

library(numDeriv)
library(R.utils)
library(Rglpk)
library(doMC)
library(combinat)
library(viridis)
library(doParallel)
library(RColorBrewer)
library(scales)
library(deSolve)
library(reshape2)
library(gridExtra)
library(igraph)
library(tidyverse)
library(here)

### 2. Set Relative Directories to execute script ####

# Assuming you have downloaded the GitHub Repo and opened the project
# Using here to source relative directories 
here()

here::i_am("Code/Simulated Experiments/1_simulation_code.R")

source(here("Code", "Simulated Experiments", "0_functions_replicator.R"))

out_folder <- here("Data", "Simulation_results")

# specify the directory to save results to
dir.create(out_folder, showWarnings = FALSE)

offset <- 1000000

# get the offset in case we've already written some files

have_files <- list.files(out_folder)
if(length(have_files) > 0){
	offset <- max(as.numeric(gsub("\\.csv", "", gsub("file_", "", have_files))))+10000000
}else{
	offset <- 0
}


### 3. Initialize Simulation run ####

# number of threads
nthreads <- 12
registerDoParallel(nthreads)

# threashold to prune the simulations and dynamics. 
THRESH_sim <<- THRESH_prune <<- THRESH <<- 1e-7 

# set the range for memory (determines Q)
pseq <- round(seq(0.25, 0.99, length=25), 3)

# set the range for the relatedness of phenotypes (determines H)
rseq <- round(seq(0, 1, length=25), 3)


# specify the two extremes -- C1 is perfectly transitive, C2 perfectly intransitive
C1 <- list(rbind(c(0.5,1,1), c(0,0.5,1), c(0,0,0.5)),
		   rbind(c(0.5, 1,1,1,1), c(0,0.5, 1, 1, 1), c(0,0,0.5, 1,1), c(0,0,0,0.5,1), c(0,0,0,0,0.5)),
		   rbind(c(0.5, 1,1,1,1,1,1),c(0,0.5,1,1,1,1,1),c(0,0,0.5,1,1,1,1),c(0,0,0,0.5,1,1,1),c(0,0,0,0,0.5,1,1),c(0,0,0,0,0,0.5,1),c(0,0,0,0,0,0,0.5)))

C2 <- list(rbind(c(0.5,0,1), c(1, 0.5, 0), c(0, 1, 0.5)), 
		   rbind(c(0.5, 1, 0, 1, 0), c(0, 0.5, 1, 0, 1), c(1, 0, 0.5, 1, 0), c(0,1,0, 0.5, 1), c(1, 0, 1, 0, 0.5)), 
		   rbind(c(0.5,1,0,0,1,1,0), c(0,0.5, 1, 1, 0, 0, 1), c(1, 0, 0.5, 1,0,0,1), c(1,0,0,0.5,1,0,1), c(0,1,1,0,0.5,1,0), c(0,1,1,1,0,0.5,0), c(1,0,0,0,1,1,0.5)))


# max number of fitting attempts and timeoue
max_count <- 200
max_timeout <- 60*2
kind <- j <- i <- z <- 1

nseq <- c(3,5) #,10,15)
mseq <- c(3,5,7)

# get the sequnce of fitting parameters
seqmat <- expand.grid(n = nseq, m = mseq, r = rseq, p = pseq, iter = 1:500) %>% data.frame() %>% as_tibble() %>% slice(sample(1:nrow(.), nrow(.), replace = FALSE))

# number of fits per thread
neach <- 25

# creat the fit squence
chunk_seq <- seq(1, nrow(seqmat), by = neach)



### 4. Run the Simulation ####

# simulate in parallel, with each thread running 100 sims at a time
res <- foreach(kind = 1:length(chunk_seq), .combine = c, .inorder = FALSE)%dopar%{
	
	# the dataframe to store the results
	outcome_df<-NULL
	
	start <- chunk_seq[kind]
	end <- ifelse(kind == length(chunk_seq), nrow(seqmat), chunk_seq[kind+1]-1)
	
	for(j in start:end){
		# print(paste0(j," of ",start,"-",end))
		my_res <- withTimeout(run_single_sim(j, seqmat, C1, C2, max_count), timeout=Inf, elapsed = max_timeout, onTimeout = "silent")
		if(!is.null(my_res)){
			outcome_df <- outcome_df %>% bind_rows(my_res)
		}
	}

	write_csv(outcome_df, paste0(out_folder,"/file_", kind + offset, ".csv"))

	return(NA)
}
