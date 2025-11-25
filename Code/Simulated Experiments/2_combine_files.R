
##### This files combines all of the individual files output by 1_simulation_code.R

### 1. Load required packages ####

# Install packages before loading if needed

library(tidyverse)
library(here)

### 2. Set Relative Directories to execute script ####

here()

here::i_am("Code/Simulated Experiments/2_combine_files.R")

# specify the folder to combine
out_folder <- here("Data", "Simulation_results")

### 3. Combine the files ####

# now list all the files and merge
my_files <- list.files(out_folder)
all_files <- tibble()
for(i in 1:length(my_files)){
	print(c(i, length(my_files)))
	all_files <- all_files %>% bind_rows(read_csv(paste0(out_folder, "/", my_files[i]), show_col_types = FALSE))
}

# write the output file
write_csv(all_files, paste0(out_folder, ".csv"))

# delete the temp directory and files
unlink(out_folder, recursive = TRUE)