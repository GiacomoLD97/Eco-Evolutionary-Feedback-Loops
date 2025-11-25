######################################################################################
##### This files combines all of the individual files output by 1_simulation_code.R

rm(list=ls())
gc()

library(tidyverse)

setwd("~/Git/ecoevo_intran/")


# specify the folder to combine
out_folder <- "data/simulation_results"


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