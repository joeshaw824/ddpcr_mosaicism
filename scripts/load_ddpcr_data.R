################################################################################
## Load ddPCR data
################################################################################

library(tidyverse)
library(readxl)

#########################
# Read in ddPCR data 
#########################

data_filepath <- "//fsdept/deptdata$/Regional Genetics Service/Validation Documents/Mosaic/ddPCR/ddpcr_mosaicism/data/"

ddpcr_files <- list.files(path = data_filepath)

#Empty data frame
ddpcr_data <- data.frame()

# Read and collate each worksheet csv
for (dataFile in ddpcr_files){
  tmp_dat <- readr::read_csv(paste0(data_filepath,dataFile), col_names = TRUE,
                             show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    # Each data file is the worksheet number. Remove ".csv" from filename
    mutate(worksheet = substr(as.character(dataFile), 1, 7),
           # Add a unique identifier for each well.
           worksheet_well_sample = paste(worksheet, well, sample, 
                                         sep = "_"))
  ddpcr_data <-rbind(ddpcr_data, tmp_dat)
  rm(tmp_dat)
}

ddpcr_data <- ddpcr_data %>%
  # Change concentration from type chr to type numeric. Rename to be explicit.
  mutate(copies_per_ul = as.numeric(ifelse(concentration == "No Call", NA, concentration))) %>%
  select(-concentration)

#########################