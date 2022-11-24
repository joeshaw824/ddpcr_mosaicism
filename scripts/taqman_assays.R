################################################################################
## ThermoFisher TaqMan assay manufacturing information
## Joseph.Shaw@gosh.nhs.uk / joseph.shaw3@nhs.net
################################################################################

#########################
# Set working directory
#########################

library(tidyverse)

setwd("//fsdept/deptdata$/Regional Genetics Service/Validation Documents/Mosaic/ddPCR/")

#########################
# Read in vial IDs
#########################

# I barcode scanned these vial IDs from the Taqman assay tray in the lab.

vial_ids <- read_csv("ddpcr_mosaicism/resources/taqman_id_vial_ids.csv") %>%
  mutate(vial_id = gsub("\\", "", vial_id, fixed = TRUE))

#########################
# Read in assay info
#########################

read_taqman <- function(filename, filepath) {
  
  taqman_file <-  read_delim(paste0(filepath, filename), 
           delim = "\t",
           skip = 4) %>%
  janitor::clean_names() %>%
  filter(customer_name != "*** End Assay Information\r")
  
  return(taqman_file)
  
}

taqman_files <- list.files(path = "ddpcr_mosaicism/resources/assay_manufacturing/")

collated_files <- data.frame()

# Read and collate each worksheet csv
for (filename in taqman_files){
  tmp_dat <- read_taqman(filename,
                     "ddpcr_mosaicism/resources/assay_manufacturing/")
  
  collated_files <-rbind(collated_files, tmp_dat)
  rm(tmp_dat)
}

collated_files_edit <- collated_files %>%
  # Remove duplicates
  filter(!duplicated(assay_id)) %>%
  # Remove whitespace from assay ID
  mutate(assay_id = gsub(" ", "", assay_id))

write.csv(collated_files_edit, "ddpcr_mosaicism/resources/collated_assay_information.csv",
          row.names = FALSE)

#########################