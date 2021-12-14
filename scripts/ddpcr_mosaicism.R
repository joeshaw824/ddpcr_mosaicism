################################################################################
## ddPCR for confirmation of mosaicism
## December 2021
## Joseph.Shaw@gosh.nhs.uk
################################################################################

#########################
# Set working directory
#########################

library(tidyverse)
library(readxl)

setwd("//fsdept/deptdata$/Regional Genetics Service/Validation Documents/Mosaic/ddPCR/")

#########################
# Get resources
#########################

ddPCR_confirmations <- read_excel("ddPCR_designs_confirmations.xlsx",
                                  sheet = "ddPCR_confirmation_list") %>%
  janitor::clean_names()

mosaicism_targets <- read_excel("ddPCR_designs_confirmations.xlsx",
                                sheet = "targets")

analysis_wells <- read_excel("ddPCR_designs_confirmations.xlsx",
                             sheet = "sample_analysis_wells") %>%
  mutate(worksheet_well_sample = paste(worksheet, well, sample, 
                                       sep = "_"))
#########################
# Read in ddPCR data 
#########################

ddpcr_files <- list.files(path = "ddpcr_mosaicism/data")

#Empty data frame
ddpcr_mosaic_data <- data.frame()

# Read and collate each worksheet csv
for (dataFile in ddpcr_files){
  tmp_dat <- readr::read_csv(paste0("ddpcr_mosaicism/data/",dataFile), col_names = TRUE) %>%
    janitor::clean_names() %>%
    # Each data file is the worksheet number. Remove ".csv" from filename
    mutate(worksheet = substr(as.character(dataFile), 1, 7),
           # Add a unique identifier for each well.
           worksheet_well_sample = paste(worksheet, well, sample, 
                                         sep = "_")) %>%
    dplyr::rename(droplets = accepted_droplets) %>%
    # Add on target (reference or variant)
    left_join(mosaicism_targets %>%
                select(target, target_category, assay), 
              by = "target")
  ddpcr_mosaic_data <-rbind(ddpcr_mosaic_data, tmp_dat)
  rm(tmp_dat)
}

# Get one row per sample well.

mosaic_data_wider <- ddpcr_mosaic_data %>% 
  pivot_wider(id_cols = c(worksheet_well_sample, worksheet, well, sample,
                          assay),
              names_from = target_category,
              values_from = c(droplets, positives, 
                              # Channel 1 +, channel 2 negative (FAM+, VIC-)
                              ch1_ch2_2,
                              concentration,
                              copies_per20u_l_well,
                              fractional_abundance, 
                              poisson_fractional_abundance_max,
                              poisson_fractional_abundance_min),
              # Use names_glue to keep new columns names with naming
              # convention
              names_glue = "{target_category}_{.value}") %>%
  # Remove columns with duplicated values
  select(-c("reference_fractional_abundance", 
            "reference_poisson_fractional_abundance_max",
            "reference_poisson_fractional_abundance_min",
            "reference_ch1_ch2_2")) %>%
  dplyr::rename(fam_positives = "variant_ch1_ch2_2")

#########################
# Get NGS mosaicism percentages 
#########################

ngs_results <- ddPCR_confirmations %>%
  mutate(ngs_percent = pmax(mosaic_miner_vaf, rc_vaf,
                            na.rm = TRUE)*100) %>%
  select(sample, mosaic_miner_vaf, rc_vaf, ngs_percent)

ngs_vs_ddpcr <- mosaic_data_wider %>%
  filter(sample %in% ngs_results$sample) %>%
  filter(worksheet_well_sample %in% analysis_wells$worksheet_well_sample) %>%
  left_join(ngs_results %>%
              select(sample, ngs_percent),
            by = "sample")

# Plot results
ggplot(ngs_vs_ddpcr, aes(x = ngs_percent,
                         y = variant_fractional_abundance)) +
  geom_errorbar(aes(ymin = variant_poisson_fractional_abundance_min,
                    ymax = variant_poisson_fractional_abundance_max),
                alpha = 0.2) +
  geom_point(size = 2, pch = 21, fill = "white") +
  scale_x_continuous(limits = c(0, 11),
                     breaks = seq(from = 0, to = 11, by = 1)) +
  scale_y_continuous(limits = c(0, 11),
                     breaks = seq(from = 0, to = 11, by = 1)) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  labs(x = "ddPCR (%)", y = "NGS (%)",
       title = "ddPCR vs NGS mosaicism results")

#########################
# Positive droplets in normal controls 
#########################

normal_controls <- c("C236", "20RG-209G0042", "21RG-074G0164",
                     "20RG-337G0031", "20RG-153G0082")

mosaic_data_wider %>%
  # Get single well data
  filter(substr(well, 1, 1) != "M" &
           sample %in% normal_controls) %>%
  mutate(sample_type = "normal_control") %>%
  rbind(
    mosaic_data_wider %>%
      filter(worksheet_well_sample %in% analysis_wells$worksheet_well_sample) %>%
      mutate(sample_type = "patient_sample")) %>%
  ggplot(aes(x = reorder(worksheet_well_sample, variant_positives),
             y = fam_positives, 
             colour = sample_type)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 300),
                     breaks = c(0, 10, 50, 100, 200, 300)) +
  labs(x = "",
       y = "FAM+ droplets",
       title = "FAM positive droplets in samples")
  #geom_hline(yintercept = 10, linetype = "dashed") +
  #geom_hline(yintercept = 50, linetype = "dashed")

mosaic_data_wider %>%
  # Get single well data
  filter(substr(well, 1, 1) != "M" &
           sample %in% normal_controls) %>%
  ggplot(aes(y = variant_positives)) +
  geom_boxplot() +
  ylim(0, 20) +
  labs(x = "",
       y = "Positive droplets (variant)",
       title = "Spread of results for normal controls")

#########################
# Concentration of normal controls 
#########################

# I was aiming for 25ng input to each 22ul well. Only 20ul is converted
# into droplets. Assuming a concentration of 1.14ng/ul, that's 22.7ng
# being measured by the reader.
# Predicted number of copies detected:
predicted_copies <- (22.7*1000)/3.3

mosaic_data_wider %>%
  filter(substr(well, 1, 1) != "M" &
           # Exclude worksheets which did not have 25ng input
           !worksheet %in% c("21-0374", "21-2298") &
           sample == "20RG-209G0042" &
           reference_concentration != "No Call") %>%
  ggplot(aes(x = reorder(worksheet_well_sample, reference_copies_per20u_l_well),
              y = reference_copies_per20u_l_well)) +
  geom_point(size = 2, pch =21) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  ylim(0, 15000) +
  geom_hline(yintercept = predicted_copies,
             linetype = "dashed")

test <- mosaic_data_wider %>%
  filter(substr(well, 1, 1) != "M" &
           # Exclude worksheets which did not have 25ng input
           !worksheet %in% c("21-0374", "21-2298") &
           sample == "20RG-209G0042" &
           reference_concentration != "No Call") %>%
  select(worksheet, well, sample, assay,
         reference_concentration)

#########################