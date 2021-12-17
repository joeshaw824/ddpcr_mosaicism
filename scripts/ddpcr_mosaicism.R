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

mosaicism_targets <- read_csv("ddpcr_mosaicism/resources/mosaicism_targets.csv")

analysis_wells <- read_csv(
  "ddpcr_mosaicism/resources/mosaicism_analysis_wells.csv") %>%
  mutate(worksheet_well_sample = paste(worksheet, well, sample, 
                                       sep = "_"),
         identity = factor(identity, levels = c("NTC",
                                                "normal",
                                                "patient",
                                                "positive_control")))

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

mosaic_data_wider <- ddpcr_mosaic_data %>% 
  pivot_wider(id_cols = c(worksheet_well_sample, worksheet, well, sample,
                          assay),
              names_from = target_category,
              values_from = c(droplets, positives, 
                              # Channel 1 +, channel 2 - (FAM+, VIC-)
                              ch1_ch2_2,
                              # Channel 1 +, channel 2 + (FAM+, VIC+)
                              ch1_ch2,
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
            "reference_ch1_ch2_2",
            "reference_ch1_ch2")) %>%
  dplyr::rename(fam_positives = "variant_ch1_ch2_2",
                double_positives = "variant_ch1_ch2")

#########################
# Get NGS mosaicism percentages 
#########################

ngs_results <- ddPCR_confirmations %>%
  mutate(ngs = ifelse(is.na(mosaic_miner_vaf), rc_vaf,
                              mosaic_miner_vaf),
         ngs_percent = ngs*100,
         mosaic_miner_percent = mosaic_miner_vaf*100) %>%
  select(sample, mosaic_miner_vaf, rc_vaf, ngs_percent,
         mosaic_miner_percent)


ngs_vs_ddpcr <- mosaic_data_wider %>%
  filter(sample %in% ngs_results$sample) %>%
  filter(worksheet_well_sample %in% analysis_wells$worksheet_well_sample) %>%
  left_join(ngs_results %>%
              select(sample, ngs_percent, mosaic_miner_percent),
            by = "sample")

# Plot all results
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
  theme(panel.grid = element_blank()) +
  labs(x = "NGS read counter (%)", y = "ddPCR (%)",
       title = "ddPCR vs NGS read counter results")

#########################
# FAM positive droplets detected 
#########################

# Number of single positive droplets in samples
mosaic_data_wider %>%
      filter(worksheet_well_sample %in% 
               analysis_wells$worksheet_well_sample) %>%
  left_join(analysis_wells %>%
              select(worksheet_well_sample, identity),
            by = "worksheet_well_sample") %>%
  arrange(identity, fam_positives) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  ggplot(aes(x = worksheet_well_sample,
             y = fam_positives,
             colour = identity)) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "FAM+VIC- droplets",
       title = "FAM positive only droplets in samples") +
  ylim(-10, 400)

# Number of double positive droplets in samples
mosaic_data_wider %>%
  filter(worksheet_well_sample %in% 
           analysis_wells$worksheet_well_sample) %>%
  left_join(analysis_wells %>%
              select(worksheet_well_sample, identity),
            by = "worksheet_well_sample") %>%
  arrange(identity, variant_positives) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  ggplot(aes(x = worksheet_well_sample,
             y = variant_positives,
             colour = identity)) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "FAM+ droplets",
       title = "Positive droplets in samples")

#########################