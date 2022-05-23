################################################################################
## ddPCR for confirmation of mosaicism
## December 2021
## Joseph.Shaw@gosh.nhs.uk
################################################################################

# Notes
# Use base R trimws to trim whitespace from the Experiment name

#########################
# Set working directory
#########################

library(tidyverse)
library(readxl)
library(ggpubr)

setwd("//fsdept/deptdata$/Regional Genetics Service/Validation Documents/Mosaic/ddPCR/")

#########################
# Get resources
#########################

mosaicism_targets <- read_csv("ddpcr_mosaicism/resources/mosaicism_targets.csv")

mosaicism_assays <- mosaicism_targets %>%
  pivot_wider(id_cols = assay,
              names_from = target_category,
              values_from = c(target, fluorophore))

analysis_wells <- read_csv(
  "ddpcr_mosaicism/resources/mosaicism_analysis_wells.csv") %>%
  mutate(worksheet_well_sample = paste(worksheet, well, sample, 
                                       sep = "_"),
         identity = factor(identity, levels = c("NTC",
                                                "normal",
                                                "patient")))
wells_for_ngs_comparison <- read_csv(
  "ddpcr_mosaicism/resources/ddpcr_wells_ngs_comparison.csv") %>%
  mutate(worksheet_well_sample = paste(worksheet, well, sample, 
                                sep = "_"))

assay_information <- read_excel(
  path = "I:/Genetics/DNA Lab/databases/Specialist_Services/Skin/ddPCR_designs_confirmations.xlsx",
  sheet = "sequences") %>%
  janitor::clean_names()

#########################
# Read in ddPCR mosaic data 
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
    # Add on target (reference or variant)
    left_join(mosaicism_targets %>%
                select(target, target_category, assay), 
              by = "target")
  ddpcr_mosaic_data <-rbind(ddpcr_mosaic_data, tmp_dat)
  rm(tmp_dat)
}

#########################
# Modify ddPCR mosaic data 
#########################

mosaic_data_mod <- ddpcr_mosaic_data %>%
  # Change concentration from type chr to type numeric. Rename to be explicit.
  mutate(copies_per_ul = as.numeric(ifelse(concentration == "No Call", NA, concentration))) %>%
  select(-concentration) %>%
  dplyr::rename(total_droplets = accepted_droplets,
                # "total_conf_max/min" refer to the "concentration" field.
                # "total_conc_max" is the maximum value of the concentration with total error.
                # I renamed these fields to be consistent with the naming convention in
                # other columns.
                total_conc_max = total_conf_max,
                total_conc_min = total_conf_min,
                poisson_conc_max = poisson_conf_max,
                poisson_conc_min = poisson_conf_min,
                copies_per_20ul_well = copies_per20u_l_well)

mosaic_data_wider <- mosaic_data_mod %>% 
  filter(!is.na(target_category)) %>%
  pivot_wider(id_cols = c(worksheet_well_sample, worksheet, well, sample,
                          assay, experiment),
              names_from = target_category,
              # Include all values in case we need them later
              values_from = -c(worksheet_well_sample, worksheet, well, sample,
                               assay, experiment, target_category),
              # Use names_glue to keep new columns names with naming
              # convention
              names_glue = "{target_category}_{.value}") %>%
  mutate(assay_gene = sub("_.*", "", assay))

#########################
# Quality monitoring 
#########################

droplets_qc_plot <- mosaic_data_wider %>%
  filter(substr(well, 1, 1) != "M") %>%
  ggplot(aes(x = worksheet, y = variant_total_droplets)) +
  geom_jitter() +
  labs(y = "Total Droplets", x = "", title = "Droplet generation QC plot") +
  geom_hline(yintercept = 10000, colour = "red", linetype = "dashed") +
  theme_bw()

# Add timestamp
ggsave(plot = droplets_qc_plot, 
       filename = "droplets_qc_plot.tiff",
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff')

ws_data <- mosaic_data_wider %>%
  filter(substr(well, 1, 1) != "M") %>%
  filter(worksheet == "22-1678")

#########################
# NGS vs ddPCR 
#########################

ngs_results <- read_csv("ddpcr_mosaicism/resources/ngs_results.csv") %>%
  janitor::clean_names() %>%
  mutate(ngs = ifelse(is.na(mosaic_miner_vaf), rc_vaf,
                              mosaic_miner_vaf),
         ngs_percent = ngs*100,
         mosaic_miner_percent = mosaic_miner_vaf*100) %>%
  select(sample, mosaic_miner_vaf, rc_vaf, ngs_percent,
         mosaic_miner_percent)

ngs_vs_ddpcr <- mosaic_data_wider %>%
  # When a sample has more than 1 ddPCR replicate, the merged well value is 
  # used. Some samples only had enough for 1 ddPCR well, so merged values 
  # were not calculated in Quantasoft.
  filter(worksheet_well_sample %in% 
           wells_for_ngs_comparison$worksheet_well_sample) %>%
  left_join(ngs_results %>%
              select(sample, ngs_percent, mosaic_miner_percent),
            by = "sample")

sample_number <- nrow(ngs_vs_ddpcr)

# Plot all results
ddpcr_ngs_plot <- ggplot(ngs_vs_ddpcr, aes(x = ngs_percent,
                         y = variant_fractional_abundance)) +
  geom_errorbar(aes(ymin = variant_poisson_fractional_abundance_min,
                    ymax = variant_poisson_fractional_abundance_max),
                alpha = 0.2) +
  geom_point(size = 2, pch = 21, fill = "white") +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "NGS read counter variant fraction (%)", y = "ddPCR variant fraction (%)",
       title = paste("ddPCR and NGS results for", sample_number, "mosaic patient samples")) +
  ggpubr::stat_cor(method = "pearson", label.x = 8, label.y = 2) +
  scale_x_continuous(limits = c(0, 19),
  breaks = seq(from = 0, to = 19, by = 1)) +
  scale_y_continuous(limits = c(0, 19),
  breaks = seq(from = 0, to = 19, by = 1))

ggsave(plot = ddpcr_ngs_plot, 
       filename = "ddpcr_vs_ngs.tiff",
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff')

#########################
# Service numbers 
#########################

# Samples tested
nrow(ngs_vs_ddpcr)

# Assays tested
length(unique(mosaic_data_wider$assay))

# Genes
length(unique(mosaic_data_wider$assay_gene))

#########################
# Selecting analysis wells
#########################

# Analysis wells: any well where an assay was under optimal/standard conditions.
# All wells on a single temp worksheet are included. For optimisation worksheets,
# only wells as the optimum annealing temperature are included.

mosaic_analysis_data <- mosaic_data_wider %>%
  filter(worksheet_well_sample %in% 
           analysis_wells$worksheet_well_sample) %>%
  left_join(analysis_wells %>%
              select(worksheet_well_sample, identity),
            by = "worksheet_well_sample") 

# Excluding samples
# %>% filter(!sample %in% c("21RG-278G0061", "21RG-343G0152",
# "22RG-004G0015"))


#########################
# Variant concentrations
#########################

variant_conc_plot <- mosaic_analysis_data %>%
  arrange(identity, variant_copies_per_ul) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  ggplot(aes(x = worksheet_well_sample,
             y = variant_copies_per_ul)) +
  scale_fill_manual(values = c("#FFFFFF", "#999999", "#333333")) +
  geom_point(pch = 21, aes(fill = identity),
             colour = "black",
             size = 2) +
  geom_errorbar(aes(ymin = variant_poisson_conc_min, 
                    ymax = variant_poisson_conc_max)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "Variant concentration (copies/ul)",
       title = "Variant concentation in ddPCR mosaic data") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 100),
                     breaks = c(0, 0.4, 1, 10, 100)) +
  geom_hline(yintercept = 0.4, linetype = "dashed")


ggsave(plot = variant_conc_plot, 
       filename = "variant_conc_plot.tiff",
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff')

#########################
# Variant fractional abundances
#########################

variant_fraction_plot <- mosaic_analysis_data %>%
  mutate(variant_fractional_abundance = ifelse(is.na(variant_fractional_abundance),
                                               0, variant_fractional_abundance)) %>%
  arrange(identity, variant_fractional_abundance) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  filter(identity != "NTC") %>%
  ggplot(aes(x = worksheet_well_sample,
             y = variant_fractional_abundance)) +
  scale_fill_manual(values = c("#FFFFFF", "#999999", "#333333")) +
  geom_point(pch = 21, aes(fill = identity),
             colour = "black",
             size = 2) +
  geom_errorbar(aes(ymin = variant_poisson_fractional_abundance_min, 
                    ymax = variant_poisson_fractional_abundance_max)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "Variant percent (%)",
       title = "Variant percentage in ddPCR mosaic data") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 10),
                     breaks = c(0, 0.1, 1, 5, 10)) +
  geom_hline(yintercept = 0.1, linetype = "dashed")

ggsave(plot = variant_fraction_plot, 
       filename = "variant_fraction_plot.tiff",
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff')

#########################
# Copies per well
#########################

variant_copies_plot <- mosaic_analysis_data %>%
  arrange(identity, variant_copies_per_20ul_well) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  ggplot(aes(x = worksheet_well_sample,
             y = variant_copies_per_20ul_well)) +
  scale_fill_manual(values = c("#FFFFFF", "#999999", "#333333")) +
  geom_point(pch = 21, aes(fill = identity),
             colour = "black",
             size = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "Variant molecules per well",
       title = "Variant molecules in ddPCR mosaic data") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 6000),
                     breaks = c(0, 10, 100, 1000, 6000)) +
  geom_hline(yintercept = 10, linetype = "dashed")

ggsave(plot = variant_copies_plot, 
       filename = "variant_copies_plot.tiff",
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff')

# Reference copies per well
# 25ng gDNA = 25000pg
# 3.3pg per haploid genome
# 25000 / 3.3 = 7575 copies per well

reference_copies_plot <- mosaic_analysis_data %>%
  ggplot(aes(x = identity, y = reference_copies_per_20ul_well))+
  geom_boxplot() +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "Reference molecules per well",
       title = "Reference molecules in ddPCR mosaic data") +
  geom_hline(yintercept = 7575, linetype = "dashed")

ggsave(plot = reference_copies_plot, 
       filename = "reference_copies_plot.tiff",
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff')

#########################
# All FAM+ droplets
#########################

fam_positive_plot <- mosaic_analysis_data %>%
  arrange(identity, variant_positives) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  ggplot(aes(x = worksheet_well_sample,
             y = variant_positives)) +
  scale_fill_manual(values = c("#FFFFFF", "#999999", "#333333")) +
  geom_point(pch = 21, aes(fill = identity),
             colour = "black",
             size = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "All FAM+ droplets",
       title = "FAM+ droplets in ddPCR wells") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 10000),
                     breaks = c(0, 10, 100, 1000, 10000)) +
  geom_hline(yintercept = 10, linetype = "dashed")

ggsave(plot = fam_positive_plot, 
       filename = "fam_positive_plot.tiff",
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff')

#########################
# FAM+VIC- droplets 
#########################

fam_only_plot <- mosaic_analysis_data %>%
  # "variant_ch1_ch2_2" = Ch1+Ch2-
  arrange(identity, variant_ch1_ch2_2) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  ggplot(aes(x = worksheet_well_sample,
                                 y = variant_ch1_ch2_2)) +
  scale_fill_manual(values = c("#FFFFFF", "#999999", "#333333")) +
  geom_point(pch = 21, aes(fill = identity),
             colour = "black",
             size = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "FAM+VIC- droplets",
       title = "FAM+VIC- droplets in ddPCR wells") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 10000),
                     breaks = c(0, 10, 100, 1000, 10000)) +
  geom_hline(yintercept = 10, linetype = "dashed")

ggsave(plot = fam_only_plot, 
       filename = "fam_only_plot.tiff",
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff')

#########################
# Values for validation document "specificity" table
#########################

# Updated well counts
mosaic_analysis_data %>%
  group_by(identity) %>%
  summarise(count = n())

#########################

## Patient cases with low droplet counts

data_inspection <- mosaic_analysis_data %>%
  filter(identity == "patient") %>%
  arrange(variant_copies_per_ul) %>%
  select(worksheet, well, sample,
          assay, variant_copies_per_ul, variant_copies_per_20ul_well,
          variant_fractional_abundance, variant_positives,
          variant_ch1_ch2_2)



