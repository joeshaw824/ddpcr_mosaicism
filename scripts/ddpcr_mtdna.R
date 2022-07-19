################################################################################
## ddPCR for mitochondrial depletion
## Joseph.Shaw@gosh.nhs.uk / joseph.shaw3@nhs.net
################################################################################

#########################
# Load libraries
#########################

library(tidyverse)
library(readxl)

setwd("//fsdept/deptdata$/Regional Genetics Service/Validation Documents/Mosaic/ddPCR/ddpcr_mosaicism")

#########################
# Load data
#########################

source("scripts/load_ddpcr_data.R")

#########################
# Isolate mtDNA data
#########################

mtdna_ws <- c("22-2608", "22-2399", "22-2113", "22-1854",
              "21-1487", "21-1485", "21-1462", "21-1405",
              "21-1353", "21-1256", "21-0961", "20-0940",
              "21-0767", "21-0716", "21-0598", "21-0485",
              "20-0568", "19-5294", "21-1525")

mtdna_data <- ddpcr_data %>%
  filter(worksheet %in% mtdna_ws) %>%
  filter(accepted_droplets > 10000) %>%
  mutate(worksheet_well = paste(worksheet, well, sep = "_"))

#########################
# Tidy up sample names
#########################

clean_names <- mtdna_data %>%
  # Epic sample names are 
  mutate(specimen_id = case_when(
    # Had to clean these ones manually
    sample == "20RG-339G056" ~"20RG-339G0056",
    sample == "20R-339G0056 5ng/ul" ~"20RG-339G0056",
    sample == "21RG048-G0052 ND1" ~"21RG-048G0052",
    sample == "21RG048-G0058 ND1" ~"21RG-048G0058",
    sample == "21RG048-G0053 ND1" ~"21RG-048G0053",
    sample == "21RG048-G0055 ND1" ~"21RG-048G0055",
    sample == "21RG048-G0056 ND1" ~"21RG-048G0056",
    sample == "21RG--014G0029 ND1" ~"21RG-014G0029",
    sample == "21RG--110G0087 ND1" ~"21RG-110G0087",
    sample == "21RG--005G0061 ND1" ~ "21RG-005G0061",
    sample == "21RG--118G0184 ND1" ~"21RG-118G0184",
    # Checked against worksheet
    sample == "21RG-05G0061 B2M" ~"21RG-005G0061",
    sample == "21RG-1104G0081 B2M" ~"21RG-110G0081",
    sample == "21RG-1104G0092 B2M" ~"21RG-110G0092",
    sample == "21RG-1104G0087 B2M" ~"21RG-110G0087",
    sample == "21RG-1104G0089 B2M" ~"21RG-110G0089",
    sample == "21RG-1104G0093 B2M" ~"21RG-110G0093",
    sample == "21RG-1104G0098 B2M" ~"21RG-110G0098",
    # The rest should have regular Epic specimen numbers
    TRUE ~str_extract(sample, pattern = "..RG-...G....")))

# 68 samples tested in total

#########################
# Tidy up target names
#########################

ND1_variants <- c("nd1", "ND1", "ND1 (HEX)")

B2M_variants <- c("B2M", "b2m", "B2M (FAM)", "B2M FAM", "B2M NTC")

mtdna_cleaned <- clean_names %>%
  mutate(target_clean = case_when(
    target %in% ND1_variants ~"ND1",
    target %in% B2M_variants ~"B2M",
    TRUE ~target))

#########################
# Investigate variation
#########################

make_mtdna_plots <- function(specimen) {
  
  b2m_plot <- mtdna_cleaned %>%
    filter(specimen_id == specimen & target_clean == "B2M") %>%
    ggplot(aes(x = worksheet_well, y = copies_per_ul)) +
    geom_point() +
    geom_errorbar(aes(ymin = poisson_conf_min, ymax = poisson_conf_max)) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste0(specimen, ": B2M Concentration (copies/ul)"), x = "")
  
  nd4_plot <- mtdna_cleaned %>%
    filter(specimen_id == specimen & target_clean == "ND4") %>%
    ggplot(aes(x = worksheet_well, y = copies_per_ul)) +
    geom_point() +
    geom_errorbar(aes(ymin = poisson_conf_min, ymax = poisson_conf_max)) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste0(specimen,": ND4 Concentration (copies/ul)"), x = "")
  
  nd1_plot <- mtdna_cleaned %>%
    filter(specimen_id == specimen & target_clean == "ND1") %>%
    ggplot(aes(x = worksheet_well, y = copies_per_ul)) +
    geom_point() +
    geom_errorbar(aes(ymin = poisson_conf_min, ymax = poisson_conf_max)) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste0(specimen,": ND1 Concentration (copies/ul)"), x = "")
  
  mtdna_plots <- ggpubr::ggarrange(b2m_plot, nd4_plot,
                                   nd1_plot,
                                   ncol = 2, nrow = 2, align = "v")
  ggsave(plot = mtdna_plots, 
         filename = paste0(specimen, "_plots_",
                           format(Sys.time(), "%Y%m%d_%H%M%S"),
                           ".tiff"),
         path = "plots/mtdna_plots/", device='tiff', dpi=100,
         units = "in",
         width = 12.5,
         height = 7)
  
}

repeated_samples <- c("21RG-118G0019", "21RG-117G0167", "21RG-119G0052",
                      "21RG-048G0053", "21RG-118G0184", "21RG-119G0051",
                      "21RG-110G0081", "21RG-118G0011")

for (specimen in repeated_samples) {
  
  make_mtdna_plots(specimen)
  
}
 
#########################