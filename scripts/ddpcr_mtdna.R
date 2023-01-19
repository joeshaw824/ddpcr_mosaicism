################################################################################
## ddPCR for mitochondrial depletion
## Joseph.Shaw@gosh.nhs.uk / joseph.shaw3@nhs.net
################################################################################

#########################
# Load libraries
#########################

library(tidyverse)
library(readxl)
library(lubridate)

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
  # filter(accepted_droplets > 10000) %>%
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
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1)) +
    labs(title = paste0(specimen, ": B2M Concentration (copies/ul)"), x = "")
  
  nd4_plot <- mtdna_cleaned %>%
    filter(specimen_id == specimen & target_clean == "ND4") %>%
    ggplot(aes(x = worksheet_well, y = copies_per_ul)) +
    geom_point() +
    geom_errorbar(aes(ymin = poisson_conf_min, ymax = poisson_conf_max)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1)) +
    labs(title = paste0(specimen,": ND4 Concentration (copies/ul)"), x = "")
  
  nd1_plot <- mtdna_cleaned %>%
    filter(specimen_id == specimen & target_clean == "ND1") %>%
    ggplot(aes(x = worksheet_well, y = copies_per_ul)) +
    geom_point() +
    geom_errorbar(aes(ymin = poisson_conf_min, ymax = poisson_conf_max)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1)) +
    labs(title = paste0(specimen,": ND1 Concentration (copies/ul)"), x = "")
  
  mtdna_plots <- ggpubr::ggarrange(b2m_plot, nd1_plot, nd4_plot,
                                   ncol = 1, nrow = 3, align = "v")
  ggsave(plot = mtdna_plots, 
         filename = paste0(specimen, "_plots_",
                           format(Sys.time(), "%Y%m%d_%H%M%S"),
                           ".tiff"),
         path = "plots/mtdna_plots/", device='tiff', dpi=100,
         units = "in",
         width = 7,
         height = 12.5)
  
}

repeated_samples <- c("21RG-048G0053", "21RG-118G0184", "21RG-119G0051",
                      "21RG-110G0081", "21RG-118G0011", "21RG-118G0019",
                      "21RG-117G0167", "21RG-119G0052", "21RG-118G0185",
                      "21RG-119G0059", "21RG-162G0048", "21RG-175G0105",
                      "20RG-289G0083", "21RG-130G0082", "21RG-267G0027")

for (specimen in repeated_samples) {
  
  make_mtdna_plots(specimen)
  
}

#########################
# Ammount of DNA input
#########################

mtdna_cleaned %>%
  filter(target_clean == "B2M" & worksheet %in% c("22-1854", "22-2113", "22-2399", "22-2608") &
           !is.na(specimen_id)) %>%
  # Haploid nuclear human genome is 3.3pg
  mutate(ng_input = (copies_per20u_l_well*3.3) / 1000) %>%
  ggplot(aes(x = worksheet_well, y = ng_input)) +
  geom_point() +
  facet_wrap(~specimen_id) +
  theme(axis.text.x = element_blank())
  
#########################
# Sample dilution factors
#########################

# 2 samples on 22-2399 (21RG-118G0011 and 21RG-162G0048) have a dilution factor of 250, not 500

dilution_factors <- read.csv("resources/mtdna_dilution_factors.csv")

#########################
# Patient ages
#########################

# To add

# Best option is to get the information from Winpath as they all have "J" numbers

patient_ages <- read_excel(path = "I:/Regional Genetics Service/Validation Documents/Disease-specific validation/Miscellaneous/Mitochondrial depletion via ddPCR/Run Data & Worksheets/FINAL_ALL_COMB_DATA_model.xlsx",
                           sheet = "allCOMB_data_YP") %>%
  janitor::clean_names()

non_date_values <- c(".", NA, "41771", "218ng/ul", "205ng/ul", "192ng/ul", "626 days", "13 days", "3364 days")

unique(patient_ages$biop)
"3/02/17"

patient_ages_clean <- patient_ages %>%
  filter(!bd %in% non_date_values & !biop %in% non_date_values) %>%
  mutate(clean_biopsy_date  = case_when(
    biop == "3/02/17" ~"03/02/2017",
    biop == "13/02/12" ~"13/02/2012",
    biop == "17/07/19" ~"17/07/2019",
    TRUE ~biop)) %>%
  dplyr::rename(birth_date = bd) %>%
  mutate(birth_date2 = as.Date(birth_date, format = "%d/%m/%Y"),
         biopsy_date2 = as.Date(clean_biopsy_date, format = "%d/%m/%Y"),
         age_days = biopsy_date2 - birth_date2) %>%
  select(specimin_id, birth_date2, biopsy_date2, age_days)

patient_age_comparison <- patient_ages_clean %>%
  left_join(patient_ages %>%
              select(specimin_id, days),
            by = "specimin_id")

#########################
# Definitions
#########################

# b2m_adjusted_nd1_cn: the dilution-adjusted concentration of ND1 mitochondrial DNA
# mt_cn_nd1: the number of mitochondrial genomes per diploid nuclear genome

#########################
# Perform mtDNA versus gDNA calculations
#########################

recent_ws <- c("22-1854", "22-2113", "22-2399", "22-2608")

recent_data <- mtdna_cleaned %>%
  filter(worksheet %in% recent_ws & substr(well, 1, 1) == "M") %>%
  select(worksheet, well, specimen_id, target_clean, copies_per_ul, 
         poisson_conf_min, poisson_conf_max) %>%
  # remove NTC wells
  filter(!is.na(specimen_id) & !is.na(target_clean)) 

recent_data_wider <- recent_data %>%
  # Pivot wider onto one line
  pivot_wider(id_cols = c(worksheet, specimen_id),
              names_from = target_clean,
              values_from = c(copies_per_ul, poisson_conf_min, poisson_conf_max)) %>%
  # Add dilution factors
  left_join(dilution_factors, by = c("worksheet", "specimen_id"))

recent_data_calc <- recent_data_wider %>%
  # Calculations according to Yogen's spreadsheet
  mutate(dilution_factor = b2m_input_ng/nd1_nd4_input_ng,
         # I've kept the names of these variables consistent with Yogen's spreadsheet
         b2m_adjusted_nd1_cn_joe = copies_per_ul_ND1 * dilution_factor,
         b2m_adjusted_nd4_cn_joe = copies_per_ul_ND4 * dilution_factor,
         # Calculation from Wachsmuth et al PMID: 26978189
         mt_cn_nd1_joe = (2*b2m_adjusted_nd1_cn_joe) / copies_per_ul_B2M,
         mt_cn_nd4_joe = (2*b2m_adjusted_nd4_cn_joe) / copies_per_ul_B2M)

#########################
# Load Yogen's results
#########################

consistent_columns <- c("worksheet", "sample", "specimin_id", "b2m_cn_merged", "nd1_cn_merged", "nd4_cn_merged",
                        "b2m_adjusted_nd1_cn", "b2m_adjusted_nd4_cn", "mt_cn_nd1", "mt_cn_nd4")

read_yogen_ws <- function(worksheet_name) {
  
  filepath <- "W:/MolecularGenetics/Neurogenetics/mtDNA/ddPCR/"
  
  output <- read_excel(path = paste0(filepath, worksheet_name, "/", worksheet_name, ".xlsm"),
                       sheet = "mtDNA % calc", skip = 1) %>%
    janitor::clean_names() %>%
    mutate(worksheet = worksheet_name) %>%
    select(all_of(consistent_columns))
  
}

ws_22_1854 <- read_yogen_ws("22-1854")

ws_22_2399 <- read_yogen_ws("22-2399")

ws_22_2608 <- read_yogen_ws("22-2608")


# Have to do individually because of non-standard name
filepath <- "W:/MolecularGenetics/Neurogenetics/mtDNA/ddPCR/"

ws_22_2113 <- read_excel(path = paste0(filepath, "22-2113/22-2113_yp.xlsm"),
                         sheet = "mtDNA % calc", skip = 1) %>%
  janitor::clean_names() %>%
  mutate(worksheet = "22-2113") %>%
  select(all_of(consistent_columns))

yogen_results <- rbind(ws_22_1854, ws_22_2399, ws_22_2608, ws_22_2113) %>%
  mutate(specimen_id = str_extract(specimin_id, pattern = "..RG-...G....")) %>%
  filter(!is.na(specimen_id))

#########################
# Check against Yogen's results
#########################

result_comparison <- recent_data_calc %>%
  full_join(yogen_results, by = c("worksheet", "specimen_id")) %>%
  mutate(nd1_adjusted_diff = abs(b2m_adjusted_nd1_cn - b2m_adjusted_nd1_cn_joe),
         nd4_adjusted_diff = abs(b2m_adjusted_nd4_cn - b2m_adjusted_nd4_cn_joe),
         nd1_cn_diff = round(abs(mt_cn_nd1 - mt_cn_nd1_joe), 0),
         nd4_cn_diff = round(abs(mt_cn_nd4 - mt_cn_nd4_joe), 0)) %>%
  select(worksheet, specimen_id, 
         nd1_adjusted_diff, b2m_adjusted_nd1_cn, b2m_adjusted_nd1_cn_joe,
         nd4_adjusted_diff, b2m_adjusted_nd4_cn, b2m_adjusted_nd4_cn_joe,
         nd1_cn_diff, mt_cn_nd1, mt_cn_nd1_joe,
         nd4_cn_diff, mt_cn_nd4, mt_cn_nd4_joe)

# There is a discrepancy for 21RG-110G0081 ND1 calculations on 22-2113, because Yogen didn't use merged values
# for this sample.

# Check results graphically
ggplot(result_comparison, aes(x = b2m_adjusted_nd1_cn, b2m_adjusted_nd1_cn_joe)) +
  geom_point() +
  labs(title = "Adjusted ND1 calculation comparison")

ggplot(result_comparison, aes(x = b2m_adjusted_nd4_cn, b2m_adjusted_nd4_cn_joe)) +
  geom_point()+
  labs(title = "Adjusted ND4 calculation comparison")

ggplot(result_comparison, aes(x = mt_cn_nd1, mt_cn_nd1_joe)) +
  geom_point()+
  labs(title = "Mitochondrial genomes per diploid nuclear genome - comparison - ND1")

ggplot(result_comparison, aes(x = mt_cn_nd4, mt_cn_nd4_joe)) +
  geom_point()+
  labs(title = "Mitochondrial genomes per diploid nuclear genome - comparison - ND4")

#########################
# Generate plots
#########################

# Generate plots of Yogen's results

recent_samples <- unique(recent_data$specimen_id)

make_cn_plot <- function(specimen) {
  
  test_df <- result_comparison %>%
    filter(specimen_id == specimen) %>%
    select(worksheet, specimen_id, mt_cn_nd1, mt_cn_nd4) %>%
    pivot_longer(cols = c(mt_cn_nd1, mt_cn_nd4),
                 names_to = "target_conc",
                 values_to = "concentration") %>%
    mutate(target = substr(target_conc, 7, 9),
           worksheet_target = paste0(worksheet, sep = "_", target))
  
  y_upper <- max(test_df$concentration)
  
  mtdna_cn_plot <- ggplot(test_df, aes(x = worksheet_target, y = concentration, colour = target)) +
    geom_point(size = 2) +
    ylim(0, y_upper+100) +
    labs(x = "Worksheet", y = "Mitochondrial genome copies per diploid nuclear genome",
         title = paste0(specimen, ": Mitochondrial genome copies per diploid nuclear genome (mt_cn value)"))
  
  ggsave(plot = mtdna_cn_plot, 
         filename = paste0(specimen, "_mtdna_cn_plot_",
                           format(Sys.time(), "%Y%m%d_%H%M%S"),
                           ".tiff"),
         path = "plots/mtdna_plots/cn_plots/", device='tiff', dpi=100,
         units = "in",
         width = 10,
         height = 7)
  
}

for (specimen in recent_samples) {
  
  make_cn_plot(specimen)
  
}

#########################