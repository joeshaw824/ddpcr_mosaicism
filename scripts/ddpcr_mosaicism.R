################################################################################
## ddPCR for confirmation of mosaicism
## Joseph.Shaw@gosh.nhs.uk / joseph.shaw3@nhs.net
################################################################################

setwd("//fsdept/deptdata$/Regional Genetics Service/Validation Documents/Mosaic/ddPCR/ddpcr_mosaicism/")

# Load data
source("scripts/load_ddpcr_data.R")

# 19/08/2022 The idea is to move away from maintaining resource csv files (they are too time-consuming
# and error-prone) and keep the majority of filtering steps visible in this script.

# 22/08/2022 Variant fractional abundance is a relative measure, and therefore not so useful when 
# looking at the dataset as a whole. The better option is values which account for the different 
# number of partitions (droplets) between wells: variant concentration and variant total molecules.

#########################
# Get resources
#########################

mosaicism_targets <- read_excel("resources/mosaic_targets_v2.xlsx",
                                     trim_ws = TRUE)

# Check for duplicate assay ids
if(anyDuplicated(mosaicism_targets$assay_id) > 0){
  print("Error: There are duplicate assay IDs in the target file")
  # Remove file so script breaks easily
  rm(mosaicism_targets)
}

# Add check here to make sure all the assays that have been ordered are in the targets file
# But don't have the script break if not.
tom_spreadsheet <- read_excel("I:/Genetics/DNA Lab/databases/Specialist_Services/Skin/ddPCR_designs_confirmations.xlsx",
                              sheet = "sequences") %>%
  janitor::clean_names()

setdiff(tom_spreadsheet$assay_id, mosaicism_targets$assay_id)

# ANWDA9U is the alternative version of GNAS c.601C>T which was ordered in error.

# Anonymised control gDNAs (samples referred for cystic fibrosis screening which were negative)
control_ids <- c("21RG-333G0100",	"C257", "21RG-333G0101",	"C258", "21RG-333G0107",	"C261", 
                 "21RG-333G0112",	"C264", "22RG-110G0026",	"C270", "22RG-110G0030",	"C272")

single_temp_worksheets <- c("21-4327", "21-4459", "22-0271", "22-1066","22-1678", "22-2490", "22-2630",
                            "22-2987")

gradient_temp_worksheets <- c("21-0374", "21-2298", "21-3894", "21-4270", "21-4391", "21-4435", 
                              "22-0227", "22-0873", "22-1395", "22-1880", "22-2099", "22-2397", "22-2889",
                              "22-3294")

# Check for worksheets on both lists
if(length(base::intersect(single_temp_worksheets, gradient_temp_worksheets)) > 0) {
  print("There are worksheets listed as both single_temp and gradient_temp")
  rm(single_temp_worksheets)
  rm(gradient_temp_worksheets)
}

# Check for worksheets absent from dataset
mosaic_worksheets <- c(single_temp_worksheets, gradient_temp_worksheets)

if(length(setdiff(mosaic_worksheets, ddpcr_data$worksheet)) > 0) {
  print("Not all mosaic worksheets have been loaded")
}

# Check for duplicate rows
if(anyDuplicated(ddpcr_data) > 0) {
  print("There are duplicate rows in the dataset")
}

failed_assays <- c(# RHOA_c.514GA
                    "ANT2RMD", 
                    # ACTB_c.439CT
                    "ANPR3GH", 
                    # HRAS_c.37GC
                    "ANXG9C4", 
                    # PIK3CA_c.1034AT
                    "ANPR4AR", 
                    # NF1_c7863ins
                    "ANT2RFN",
                    # MAP2K1_c.168GT
                    "ANRWXVN",
                    # NF1 indel design by Thermo team
                    "NF1_c.7863_7864ins")

assay_information <- read_excel(
  path = "I:/Genetics/DNA Lab/databases/Specialist_Services/Skin/ddPCR_designs_confirmations.xlsx",
  sheet = "sequences") %>%
  janitor::clean_names()

# Use Tom's spreadsheet to avoid duplication of effort
ngs_results <- read_excel(path = "I:/Genetics/DNA Lab/databases/Specialist_Services/Skin/ddPCR_designs_confirmations.xlsx",
                          sheet = "ddPCR_confirmation_list") %>%
  janitor::clean_names() %>%
  mutate(mosaic_miner_vaf = as.numeric(mosaic_miner_vaf),
         rc_vaf = as.numeric(rc_vaf),
         
         ngs = ifelse(is.na(mosaic_miner_vaf), rc_vaf,
                      mosaic_miner_vaf),
         ngs_percent = ngs*100,
         mosaic_miner_percent = mosaic_miner_vaf*100)

if(anyDuplicated(ngs_results$sample) > 0){
  print("Error: There are duplicate samples in the NGS results")
  rm(ngs_results)
}

# Samples not reported as having mosaicism confirmed

not_confirmed <- c(# Matched blood from a skin sample (21RG-326G0125) with a confirmed PIK3CA variant.
  "22-1678_B03_21RG-343G0152", "22-1678_C03_21RG-343G0152",
  
  # GNAS c.602GA - not confirmed on ddPCR or NIPD panel
  "22-2490_B03_22RG-097G0042", "22-2490_A03_22RG-097G0042", 
  
  # Paired blood sample for 22RG-097G0042 - not confirmed on NIPD
  "22-2490_C03_19RG-178G0080", "22-2490_D03_19RG-178G0080", 
  
  # Blood and skin for same patient - 22RG-004G0015 
  # had a confirmed NF1 SNV at 1%
  "22-0271_H04_22RG-004G0015", "22-0271_G04_22RG-004G0015",
  "22-0271_B05_20RG-209G0042", "22-0271_A05_20RG-209G0042",
  "21-4459_F03_20RG-209G0042", "21-4459_E03_20RG-209G0042",
  
  # BRAF variant not detected
  "22-2630_A03_21RG-138G0065", "22-2630_B03_21RG-138G0065",
  
  # GNAS c602 not confirmed
  "22-2630_A07_22RG-097G0042", "22-2630_B07_22RG-097G0042")

#########################
# Patient demographics
#########################

patient_info <- read_excel("resources/Mosaic_sequencing_referrals_20220915_0805.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::rename(specimen_id = test_specimen_id) %>%
  filter(!base::duplicated(specimen_id))

#########################
# Restructure target file
#########################

# There is probably a cleaner way to do this with pivot_wider but I 
# couldn't work it out.

targets_ch1 <- mosaicism_targets %>%
  select(assay_id, assay_name, channel1_target,
         channel1_category, channel1_fluorophore) %>%
  dplyr::rename(category = channel1_category,
                fluorophore = channel1_fluorophore,
                target = channel1_target) %>%
  mutate(channel = 1) %>%
  select(assay_id, assay_name, channel, target, category, fluorophore)

targets_ch2 <- mosaicism_targets %>%
  select(assay_id, assay_name, channel2_target,
         channel2_category, channel2_fluorophore) %>%
  dplyr::rename(category = channel2_category,
                fluorophore = channel2_fluorophore,
                target = channel2_target) %>%
  mutate(channel = 2) %>%
  select(assay_id, assay_name, channel, target, category, fluorophore)

targets_rearranged <- rbind(targets_ch1, targets_ch2)

#########################
# Targets with assay ids
#########################

# Probably the best way to join up the old worksheets is with a csv which states the 
# worksheet, assay_id and assay_name.
# Not that neat but it will get all the data in the same format.

ws_assay_ids <- read_csv("resources/worksheet_assay_ids.csv") %>%
  select(worksheet, assay_id) %>%
  inner_join(targets_rearranged, by = "assay_id",
             keep = FALSE)

#########################
# Add targets to previous data 
#########################

# These worksheets were run without the "assay id" value in the "experiment" column of the Quantasoft input.

previous_worksheets <- c("21-0374", "21-2298", "21-3894", "21-4270", "21-4327", "21-4391", "21-4435", 
                    "21-4459", "22-0227", "22-0271", "22-0873", "22-1066", "22-1395")

previous_data <- ddpcr_data %>%
  filter(worksheet %in% previous_worksheets) %>%
  left_join(ws_assay_ids,
            by = c("worksheet", "target")) %>%
  # Change column order to allow rbind in next step
  dplyr::relocate(assay_id, .before = sample) %>%
  select(-experiment)

#########################
# Add targets to recent data 
#########################

new_data <- ddpcr_data %>%
  filter(worksheet %in% mosaic_worksheets & !worksheet %in% previous_worksheets) %>%
  dplyr::rename(assay_id = experiment) %>%
  left_join(targets_rearranged, 
            by = c("assay_id", "target")) %>%
  # Remove non-mosaic (mtDNA) data
  filter(!is.na(assay_name))

data_with_targets <- rbind(previous_data, new_data)

#########################
# Modify ddPCR mosaic data 
#########################

mosaic_data_mod <- data_with_targets %>%
  dplyr::rename(total_droplets = accepted_droplets,
                # "total_conf_max/min" refer to the "concentration" field.
                # "total_conc_max" is the maximum value of the concentration with total error.
                # I renamed these fields to be consistent with the naming convention in
                # other columns.
                total_conc_max = total_conf_max,
                total_conc_min = total_conf_min,
                poisson_conc_max = poisson_conf_max,
                poisson_conc_min = poisson_conf_min,
                copies_per_20ul_well = copies_per20u_l_well) %>%
  dplyr::rename(target_category = category)

mosaic_ddpcr_db <- mosaic_data_mod %>% 
  filter(!is.na(target_category)) %>%
  pivot_wider(id_cols = c(worksheet_well_sample, worksheet, well, sample,
                          assay_name, assay_id),
              names_from = target_category,
              # Include all values in case we need them later
              values_from = -c(worksheet_well_sample, worksheet, well, sample,
                               assay_name, assay_id),
              # Use names_glue to keep new columns names with naming
              # convention
              names_glue = "{target_category}_{.value}") %>%
  mutate(assay_gene = sub("_.*", "", assay_name),
         assay_gene = sub("c.*", "", assay_gene),
         assay_id_name = paste0(assay_name, " ", assay_id),
         
         # Add identity
         identity = case_when(
           sample == "NTC" ~"NTC",
           sample %in% control_ids ~"normal",
           
           # Some patient samples used as controls for other assays
           sample %in% c("20RG-209G0042", "21RG-074G0164", "20RG-337G0031") & worksheet == "21-3894" ~"normal",
           
           sample %in% c("20RG-209G0042", "21RG-074G0164", "20RG-337G0031", "20RG-153G0082") & worksheet == "21-4270"  ~"normal",
           
           sample %in% c("20RG-209G0042", "21RG-074G0164", "20RG-337G0031", "20RG-153G0082") & worksheet == "21-4327"  ~"normal",
           
           sample %in% c("20RG-209G0042", "21RG-074G0164", "20RG-337G0031", "20RG-153G0082") & worksheet == "21-4391"  ~"normal",
           
           sample == "20RG-209G0042" & worksheet == "21-2298"  ~"normal",
           
           sample == "20RG-260G0066" & worksheet == "21-0374" ~"normal",
           
           TRUE ~"patient"),
         
         # Add factor levels for plots
         identity = factor(identity, levels = c("NTC", "normal", "patient")),
         
         status = case_when(
           identity == "NTC" ~"NTC",
           identity == "normal" ~"normal",
           identity == "patient" & !worksheet_well_sample %in% not_confirmed ~"patient - mosaicism confirmed",
           identity == "patient" & worksheet_well_sample %in% not_confirmed ~"patient - not confirmed"),
         
         status = factor(status, levels = c("NTC", "normal", "patient - not confirmed",
                                            "patient - mosaicism confirmed"))) %>%
  filter(sample != "G_block") %>%
  filter(!assay_id %in% failed_assays)

controls_to_check <- mosaic_ddpcr_db %>%
  filter(identity == "normal" & !sample %in% control_ids) %>%
  filter(!sample %in% c("20RG-209G0042", "21RG-074G0164", "20RG-337G0031", "20RG-153G0082",
                        "20RG-260G0066")) %>%
  select(worksheet_well_sample, assay_name, identity)

write.csv(mosaic_ddpcr_db,
          file = paste0("database/mosaic_ddpcr_db_", 
                        format(Sys.time(), "%Y%m%d"), ".csv"),
          row.names = FALSE)

#########################
# Selecting analysis wells
#########################

# Analysis wells: any well where an assay was under standard conditions of 59 degrees annealing temperature.
# All wells on a single temp worksheet are included. For optimisation worksheets,
# only wells at 59 degrees (row E) are included.

mosaic_analysis_data <- mosaic_ddpcr_db %>%
  filter(
    (worksheet %in% single_temp_worksheets & substr(well, 1, 1) != "M") | 
      (worksheet %in% gradient_temp_worksheets & worksheet != "21-4391" & substr(well, 1, 1) == "E") |
      (worksheet == "21-4391" & substr(well, 1, 1) == "D"))

#########################
# NGS vs ddPCR 
#########################

ngs_vs_ddpcr <- mosaic_ddpcr_db %>%
  filter(identity == "patient") %>%
  filter((worksheet %in% single_temp_worksheets & substr(well, 1, 1) == "M") |
           # Samples with only one well tested (so no merged value)
           worksheet_well_sample %in% c("21-4327_C03_21RG-188G0052", "22-1678_A01_21RG-264G0082") |
           (worksheet %in% gradient_temp_worksheets 
            & worksheet != "21-4391"
            & substr(well, 1, 1) == "E") |
           # Worksheet 21-4391 was inverted on the PCR block
           worksheet == "21-4391" & substr(well, 1, 1) == "D") %>%
  # Sample tested twice for 2 variants - GNAS c601CT not detected on NGS
  filter(worksheet_well_sample != "22-0271_M12_22RG-004G0015") %>%
  inner_join(ngs_results,
              by = "sample") %>%
  mutate(variant_check = ifelse(assay_name == variant_assay, TRUE, FALSE))

# Check data is for the same variant
ngs_vs_ddpcr %>%
  select(worksheet_well_sample, assay_name, variant_assay, variant_check)

sample_number <- nrow(ngs_vs_ddpcr)

single_temp_all <- mosaic_ddpcr_db %>%
  filter(worksheet %in% single_temp_worksheets & identity == "patient")

setdiff(single_temp_all$sample, ngs_vs_ddpcr$sample)
  
# "22RG-181G0203" - paired buccal sample - not tested by NGS
# "10G07516" - positive GNAS control
# "21RG-343G0152" - paired blood sample
# "19RG-178G0080" - paired bood sample
# 22RG-123G0141 - paired blood sample
# 20RG-209G0042 - paired blood
# 21RG-131G0259 - BRAF variant

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
       title = paste("ddPCR and NGS results for", sample_number, " patient samples")) +
  ggpubr::stat_cor(method = "pearson", label.x = 8, label.y = 2) +
  scale_x_continuous(limits = c(0, 20),
  breaks = seq(from = 0, to = 20, by = 1)) +
  scale_y_continuous(limits = c(0, 20),
  breaks = seq(from = 0, to = 20, by = 1))

ggsave(plot = ddpcr_ngs_plot, 
       filename = paste0("ddpcr_vs_ngs_", format(Sys.time(), "%Y%m%d"), ".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 20)

#########################
# Plot themes
#########################

plot_theme <- theme(legend.position = "none",
                    panel.grid = element_blank(),
                    legend.title = element_blank())
  
plot_fill <- scale_fill_manual(values = c(
  # white
  "#FFFFFF", 
  # grey
  "#999999", 
  # light blue
  "#99CCFF", 
  # red
  "#FF3333"))
  
plot_jitter <- geom_jitter(pch = 21, aes(fill = status), colour = "black", size = 2)

plot_point <- geom_point(pch = 21, aes(fill = status), colour = "black", size = 2)

#########################
# Variant concentrations
#########################

variant_conc_plot <- mosaic_analysis_data %>%
  ggplot(aes(x = status,
             y = variant_copies_per_ul)) +
  plot_fill +
  plot_jitter +
  theme_bw() +
  plot_theme +
  labs(x = "",
       y = "Variant concentration (copies/ul)",
       title = "Variant concentation in ddPCR mosaic data") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 1, 10, 100)) +
  geom_hline(yintercept = 1, linetype = "dashed")

ggsave(plot = variant_conc_plot, 
       filename = paste0("variant_conc_plot_error_bars", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/",
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 15)

#########################
# Variant fractional abundances
#########################

variant_fraction_plot <- mosaic_analysis_data %>%
  mutate(variant_fractional_abundance = ifelse(is.na(variant_fractional_abundance),
                                               0, variant_fractional_abundance)) %>%
  filter(variant_fractional_abundance < 5) %>%
  ggplot(aes(x = status,
             y = variant_fractional_abundance)) +
  plot_fill +
  plot_jitter +
  theme_bw() +
  plot_theme +
  labs(x = "",
       y = "Variant percent (%)",
       title = "Variant percentages under 5% in ddPCR mosaic data") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 5),
                     breaks = c(0, 0.3, 1, 5)) +
  geom_hline(yintercept = 0.2, linetype = "dashed") 

ggsave(plot = variant_fraction_plot, 
       filename = paste0("variant_fraction_plot_error_bars", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 15)

#########################
# Variant copies per well
#########################

well_number <- nrow(mosaic_analysis_data)
assay_number <- length(unique(mosaic_analysis_data$assay_id))

variant_copies_plot <- mosaic_analysis_data %>%
  ggplot(aes(x = status,
             y = variant_copies_per_20ul_well)) +
  plot_fill +
  plot_jitter +
  theme_bw() +
  plot_theme +
  labs(x = "",
       y = "Variant molecules per well",
       title = paste0("Variant molecules in ddPCR mosaic data (", well_number, " wells, ", 
                      assay_number, " assays)")) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 6000)) +
  geom_hline(yintercept = 20, linetype = "dashed")

ggsave(plot = variant_copies_plot, 
       filename = paste0("variant_copies_plot_", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 15)

#########################
# Reference copies per well
#########################

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
       filename = paste0("reference_copies_plot_", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 15,
       height = 15)

#########################
# All FAM+ droplets
#########################

fam_positive_plot <- mosaic_analysis_data %>%
  ggplot(aes(x = status,
             y = variant_positives)) +
  plot_fill +
  plot_jitter +
  theme_bw() +
  plot_theme +
  labs(x = "",
       y = "All FAM+ droplets",
       title = "FAM+ droplets in ddPCR wells") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     #limits = c(0, 10000),
                     breaks = c(0, 10, 100, 1000, 10000)) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed")

ggsave(plot = fam_positive_plot, 
       filename = paste0("fam_positive_plot_", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 15)

#########################
# FAM+VIC- droplets 
#########################

fam_only_plot <- mosaic_analysis_data %>%
  ggplot(aes(x = status, y = variant_ch1_ch2_2)) +
  plot_fill +
  plot_jitter +
  theme_bw() +
  plot_theme +
  labs(x = "",
       y = "FAM+VIC- droplets",
       title = "FAM+VIC- droplets in ddPCR wells") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     #limits = c(0, 10000),
                     breaks = c(0, 10, 100, 1000, 10000))

ggsave(plot = fam_only_plot, 
       filename = paste0("fam_only_plot_", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 15)

#########################
# Service numbers 
#########################

mosaic_analysis_data %>%
  group_by(status) %>%
  summarise(total = n())

# Assays tested
length(unique(mosaic_analysis_data$assay_id))

# Genes
length(unique(mosaic_analysis_data$assay_gene))

# Worksheets
length(unique(mosaic_analysis_data$worksheet))

# Patients tested
patients_only <- mosaic_analysis_data %>%
  filter(identity == "patient")

length(unique(patients_only$sample))

mosaicism_targets %>%
  mutate(assay_gene = sub("_.*", "", channel1_target),
         assay_gene = sub("c.*", "", assay_gene)) %>%
  group_by(assay_gene) %>%
  summarise(total = n()) %>%
  arrange(desc(total))

patients_with_data <- patient_info %>%
  left_join(mosaic_analysis_data %>%
              filter(!base::duplicated(sample)) %>%
              select(sample, assay_name, variant_fractional_abundance) %>%
              dplyr::rename(specimen_id = sample),
            by = "specimen_id") %>%
  filter(!is.na(assay_name)) %>%
  mutate(age_years = as.numeric(gsub("-year old", "", age)))

# Children under 10 tested
patients_with_data %>%
  filter(age_years <+10) %>%
  select(specimen_id, age_years)

#########################
# Patients per assay 
#########################

patients_per_assay <- mosaic_analysis_data %>%
  filter(identity == "patient") %>%
  group_by(assay_id, assay_name) %>%
  summarise(patients_tested = length(unique(sample))) %>%
  arrange(desc(patients_tested))

write.csv(patients_per_assay, "ddpcr_mosaicism/database/patients_per_assay.csv",
          row.names = FALSE)

#########################
# All assays facet wrap
#########################

facet_plot_theme <- theme(axis.text.x = element_blank(),
                          legend.position = "bottom")

variant_conc_facet <- variant_conc_plot +
  facet_wrap(~assay_id_name) +
  facet_plot_theme

# Some data missing - probably due to NA values
variant_fraction_facet <- variant_fraction_plot +
  facet_wrap(~assay_id_name) +
  facet_plot_theme

variant_copies_facet <- variant_copies_plot +
  facet_wrap(~assay_id_name) +
  facet_plot_theme

positives_facet_plot <- fam_positive_plot +
  facet_wrap(~assay_id_name) +
  facet_plot_theme

#########################
# Do the error bars overlap?
#########################

# Very slightly at the end of the normal range

mosaic_analysis_data %>%
  filter(variant_copies_per_ul > 0.2 & variant_copies_per_ul < 3 &
           status != "patient - not confirmed") %>%
  arrange(status, variant_copies_per_ul) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  ggplot(aes(x = worksheet_well_sample,
             y = variant_copies_per_ul))  +
  scale_fill_manual(values = c(
    # grey
    "#999999",
    # red
    "#FF3333")) +
  plot_point +
  theme_bw() +
  plot_theme +
  geom_errorbar(aes(ymin = variant_poisson_conc_min, ymax = variant_poisson_conc_max)) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(x = "",
       y = "Variant concentration (copies/ul)",
       title = "Variant concentation between normal controls and confirmed mosaic samples") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 1000),
                     breaks = c(0, 1, 10, 100, 1000)) +
  geom_hline(yintercept = 1, linetype = "dashed")

#########################
# Potential G-block contamination
#########################

# Plot for patient with NF1 variant and background GNAS c601CT
mosaic_analysis_data %>%
  filter(sample %in% c(
    # Skin
    "22RG-004G0015", 
    # Blood
    "20RG-209G0042")) %>%
  ggplot(aes(x = status,
             y = variant_positives)) +
  plot_fill +
  plot_jitter +
  theme_bw() +
  plot_theme +
  labs(x = "",
       y = "All FAM+ droplets",
       title = "FAM+ droplets in ddPCR wells") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     #limits = c(0, 10000),
                     breaks = c(0, 10, 100, 1000, 10000)) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  facet_wrap(~assay_id_name) +
  facet_plot_theme
  
#########################