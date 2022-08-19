################################################################################
## ddPCR for confirmation of mosaicism
## Joseph.Shaw@gosh.nhs.uk / joseph.shaw3@nhs.net
################################################################################

setwd("//fsdept/deptdata$/Regional Genetics Service/Validation Documents/Mosaic/ddPCR/ddpcr_mosaicism/")

# Load data
source("scripts/load_ddpcr_data.R")

# 19/08/2022 The idea is to move away from maintaining resource csv files (they are too time-consuming
# and error-prone) and keep the majority of filtering steps visible in this script.

## TO DO: go through data and identify failed assays, patients used as normal controls etc
# Probably best to do in one go when I have lots of time.
# Plot jitter plots

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

controls <- read_csv("resources/ddpcr_mosaic_controls.csv")

control_ids <- c(controls$specimen_id, controls$control_code)

mosaic_worksheets <- read_csv("resources/mosaic_ddpcr_worksheets.csv")

single_temp_worksheets <- mosaic_worksheets %>%
  filter(category == "single_temp")

gradient_temp_worksheets <- mosaic_worksheets %>%
  filter(category == "gradient_temp")

failed_assays <- c("ANT2RMD", "ANPR3GH", "ANXG9C4", "ANPR4AR", "ANT2RFN")

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

#########################
# old resources
#########################

analysis_wells <- read_csv(
  "resources/mosaicism_analysis_wells.csv") %>%
  mutate(worksheet_well_sample = paste(worksheet, well, sample, 
                                       sep = "_"),
         identity = factor(identity, levels = c("NTC",
                                                "normal",
                                                "patient")))

wells_for_ngs_comparison <- read_csv(
  "resources/ddpcr_wells_ngs_comparison.csv") %>%
  mutate(worksheet_well_sample = paste(worksheet, well, sample, 
                                sep = "_"))

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
  filter(!worksheet %in% previous_worksheets) %>%
  dplyr::rename(assay_id = experiment) %>%
  left_join(targets_rearranged, 
            by = c("assay_id", "target")) %>%
  # Remove non-mosaic (mtDNA) data
  filter(!is.na(assay_name))

data_with_targets <- rbind(previous_data, new_data)

#########################
# Samples not reported 
#########################

not_reported <- c(# matched blood from a skin sample (21RG-326G0125) with a confirmed PIK3CA variant.
                  "22-1678_B03_21RG-343G0152", "22-1678_C03_21RG-343G0152",
                  
                  # GNAS c.602GA - not confirmed on ddPCR or NIPD panel
                  "22-2490_B03_22RG-097G0042", "22-2490_A03_22RG-097G0042", 
                  
                  # Paired blood sample for 22RG-097G0042 - not confirmed on NIPD
                  "22-2490_C03_19RG-178G0080", "22-2490_D03_19RG-178G0080", 
                  
                  # Blood and skin for same patient - 22RG-004G0015 
                  # had a confirmed NF1 SNV at 1%
                  "22-0271_H04_22RG-004G0015", "22-0271_G04_22RG-004G0015",
                  "22-0271_B05_20RG-209G0042", "22-0271_A05_20RG-209G0042",
                  
                  # BRAF variant not detected
                  "22-2630_A03_21RG-138G0065", "22-2630_B03_21RG-138G0065",
                  
                  # GNAS c601 not detected
                  "21-4459_F03_20RG-209G0042", "21-4459_E03_20RG-209G0042",
                  
                  # GNAS c602 not confirmed
                  "22-2630_A07_22RG-097G0042", "22-2630_B07_22RG-097G0042")

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
         identity = factor(identity, levels = c("NTC",
                                                "normal",
                                                "patient")))

write.csv(mosaic_ddpcr_db,
          file = paste0("database/mosaic_ddpcr_db_", 
                        format(Sys.time(), "%Y%m%d"), ".csv"),
          row.names = FALSE)

#########################
# NGS vs ddPCR 
#########################

# Use sample_worksheet as a filter

ngs_vs_ddpcr <- mosaic_ddpcr_db %>%
  # Add section to catch samples which only had 1 well on a single_temp worksheet
  filter((worksheet %in% single_temp_worksheets$worksheet & substr(well, 1, 1) == "M") |
           (worksheet %in% gradient_temp_worksheets$worksheet & substr(well, 1, 1) == "E")) %>%
  inner_join(ngs_results,
              by = "sample") %>%
  mutate(variant_check = ifelse(assay_name == variant_assay, TRUE, FALSE))

# Check data is for the same variant
ngs_vs_ddpcr %>%
  filter(variant_check == FALSE)

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
# Selecting analysis wells
#########################

# Analysis wells: any well where an assay was under standard conditions of 59 degrees annealing temperature.
# All wells on a single temp worksheet are included. For optimisation worksheets,
# only wells at 59 degrees (row E) are included.

mosaic_analysis_data <- mosaic_ddpcr_db %>%
  filter((worksheet %in% single_temp_worksheets$worksheet & substr(well, 1, 1) != "M" | 
         (worksheet %in% gradient_temp_worksheets$worksheet & substr(well, 1, 1) == "E"))
         & !assay_id %in% failed_assays)

#########################
# Variant concentrations
#########################

mosaic_analysis_data %>%
  filter(identity == "patient") %>%
  filter(!worksheet_well_sample %in% not_reported) %>%
  arrange(variant_copies_per_ul) %>%
  select(worksheet_well_sample, variant_copies_per_ul)


variant_conc_plot <- mosaic_analysis_data %>%
  filter(!worksheet_well_sample %in% not_reported) %>%
  arrange(identity, variant_copies_per_ul) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample))) %>%
  ggplot(aes(x = worksheet_well_sample,
             y = variant_copies_per_ul)) +
  scale_fill_manual(values = c("#FFFFFF", "#999999", "#333333")) +
  geom_point(pch = 21, aes(fill = identity),
             colour = "black",
             size = 2) +
  #geom_errorbar(aes(ymin = variant_poisson_conc_min, 
                    #ymax = variant_poisson_conc_max)) +
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
                     breaks = c(0, 1, 10, 100)) +
  # geom_hline(yintercept = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  annotate(geom = "text", x = 100, y = 1.5,
         label = "Variants detected above 1 copy per ul were reported")

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
  #geom_errorbar(aes(ymin = variant_poisson_fractional_abundance_min, 
                    #ymax = variant_poisson_fractional_abundance_max)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_blank()) +
  labs(x = "",
       y = "Variant percent (%)",
       title = "Variant percentages under 5% in ddPCR mosaic data") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 5),
                     breaks = c(0, 0.3, 1, 5)) +
  #geom_hline(yintercept = 0.1, linetype = "dashed") +
  geom_hline(yintercept = 0.3, linetype = "dashed") +
  annotate(geom = "text", x = 80, y = 0.5,
           label = "Variants detected above 0.3% were reported")

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
       title = paste0("Variant molecules in ddPCR mosaic data (", well_number, " wells, ", 
                      assay_number, " assays)")) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 6000),
                     breaks = c(0, 10, 100, 1000, 6000)) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  annotate(geom = "text", x = 120, y = 25,
           label = "Variants with over 20 molecules detected were reported") 

ggsave(plot = variant_copies_plot, 
       filename = paste0("variant_copies_plot_", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 15)

mosaic_analysis_data %>%
  ggplot(aes(x = identity,
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
       title = paste0("Variant molecules in ddPCR mosaic data (", well_number, " wells, ", 
                      assay_number, " assays)")) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 6000),
                     breaks = c(0, 10, 100, 1000, 6000)) +
  facet_wrap(~assay_name)


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
  geom_hline(yintercept = 20, linetype = "dashed")

ggsave(plot = fam_positive_plot, 
       filename = paste0("fam_positive_plot_", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 15)

single_assay_plot <- mosaic_analysis_data %>%
  filter(assay_name == "GNAQ_c.548GA") %>%
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
       title = "GNAQ_c.548GA: FAM+ droplets in ddPCR wells") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 10000),
                     breaks = c(0, 10, 100, 1000, 10000))

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
                     breaks = c(0, 10, 100, 1000, 10000))

ggsave(plot = fam_only_plot, 
       filename = paste0("fam_only_plot_", format(Sys.time(), "%Y%m%d"),".tiff"),
       path = "ddpcr_mosaicism/plots/", 
       device= 'tiff',
       units = "cm",
       width = 20,
       height = 15)

#########################
# BSGM 2022 abstract - Samples reported as "mosaicism confirmed"
#########################

check_this <- mosaic_analysis_data %>%
  filter(identity == "patient") %>%
  arrange(variant_copies_per_20ul_well) %>%
  select(worksheet_well_sample, assay_name, variant_copies_per_20ul_well)



reported_patient_data <- mosaic_analysis_data %>%
  arrange(identity, variant_copies_per_20ul_well) %>%
  filter(!worksheet_well_sample %in% not_reported ) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample)))

reported_patients <- reported_patient_data %>%
  filter(identity == "patient")

length(unique(reported_patients$sample))

ggplot(reported_patient_data, aes(x = worksheet_well_sample,
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
       title = paste0("Variant molecules in ddPCR mosaic data (", well_number, " wells, ", 
                      assay_number, " assays)")) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 6000),
                     breaks = c(0, 10, 100, 1000, 6000)) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  annotate(geom = "text", x = 120, y = 25,
           label = "Variants with over 20 molecules detected were reported")


normals_only <- mosaic_analysis_data %>%
  filter(identity == "normal") 

max(normals_only$variant_copies_per_20ul_well)
median(normals_only$variant_copies_per_20ul_well)

#########################
# Service numbers 
#########################

reportedpatientsonly <- mosaic_analysis_data %>%
  filter(identity == "patient") %>%
  filter(variant_fractional_abundance > 0.29)

patients_only <- mosaic_analysis_data %>%
  filter(identity == "patient") 

# Samples tested
length(unique(patients_only$sample))

# Assays tested
length(unique(patientsonly$assay_id))

# Genes
length(unique(patientsonly$assay_gene))

# Worksheets
length(unique(mosaic_ddpcr_db$worksheet))

# Updated well counts
mosaic_analysis_data %>%
  #filter(assay_name == "GNAQ_c.548GA") %>%
  group_by(identity) %>%
  summarise(count = n())

mosaicism_targets %>%
  mutate(assay_gene = sub("_.*", "", channel1_target),
         assay_gene = sub("c.*", "", assay_gene)) %>%
  group_by(assay_gene) %>%
  summarise(total = n()) %>%
  arrange(desc(total))

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