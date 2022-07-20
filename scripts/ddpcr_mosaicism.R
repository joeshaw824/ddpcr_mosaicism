################################################################################
## ddPCR for confirmation of mosaicism
## Joseph.Shaw@gosh.nhs.uk / joseph.shaw3@nhs.net
################################################################################

# Load data
source("scripts/load_ddpcr_data.R")

#########################
# Get resources
#########################

mosaicism_targets <- read_excel("ddpcr_mosaicism/resources/mosaic_targets_v2.xlsx")
# Add check to make sure assay_id is unique

#mosaicism_assays <- mosaicism_targets %>%
  #pivot_wider(id_cols = assay,
              #names_from = target_category,
              #values_from = c(target, fluorophore))

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

ws_assay_ids <- read_csv("ddpcr_mosaicism/resources/worksheet_assay_ids.csv") %>%
  select(worksheet, assay_id) %>%
  inner_join(targets_rearranged, by = "assay_id",
             keep = FALSE)

#########################
# Add targets to previous data 
#########################

# These worksheets were run without the "assay id" value in the "experiment" column of the Quantasoft input.

previous_worksheets <- c("21-0374", "21-2298", "21-3894", "21-4270", "21-4327", "21-4391", "21-4435", 
                    "21-4459", "22-0227", "22-0271", "22-0873", "22-1066", "22-1395")

previous_data <- ddpcr_mosaic_data %>%
  filter(worksheet %in% previous_worksheets) %>%
  left_join(ws_assay_ids,
            by = c("worksheet", "target")) %>%
  # Change column order to allow rbind in next step
  dplyr::relocate(assay_id, .before = sample) %>%
  select(-experiment)

#########################
# Add targets to recent data 
#########################

new_data <- ddpcr_mosaic_data %>%
  filter(!worksheet %in% previous_worksheets) %>%
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
         assay_gene = sub("c.*", "", assay_gene))

write.csv(mosaic_ddpcr_db,
          file = paste0("ddpcr_mosaicism/database/mosaic_ddpcr_db_", 
                        format(Sys.time(), "%Y%m%d"), ".csv"),
          row.names = FALSE)

#########################
# NGS vs ddPCR 
#########################

# Change so that this step does not rely on the manually-curated csv. Exclude all the
# wells from samples with merged wells.

ngs_results <- read_csv("ddpcr_mosaicism/resources/ngs_results.csv") %>%
  janitor::clean_names() %>%
  mutate(ngs = ifelse(is.na(mosaic_miner_vaf), rc_vaf,
                              mosaic_miner_vaf),
         ngs_percent = ngs*100,
         mosaic_miner_percent = mosaic_miner_vaf*100) %>%
  select(sample, mosaic_miner_vaf, rc_vaf, ngs_percent,
         mosaic_miner_percent)

ngs_vs_ddpcr <- mosaic_ddpcr_db %>%
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

# Analysis wells: any well where an assay was under optimal/standard conditions.
# All wells on a single temp worksheet are included. For optimisation worksheets,
# only wells as the optimum annealing temperature are included.

mosaic_analysis_data <- mosaic_ddpcr_db %>%
  filter(worksheet_well_sample %in% 
           analysis_wells$worksheet_well_sample) %>%
  left_join(analysis_wells %>%
              select(worksheet_well_sample, identity),
            by = "worksheet_well_sample") 

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
  geom_hline(yintercept = 0.4, linetype = "dashed") +
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
       title = "Variant percentages under 5% in ddPCR mosaic data") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 5),
                     breaks = c(0, 0.1, 0.3, 1, 5)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
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

# "22-0227_H07_21RG-278G0061" - NF1 insertion - potential failed assay

# "22-1678_B03_21RG-343G0152" and "22-1678_C03_21RG-343G0152" - matched blood from a skin 

# sample (21RG-326G0125) with a confirmed PIK3CA variant.

# "22-2490_B03_22RG-097G0042" and "22-2490_A03_22RG-097G0042" - GNAS c.602GA - not confirmed on ddPCR or NIPD panel

# "22-2490_C03_19RG-178G0080" and "22-2490_D03_19RG-178G0080" - paired blood sample for 22RG-097G0042 - not confirmed on NIPD

# "22-0271_H04_22RG-004G0015" and "22-0271_G04_22RG-004G0015" and "22-0271_B05_20RG-209G0042" and
# "22-0271_A05_20RG-209G0042" - blood and skin for same patient - 22RG-004G0015 
# had a confirmed NF1 SNV at 1%

not_reported <- c("22-0227_H07_21RG-278G0061", "22-1678_B03_21RG-343G0152", "22-1678_C03_21RG-343G0152",
                  "22-2490_B03_22RG-097G0042", "22-2490_A03_22RG-097G0042", "22-2490_C03_19RG-178G0080",
                  "22-2490_D03_19RG-178G0080", "22-0271_H04_22RG-004G0015", "22-0271_G04_22RG-004G0015",
                  "22-0271_B05_20RG-209G0042", "22-0271_A05_20RG-209G0042")

reported_patient_data <- mosaic_analysis_data %>%
  arrange(identity, variant_copies_per_20ul_well) %>%
  filter(!worksheet_well_sample %in% not_reported ) %>%
  mutate(worksheet_well_sample = factor(worksheet_well_sample,
                                        levels = c(worksheet_well_sample)))

length(unique(reported_patient_data$assay_id))

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

# Samples tested
length(unique(patientsonly$sample))

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


max(reportedpatientsonly$variant_fractional_abundance, na.rm = TRUE)
min(reportedpatientsonly$variant_fractional_abundance, na.rm = TRUE)
median(reportedpatientsonly$variant_fractional_abundance, na.rm = TRUE)




mosaicism_targets %>%
  mutate(assay_gene = sub("_.*", "", channel1_target),
         assay_gene = sub("c.*", "", assay_gene)) %>%
  group_by(assay_gene) %>%
  summarise(total = n()) %>%
  arrange(desc(total))

patients_per_assay <- mosaic_analysis_data %>%
  filter(!base::duplicated(sample)) %>%
  group_by(assay_id) %>%
  summarise(patients_per_assay = n()) %>%
  arrange(desc(patients_per_assay)) %>%
  left_join(mosaicism_targets %>%
              select(assay_id, assay_name),
            by = "assay_id") %>%
  select(assay_id, assay_name, patients_per_assay)

nrow(patients_per_assay %>%
       filter(patients_per_assay == 1))


normals_only <- mosaic_analysis_data %>%
  filter(identity =="normal")

max(normals_only$variant_fractional_abundance, na.rm = TRUE)

colnames(normals_only)

ggplot(normals_only, aes(x = identity, y = variant_copies_per_20ul_well)) +
  geom_jitter() 



patients_per_assay <- mosaic_analysis_data %>%
  filter(identity == "patient") %>%
  group_by(assay_id, assay_name) %>%
  summarise(patients_tested = length(unique(sample))) %>%
  arrange(desc(patients_tested))

write.csv(patients_per_assay, "ddpcr_mosaicism/database/patients_per_assay.csv",
          row.names = FALSE)

mosaic_analysis_data %>%
  filter(assay_name == "GNAQ_c.548GA") %>%
  ggplot(aes(x = reference_threshold, y = variant_threshold)) +
  geom_point()


# Row E is 59 degrees
# How many assays wouldn't work at 59 degrees? I.e. if we just ran at 59 to start with?
# I.e. didn't do optimisation worksheets?

