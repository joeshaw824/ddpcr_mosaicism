################################################################################
## Quality monitoring
################################################################################

setwd("//fsdept/deptdata$/Regional Genetics Service/Validation Documents/Mosaic/ddPCR/ddpcr_mosaicism")

source("scripts/load_ddpcr_data.R")

#########################
# Calculations
#########################

ws_2022 <- unique(grep(pattern = "^22", x = ddpcr_data$worksheet, value = TRUE))

failed_wells <- nrow(ddpcr_data %>%
       filter(substr(well, 1, 1) != "M" &
                accepted_droplets < 10000))

total_wells <- nrow(ddpcr_data %>%
                      filter(substr(well, 1, 1) != "M"))

fail_rate <- (failed_wells / total_wells) * 100

#########################
# Plot
#########################

droplets_qc_plot <- ddpcr_data %>%
  filter(substr(well, 1, 1) != "M" &
           worksheet %in% ws_2022) %>%
  # Remove duplicate rows for channel 1 and channel 2
  distinct(worksheet_well_sample, .keep_all = TRUE) %>%
  ggplot(aes(x = worksheet, y = accepted_droplets)) +
  geom_jitter() +
  labs(y = "Total Droplets", x = "", 
       title = paste0("Droplet generation QC plot for 2022: ",
                      round(fail_rate, 1), "% fail rate")) +
  geom_hline(yintercept = 10000, colour = "red", linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = droplets_qc_plot, 
       filename = paste0("droplets_qc_plot_", format(Sys.time(), "%Y%m%d"), ".tiff"),
       path = "ddpcr_mosaicism/plots", 
       device= 'tiff')

#########################