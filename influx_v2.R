library(tidyverse)
library(dplyr)
library(scales)
library(viridis)
library(ggplot2)
library(cowplot)

# Set the working directory
path_to_git_repository <- "C:/Users/SNEHA/OneDrive/Documents/seaflow/fsc-poc-calibration"
setwd(path_to_git_repository)

# Read data
poc <- read_csv("poc-data.csv", show_col_types = FALSE)
cultures <- read_csv("influx-cultures.csv", show_col_types = FALSE)

# Rename the column
poc <- poc %>% rename("Sample.ID" = "Sample ID")

# Add Sample ID categories to cultures data
cultures <- cultures %>%
  mutate(Sample.ID = rep(
    c("EHUX", "LICMO", "Micromonas pusilla", "Navicula transitans", 
      "Phaeodactylum tricornutum", "Thalassiosira pseudonana (1135)", 
      "Thalassiosira pseudonana (3367)", "TW 3365", 
      "Prochlorococcus (1314)", "Synechococcus (WH7803)", 
      "Prochlorococcus (AS9601)", "Prochlorococcus (MED4)", 
      "Prochlorococcus (NATL12A)", "Synechococcus (WH8102)"),
    times = c(2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2)  # Adjusted replicates to sum to 30
  ))

#summary calculation
poc_summary <- poc %>%
  group_by(Sample.ID) %>%
  summarize(
    C_ug_ml_mean = mean(`C (ug/ml)`, na.rm = TRUE),
    C_ug_ml_sd = sd(`C (ug/ml)`, na.rm = TRUE),
    N_ug_ml_mean = mean(`N (ug/ml)`, na.rm = TRUE),
    N_ug_ml_sd = sd(`N (ug/ml)`, na.rm = TRUE)
  )

# Calculate mean and standard deviation for cultures data
cultures_summary <- cultures %>%
  group_by(Sample.ID) %>%
  summarize(
    abundance_cells_mL_mean = mean(abundance_cells.mL, na.rm = TRUE),
    abundance_cells_mL_sd = sd(abundance_cells.mL, na.rm = TRUE),
    norm_fsc_mean = mean(norm.fsc, na.rm = TRUE),
    norm_fsc_sd = sd(norm.fsc, na.rm = TRUE),
    norm_chl_mean = mean(norm.chl, na.rm = TRUE),
    norm_chl_sd = sd(norm.chl, na.rm = TRUE)
  )

# Merge summaries for POC and cultures data
merged_data <- inner_join(cultures_summary, poc_summary, by ="Sample.ID")

# Calculate cell quotas and standard deviation for pgC and pgN per cell
merged_data <- merged_data %>%
  mutate(
    pgC_cell = 10^6 * C_ug_ml_mean / abundance_cells_mL_mean,
    pgN_cell = 10^6 * N_ug_ml_mean / abundance_cells_mL_mean,
    pgC_cell_sd = pgC_cell * sqrt((C_ug_ml_sd / C_ug_ml_mean)^2 + 
                                    (abundance_cells_mL_sd / abundance_cells_mL_mean)^2),
    pgN_cell_sd = pgN_cell * sqrt((N_ug_ml_sd / N_ug_ml_mean)^2 + 
                                    (abundance_cells_mL_sd / abundance_cells_mL_mean)^2)
  )

# Save the output with selected columns
merged_data %>%
  select(Sample.ID, norm_fsc_mean, norm_fsc_sd, norm_chl_mean, norm_chl_sd, 
         abundance_cells_mL_mean, abundance_cells_mL_sd, pgC_cell, pgN_cell, 
         pgC_cell_sd, pgN_cell_sd) %>%
  write_csv("Influx-Qc-cultures.csv")

merged_data %>%
  select(Sample.ID, abundance_cells_mL_mean, abundance_cells_mL_sd, 
         pgC_cell, pgN_cell, pgC_cell_sd, pgN_cell_sd) %>%
  write_csv("Qc-cultures.csv")



# Define file paths and set working directory
path.to.git.repository <- "~/Documents/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)

# Read data
mie <- read_csv("calibrated-mieINFLUX.csv")
merge <- read_csv("Influx-Qc-cultures.csv") %>%
  filter(Sample.ID != "Phaeodactylum tricornutum") %>% # remove non-spherical cells
  arrange('norm_fsc_mean')



# Create a list to hold the plots
plots <- list()

# Loop through each instrument and create plots
for(inst in c("Leo", "Penny")) {
  
  # Generate the plot
  p <- ggplot(merge, aes(x = norm_fsc_mean, y = pgC_cell)) +
    geom_errorbar(aes(ymin = pgC_cell - pgC_cell_sd, ymax = pgC_cell + pgC_cell_sd), 
                  color = 'darkgrey', width = 0.05) +
    geom_errorbarh(aes(xmin = norm_fsc_mean - norm_fsc_sd, xmax = norm_fsc_mean + norm_fsc_sd), 
                   color = 'darkgrey', height = 0.05) +
    # Add the Mie-based model lines (red line for mid, dashed lines for upr and lwr)
    geom_line(data = mie, aes(x = scatter, y = .data[[paste0("Qc_", inst, "_mid")]], 
                              color = "Mie-based model (n = 1.38 ± 0.3)"), size = 1) +
    geom_line(data = mie, aes(x = scatter, y = .data[[paste0("Qc_", inst, "_upr")]]), 
              color = 'grey', linetype = "dashed", size = 1, show.legend = FALSE) +
    geom_line(data = mie, aes(x = scatter, y = .data[[paste0("Qc_", inst, "_lwr")]]), 
              color = 'grey', linetype = "dashed", size = 1, show.legend = FALSE) +
    geom_point(aes(fill = Sample.ID), size = 6, shape = 21, color = "black") +
    # Adjust color fill for Sample.ID and add the red line to legend
    scale_fill_viridis_d(option = "D", alpha = 0.5) +
    scale_color_manual(values = c("Mie-based model (n = 1.38 ± 0.3)" = "red3")) + # Custom legend label
    scale_x_log10(limits = c(0.002, 10), 
                  breaks = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    scale_y_log10(limits = c(0.005, 100), 
                  breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 1000)) +
    labs(title = inst, 
         x = "Normalized scatter (dimensionless)", 
         y = expression(paste("Qc (pgC cell"^{-1},")"))) +
    theme_classic() +  # Ensure white background
    theme(
      legend.position = c(0.15, 0.7),  # Top-left corner inside the box
      legend.background = element_blank(),  # Remove legend background box
      legend.title = element_blank(), # Remove "Sample.ID" title
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, size = 1) # Box settings around the plot
    )
  
  # Save the plot in the list
  plots[[inst]] <- p
}

# Combine the plots side by side without labels A and B
combined_plot <- plot_grid(plots$Leo, plots$Penny, ncol = 2)

# Save the combined plot to a PNG file with white background
ggsave("INFLUX_Qc-scatter.png", combined_plot, 
       width = 18, height = 8, dpi = 300, bg = "white")  # Explicitly set bg to white.
