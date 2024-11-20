# Load necessary libraries
library(tidyverse)    # For data manipulation
library(scales)       # For color transparency
library(viridis)      # For color palette
library(gridExtra)    # To arrange multiple plots in a single output
library(purrr)
library(dplyr)

path.to.git.repository <- "~/Documents/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)
poc <- read.csv("Qc-cultures.csv")


# Define the instruments to loop over
instruments <- c(740, 751)

# Loop through each instrument
for (inst in instruments) {
  
  # Load culture data for each instrument with explicit column specification
  cultures <- read_csv(paste0(inst, "-cultures.csv"), show_col_types = FALSE)
  
  # View column specification to ensure correct data types
  print(spec(cultures))
  
  # Assign 'Sample.ID' labels based on the instrument number
  cultures <- cultures %>%
    mutate(Sample.ID = case_when(
      inst == 740 ~ rep(c("Thalassiosira weissflogii", "Navicula transitans", 
                          "Thalassiosira pseudonana (1135)", "Thalassiosira pseudonana (3367)", 
                          "Phaedactylum tricornutum", "Micromonas pusilla", 
                          "Prochlorococcus (AS9601)", "Prochlorococcus (1314)", 
                          "Prochlorococcus (MED4)", "Prochlorococcus (NATL12A)", 
                          "Synechococcus (WH8102)", "Synechococcus (WH7803)"), each = 2),
      inst == 751 ~ rep(c("Thalassiosira weissflogii", "Navicula transitans", 
                          "Thalassiosira pseudonana (1135)", "Thalassiosira pseudonana (3367)", 
                          "Phaedactylum tricornutum", "Micromonas pusilla", 
                          "Prochlorococcus (MED4)", "Prochlorococcus (AS9601)", 
                          "Prochlorococcus (1314)", "Prochlorococcus (NATL12A)", 
                          "Synechococcus (WH8102)", "Synechococcus (WH7803)"), each = 2)
    ))
  
  # Calculate the standard deviation for each numeric group in 'cultures'
  cultures_sd <- cultures %>%
    group_by(Sample.ID) %>%
    summarize(across(where(is.numeric), ~ sd(.x, na.rm = TRUE)), .groups = 'drop') %>%
    rename_with(~ paste0(.x, ".sd"), -Sample.ID)
  
  # Calculate the mean for each numeric group in 'cultures'
  cultures_mean <- cultures %>%
    group_by(Sample.ID) %>%
    summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
  
  # Add the standard deviation values as new columns in 'cultures_mean'
  cultures_mean <- cultures_mean %>%
    inner_join(cultures_sd, by = "Sample.ID") %>%
    rename(norm.fsc.sd = norm.fsc.sd, norm.chl.sd = norm.chl.sd)
  
  ### MERGE POC with Cell Abundance
  # Merge with 'poc' data on 'Sample.ID' to combine related measurements
  merge_data <- poc %>%
    inner_join(cultures_mean, by = "Sample.ID")
  
  # Write the final merged data to a CSV file for the current instrument
  merge_data %>%
    select(Sample.ID, norm.fsc, norm.fsc.sd, norm.chl, norm.chl.sd, 
           abundance_cells_mL_mean, abundance_cells_mL_sd, pgC_cell, pgN_cell, 
           pgC_cell_sd, pgN_cell_sd) %>%
    write_csv(paste0(inst, "-Qc-cultures.csv"))
}



# Set working directory (adjust this to your actual file path)
path_to_git_repository <- "~/Documents/Codes/fsc-poc-calibration"
setwd(path_to_git_repository)

# Load the Mie theory data from CSV file
mie <- read_csv("calibrated-mie.csv")

# Define output PNG for side-by-side plots, double the width
png("Qc-scatter.png", width = 18, height = 8, units = "in", res = 300)

# Initialize list to store plots
plots <- list()

# Loop through each instrument (740 and 751)
for (inst in c(740, 751)) {
  
  # Load culture data for each instrument
  merge <- read_csv(paste0(inst, "-Qc-cultures.csv"))
  
  merge2 <- merge %>%
    dplyr::filter(Sample.ID != "Phaeodactylum tricornutum") %>%
    arrange(norm.fsc)
  
  # Create a ggplot for the current instrument
  p <- ggplot(merge2, aes(x = norm.fsc, y = pgC_cell)) +
    
    # Add error bars for pgC.cell with standard deviations
    geom_errorbar(aes(ymin = pgC_cell - pgC_cell_sd, ymax = pgC_cell + pgC_cell_sd), 
                  width = 0, color = "grey", size = 1) +
    geom_errorbarh(aes(xmin = norm.fsc - norm.fsc.sd, xmax = norm.fsc + norm.fsc.sd),
                   height = 0, color = "grey", size = 1) +
    
    # Add lines for Mie-based model with upper, mid, and lower bounds
    geom_line(data = mie, aes(x = scatter, y = !!sym(paste0("Qc_", inst, "_mid")), 
                              color = "Mie-based model (n = 1.38 ± 0.3)"), 
              size = 1.2) +
    geom_line(data = mie, aes(x = scatter, y = !!sym(paste0("Qc_", inst, "_upr"))), 
              color = "grey", size = 1, show.legend = FALSE) +
    geom_line(data = mie, aes(x = scatter, y = !!sym(paste0("Qc_", inst, "_lwr"))), 
              color = "grey", size = 1, show.legend = FALSE) +
    
    # Add larger points for each culture sample with semi-transparent colors
    geom_point(aes(fill = Sample.ID), shape = 21, size = 8, 
               color = "black", alpha = 0.7) +  # Larger points and slightly more opaque
    scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.9) +
    
    # Custom legend entry for the Mie-based model
    scale_color_manual(values = c("Mie-based model (n = 1.38 ± 0.3)" = "red3")) +
    
    # Set log scales for x and y axes
    scale_x_log10(limits = c(0.002, 10), 
                  breaks = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
                  labels = scales::label_number()) +
    scale_y_log10(limits = c(0.005, 100),
                  breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100),
                  labels = scales::label_number()) +
    
    # Add labels and title
    labs(x = "Normalized scatter (dimensionless)", 
         y = expression(paste("Qc (pgC cell"^{-1},")")),
         title = paste("Instrument:", inst)) +
    
    # Customize theme for clarity and format
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = c(0.2, 0.7),  # Position legend in the top-left corner
      legend.background = element_blank(),  # No legend background box
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  
  # Save the plot in the list
  plots[[as.character(inst)]] <- p
}

# Arrange plots side-by-side in the PNG output
grid.arrange(grobs = plots, ncol = 2)

# Close the PNG device
dev.off()

###CORREALTION


# Initialize list to store plots
plots <- list()

generate_plot <- function(inst, merge2, mie) {
  # Generate column name dynamically
  col_name <- paste0("Qc_", inst, "_mid")
  
  # Match scatter intervals to mie scatter values
  id <- findInterval(merge2$norm.fsc, mie$scatter)
  id[is.na(id)] <- 1  # Replace NA indices with a valid default
  
  # Add predicted values to merge2 using the mie table
  merge2 <- merge2 %>%
    mutate(predicted = mie[id, col_name] %>% pull())  # Ensure numeric type
  
  # Linear regression
  reg <- lm(pgC_cell ~ predicted, data = merge2)
  print(paste("Regression summary for instrument", inst))
  print(summary(reg))
  
  # Generate the ggplot
  p <- merge2 %>%
    ggplot(aes(x = norm.fsc, y = pgC_cell)) +
    
    # Add vertical error bars
    geom_errorbar(aes(ymin = pgC_cell - pgC_cell_sd, ymax = pgC_cell + pgC_cell_sd), 
                  color = "grey", width = 0) +
    
    # Add horizontal error bars
    geom_errorbarh(aes(xmin = norm.fsc - norm.fsc.sd, xmax = norm.fsc + norm.fsc.sd), 
                   color = "grey", height = 0) +
    
    # Add larger points with transparency
    geom_point(aes(fill = Sample.ID), size = 8, shape = 21, alpha = 0.7) +
    
    # Set color scale for points with title for the legend
    scale_fill_viridis_d(name = "Sample ID") +
    
    # Log scale for x and y axes with limits and breaks
    scale_x_log10(limits = c(0.002, 10), 
                  breaks = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    scale_y_log10(limits = c(0.005, 100), 
                  breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)) +
    
    # Axis labels and plot title
    labs(x = "Normalized scatter (dimensionless)", 
         y = expression(paste("Qc (pgC cell"^{-1}, ")")),
         title = paste("Instrument", inst)) +
    
    # Add red reference line for the Mie-based model
    geom_abline(intercept = 0.5239, slope = 1.0462, color = "red3", size = 1) +
    
    # Apply classic theme to add a border around the plot
    theme_classic(base_size = 14) +
    
    # Customize theme elements for legend, plot border, and text size
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
      legend.position = c(0.2,0.7),                                          # Position legend on the right
      legend.background = element_blank(),                                # No background box for the legend
      legend.title = element_text(size = 16),                             # Increase legend title font size
      legend.text = element_text(size = 14),                              # Increase legend text font size
      plot.title = element_text(hjust = 0.5, size = 16),                  # Centered and larger title
      axis.text = element_text(size = 12),                                # Axis text size
      axis.title = element_text(size = 14),                               # Axis title size
      panel.grid.major = element_line(color = "grey", size = 0.5),        # Add major grid lines
      panel.grid.minor = element_line(color = "lightgrey", size = 0.25)   # Add minor grid lines
    ) +
    
    # Adjust legend to control point size in the legend
    guides(fill = guide_legend(title = "Sample ID", override.aes = list(size = 5)))
  
  # Return the plot object
  return(p)
}

# Generate plots for both instruments
plots <- map(c("740", "751"), ~generate_plot(.x, merge2, mie))

# Save all plots in a single PNG file using gridExtra
png("Qc-scatter-combined.png", width = 24, height = 10, units = "in", res = 300)
grid.arrange(grobs = plots, ncol = 2)  # Arrange plots side by side
dev.off()

# Log-transformed regression for both instruments
for (inst in instruments) {
  log_reg <- merge2 %>%
    dplyr::filter(pgC_cell > 0, norm.fsc > 0) %>%  # Remove invalid values for log
    mutate(log_pgC.cell = log10(pgC_cell), log_norm.fsc = log10(norm.fsc)) %>%
    lm(log_pgC.cell ~ log_norm.fsc, data = .)
  
  print(paste("Log-transformed regression summary for instrument", inst))
  print(summary(log_reg))
}
