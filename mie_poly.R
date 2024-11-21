library(tidyverse)

# Check the structure of the data
mie <- read_csv("calibrated-mie.csv")

mie_subset <- mie %>%
  select(scatter, diam_740_mid)

# --- Linear Fitting ---

# Polynomial fitting on original data
model2 <- lm(diam_740_mid ~ poly(scatter, 2), data = mie_subset)
model3 <- lm(diam_740_mid ~ poly(scatter, 3), data = mie_subset)
model4 <- lm(diam_740_mid ~ poly(scatter, 4), data = mie_subset)

# View the summaries
summary(model2)
summary(model3)
summary(model4)

# Compare using AIC
AIC(model2, model3, model4)

# Generate a sequence of scatter values
scatter_values <- seq(min(mie_subset$scatter), max(mie_subset$scatter), length.out = 1000)

# Predict diam_740_mid for the new scatter values using the best model
predictions <- predict(model3, newdata = data.frame(scatter = scatter_values))

# Create a data frame for plotting
prediction_df <- data.frame(scatter = scatter_values, diam_pred = predictions)




# --- Log Transform Fitting ---

# Log transform the data
mie_subset <- mie_subset %>%
  mutate(log_scatter = log10(scatter), 
         log_diam_740_mid = log10(diam_740_mid))

# Polynomial fitting on log-transformed data
model2_log <- lm(log_diam_740_mid ~ poly(log_scatter, 2), data = mie_subset)
model3_log <- lm(log_diam_740_mid ~ poly(log_scatter, 3), data = mie_subset)
model4_log <- lm(log_diam_740_mid ~ poly(log_scatter, 4), data = mie_subset)

# View the summaries
summary(model2_log)
summary(model3_log)
summary(model4_log)

# Compare using AIC
AIC(model2_log, model3_log, model4_log)
# Generate a sequence of log_scatter values
log_scatter_values <- seq(min(mie_subset$log_scatter), 
                          max(mie_subset$log_scatter), 
                          length.out = 1000)

# Predict log_diam_740_mid for the new log_scatter values (using model3_log as example)
log_predictions <- predict(model3_log, 
                           newdata = data.frame(log_scatter = log_scatter_values))

# Back-transform to original scale
predictions_log <- 10^log_predictions

# Create a data frame for plotting
prediction_df_log <- data.frame(scatter = 10^log_scatter_values, 
                                diam_pred = predictions_log)


# --- Plotting ---
data_plotting <- bind_rows(mie_subset %>% select(scatter, diameter = diam_740_mid) %>% mutate(type = "Original"), 
                           prediction_df %>% rename(diameter = diam_pred) %>% mutate(type = "Linear"), 
                           prediction_df_log %>% rename(diameter = diam_pred) %>% mutate(type = "Log"))
# Linear scale
a <- data_plotting %>% ggplot(aes(x = scatter, y = diameter)) +
  geom_line(aes(col = type)) +
  labs(title = "Linear scale")

# Log scale
b <- data_plotting %>% ggplot(aes(x = scatter, y = diameter)) +
  geom_line(aes(col = type)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Log scale")
  
ggpubr::ggarrange(a, b, ncol = 2, common.legend = TRUE)


