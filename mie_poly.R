library(dplyr)
library(ggplot2)

# Assuming 'mie' is your dataset
# Check the structure of the data
str(mie)

mie_subset <- mie %>%
  select(scatter, diam_740_mid)

# Fit polynomial models with scatter as the predictor
model1 <- lm(diam_740_mid ~ poly(scatter, 1), data = mie_subset) # Linear
model2 <- lm(diam_740_mid ~ poly(scatter, 2), data = mie_subset) # Quadratic
model3 <- lm(diam_740_mid ~ poly(scatter, 3), data = mie_subset) # Cubic

# View the summaries
summary(model1)
summary(model2)
summary(model3)

# Compare using AIC
AIC(model1, model2, model3)

# Generate a sequence of scatter values
scatter_values <- seq(min(mie_subset$scatter), max(mie_subset$scatter), length.out = 100)

# Predict diam_740_mid for the new scatter values using the best model
predictions <- predict(model3, newdata = data.frame(scatter = scatter_values))

# Create a data frame for plotting
prediction_df <- data.frame(scatter = scatter_values, diam_pred = predictions)

ggplot(mie_subset, aes(x = scatter, y = diam_740_mid)) +
  geom_point(color = "blue", alpha = 0.1) +
  geom_line(data = prediction_df, aes(x = scatter, y = diam_pred), color = "red", linewidth = 1) +
  labs(title = "Polynomial Regression of Diameter vs. Scatter",
       x = "scatter", y = "diam_740_mid")

