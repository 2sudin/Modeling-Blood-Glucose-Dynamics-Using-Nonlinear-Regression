# Modeling Blood Glucose Dynamics Using Nonlinear Regression

# Package Installation
install.packages(c("tidyverse", "zoo", "GGally", "reshape2"))
install.packages(c("moments", "knitr"))
install.packages("e1071", "corrplot")
install.packages("nortest")

# Library Loading
library(dplyr)
library(tidyverse)
library(zoo)
library(GGally)
library(reshape2)
library(moments)
library(knitr)
library(corrplot)
library(e1071)
library(nortest)

# Reading the dataset
data <- read.csv("D:/M Data Science/sem 1/001 Statistical methods [St. Methods]/Assignments/My Assignment/dataset_MSC7.csv")

# Preliminary data analysis 
nrow(data)        # Check rows
ncol(data)        # Check columns
dim(data)         # Check dimensions
head(data)        # Display first few rows
str(data)         # View structure
sapply(data, class) # Check data types

# Identifying numeric columns
numeric_cols <- sapply(data, is.numeric)
numeric_data <- data[, numeric_cols]

# Generating a descriptive statistics table
stats_table <- data.frame(
  variable = names(numeric_data),
  mean = sapply(numeric_data, mean, na.rm = TRUE),
  median = sapply(numeric_data, median, na.rm = TRUE),
  min = sapply(numeric_data, min, na.rm = TRUE),
  max = sapply(numeric_data, max, na.rm = TRUE),
  sd = sapply(numeric_data, sd, na.rm = TRUE),
  q25 = sapply(numeric_data, function(x) quantile(x, 0.25, na.rm = TRUE)),
  q75 = sapply(numeric_data, function(x) quantile(x, 0.75, na.rm = TRUE)),
  skew = sapply(numeric_data, skewness, na.rm = TRUE),
  kurt = sapply(numeric_data, kurtosis, na.rm = TRUE)
)

# Format output nicely
kable(stats_table, digits = 3, caption = "Descriptive Statistics for Numeric Variables:")

# Plot of hr mean before interpolation
df$time <- seq_len(nrow(df))
plot(df$time, df$hr_mean, type = "l", col = "blue", 
     xlab = "Time", ylab = "hr_mean", main = "hr")

# Missing values handling
sum(is.na(data))
colSums(is.na(data))

# Fill missing hr_mean values using linear interpolation
data$hr_mean <- na.approx(data$hr_mean, na.rm = FALSE)

# Fill remaining missing values at start/end
data$hr_mean <- na.locf(data$hr_mean, fromLast = TRUE)
data$hr_mean <- na.locf(data$hr_mean)

# Duplicates check
sum(duplicated(data))

# Recalculate stats for column 4
mean(data[[4]], na.rm = TRUE)
median(data[[4]], na.rm = TRUE)
summary(data[[4]])

# Renaming columns
colnames(data)[7] <- "Y"
names(data)[names(data) %in% c("bg_mean","insulin_sum","carbs_sum","hr_mean","steps_sum","cals_sum")] <- 
  c("x1","x2","x3","x4","x5","x6")

# Add time column (1 unit = 1 hour)
data$time <- 1:nrow(data)

# Reorder columns: Y at the end
data <- data[, c("time", "x1", "x2", "x3", "x4", "x5", "x6", "Y")]
sapply(data, class)

# Select variables for multivariate time series
ts_data <- dplyr::select(data, time, x1, x2, x3, x4, x5, x6)

# Convert to long format
ts_long <- melt(ts_data, id.vars = "time")

# Create time series object
ts_data_obj <- ts(ts_data %>% select(-time), frequency = 12 * 24)

# Variable mapping for labels and colors
var_labels <- c(x1 = "bg_mean", x2 = "insulin_sum", x3 = "carbs_sum", 
                x4 = "hr_mean", x5 = "steps_sum", x6 = "cals_sum")

var_colors <- c(x1 = "#F4A3A3", x2 = "#FFD966", x3 = "#8FD694", 
                x4 = "#87CEEB", x5 = "#4A90E2", x6 = "#F4A6C1")

# Faceted time series plot
ggplot(ts_long, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~variable, scales = "free_y", ncol = 1, 
             labeller = labeller(variable = function(x) paste(x, "-", var_labels[x]))) +
  scale_color_manual(values = var_colors) +
  labs(title = "Time Series of Input Variables", x = "Time", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# Plot Y against time
ggplot(data, aes(x = time, y = Y)) +
  geom_line(color = "steelblue", linewidth = 0.1) +
  labs(title = "Time Series of Blood Glucose (Y)", x = "Time", y = "Blood Glucose (Y)") +
  theme_minimal()

# 100-point rolling average
ts_long_roll <- ts_long %>%
  group_by(variable) %>%
  arrange(time) %>%
  mutate(roll_avg = rollmean(value, k = 100, fill = NA, align = "right")) %>%
  ungroup()

# Plot original values with rolling average
ggplot(ts_long_roll, aes(x = time, color = variable)) +
  geom_line(aes(y = value), linewidth = 0.35, alpha = 0.4) +
  geom_line(aes(y = roll_avg), linewidth = 0.6, alpha = 0.9) +
  facet_wrap(~variable, scales = "free_y", ncol = 1,
             labeller = labeller(variable = function(x) paste(x, ":", var_labels[x]))) +
  scale_color_manual(values = var_colors) +
  labs(title = "Time Series of Input Variables with 100-point Rolling Average", x = "Time", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# Rolling average for Y
data$Y_roll <- rollmean(data$Y, k = 100, fill = NA, align = "right")

ggplot(data, aes(x = time)) +
  geom_line(aes(y = Y), color = "steelblue", alpha = 0.5) +
  geom_line(aes(y = Y_roll), color = "red", size = 1) +
  labs(title = "Time Series of Blood Glucose (Y) with 100-point Rolling Average", x = "Time", y = "Blood Glucose (Y)") +
  theme_minimal()

# Simplified Trends (200-row Sections)
ts_long <- ts_long %>% mutate(section = ((row_number() - 1) %/% 200) + 1)
ts_summary <- ts_long %>%
  group_by(variable, section) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), mean_time = mean(time, na.rm = TRUE)) %>%
  ungroup()

ggplot(ts_summary, aes(x = mean_time, y = mean_value, color = variable)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = var_colors) +
  facet_wrap(~variable, scales = "free_y", ncol = 1,
             labeller = labeller(variable = function(x) paste(x, ":", var_labels[x]))) +
  labs(title = "Simplified Time Series Trends (200-row Sections)", x = "Time", y = "Mean Value") +
  theme_minimal() +
  theme(legend.position = "none")

# Simplified Y
data <- data %>% mutate(section = ((row_number() - 1) %/% 200) + 1)
y_summary <- data %>%
  group_by(section) %>%
  summarise(mean_Y = mean(Y, na.rm = TRUE), mean_time = mean(time, na.rm = TRUE)) %>%
  ungroup()

ggplot(y_summary, aes(x = mean_time, y = mean_Y)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "Simplified Time Series of Y (200-row Sections)", x = "Time", y = "Mean Blood Glucose (Y)") +
  theme_minimal()

# Distribution plots for inputs (x1-x6)
plot_data <- data %>%
  select(x1, x2, x3, x4, x5, x6) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "value")

skew_df <- plot_data %>%
  group_by(Variable) %>%
  summarise(skewness = round(skewness(value, na.rm = TRUE), 3))

plot_data <- plot_data %>% left_join(skew_df, by = "Variable")

ggplot(plot_data, aes(x = value, fill = Variable)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.7, color = "black") +
  geom_density(alpha = 0.3, linewidth = 1) +
  facet_wrap(~Variable, scales = "free", ncol = 3) +
  geom_text(data = skew_df, aes(x = Inf, y = Inf, label = paste("Skewness:", skewness)),
            hjust = 1.1, vjust = 1.5, size = 3, inherit.aes = FALSE) +
  labs(title = "Distribution of Input Variables", x = "Value", y = "Density", fill = "Variable") +
  theme_minimal() +
  theme(legend.position = "right", strip.text = element_text(face = "bold"))

# Distribution of Y
skew_fun <- function(x) {
  x <- x[!is.na(x)]; n <- length(x); m <- mean(x); s <- sd(x)
  sum((x - m)^3) / ((n - 1) * s^3)
}

y_skew <- round(skew_fun(data$Y), 3)

ggplot(data, aes(x = Y)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_density(color = "red", linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = paste("Skewness:", y_skew), hjust = 1.1, vjust = 1.5, size = 4) +
  labs(title = "Distribution of Output Variable (Y)", x = "Blood Glucose(Y)", y = "Density") +
  theme_minimal()

# Correlations
corr_df <- data %>%
  select(x1, x2, x3, x4, x5, x6, Y) %>%
  summarise(across(x1:x6, ~ cor(.x, Y, use = "complete.obs"))) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "correlation")

ggplot(corr_df, aes(x = variable, y = correlation, fill = variable)) +
  geom_col(color = "black") +
  scale_fill_manual(values = var_colors) +
  geom_text(aes(label = round(correlation, 2)), vjust = -0.5) +
  labs(title = "Correlation of Variables x1,x2,x3,x4,x5,x6 with Y", x = "Variable", y = "Correlation with Y") +
  theme_minimal() +
  theme(legend.position = "none")

# Correlation Matrix
corr_data <- data %>% select(x1, x2, x3, x4, x5, x6, Y)
corr_matrix <- cor(corr_data, use = "complete.obs")
corr_long <- melt(corr_matrix, varnames = c("Var1", "Var2"), value.name = "Correlation")

ggplot(corr_long, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "Correlation Matrix of x1–x6 with Y", x = "", y = "")

# Scatter Plots
scatter_long <- data %>%
  select(x1:x6, Y) %>%
  pivot_longer(cols = x1:x6, names_to = "variable", values_to = "value")

ggplot(scatter_long, aes(x = value, y = Y, color = variable)) +
  geom_point(alpha = 0.15, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  scale_color_manual(values = var_colors) +
  facet_wrap(~variable, scales = "free_x") +
  labs(title = "Scatter Plots of x1 to x6 vs Y", x = "Input Variable", y = "Y (Blood Glucose)") +
  theme_minimal()

# Regression Modeling
y <- data$Y
X1 <- model.matrix(~ I(x1^3) + I(x2^2) + I(x3^2) + x4 + x5 + x6, data)
X2 <- model.matrix(~ I(x1^2) + I(x2^2) + I(x3^3) + x4 + x5 + x6, data)
X3 <- model.matrix(~ x1 + x2 + x3 + I(x4^2) + x5 + I(x6^2), data)
X4 <- model.matrix(~ I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2) + I(x6^2), data)
X5 <- model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + I(x1*x2) + I(x3*x4) + I(x2*x6), data)

# Least Squares estimator function
ls_est <- function(X, Y) { qr.solve(X, Y) }

# Estimate thetas
theta_1 <- ls_est(X1, y)
theta_2 <- ls_est(X2, y)
theta_3 <- ls_est(X3, y)
theta_4 <- ls_est(X4, y)
theta_5 <- ls_est(X5, y)

options(scipen = 999)
theta_1; theta_2; theta_3; theta_4; theta_5

# RSS calculations
rss <- function(X, Y, theta) { sum((Y - X %*% theta)^2) }
RSS1 <- rss(X1, y, theta_1)
RSS2 <- rss(X2, y, theta_2)
RSS3 <- rss(X3, y, theta_3)
RSS4 <- rss(X4, y, theta_4)
RSS5 <- rss(X5, y, theta_5)

# Log-likelihood calculations
m1 <- lm(Y ~ I(x1^3) + I(x2^2) + I(x3^2) + x4 + x5 + x6, data)
m2 <- lm(Y ~ I(x1^2) + I(x2^2) + I(x3^3) + x4 + x5 + x6, data)
m3 <- lm(Y ~ x1 + x2 + x3 + I(x4^2) + x5 + I(x6^2), data)
m4 <- lm(Y ~ I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2) + I(x6^2), data)
m5 <- lm(Y ~ x1 + x2 + x3 + x4 + x5 + x6 + I(x1*x2) + I(x3*x4) + I(x2*x6), data)

LogLikelihood_var <- function(m) {
  n <- length(residuals(m))
  RSS <- sum(residuals(m)^2)
  s2 <- RSS / (n - 1)
  ll <- (n/2)*log(2*pi) - (n/2)*log(s2) - RSS/(2*s2)
  c(Variance = s2, LogLikelihood = ll)
}

result <- data.frame(Model = paste0("Model ", 1:5), t(sapply(list(m1, m2, m3, m4, m5), LogLikelihood_var)))
print(result)

# AIC and BIC
AIC_BIC <- function(m) {
  n <- length(residuals(m)); k <- length(coef(m))
  ll <- LogLikelihood_var(m)["LogLikelihood"]
  c(AIC = 2*k - 2*ll, BIC = k*log(n) - 2*ll)
}

# Q-Q plots for residuals
model_colors <- as.vector(var_colors[1:5])
models <- list(m1, m2, m3, m4, m5)
par(mfrow = c(2, 3), mar = c(3, 3, 2, 1))

for (i in 1:5) {
  r <- residuals(models[[i]])
  r <- r[!is.na(r)]
  qqnorm(r, main = paste("Q-Q Plot: Model", i), pch = 1, col = model_colors[i])
  qqline(r, col = "red", lwd = 2)
}
par(mfrow = c(1, 1))

# Histograms of residuals
par(mfrow = c(2, 3))
for (i in 1:5) {
  r <- residuals(models[[i]])
  hist(r, breaks = 30, col = model_colors[i], border = "white", main = paste("Histogram: Model", i))
  abline(v = 0, col = "red", lwd = 2)
}
par(mfrow = c(1, 1))

# Normality Tests
normality_results <- data.frame(
  Model = paste0("Model ", 1:5),
  Anderson_Darling_p = sapply(models, function(m) ad.test(residuals(m))$p.value)
)
print(normality_results)

# Final Results Table
final_results <- data.frame(
  Model = paste0("Model ", 1:5),
  t(sapply(models, function(m) {
    n <- length(residuals(m)); k <- length(coef(m))
    ll_out <- LogLikelihood_var(m)
    c(Variance = ll_out["Variance"], LogLikelihood = ll_out["LogLikelihood"],
      AIC = 2*k - 2*ll_out["LogLikelihood"], BIC = k*log(n) - 2*ll_out["LogLikelihood"])
  }))
)
print(final_results)

# Split Data (70/30)
set.seed(123)
n_val <- nrow(data)
train_index <- sample(seq_len(n_val), size = 0.7 * n_val)
train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

# Model 5 on Training Data
model5_train <- lm(Y ~ x1 + x2 + x3 + x4 + x5 + x6 + I(x1 * x2) + I(x3 * x4) + I(x2 * x6), data = train_data)
summary(model5_train)

# Predictions and Confidence Intervals
pred_ci <- predict(model5_train, newdata = test_data, interval = "prediction", level = 0.95)

plot_data_pred <- data.frame(
  y_true = test_data$Y,
  y_pred = pred_ci[, "fit"],
  lwr    = pred_ci[, "lwr"],
  upr    = pred_ci[, "upr"]
)

# Plot Observed vs Predicted
plot(plot_data_pred$y_true, plot_data_pred$y_pred, pch = 16, xlab = "Observed Y", ylab = "Predicted Y", 
     main = "Model 5 Predictions with 95% PIs")
arrows(x0 = plot_data_pred$y_true, y0 = plot_data_pred$lwr, x1 = plot_data_pred$y_true, y1 = plot_data_pred$upr, 
       angle = 90, code = 3, length = 0.05, col = "gray")
abline(0, 1, col = "red", lwd = 2)

# Training vs Testing Distributions
train_data$set <- "Training"; test_data$set <- "Testing"
combined_data <- rbind(train_data, test_data)

hist(train_data$Y, breaks = 30, col = rgb(0, 0, 1, 0.5), xlab = "Y", main = "Distribution of Y: Training vs Testing")
hist(test_data$Y, breaks = 30, col = rgb(1, 0, 0, 0.5), add = TRUE)
legend("topright", legend = c("Training", "Testing"), fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))

# Density and Boxplots
plot(density(train_data$Y), col = "blue", lwd = 2, main = "Density Distribution of Y")
lines(density(test_data$Y), col = "red", lwd = 2)
boxplot(Y ~ set, data = combined_data, col = c("lightblue", "lightpink"), main = "Boxplot of Y", ylab = "Y")

# ABC Section
theta_5_no_intercept <- theta_5[-1]
theta_5_abs <- abs(theta_5_no_intercept)
top_2_indices <- order(theta_5_abs, decreasing = TRUE)[1:2]

param_1_index <- top_2_indices[1] + 1
param_2_index <- top_2_indices[2] + 1
param_1_value <- theta_5[param_1_index]
param_2_value <- theta_5[param_2_index]

# Priors
range_factor <- 0.3
prior_1_min <- param_1_value - abs(param_1_value) * range_factor
prior_1_max <- param_1_value + abs(param_1_value) * range_factor
prior_2_min <- param_2_value - abs(param_2_value) * range_factor
prior_2_max <- param_2_value + abs(param_2_value) * range_factor

# Rejection ABC
set.seed(456)
n_samples <- 50000
epsilon <- 0.5
n_accepted <- 0
max_accepted <- 2000
accepted_params <- matrix(NA, nrow = max_accepted, ncol = 2)
colnames(accepted_params) <- c(paste0("beta_", top_2_indices[1]), paste0("beta_", top_2_indices[2]))

Y_obs <- data$Y
theta_fixed <- theta_5

for (i in 1:n_samples) {
  if (n_accepted >= max_accepted) break
  theta_candidate <- theta_fixed
  theta_candidate[param_1_index] <- runif(1, prior_1_min, prior_1_max)
  theta_candidate[param_2_index] <- runif(1, prior_2_min, prior_2_max)
  
  Y_sim <- X5 %*% theta_candidate
  distance <- mean(abs(Y_obs - Y_sim))
  
  if (distance < epsilon) {
    n_accepted <- n_accepted + 1
    accepted_params[n_accepted, ] <- c(theta_candidate[param_1_index], theta_candidate[param_2_index])
  }
}

accepted_params <- accepted_params[1:n_accepted, ]

# Plot Posterior Distributions
par(mfrow = c(2, 2))
hist(accepted_params[, 1], breaks = 30, col = "lightblue", main = "Marginal Posterior: Beta 1", freq = FALSE)
lines(density(accepted_params[, 1]), col = "blue", lwd = 2)
abline(v = param_1_value, col = "red", lty = 2)

hist(accepted_params[, 2], breaks = 30, col = "lightgreen", main = "Marginal Posterior: Beta 2", freq = FALSE)
lines(density(accepted_params[, 2]), col = "darkgreen", lwd = 2)
abline(v = param_2_value, col = "red", lty = 2)
par(mfrow = c(1, 1))

# Joint posterior (2D scatter plot)
plot(
  accepted_params[, 1],
  accepted_params[, 2],
  pch  = 16,
  cex  = 0.5,
  col  = rgb(0, 0, 1, 0.3),
  xlab = colnames(accepted_params)[1],
  ylab = colnames(accepted_params)[2],
  main = "Joint Posterior Distribution"
)
points(
  param_1_value, 
  param_2_value, 
  pch = 3, 
  col = "red", 
  cex = 2, 
  lwd = 3
)
legend(
  "topright",
  legend = c("Posterior samples", "LS Estimate"),
  pch    = c(16, 3),
  col    = c(rgb(0, 0, 1, 0.5), "red"),
  cex    = 0.8
)

# Joint posterior (2D density contour)
library(MASS)
kde <- kde2d(accepted_params[, 1], accepted_params[, 2], n = 50)
contour(
  kde,
  xlab = colnames(accepted_params)[1],
  ylab = colnames(accepted_params)[2],
  main = "Joint Posterior Density Contours",
  col  = "blue",
  lwd  = 1.5
)
points(
  param_1_value, 
  param_2_value, 
  pch = 3, 
  col = "red", 
  cex = 2, 
  lwd = 3
)

par(mfrow = c(1, 1))

# Step 5: Summary statistics of posterior distributions
cat("\n=== Posterior Summary Statistics ===\n")
cat("\nParameter 1 (", colnames(accepted_params)[1], "):\n")
cat("  LS Estimate:      ", round(param_1_value, 6), "\n")
cat("  Posterior Mean:   ", round(mean(accepted_params[, 1]), 6), "\n")
cat("  Posterior Median: ", round(median(accepted_params[, 1]), 6), "\n")
cat("  Posterior SD:     ", round(sd(accepted_params[, 1]), 6), "\n")
cat("  95% Credible Interval: [", 
    round(quantile(accepted_params[, 1], 0.025), 6), ",", 
    round(quantile(accepted_params[, 1], 0.975), 6), "]\n")

cat("\nParameter 2 (", colnames(accepted_params)[2], "):\n")
cat("  LS Estimate:      ", round(param_2_value, 6), "\n")
cat("  Posterior Mean:   ", round(mean(accepted_params[, 2]), 6), "\n")
cat("  Posterior Median: ", round(median(accepted_params[, 2]), 6), "\n")
cat("  Posterior SD:     ", round(sd(accepted_params[, 2]), 6), "\n")
cat("  95% Credible Interval: [", 
    round(quantile(accepted_params[, 2], 0.025), 6), ",", 
    round(quantile(accepted_params[, 2], 0.975), 6), "]\n")

# Correlation between the two parameters in posterior
posterior_cor <- cor(accepted_params[, 1], accepted_params[, 2])
cat("\nPosterior Correlation between parameters:", round(posterior_cor, 4), "\n")

# Optional: Save accepted parameters for further analysis
write.csv(accepted_params, "abc_posterior_samples.csv", row.names = FALSE)
cat("\nPosterior samples saved to 'abc_posterior_samples.csv'\n")