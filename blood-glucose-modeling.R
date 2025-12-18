# Blood Glucose Dynamics: Nonlinear Regression Modeling

# Setup and Dependencies
required_packages <- c("tidyverse", "zoo", "GGally", "reshape2", "moments", "knitr", "e1071", "corrplot", "MASS")
library(tidyverse)
library(zoo)
library(GGally)
library(reshape2)
library(moments)
library(knitr)
library(corrplot)
library(e1071)
library(MASS)

# Load dataset
data_path <- "D:/M Data Science/sem 1/001 Statistical methods [St. Methods]/Assignments/My Assignment/dataset_MSC7.csv"
data <- read.csv(data_path)

# Task 1: Preliminary Data Analysis & Preprocessing

print("Dataset Dimensions:")
print(dim(data))

# Quick check on data structure
str(data)
sapply(data, class)

# Descriptive stats for numeric variables
numeric_cols <- sapply(data, is.numeric)
numeric_data <- data[, numeric_cols]

stats_table <- data.frame(
  variable = names(numeric_data),
  mean   = sapply(numeric_data, mean, na.rm = TRUE),
  median = sapply(numeric_data, median, na.rm = TRUE),
  min    = sapply(numeric_data, min, na.rm = TRUE),
  max    = sapply(numeric_data, max, na.rm = TRUE),
  sd     = sapply(numeric_data, sd, na.rm = TRUE),
  q25    = sapply(numeric_data, function(x) quantile(x, 0.25, na.rm = TRUE)),
  q75    = sapply(numeric_data, function(x) quantile(x, 0.75, na.rm = TRUE)),
  skew   = sapply(numeric_data, skewness, na.rm = TRUE),
  kurt   = sapply(numeric_data, kurtosis, na.rm = TRUE)
)

kable(stats_table, digits = 3, caption = "Summary Statistics")

# Missing Values Treatment
data$time <- seq_len(nrow(data))
print("NAs per column:")
print(colSums(is.na(data)))

# Linear interpolation for heart rate
data$hr_mean <- na.approx(data$hr_mean, na.rm = FALSE)
data$hr_mean <- na.locf(data$hr_mean, fromLast = TRUE)
data$hr_mean <- na.locf(data$hr_mean)

# Variable Renaming
colnames(data)[7] <- "Y" 
names(data)[names(data) %in% c("bg_mean","insulin_sum","carbs_sum","hr_mean","steps_sum","cals_sum")] <- 
  c("x1","x2","x3","x4","x5","x6")

# Reordering columns for standard convention
data <- data[, c("time", "x1", "x2", "x3", "x4", "x5", "x6", "Y")]