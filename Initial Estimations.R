rm(list=ls())

# Load required libraries
library(plyr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(haven)
library(labelled)
library(data.table)
library(plm)
library(xlsx)
library(readxl)
library(panelr)
library(foreign)
library(stargazer)
library(rdrobust)
library(rdd)
library(rms)
library(logistf)
library(data.table)
library(biglm)
library(lmtest)

setwd("~/Library/CloudStorage/Box-Box/Work zone [TX]")

## Import the workzone df
path = "data/complete_dataset_6_weeks.csv"
#workzone_df <- fread(path)
workzone_df <- read.csv(path)

treated_workzones <- workzone_df %>%
  filter(Treatment == 1) 

workzone_df$workzone_month <- factor(workzone_df$workzone_month)
workzone_df$workzone_year <- factor(workzone_df$workzone_year)


## Basic Logistic Model
num_df <- workzone_df[, sapply(workzone_df, is.numeric), drop = FALSE]

# Compute summary statistics
desc_table <- data.frame(
  Variable = names(num_df),
  N = sapply(num_df, function(x) sum(!is.na(x))),
  Mean = sapply(num_df, function(x) round(mean(x, na.rm = TRUE), 3)),
  SD = sapply(num_df, function(x) round(sd(x, na.rm = TRUE), 3)),
  Min = sapply(num_df, function(x) round(min(x, na.rm = TRUE), 3)),
  Q25 = sapply(num_df, function(x) round(quantile(x, 0.25, na.rm = TRUE), 3)),
  Q75 = sapply(num_df, function(x) round(quantile(x, 0.75, na.rm = TRUE), 3)),
  Max = sapply(num_df, function(x) round(max(x, na.rm = TRUE), 3))
)

# Save to CSV
write.csv(
  desc_table, 
  paste('results/descriptive_statistics_6_weeks_model.csv', sep = ''), 
  row.names = FALSE
)


## This version of the code does not account for NHS flag
initial_logistic_model <- crash_flag ~ Treatment + r + log_duration + 
  log_length + workzone_week_day + datetime + log_aadt + 
  log_intersections + NUM_LANES + speed_limit + HourlyWindSpeed + 
  HourlyDryBulbTemperature + HourlyPrecipitation + processed_avg_speed + 
  workzone_month + workzone_year


## 6 weeks for the base model
fitted_model <- glm(initial_logistic_model, data = workzone_df, family = binomial())
 
# Fit the null model (intercept only)
null_model <- glm(crash_flag ~ 1, data = workzone_df, family = binomial())

# Log-likelihoods
ll_full <- as.numeric(logLik(fitted_model))
ll_null <- as.numeric(logLik(null_model))

# McFadden's R²
pseud_r2 <- round(1 - (ll_full / ll_null), 3)

chi2 <- -2 * (ll_null - ll_full)
df_chi2 <- length(coef(fitted_model)) - length(coef(null_model))
p_val <- pchisq(chi2, df = df_chi2, lower.tail = FALSE)
sig_star <- ifelse(p_val < 0.001, "***",
                   ifelse(p_val < 0.01, "**",
                          ifelse(p_val < 0.05, "*",
                                 ifelse(p_val < 0.1, ".", ""))))
chi2_text <- paste0(round(chi2, 3),
                    sig_star,
                    " (df = ", df_chi2, ")"
                    )

# Create a summary table
stargazer(
  fitted_model, 
  type = "latex", 
  title = "Results from the Firth logistic regression model.",
  dep.var.labels = "Crash Occurrence",
  covariate.labels = c("Work zone presence", 
                       "Week",
                       "Duration (second; log)", 
                       "Length (meter; log)",
                       "Weekday of week",
                       "Daytime of day",
                       "AADT (vehicles per day; log)",
                       #"NHS major roads",
                       "Number of intersections (log)",
                       "Lane counts = 1",
                       "Speed limit",
                       "Average wind speed (mph)",
                       "Average temperature (F)",
                       "Average precipitation (inch)",
                       "Actual speed processed (mph)"),
  digits = 3,
  omit = 'workzone_(month|year)[0-9]+',
  align = TRUE,
  keep.stat = c('n'),
  add.lines = list(
    c("R2", pseud_r2),
    c("Chi2", chi2_text)
  ),
  out = 'results/6_week_logistic_model.tex'
)
  

################################################################################
# Estimations for different AADTs (TXDOT quartiles)
################################################################################

min_aadt = 0
ecuaciones = list()
p_r2_list = list()
chi_2_list = list()

for (aadt in (c(11580, 93060, 280662))){
  
  nombre_eq = paste('eq', (exp(min_aadt)), aadt, sep = '_')
  
  temp_treated_workzones <- treated_workzones %>%
    filter((log_aadt > min_aadt) & (log_aadt <= log(aadt))) %>%
    select(GLOBALID) %>%
    distinct() %>%
    pull()
  
  temp_workzone <- workzone_df %>%
    filter(GLOBALID %in% temp_treated_workzones)

  num_df <- temp_workzone[, sapply(temp_workzone, is.numeric), drop = FALSE]
  
  # Compute summary statistics
  desc_table <- data.frame(
    Variable = names(num_df),
    N = sapply(num_df, function(x) sum(!is.na(x))),
    Mean = sapply(num_df, function(x) round(mean(x, na.rm = TRUE), 3)),
    SD = sapply(num_df, function(x) round(sd(x, na.rm = TRUE), 3)),
    Min = sapply(num_df, function(x) round(min(x, na.rm = TRUE), 3)),
    Q25 = sapply(num_df, function(x) round(quantile(x, 0.25, na.rm = TRUE), 3)),
    Q75 = sapply(num_df, function(x) round(quantile(x, 0.75, na.rm = TRUE), 3)),
    Max = sapply(num_df, function(x) round(max(x, na.rm = TRUE), 3))
  )
  
  # Save to CSV
  write.csv(
    desc_table, 
    paste('results/descriptive_statistics_', nombre_eq, '.csv', sep = ''), 
    row.names = FALSE
  )
  
  fitted_model <- glm(initial_logistic_model, data = temp_workzone, family = binomial())
  
  ecuaciones[[nombre_eq]] <- fitted_model
  
  null_model <- glm(crash_flag ~ 1, data = temp_workzone, family = binomial())
  
  ll_full <- as.numeric(logLik(fitted_model))
  ll_null <- as.numeric(logLik(null_model))
  
  pseud_r2 <- round(1 - (ll_full / ll_null), 3)
  
  chi2 <- -2 * (ll_null - ll_full)
  df_chi2 <- length(coef(fitted_model)) - length(coef(null_model))
  p_val <- pchisq(chi2, df = df_chi2, lower.tail = FALSE)
  sig_star <- ifelse(p_val < 0.001, "***",
                     ifelse(p_val < 0.01, "**",
                            ifelse(p_val < 0.05, "*",
                                   ifelse(p_val < 0.1, ".", ""))))
  chi2_text <- paste0(round(chi2, 3),
                      sig_star,
                      " (df = ", df_chi2, ")"
  )
  
  # Store for later
  p_r2_list <- c(p_r2_list, pseud_r2)
  chi_2_list <- c(chi_2_list, chi2_text)
  
  min_aadt = log(aadt)  
}

stargazer(
  ecuaciones, 
  style = 'default', 
  type = 'latex', 
  title = "Results for strata AADT",
  dep.var.labels = "Crash Occurrence",
  covariate.labels = c("Work zone presence", 
                       "Week",
                       "Duration (second; log)", 
                       "Length (meter; log)",
                       "Weekday of week",
                       "Daytime of day",
                       "AADT (vehicles per day; log)",
                       #"NHS major roads",
                       "Number of intersections (log)",
                       "Lane counts = 1",
                       "Speed limit",
                       "Average wind speed (mph)",
                       "Average temperature (F)",
                       "Average precipitation (inch)",
                       "Actual speed processed (mph)"),
  digits = 3,
  align = TRUE,
  keep.stat = c('n'),
  add.lines = list(
    c("R2", as.character(p_r2_list)),
    c("Chi2", as.character(chi_2_list))
  ),
  omit = 'workzone_(month|year)[0-9]+',
  out = 'results/multiple_aadt_model_quartiles.tex'
)


################################################################################
# Estimations for different Lenghts (TXDOT mean)
################################################################################

ecuaciones = list()
p_r2_list = list()
chi_2_list = list()
min_length = 0
for (max_length in (c(5400,175000))){
  
  nombre_eq = paste('eq', 'lenghts', max_length, sep = '_')
  
  temp_treated_workzones <- treated_workzones %>%
    filter((log_length >= min_length) & (log_length < log(max_length))) %>%
    select(GLOBALID) %>%
    distinct() %>%
    pull()
  
  temp_workzone <- workzone_df %>%
    filter(GLOBALID %in% temp_treated_workzones)
  
  num_df <- temp_workzone[, sapply(temp_workzone, is.numeric), drop = FALSE]
  
  # Compute summary statistics
  desc_table <- data.frame(
    Variable = names(num_df),
    N = sapply(num_df, function(x) sum(!is.na(x))),
    Mean = sapply(num_df, function(x) round(mean(x, na.rm = TRUE), 3)),
    SD = sapply(num_df, function(x) round(sd(x, na.rm = TRUE), 3)),
    Min = sapply(num_df, function(x) round(min(x, na.rm = TRUE), 3)),
    Q25 = sapply(num_df, function(x) round(quantile(x, 0.25, na.rm = TRUE), 3)),
    Q75 = sapply(num_df, function(x) round(quantile(x, 0.75, na.rm = TRUE), 3)),
    Max = sapply(num_df, function(x) round(max(x, na.rm = TRUE), 3))
  )
  
  # Save to CSV
  write.csv(
    desc_table, 
    paste('results/descriptive_statistics_', nombre_eq, '.csv', sep = ''), 
    row.names = FALSE
  )
  
  fitted_model <- glm(initial_logistic_model, data = temp_workzone, family = binomial())
  
  ecuaciones[[nombre_eq]] <- fitted_model
  
  null_model <- glm(crash_flag ~ 1, data = temp_workzone, family = binomial())
  
  ll_full <- as.numeric(logLik(fitted_model))
  ll_null <- as.numeric(logLik(null_model))
  
  pseud_r2 <- round(1 - (ll_full / ll_null), 3)
  
  chi2 <- -2 * (ll_null - ll_full)
  df_chi2 <- length(coef(fitted_model)) - length(coef(null_model))
  p_val <- pchisq(chi2, df = df_chi2, lower.tail = FALSE)
  sig_star <- ifelse(p_val < 0.001, "***",
                     ifelse(p_val < 0.01, "**",
                            ifelse(p_val < 0.05, "*",
                                   ifelse(p_val < 0.1, ".", ""))))
  chi2_text <- paste0(round(chi2, 3),
                      sig_star,
                      " (df = ", df_chi2, ")"
  )
  
  # Store for later
  p_r2_list <- c(p_r2_list, pseud_r2)
  chi_2_list <- c(chi_2_list, chi2_text)
  
  min_length = log(max_length)
}

stargazer(
  ecuaciones, 
  style = 'default', 
  type = 'latex', 
  title = "Results for strata AADT",
  dep.var.labels = "Crash Occurrence",
  covariate.labels = c("Work zone presence", 
                       "Week",
                       "Duration (second; log)", 
                       "Length (meter; log)",
                       "Weekday of week",
                       "Daytime of day",
                       "AADT (vehicles per day; log)",
                       #"NHS major roads",
                       "Number of intersections (log)",
                       "Lane counts = 1",
                       "Speed limit",
                       "Average wind speed (mph)",
                       "Average temperature (F)",
                       "Average precipitation (inch)",
                       "Actual speed processed (mph)"),
  digits = 3,
  align = TRUE,
  keep.stat = c('n'),
  add.lines = list(
    c("R2", as.character(p_r2_list)),
    c("Chi2", as.character(chi_2_list))
  ),
  omit = 'workzone_(month|year)[0-9]+',
  out = 'results/multiple_length_results_mean.tex'
)


################################################################################
# Estimations for different daytime conditions
################################################################################

model <- crash_flag ~ Treatment + r + log_duration + 
  log_length + workzone_week_day + log_aadt + 
  log_intersections + NUM_LANES + speed_limit + HourlyWindSpeed + 
  HourlyDryBulbTemperature + HourlyPrecipitation + processed_avg_speed + 
  workzone_month + workzone_year

ecuaciones = list()
p_r2_list = list()
chi_2_list = list()

for (daytime_flag in (c(0,1))){
  
  nombre_eq = paste('eq', 'daytime', daytime_flag, sep = '_')

  temp_workzone <- workzone_df %>%
    filter(datetime == daytime_flag)
  
  num_df <- temp_workzone[, sapply(temp_workzone, is.numeric), drop = FALSE]
  
  # Compute summary statistics
  desc_table <- data.frame(
    Variable = names(num_df),
    N = sapply(num_df, function(x) sum(!is.na(x))),
    Mean = sapply(num_df, function(x) round(mean(x, na.rm = TRUE), 3)),
    SD = sapply(num_df, function(x) round(sd(x, na.rm = TRUE), 3)),
    Min = sapply(num_df, function(x) round(min(x, na.rm = TRUE), 3)),
    Q25 = sapply(num_df, function(x) round(quantile(x, 0.25, na.rm = TRUE), 3)),
    Q75 = sapply(num_df, function(x) round(quantile(x, 0.75, na.rm = TRUE), 3)),
    Max = sapply(num_df, function(x) round(max(x, na.rm = TRUE), 3))
  )
  
  # Save to CSV
  write.csv(
    desc_table, 
    paste('results/descriptive_statistics_', nombre_eq, '.csv', sep = ''), 
    row.names = FALSE
  )

  fitted_model <- glm(model, data = temp_workzone, family = binomial())
  
  ecuaciones[[nombre_eq]] <- fitted_model
  
  null_model <- glm(crash_flag ~ 1, data = temp_workzone, family = binomial())
  
  ll_full <- as.numeric(logLik(fitted_model))
  ll_null <- as.numeric(logLik(null_model))
  
  pseud_r2 <- round(1 - (ll_full / ll_null), 3)
  
  chi2 <- -2 * (ll_null - ll_full)
  df_chi2 <- length(coef(fitted_model)) - length(coef(null_model))
  p_val <- pchisq(chi2, df = df_chi2, lower.tail = FALSE)
  sig_star <- ifelse(p_val < 0.001, "***",
                     ifelse(p_val < 0.01, "**",
                            ifelse(p_val < 0.05, "*",
                                   ifelse(p_val < 0.1, ".", ""))))
  chi2_text <- paste0(round(chi2, 3),
                      sig_star,
                      " (df = ", df_chi2, ")"
  )
  
  # Store for later
  p_r2_list <- c(p_r2_list, pseud_r2)
  chi_2_list <- c(chi_2_list, chi2_text)
  
}

stargazer(
  ecuaciones, 
  style = 'default', 
  type = 'latex', 
  title = "Results for strata AADT",
  dep.var.labels = "Crash Occurrence",
  covariate.labels = c("Work zone presence", 
                       "Week",
                       "Duration (second; log)", 
                       "Length (meter; log)",
                       "Weekday of week",
                       #"Daytime of day",
                       "AADT (vehicles per day; log)",
                       #"NHS major roads",
                       "Number of intersections (log)",
                       "Lane counts = 1",
                       "Speed limit",
                       "Average wind speed (mph)",
                       "Average temperature (F)",
                       "Average precipitation (inch)",
                       "Actual speed processed (mph)"),
  digits = 3,
  align = TRUE,
  keep.stat = c('n'),
  add.lines = list(
    c("R2", as.character(p_r2_list)),
    c("Chi2", as.character(chi_2_list))
  ),
  omit = 'workzone_(month|year)[0-9]+',
  out = 'results/multiple_daytime_results.tex'
)

################################################################################
# Estimations for different speed ranges
################################################################################

ecuaciones = list()
p_r2_list = list()
chi_2_list = list()
min_speed = -41
for (max_speed in (c(-12.1, 11, 38))){
  
  nombre_eq = paste('eq', 'speed', max_speed, sep = '_')
  
  temp_treated_workzones <- treated_workzones %>%
    filter((processed_avg_speed >= min_speed) & (processed_avg_speed < max_speed)) %>%
    select(GLOBALID) %>%
    distinct() %>%
    pull()
  
  temp_workzone <- workzone_df %>%
    filter(GLOBALID %in% temp_treated_workzones)  
  
  num_df <- temp_workzone[, sapply(temp_workzone, is.numeric), drop = FALSE]
  
  # Compute summary statistics
  desc_table <- data.frame(
    Variable = names(num_df),
    N = sapply(num_df, function(x) sum(!is.na(x))),
    Mean = sapply(num_df, function(x) round(mean(x, na.rm = TRUE), 3)),
    SD = sapply(num_df, function(x) round(sd(x, na.rm = TRUE), 3)),
    Min = sapply(num_df, function(x) round(min(x, na.rm = TRUE), 3)),
    Q25 = sapply(num_df, function(x) round(quantile(x, 0.25, na.rm = TRUE), 3)),
    Q75 = sapply(num_df, function(x) round(quantile(x, 0.75, na.rm = TRUE), 3)),
    Max = sapply(num_df, function(x) round(max(x, na.rm = TRUE), 3))
  )
  
  # Save to CSV
  write.csv(
    desc_table, 
    paste('results/descriptive_statistics_', nombre_eq, '.csv', sep = ''), 
    row.names = FALSE
  )
  
  fitted_model <- glm(initial_logistic_model, data = temp_workzone, family = binomial())
  
  ecuaciones[[nombre_eq]] <- fitted_model
  
  null_model <- glm(crash_flag ~ 1, data = temp_workzone, family = binomial())
  
  ll_full <- as.numeric(logLik(fitted_model))
  ll_null <- as.numeric(logLik(null_model))
  
  pseud_r2 <- round(1 - (ll_full / ll_null), 3)
  
  chi2 <- -2 * (ll_null - ll_full)
  df_chi2 <- length(coef(fitted_model)) - length(coef(null_model))
  p_val <- pchisq(chi2, df = df_chi2, lower.tail = FALSE)
  sig_star <- ifelse(p_val < 0.001, "***",
                     ifelse(p_val < 0.01, "**",
                            ifelse(p_val < 0.05, "*",
                                   ifelse(p_val < 0.1, ".", ""))))
  chi2_text <- paste0(round(chi2, 3),
                      sig_star,
                      " (df = ", df_chi2, ")"
  )
  
  # Store for later
  p_r2_list <- c(p_r2_list, pseud_r2)
  chi_2_list <- c(chi_2_list, chi2_text)
  
  min_speed = max_speed
}

stargazer(
  ecuaciones, 
  style = 'default', 
  type = 'latex', 
  title = "Results for strata AADT",
  dep.var.labels = "Crash Occurrence",
  covariate.labels = c("Work zone presence", 
                       "Week",
                       "Duration (second; log)", 
                       "Length (meter; log)",
                       "Weekday of week",
                       "Daytime of day",
                       "AADT (vehicles per day; log)",
                       #"NHS major roads",
                       "Number of intersections (log)",
                       "Lane counts = 1",
                       "Speed limit",
                       "Average wind speed (mph)",
                       "Average temperature (F)",
                       "Average precipitation (inch)",
                       "Actual speed processed (mph)"),
  digits = 3,
  align = TRUE,
  keep.stat = c('n'),
  add.lines = list(
    c("R2", as.character(p_r2_list)),
    c("Chi2", as.character(chi_2_list))
  ),
  omit = 'workzone_(month|year)[0-9]+',
  out = 'results/multiple_speed_ranges.tex'
)



