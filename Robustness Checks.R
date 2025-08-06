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

## Basic Logistic Model
initial_logistic_model <- crash_flag ~ Treatment + r + log_duration + 
  log_length + workzone_week_day + datetime + log_aadt + #SEC_NHS + 
  log_intersections + NUM_LANES + speed_limit + HourlyWindSpeed + 
  HourlyDryBulbTemperature + HourlyPrecipitation + processed_avg_speed + 
  workzone_month + workzone_year

quadratic_logistic_model <- crash_flag ~ Treatment + r + r_2 + log_duration + 
  log_length + workzone_week_day + datetime + log_aadt + #SEC_NHS + 
  log_intersections + NUM_LANES + speed_limit + HourlyWindSpeed + 
  HourlyDryBulbTemperature + HourlyPrecipitation + processed_avg_speed + 
  workzone_month + workzone_year


################################################################################
# Estimations for different bandwidths linear and quadratic
################################################################################

quadratic_equations = list()
ecuaciones = list()
p_r2 = list()
chi_2 = list()
chi_df = list()
chi_sig = list()

for (week_num in (c(1:10))){
  
  path = paste("data/complete_dataset_", week_num, "_weeks.csv", sep = "")
  workzone_df <- read.csv(path)
  
  workzone_df$workzone_month <- factor(workzone_df$workzone_month)
  workzone_df$workzone_year <- factor(workzone_df$workzone_year)
  workzone_df$r_2 <- (workzone_df$r)^2

  nombre_eq = paste('eq_complete_model', week_num, sep = '_')
  
  model_result <- bigglm(
    initial_logistic_model, 
    data = workzone_df, 
    family = binomial()
  )
  
  ecuaciones[nombre_eq] <- list(coeftest(model_result))
  
  
  
  if(week_num >1){
    quadratic_equations[nombre_eq] <- list(coeftest(bigglm(
      quadratic_logistic_model, 
      data = workzone_df, 
      family = binomial()
      )))
  }
  rm(workzone_df)
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
  omit = 'workzone_(month|year)[0-9]+',
  out = 'results/multiple_bandwidths_linear_model.tex'
)

stargazer(
  quadratic_equations, 
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
  omit = 'workzone_(month|year)[0-9]+',
  out = 'results/multiple_bandwidths_quadratic_model.tex'
)


################################################################################
# Estimations for different bandwidths - AIC and BIC
################################################################################

quadratic_equations = list()
aic_linear = list()
bic_linear = list()
aic_quadratic = list()
bic_quadratic = list()

for (week_num in (c(2:8))){
  
  path = paste("data/complete_dataset_", week_num, "_weeks.csv", sep = "")
  workzone_df <- read.csv(path)
  
  workzone_df$workzone_month <- factor(workzone_df$workzone_month)
  workzone_df$workzone_year <- factor(workzone_df$workzone_year)
  workzone_df$r_2 <- (workzone_df$r)^2
  

  if(week_num<8){  
    linear_model <- glm(initial_logistic_model, data = workzone_df, family = binomial())
    quadratic_model <- glm(quadratic_logistic_model, data = workzone_df, family = binomial())
    
    aic_linear[week_num] <- AIC(linear_model)
    bic_linear[week_num] <- BIC(linear_model)
  
    aic_quadratic[week_num] <- AIC(quadratic_model)
    bic_quadratic[week_num] <- BIC(quadratic_model)
  } else {
    linear_model <- bigglm(initial_logistic_model, data = workzone_df, family = binomial())
    quadratic_model <- bigglm(quadratic_logistic_model, data = workzone_df, family = binomial())
    k <- length(coef(linear_model))
    n <- nrow(workzone_df)
    
    pred_probs_linear <- predict(linear_model, newdata = workzone_df, type = "response")
    
    loglik_linear <- sum(
      workzone_df$crash_flag * log(pred_probs_linear) +
        (1 - workzone_df$crash_flag) * log(1 - pred_probs_linear)
    )
    
    aic_linear[week_num] <- AIC(linear_model)
    bic_linear[week_num] <- -2 * loglik_linear + k * log(n)
    
    pred_probs_quadratic <- predict(quadratic_model, newdata = workzone_df, type = "response")
    
    loglik_quadratic <- sum(
      workzone_df$crash_flag * log(pred_probs_quadratic) +
        (1 - workzone_df$crash_flag) * log(1 - pred_probs_quadratic)
    )
    
    aic_quadratic[week_num] <- AIC(quadratic_model)
    bic_quadratic[week_num] <- -2 * loglik_quadratic + k * log(n)
  }
  rm(workzone_df)
}

print(aic_linear)
print(aic_quadratic)

print(bic_linear)
print(bic_quadratic)

################################################################################
################################################################################
# Estimations for different bandwidths - R2 and chi2
################################################################################


p_r2 = list()
chi_2 = list()
chi_df = list()
chi_sig = list()

for (week_num in (c(2:10))){
  
  path = paste("data/complete_dataset_", week_num, "_weeks.csv", sep = "")
  workzone_df <- read.csv(path)
  
  workzone_df$workzone_month <- factor(workzone_df$workzone_month)
  workzone_df$workzone_year <- factor(workzone_df$workzone_year)
  workzone_df$r_2 <- (workzone_df$r)^2
  
  
  if(week_num<8){  
    model_result <- glm(initial_logistic_model, data = workzone_df, family = binomial())
    
    logLik_model <- as.numeric(logLik(model_result))
    logLik_null <- as.numeric(logLik(update(model_result, . ~ 1)))
    
    r2_mcfadden <- round(1 - (logLik_model / logLik_null), 3)
    lr_stat <- -2 * (logLik_null - logLik_model)
    df_diff <- attr(logLik(model_result), "df") - attr(logLik(update(model_result, . ~ 1)), "df")
    p_sig <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
    
  } else {
    big_model <- bigglm(initial_logistic_model, data = workzone_df, family = binomial())
    
    # Must provide newdata explicitly!
    pred_probs_full <- predict(big_model, newdata = workzone_df, type = "response")
    
    null_model <- bigglm(crash_flag ~ 1, data = workzone_df, family = binomial())
    
    pred_probs_null <- predict(null_model, newdata = workzone_df, type = "response")
    
    y <- workzone_df$crash_flag
    
    logLik_full <- sum(y * log(pred_probs_full) + (1 - y) * log(1 - pred_probs_full))
    logLik_null <- sum(y * log(pred_probs_null) + (1 - y) * log(1 - pred_probs_null))
    
    
    r2_mcfadden <- 1 - (logLik_full / logLik_null)
    lr_stat <- -2 * (logLik_null - logLik_full)
    df_diff <- length(coef(big_model)) - 1
    p_sig <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
    
  }
  
  p_r2[paste('eq_', week_num, sep = '')] = r2_mcfadden
  chi_2[paste('eq_', week_num, sep = '')] = lr_stat
  chi_df[paste('eq_', week_num, sep = '')] = df_diff
  chi_sig[paste('eq_', week_num, sep = '')] = p_sig
  
  rm(workzone_df)
}

summary_df <- data.frame(
  Model = names(p_r2),
  McFadden_R2 = unlist(p_r2),
  Chi2 = unlist(chi_2),
  DF = unlist(chi_df),
  sig= unlist(chi_sig)
)

stargazer(
  summary_df, 
  type = "html", 
  title = "Model Comparison: McFadden R² and Chi-squared Statistics",
  summary = FALSE,
  digits = 3,
  align = TRUE,
  out = 'results/multiple_bandwiths_results.html'
)


################################################################################
# Estimations for different time intervals
################################################################################

ecuaciones = list()

for (time_int in (c(60, 120, 180, 240))){
  
  path = paste("data/complete_dataset_", time_int, "_mins.csv", sep = "")
  workzone_df <- read.csv(path)
  nombre_eq = paste('eq_complete_model', time_int, sep = '_')
  
  model <- bigglm(initial_logistic_model, data = workzone_df, family = binomial())
  
  ecuaciones[nombre_eq] <- list(coeftest(model))
  
  rm(workzone_df)
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
  omit = 'workzone_(month|year)[0-9]+',
  out = 'results/multiple_time_interval_model.tex'
)

################################################################################
# Estimations for different time intervals - R2 and chi2
################################################################################

ecuaciones = list()

for (time_int in (c(60, 120, 180, 240))){
  
  path = paste("data/complete_dataset_", time_int, "_mins.csv", sep = "")
  workzone_df <- read.csv(path)
  nombre_eq = paste('eq_complete_model', time_int, sep = '_')
  
  model_result <- glm(initial_logistic_model, data = workzone_df, family = binomial())
  
  logLik_model <- as.numeric(logLik(model_result))
  logLik_null <- as.numeric(logLik(update(model_result, . ~ 1)))
  
  r2_mcfadden <- round(1 - (logLik_model / logLik_null), 3)
  lr_stat <- -2 * (logLik_null - logLik_model)
  df_diff <- attr(logLik(model_result), "df") - attr(logLik(update(model_result, . ~ 1)), "df")
  p_sig <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
  
  p_r2[nombre_eq] = r2_mcfadden
  chi_2[nombre_eq] = lr_stat
  chi_df[nombre_eq] = df_diff
  chi_sig[nombre_eq] = p_sig
  
  rm(workzone_df)
}
summary_df <- data.frame(
  Model = names(p_r2),
  McFadden_R2 = unlist(p_r2),
  Chi2 = unlist(chi_2),
  DF = unlist(chi_df),
  sig= unlist(chi_sig)
)

stargazer(
  summary_df, 
  type = "html", 
  title = "Model Comparison: McFadden R² and Chi-squared Statistics",
  summary = FALSE,
  digits = 3,
  align = TRUE,
  out = 'results/multiple_time_intervals_results.html'
)

################################################################################

################################################################################
