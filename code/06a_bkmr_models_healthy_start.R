#' -----------------------------------------------------------------------------
#' Date created: June 2, 2023
#' Date updated: June 23, 2020
#' Author: Sheena Martenies
#' Contact: smarte4@illinois.edu
#' 
#' Description: BKMR models between exposures and continuous BMI z-score for
#' the Healthy Start cohort
#' 
#' Note: this code relies heavily on the BKMR tutorial from Jennifer Bobb:
#' https://jenfb.github.io/bkmr/overview.html
#' 
#' Bobb, JF, Valeri L, Claus Henn B, Christiani DC, Wright RO, Mazumdar M, 
#' Godleski JJ, Coull BA. Bayesian Kernel Machine Regression for Estimating the 
#' Health Effects of Multi-Pollutant Mixtures. Biostatistics 16, no. 3 (July 1, 
#' 2015): 493–508. doi:10.1093/biostatistics/kxu058
#' ----------------------------------------------------------------------------

#' Load required libraries
library(ggplot2)
library(viridis)
library(ggthemes)
library(tidyverse)
library(lubridate)
library(writexl)
library(readxl)
library(bkmr)

set.seed(123)

#' [Box Health] directory where the data MUST live
#' These are PHI and must be protected as such!
#' Do not save any object containing these data in any other director!!
box_path <- "/Users/sheenamartenies/Library/CloudStorage/Box-Box/[Box Health - External] Echo_oif4d"

scale_this <- function(x) as.vector(scale(x))

#' -----------------------------------------------------------------------------
#' Read in HS data sets
#' Note: these are sensitive data and can ONLY live in the [Box Health] folder
#' -----------------------------------------------------------------------------

hs_data <- read_csv(paste0(box_path, "/clean_data/hs_analytic_data_clean.csv")) %>%
  mutate(cohort = "hs") %>%
  na.omit()

#' -----------------------------------------------------------------------------
#' Fit the BKMR models using the GIS based traffic measures
#' -----------------------------------------------------------------------------

scale_exp <- c("pm", "no2", "o3", "rhavg", "tavg_c", "dist_roads_m", 
               "ndvi_500m", "dist_parks_m", "park_area_500m", 
               "tree_cover_500m","impervious_500m", "pct_hh_poverty", "pct_less_hs_edu")

#' --------------------------------------
#' Outcomes at delivery
#' out_period = 0
#' --------------------------------------

df <- filter(hs_data, exp_period == "pregnancy" & out_period == 0) %>%
  na.omit()

y <- df$zbmi

Z <- select(df, pm, o3, no2, rhavg, tavg_c,
            dist_roads_m,
            ndvi_500m, dist_parks_m, park_area_500m,
            park_count_500m, tree_cover_500m, impervious_500m,
            lila, pct_hh_poverty, pct_less_hs_edu) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this) %>%
  as.matrix

X <- select(df,
            maternal_age, maternal_partner,
            ed_no_hs, ed_hs, ed_aa, ed_4yr, smokesh, gestsmoking,
            ppbmi_underweight, ppbmi_overweight, ppbmi_obese_i, ppbmi_obese_ii, ppbmi_obese_iii,
            concep_spring, concep_summer, concep_fall, male) %>%
  mutate(maternal_age = scale_this(maternal_age)) %>%
  as.matrix

#' Fit the BKMR model using the scaled variables
#' Using default priors here
#' Since BMI z-scores are already standardized, not scaling this one

fit_bkmr <- kmbayes(y = y, Z = Z, X = X, iter = 20000,
                    verbose = FALSE, varsel = TRUE)

#' Save current results
model_name <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_hs.rdata")
save(df, y, Z, X, fit_bkmr, file = here::here("results/bkmr_results", model_name))

#' Trace plots for inspection
graphics.off()
jpeg(file = here::here("figs/bkmr_plots/", paste0("BKMR_trace_", "pregnancy",
                                                  "_", 0, "_zbmi_hs.jpeg")))
par(mfrow=c(3,1))
TracePlot(fit = fit_bkmr, par = "beta", comp = 9)
TracePlot(fit = fit_bkmr, par = "sigsq.eps")
TracePlot(fit = fit_bkmr, par = "r", comp = 15)
dev.off()
par(mfrow=c(1,1))

#' Extract the PIPs
fit_pips <- ExtractPIPs(fit_bkmr) %>%
  mutate(exp_period = "pregnancy",
         out_period = 0,
         group = "gis variables")
pips_name <- paste0("BKMR_PIPs_", "pregnancy", "_", 0, "_zbmi_hs.csv")
write_csv(fit_pips, file = here::here("results/bkmr_results/pips", pips_name))

#' Summarize output
#' 1) Plot the predictor-response function
pred_univariate_25 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.25)
pred_univariate_50 <- PredictorResponseUnivar(fit = fit_bkmr)
pred_univariate_75 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.75)

ggplot(pred_univariate_25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 25th percentile")
plot_name2 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_h25_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name2),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 50th percentile")
plot_name1 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_h50_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name1),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 75th percentile")
plot_name3 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_h75_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name3),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Summary statistics for the predictor-response function
o_risk <- OverallRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs = seq(0.25, 0.75, by = 0.05),
                                      q.fixed = 0.50, method = "approx")
ggplot(o_risk, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
plot_name4 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_overall_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name4),
       device = "jpeg", height = 4, width = 6, units = "in")
overall_name <- paste0("BKMR_overall_", "pregnancy", "_", 0, "_zbmi_hs.csv")
write_csv(o_risk, file = here::here("results/bkmr_results/risks", overall_name))

#' 3) Results for individual predictors
s_risk <- SingVarRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs.diff = c(0.25, 0.75),
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")
ggplot(s_risk, aes(variable, est, ymin = est - 1.96*sd,
                          ymax = est + 1.96*sd, col = q.fixed)) +
  geom_pointrange(position = position_dodge(width = 0.75)) +
  coord_flip()
plot_name5 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_singvar_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name5),
       device = "jpeg", height = 6, width = 5, units = "in")
singvar_name <- paste0("BKMR_singvar_", "pregnancy", "_", 0, "_zbmi_hs.csv")
write_csv(s_risk, file = here::here("results/bkmr_results/risks", singvar_name))

#' Saving results
model_name <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_results_hs.rdata")
save(fit_pips, 
     pred_univariate_25, pred_univariate_50, pred_univariate_75,
     o_risk, s_risk, 
     file = here::here("results/bkmr_results", model_name))

#' --------------------------------------
#' Outcomes at 6 months
#' out_period = 6
#' --------------------------------------

df <- filter(hs_data, exp_period == "pregnancy" & out_period == 6) %>%
  na.omit()

y <- df$zbmi

Z <- select(df, pm, o3, no2, rhavg, tavg_c,
            dist_roads_m,
            ndvi_500m, dist_parks_m, park_area_500m,
            park_count_500m, tree_cover_500m, impervious_500m,
            lila, pct_hh_poverty, pct_less_hs_edu) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this) %>%
  as.matrix

X <- select(df,
            maternal_age, maternal_partner,
            ed_no_hs, ed_hs, ed_aa, ed_4yr, smokesh, gestsmoking,
            ppbmi_underweight, ppbmi_overweight, ppbmi_obese_i, ppbmi_obese_ii, ppbmi_obese_iii,
            concep_spring, concep_summer, concep_fall, male) %>%
  mutate(maternal_age = scale_this(maternal_age)) %>%
  as.matrix

#' Fit the BKMR model using the scaled variables
#' Using default priors here
#' Since BMI z-scores are already standardized, not scaling this one

fit_bkmr <- kmbayes(y = y, Z = Z, X = X, iter = 20000,
                    verbose = FALSE, varsel = TRUE)

#' Save current results
model_name <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_hs.rdata")
save(df, y, Z, X, fit_bkmr, file = here::here("results/bkmr_results", model_name))

#' Trace plots for inspection
graphics.off()
jpeg(file = here::here("figs/bkmr_plots/", paste0("BKMR_trace_", "pregnancy",
                                                  "_", 6, "_zbmi_hs.jpeg")))
par(mfrow=c(3,1))
TracePlot(fit = fit_bkmr, par = "beta", comp = 1)
TracePlot(fit = fit_bkmr, par = "sigsq.eps")
TracePlot(fit = fit_bkmr, par = "r", comp = 1)
dev.off()
par(mfrow=c(1,1))

#' Extract the PIPs
fit_pips <- ExtractPIPs(fit_bkmr) %>%
  mutate(exp_period = "pregnancy",
         out_period = 6,
         group = "gis variables")
pips_name <- paste0("BKMR_PIPs_", "pregnancy", "_", 6, "_zbmi_hs.csv")
write_csv(fit_pips, file = here::here("results/bkmr_results/pips", pips_name))

#' Summarize output
#' 1) Plot the predictor-response function
pred_univariate_25 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.25)
pred_univariate_50 <- PredictorResponseUnivar(fit = fit_bkmr)
pred_univariate_75 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.75)

ggplot(pred_univariate_25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 25th percentile")
plot_name2 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_h25_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name2),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 50th percentile")
plot_name1 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_h50_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name1),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 75th percentile")
plot_name3 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_h75_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name3),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Summary statistics for the predictor-response function
o_risk <- OverallRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs = seq(0.25, 0.75, by = 0.05),
                                      q.fixed = 0.50, method = "approx")
ggplot(o_risk, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
plot_name4 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_overall_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name4),
       device = "jpeg", height = 4, width = 6, units = "in")
overall_name <- paste0("BKMR_overall_", "pregnancy", "_", 6, "_zbmi_hs.csv")
write_csv(o_risk, file = here::here("results/bkmr_results/risks", overall_name))

#' 3) Results for individual predictors
s_risk <- SingVarRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs.diff = c(0.25, 0.75),
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")
ggplot(s_risk, aes(variable, est, ymin = est - 1.96*sd,
                          ymax = est + 1.96*sd, col = q.fixed)) +
  geom_pointrange(position = position_dodge(width = 0.75)) +
  coord_flip()
plot_name5 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_singvar_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name5),
       device = "jpeg", height = 6, width = 5, units = "in")
singvar_name <- paste0("BKMR_singvar_", "pregnancy", "_", 6, "_zbmi_hs.csv")
write_csv(s_risk, file = here::here("results/bkmr_results/risks", singvar_name))

#' Saving results
model_name <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_results_hs.rdata")
save(fit_pips, 
     pred_univariate_25, pred_univariate_50, pred_univariate_75,
     o_risk, s_risk,
     file = here::here("results/bkmr_results", model_name))

#' --------------------------------------
#' Outcomes at 12 months
#' out_period = 12
#' --------------------------------------

df <- filter(hs_data, exp_period == "pregnancy" & out_period == 12) %>%
  na.omit()

y <- df$zbmi

Z <- select(df, pm, o3, no2, rhavg, tavg_c, 
            dist_roads_m,
            ndvi_500m, dist_parks_m, park_area_500m, 
            park_count_500m, tree_cover_500m, impervious_500m,
            lila, pct_hh_poverty, pct_less_hs_edu) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this) %>%
  as.matrix

X <- select(df, 
            maternal_age, maternal_partner, 
            ed_no_hs, ed_hs, ed_aa, ed_4yr, smokesh, gestsmoking, 
            ppbmi_underweight, ppbmi_overweight, ppbmi_obese_i, ppbmi_obese_ii, ppbmi_obese_iii, 
            concep_spring, concep_summer, concep_fall, male) %>%
  mutate(maternal_age = scale_this(maternal_age)) %>%
  as.matrix

#' Fit the BKMR model using the scaled variables
#' Using default priors here
#' Since BMI z-scores are already standardized, not scaling this one

fit_bkmr <- kmbayes(y = y, Z = Z, X = X, iter = 20000,
                    verbose = FALSE, varsel = TRUE)

#' Save current results
model_name <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_hs.rdata")
save(df, y, Z, X, fit_bkmr, file = here::here("results/bkmr_results", model_name))

#' Trace plots for inspection
graphics.off()
jpeg(file = here::here("figs/bkmr_plots/", paste0("BKMR_trace_", "pregnancy", 
                                                  "_", 12, "_zbmi_hs.jpeg")))
par(mfrow=c(3,1))
TracePlot(fit = fit_bkmr, par = "beta", comp = 1)
TracePlot(fit = fit_bkmr, par = "sigsq.eps")
TracePlot(fit = fit_bkmr, par = "r", comp = 1)
dev.off()
par(mfrow=c(1,1))

#' Extract the PIPs
fit_pips <- ExtractPIPs(fit_bkmr) %>%
  mutate(exp_period = "pregnancy",
         out_period = 12, 
         group = "gis variables")
pips_name <- paste0("BKMR_PIPs_", "pregnancy", "_", 12, "_zbmi_hs.csv") 
write_csv(fit_pips, file = here::here("results/bkmr_results/pips", pips_name))

#' Summarize output
#' 1) Plot the predictor-response function
pred_univariate_25 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.25)
pred_univariate_50 <- PredictorResponseUnivar(fit = fit_bkmr)
pred_univariate_75 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.75)

ggplot(pred_univariate_25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 25th percentile")
plot_name2 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_h25_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name2),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 50th percentile")
plot_name1 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_h50_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name1),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 75th percentile")
plot_name3 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_h75_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name3),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Summary statistics for the predictor-response function
o_risk <- OverallRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs = seq(0.25, 0.75, by = 0.05),
                                      q.fixed = 0.50, method = "approx")
ggplot(o_risk, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
plot_name4 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_overall_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name4),
       device = "jpeg", height = 4, width = 6, units = "in")
overall_name <- paste0("BKMR_overall_", "pregnancy", "_", 12, "_zbmi_hs.csv") 
write_csv(o_risk, file = here::here("results/bkmr_results/risks", overall_name))

#' 3) Results for individual predictors
s_risk <- SingVarRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")
ggplot(s_risk, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
plot_name5 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_singvar_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name5),
       device = "jpeg", height = 6, width = 5, units = "in")
singvar_name <- paste0("BKMR_singvar_", "pregnancy", "_", 12, "_zbmi_hs.csv") 
write_csv(s_risk, file = here::here("results/bkmr_results/risks", singvar_name))

#' Saving results
model_name <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_results_hs.rdata")
save(fit_pips, 
     pred_univariate_25, pred_univariate_50, pred_univariate_75,
     o_risk, s_risk, 
     file = here::here("results/bkmr_results", model_name))

#' --------------------------------------
#' Outcomes at 24 months
#' out_period = 24
#' --------------------------------------

df <- filter(hs_data, exp_period == "pregnancy" & out_period == 24) %>%
  na.omit()

y <- df$zbmi

Z <- select(df, pm, o3, no2, rhavg, tavg_c, 
            dist_roads_m,
            ndvi_500m, dist_parks_m, park_area_500m, 
            park_count_500m, tree_cover_500m, impervious_500m,
            lila, pct_hh_poverty, pct_less_hs_edu) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this) %>%
  as.matrix

X <- select(df, 
            maternal_age, maternal_partner, 
            ed_no_hs, ed_hs, ed_aa, ed_4yr, smokesh, gestsmoking, 
            ppbmi_underweight, ppbmi_overweight, ppbmi_obese_i, ppbmi_obese_ii, ppbmi_obese_iii, 
            concep_spring, concep_summer, concep_fall, male) %>%
  mutate(maternal_age = scale_this(maternal_age)) %>%
  as.matrix

#' Fit the BKMR model using the scaled variables
#' Using default priors here
#' Since BMI z-scores are already standardized, not scaling this one

fit_bkmr <- kmbayes(y = y, Z = Z, X = X, iter = 20000,
                    verbose = FALSE, varsel = TRUE)

#' Save current results
model_name <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_hs.rdata")
save(df, y, Z, X, fit_bkmr, file = here::here("results/bkmr_results", model_name))

#' Trace plots for inspection
graphics.off()
jpeg(file = here::here("figs/bkmr_plots/", paste0("BKMR_trace_", "pregnancy", 
                                                  "_", 24, "_zbmi_hs.jpeg")))
par(mfrow=c(3,1))
TracePlot(fit = fit_bkmr, par = "beta", comp = 1)
TracePlot(fit = fit_bkmr, par = "sigsq.eps")
TracePlot(fit = fit_bkmr, par = "r", comp = 1)
dev.off()
par(mfrow=c(1,1))

#' Extract the PIPs
fit_pips <- ExtractPIPs(fit_bkmr) %>%
  mutate(exp_period = "pregnancy",
         out_period = 24, 
         group = "gis variables")
pips_name <- paste0("BKMR_PIPs_", "pregnancy", "_", 24, "_zbmi_hs.csv") 
write_csv(fit_pips, file = here::here("results/bkmr_results/pips", pips_name))

#' Summarize output
#' 1) Plot the predictor-response function
pred_univariate_25 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.25)
pred_univariate_50 <- PredictorResponseUnivar(fit = fit_bkmr)
pred_univariate_75 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.75)

ggplot(pred_univariate_25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 25th percentile")
plot_name2 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_h25_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name2),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 50th percentile")
plot_name1 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_h50_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name1),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 75th percentile")
plot_name3 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_h75_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name3),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Summary statistics for the predictor-response function
o_risk <- OverallRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs = seq(0.25, 0.75, by = 0.05),
                                      q.fixed = 0.50, method = "approx")
ggplot(o_risk, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
plot_name4 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_overall_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name4),
       device = "jpeg", height = 4, width = 6, units = "in")
overall_name <- paste0("BKMR_overall_", "pregnancy", "_", 24, "_zbmi_hs.csv") 
write_csv(o_risk, file = here::here("results/bkmr_results/risks", overall_name))

#' 3) Results for individual predictors
s_risk <- SingVarRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")
ggplot(s_risk, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
plot_name5 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_singvar_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name5),
       device = "jpeg", height = 6, width = 5, units = "in")
singvar_name <- paste0("BKMR_singvar_", "pregnancy", "_", 24, "_zbmi_hs.csv") 
write_csv(s_risk, file = here::here("results/bkmr_results/risks", singvar_name))

#' Saving results
model_name <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_results_hs.rdata")
save(fit_pips, 
     pred_univariate_25, pred_univariate_50, pred_univariate_75,
     o_risk, s_risk, 
     file = here::here("results/bkmr_results", model_name))

#' -----------------------------------------------------------------------------
#' Fit the BKMR models using the spatiotemporal traffic measures
#' -----------------------------------------------------------------------------

scale_exp <- c("pm", "no2", "o3", "rhavg", "tavg_c", "bc_pred", 
               "ndvi_500m", "dist_parks_m", "park_area_500m", 
               "tree_cover_500m","impervious_500m", "pct_hh_poverty", "pct_less_hs_edu")

#' --------------------------------------
#' Outcomes at delivery
#' out_period = 0
#' --------------------------------------

df <- filter(hs_data, exp_period == "pregnancy" & out_period == 0) %>%
  na.omit()

y <- df$zbmi

Z <- select(df, pm, o3, no2, rhavg, tavg_c, 
            bc_pred,
            ndvi_500m, dist_parks_m, park_area_500m, 
            park_count_500m, tree_cover_500m, impervious_500m,
            lila, pct_hh_poverty, pct_less_hs_edu) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this) %>%
  as.matrix

X <- select(df, 
            maternal_age, maternal_partner, 
            ed_no_hs, ed_hs, ed_aa, ed_4yr, smokesh, gestsmoking, 
            ppbmi_underweight, ppbmi_overweight, ppbmi_obese_i, ppbmi_obese_ii, ppbmi_obese_iii, 
            concep_spring, concep_summer, concep_fall, male) %>%
  mutate(maternal_age = scale_this(maternal_age)) %>%
  as.matrix

#' Fit the BKMR model using the scaled variables
#' Using default priors here
#' Since BMI z-scores are already standardized, not scaling this one

fit_bkmr <- kmbayes(y = y, Z = Z, X = X, iter = 20000,
                    verbose = FALSE, varsel = TRUE)

#' Save current results
model_name <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_st_hs.rdata")
save(df, y, Z, X, fit_bkmr, file = here::here("results/bkmr_results", model_name))

#' Trace plots for inspection
graphics.off()
jpeg(file = here::here("figs/bkmr_plots/", paste0("BKMR_trace_", "pregnancy", 
                                                  "_", 0, "_zbmi_st_hs.jpeg")))
par(mfrow=c(3,1))
TracePlot(fit = fit_bkmr, par = "beta", comp = 6)
TracePlot(fit = fit_bkmr, par = "sigsq.eps")
TracePlot(fit = fit_bkmr, par = "r", comp = 13)
dev.off()
par(mfrow=c(1,1))

#' Extract the PIPs
fit_pips <- ExtractPIPs(fit_bkmr) %>%
  mutate(exp_period = "pregnancy",
         out_period = 0, 
         group = "st variables")
pips_name <- paste0("BKMR_PIPs_", "pregnancy", "_", 0, "_zbmi_st_hs.csv") 
write_csv(fit_pips, file = here::here("results/bkmr_results/pips", pips_name))

#' Summarize output
#' 1) Plot the predictor-response function
pred_univariate_25 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.25)
pred_univariate_50 <- PredictorResponseUnivar(fit = fit_bkmr)
pred_univariate_75 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.75)

ggplot(pred_univariate_25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 25th percentile")
plot_name2 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_h25_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name2),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 50th percentile")
plot_name1 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_h50_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name1),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 75th percentile")
plot_name3 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_h75_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name3),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Summary statistics for the predictor-response function
o_risk <- OverallRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs = seq(0.25, 0.75, by = 0.05),
                                      q.fixed = 0.50, method = "approx")
ggplot(o_risk, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
plot_name4 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_overall_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name4),
       device = "jpeg", height = 4, width = 6, units = "in")
overall_name <- paste0("BKMR_overall_", "pregnancy", "_", 0, "_zbmi_st_hs.csv") 
write_csv(o_risk, file = here::here("results/bkmr_results/risks", overall_name))

#' 3) Results for individual predictors
s_risk <- SingVarRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")
ggplot(s_risk, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
plot_name5 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_singvar_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name5),
       device = "jpeg", height = 6, width = 5, units = "in")
singvar_name <- paste0("BKMR_singvar_", "pregnancy", "_", 0, "_zbmi_st_hs.csv") 
write_csv(s_risk, file = here::here("results/bkmr_results/risks", singvar_name))

#' Saving results
model_name <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_results_st_hs.rdata")
save(fit_pips, 
     pred_univariate_25, pred_univariate_50, pred_univariate_75,
     o_risk, s_risk, 
     file = here::here("results/bkmr_results", model_name))

#' --------------------------------------
#' Outcomes at 6 months
#' out_period = 6
#' --------------------------------------

df <- filter(hs_data, exp_period == "pregnancy" & out_period == 6) %>%
  na.omit()

y <- df$zbmi

Z <- select(df, pm, o3, no2, rhavg, tavg_c, 
            bc_pred,
            ndvi_500m, dist_parks_m, park_area_500m, 
            park_count_500m, tree_cover_500m, impervious_500m,
            lila, pct_hh_poverty, pct_less_hs_edu) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this) %>%
  as.matrix

X <- select(df, 
            maternal_age, maternal_partner, 
            ed_no_hs, ed_hs, ed_aa, ed_4yr, smokesh, gestsmoking, 
            ppbmi_underweight, ppbmi_overweight, ppbmi_obese_i, ppbmi_obese_ii, ppbmi_obese_iii, 
            concep_spring, concep_summer, concep_fall, male) %>%
  mutate(maternal_age = scale_this(maternal_age)) %>%
  as.matrix

#' Fit the BKMR model using the scaled variables
#' Using default priors here
#' Since BMI z-scores are already standardized, not scaling this one

fit_bkmr <- kmbayes(y = y, Z = Z, X = X, iter = 20000,
                    verbose = FALSE, varsel = TRUE)

#' Save current results
model_name <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_st_hs.rdata")
save(df, y, Z, X, fit_bkmr, file = here::here("results/bkmr_results", model_name))

#' Trace plots for inspection
graphics.off()
jpeg(file = here::here("figs/bkmr_plots/", paste0("BKMR_trace_", "pregnancy", 
                                                  "_", 6, "_zbmi_st_hs.jpeg")))
par(mfrow=c(3,1))
TracePlot(fit = fit_bkmr, par = "beta", comp = 1)
TracePlot(fit = fit_bkmr, par = "sigsq.eps")
TracePlot(fit = fit_bkmr, par = "r", comp = 1)
dev.off()
par(mfrow=c(1,1))

#' Extract the PIPs
fit_pips <- ExtractPIPs(fit_bkmr) %>%
  mutate(exp_period = "pregnancy",
         out_period = 6, 
         group = "st variables")
pips_name <- paste0("BKMR_PIPs_", "pregnancy", "_", 6, "_zbmi_st_hs.csv") 
write_csv(fit_pips, file = here::here("results/bkmr_results/pips", pips_name))

#' Summarize output
#' 1) Plot the predictor-response function
pred_univariate_25 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.25)
pred_univariate_50 <- PredictorResponseUnivar(fit = fit_bkmr)
pred_univariate_75 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.75)

ggplot(pred_univariate_25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 25th percentile")
plot_name2 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_h25_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name2),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 50th percentile")
plot_name1 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_h50_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name1),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 75th percentile")
plot_name3 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_h75_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name3),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Summary statistics for the predictor-response function
o_risk <- OverallRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs = seq(0.25, 0.75, by = 0.05),
                                      q.fixed = 0.50, method = "approx")
ggplot(o_risk, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
plot_name4 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_overall_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name4),
       device = "jpeg", height = 4, width = 6, units = "in")
overall_name <- paste0("BKMR_overall_", "pregnancy", "_", 6, "_zbmi_st_hs.csv") 
write_csv(o_risk, file = here::here("results/bkmr_results/risks", overall_name))

#' 3) Results for individual predictors
s_risk <- SingVarRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")
ggplot(s_risk, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
plot_name5 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_singvar_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name5),
       device = "jpeg", height = 6, width = 5, units = "in")
singvar_name <- paste0("BKMR_singvar_", "pregnancy", "_", 6, "_zbmi_st_hs.csv") 
write_csv(s_risk, file = here::here("results/bkmr_results/risks", singvar_name))

#' Saving results
model_name <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_results_st_hs.rdata")
save(fit_pips,      
     pred_univariate_25, pred_univariate_50, pred_univariate_75,
     o_risk, s_risk,
     file = here::here("results/bkmr_results", model_name))

#' --------------------------------------
#' Outcomes at 12 months
#' out_period = 12
#' --------------------------------------

df <- filter(hs_data, exp_period == "pregnancy" & out_period == 12) %>%
  na.omit()

y <- df$zbmi

Z <- select(df, pm, o3, no2, rhavg, tavg_c, 
            bc_pred,
            ndvi_500m, dist_parks_m, park_area_500m, 
            park_count_500m, tree_cover_500m, impervious_500m,
            lila, pct_hh_poverty, pct_less_hs_edu) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this) %>%
  as.matrix

X <- select(df, 
            maternal_age, maternal_partner, 
            ed_no_hs, ed_hs, ed_aa, ed_4yr, smokesh, gestsmoking, 
            ppbmi_underweight, ppbmi_overweight, ppbmi_obese_i, ppbmi_obese_ii, ppbmi_obese_iii, 
            concep_spring, concep_summer, concep_fall, male) %>%
  mutate(maternal_age = scale_this(maternal_age)) %>%
  as.matrix

#' Fit the BKMR model using the scaled variables
#' Using default priors here
#' Since BMI z-scores are already standardized, not scaling this one

fit_bkmr <- kmbayes(y = y, Z = Z, X = X, iter = 20000,
                    verbose = FALSE, varsel = TRUE)

#' Save current results
model_name <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_st_hs.rdata")
save(df, y, Z, X, fit_bkmr, file = here::here("results/bkmr_results", model_name))

#' Trace plots for inspection
graphics.off()
jpeg(file = here::here("figs/bkmr_plots/", paste0("BKMR_trace_", "pregnancy", 
                                                  "_", 12, "_zbmi_st_hs.jpeg")))
par(mfrow=c(3,1))
TracePlot(fit = fit_bkmr, par = "beta", comp = 1)
TracePlot(fit = fit_bkmr, par = "sigsq.eps")
TracePlot(fit = fit_bkmr, par = "r", comp = 1)
dev.off()
par(mfrow=c(1,1))

#' Extract the PIPs
fit_pips <- ExtractPIPs(fit_bkmr) %>%
  mutate(exp_period = "pregnancy",
         out_period = 12, 
         group = "st variables")
pips_name <- paste0("BKMR_PIPs_", "pregnancy", "_", 12, "_zbmi_st_hs.csv") 
write_csv(fit_pips, file = here::here("results/bkmr_results/pips", pips_name))

#' Summarize output
#' 1) Plot the predictor-response function
pred_univariate_25 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.25)
pred_univariate_50 <- PredictorResponseUnivar(fit = fit_bkmr)
pred_univariate_75 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.75)

ggplot(pred_univariate_25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 25th percentile")
plot_name2 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_h25_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name2),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 50th percentile")
plot_name1 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_h50_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name1),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 75th percentile")
plot_name3 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_h75_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name3),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Summary statistics for the predictor-response function
o_risk <- OverallRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs = seq(0.25, 0.75, by = 0.05),
                                      q.fixed = 0.50, method = "approx")
ggplot(o_risk, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
plot_name4 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_overall_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name4),
       device = "jpeg", height = 4, width = 6, units = "in")
overall_name <- paste0("BKMR_overall_", "pregnancy", "_", 12, "_zbmi_st_hs.csv") 
write_csv(o_risk, file = here::here("results/bkmr_results/risks", overall_name))

#' 3) Results for individual predictors
s_risk <- SingVarRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")
ggplot(s_risk, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
plot_name5 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_singvar_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name5),
       device = "jpeg", height = 6, width = 5, units = "in")
singvar_name <- paste0("BKMR_singvar_", "pregnancy", "_", 12, "_zbmi_st_hs.csv") 
write_csv(s_risk, file = here::here("results/bkmr_results/risks", singvar_name))

#' Saving results
model_name <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_results_st_hs.rdata")
save(fit_pips, 
     pred_univariate_25, pred_univariate_50, pred_univariate_75,
     o_risk, s_risk, 
     file = here::here("results/bkmr_results", model_name))

#' --------------------------------------
#' Outcomes at 24 months
#' out_period = 24
#' --------------------------------------

df <- filter(hs_data, exp_period == "pregnancy" & out_period == 24) %>%
  na.omit()

y <- df$zbmi

Z <- select(df, pm, o3, no2, rhavg, tavg_c, 
            bc_pred,
            ndvi_500m, dist_parks_m, park_area_500m, 
            park_count_500m, tree_cover_500m, impervious_500m,
            lila, pct_hh_poverty, pct_less_hs_edu) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this) %>%
  as.matrix

X <- select(df, 
            maternal_age, maternal_partner, 
            ed_no_hs, ed_hs, ed_aa, ed_4yr, smokesh, gestsmoking, 
            ppbmi_underweight, ppbmi_overweight, ppbmi_obese_i, ppbmi_obese_ii, ppbmi_obese_iii, 
            concep_spring, concep_summer, concep_fall, male) %>%
  mutate(maternal_age = scale_this(maternal_age)) %>%
  as.matrix

#' Fit the BKMR model using the scaled variables
#' Using default priors here
#' Since BMI z-scores are already standardized, not scaling this one

fit_bkmr <- kmbayes(y = y, Z = Z, X = X, iter = 20000,
                    verbose = FALSE, varsel = TRUE)

#' Save current results
model_name <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_st_hs.rdata")
save(df, y, Z, X, fit_bkmr, file = here::here("results/bkmr_results", model_name))

#' Trace plots for inspection
graphics.off()
jpeg(file = here::here("figs/bkmr_plots/", paste0("BKMR_trace_", "pregnancy", 
                                                  "_", 24, "_zbmi_st_hs.jpeg")))
par(mfrow=c(3,1))
TracePlot(fit = fit_bkmr, par = "beta", comp = 1)
TracePlot(fit = fit_bkmr, par = "sigsq.eps")
TracePlot(fit = fit_bkmr, par = "r", comp = 1)
dev.off()
par(mfrow=c(1,1))

#' Extract the PIPs
fit_pips <- ExtractPIPs(fit_bkmr) %>%
  mutate(exp_period = "pregnancy",
         out_period = 24, 
         group = "st variables")
pips_name <- paste0("BKMR_PIPs_", "pregnancy", "_", 24, "_zbmi_st_hs.csv") 
write_csv(fit_pips, file = here::here("results/bkmr_results/pips", pips_name))

#' Summarize output
#' 1) Plot the predictor-response function
pred_univariate_25 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.25)
pred_univariate_50 <- PredictorResponseUnivar(fit = fit_bkmr)
pred_univariate_75 <- PredictorResponseUnivar(fit = fit_bkmr, q.fixed = 0.75)

ggplot(pred_univariate_25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 25th percentile")
plot_name2 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_h25_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name2),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 50th percentile")
plot_name1 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_h50_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name1),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred_univariate_75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  ylab("h(z) holding other variables at 75th percentile")
plot_name3 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_h75_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name3),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Summary statistics for the predictor-response function
o_risk <- OverallRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X,
                                      qs = seq(0.25, 0.75, by = 0.05),
                                      q.fixed = 0.50, method = "approx")
ggplot(o_risk, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
plot_name4 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_overall_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name4),
       device = "jpeg", height = 4, width = 6, units = "in")
overall_name <- paste0("BKMR_overall_", "pregnancy", "_", 24, "_zbmi_st_hs.csv") 
write_csv(o_risk, file = here::here("results/bkmr_results/risks", overall_name))

#' 3) Results for individual predictors
s_risk <- SingVarRiskSummaries(fit = fit_bkmr, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")
ggplot(s_risk, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
plot_name5 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_singvar_st_hs.jpeg")
ggsave(filename = here::here("figs/bkmr_plots", plot_name5),
       device = "jpeg", height = 6, width = 5, units = "in")
singvar_name <- paste0("BKMR_singvar_", "pregnancy", "_", 24, "_zbmi_st_hs.csv") 
write_csv(s_risk, file = here::here("results/bkmr_results/risks", singvar_name))

#' Saving results
model_name <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_results_st_hs.rdata")
save(fit_pips, 
     pred_univariate_25, pred_univariate_50, pred_univariate_75,
     o_risk, s_risk, 
     file = here::here("results/bkmr_results", model_name))