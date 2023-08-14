#' -----------------------------------------------------------------------------
#' Date created: March 7, 2023
#' Author: Sheena Martenies
#' Contact: smarte4@illinois.edu
#' 
#' Description: Quantile G-computation models for "mixture" effects
#' for continuous bmi z-scores: MADRES cohort
#' 
#' First trimester exposures
#' 
#' Note: this code relies heavily on the BKMR tutorial from Alexander Keil:
#' https://cran.r-project.org/web/packages/qgcomp/vignettes/qgcomp-vignette.html
#' 
#' Alexander P. Keil, Jessie P. Buckley, Katie M. O’Brien, Kelly K. Ferguson, 
#' Shanshan Zhao, Alexandra J. White. A quantile-based g-computation approach to 
#' addressing the effects of exposure mixtures. https://doi.org/10.1289/EHP5838
#' ----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(haven)
library(readxl)
library(ggcorrplot)
library(qgcomp)

set.seed(123)

#' [Box Health] directory where the data MUST live
#' These are PHI and must be protected as such!
#' Do not save any object containing these data in any other directory!!
box_path <- "/Users/sheenamartenies/Library/CloudStorage/Box-Box/[Box Health - External] Echo_oif4d"

scale_this <- function(x) as.vector(scale(x))

#' -----------------------------------------------------------------------------
#' Read in MADRES data set
#' Note: these are sensitive data and can ONLY live in the [Box Health] folder
#' -----------------------------------------------------------------------------

md_data <- read_csv(paste0(box_path, "/clean_data/md_analytic_data_clean.csv")) %>%
  mutate(cohort = "md") %>%
  na.omit()
glimpse(md_data)

gcomp_psi_results <- tibble()
gcomp_weights_results <- tibble()

#' -----------------------------------------------------------------------------
#' Fit the g-computation models using the GIS based traffic measures
#' -----------------------------------------------------------------------------

scale_exp <- c("pm", "no2", "o3", "rhavg", "tavg_c", "dist_roads_m",
               "ndvi_500m", "dist_parks_m", "park_area_500m", "tree_cover_500m",
               "impervious_500m", "pct_hh_poverty", "pct_less_hs_edu")

#' Exposures
exps <- c("pm", "no2", "o3", "rhavg", "tavg_c", "dist_roads_m",
          "ndvi_500m", "dist_parks_m", "park_area_500m", "park_count_500m",
          "tree_cover_500m", "impervious_500m", "lila", 
          "pct_hh_poverty", "pct_less_hs_edu")

#'Covariates
covars <- c("maternal_age", "maternal_partner", 
            "ed_no_hs", "ed_hs", "ed_aa", "ed_4yr", "smokesh", "gestsmoking", 
            "ppbmi_underweight", "ppbmi_overweight", "ppbmi_obese_i", "ppbmi_obese_ii", "ppbmi_obese_iii", 
            "concep_spring", "concep_summer", "concep_fall", "male")

#' --------------------------------------
#' Outcomes at delivery
#' out_period = 0
#' --------------------------------------

df_0 <- md_data %>%
  filter(exp_period == "t1" & out_period == 0) %>%
  select(all_of(exps), all_of(covars), zbmi) %>%
  na.omit() %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Fit the quantile g computation models adjusted for all covariates
fit.qgcomp.c_0 <- qgcomp(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + dist_roads_m + 
                           ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                           tree_cover_500m + impervious_500m + lila + 
                           pct_hh_poverty + pct_less_hs_edu +
                           maternal_age + maternal_partner +
                           ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                           ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                           concep_spring + concep_summer + concep_fall + male,
                         data = df_0, expnms = exps, family = gaussian(),
                         bayes = T)
fit.qgcomp.b_0 <- qgcomp.boot(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + dist_roads_m + 
                                ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                                tree_cover_500m + impervious_500m + lila + 
                                pct_hh_poverty + pct_less_hs_edu +
                                maternal_age + maternal_partner +
                                ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                                ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                                concep_spring + concep_summer + concep_fall + male, 
                              data = df_0, expnms = exps, family = gaussian(),
                              bayes = T, seed = 123, B = 1000)

#' Summarize results
temp <- tibble(exposure_period = "t1",
               time_point = 0,
               outcome = "zbmi",
               group = "gis variables",
               model = c("qgcomp", "qgcomp.boot"),
               psi = c(fit.qgcomp.c_0$psi, fit.qgcomp.b_0$psi),
               psi_ll = c(fit.qgcomp.c_0$ci[1], fit.qgcomp.b_0$ci[1]),
               psi_ul = c(fit.qgcomp.c_0$ci[2], fit.qgcomp.b_0$ci[2]))
gcomp_psi_results <- bind_rows(gcomp_psi_results, temp)

temp2 <- tibble(exposure_period = "t1",
                time_point = 0,
                outcome = "zbmi",
                group = "gis variables",
                exposure = names(fit.qgcomp.c_0$fit$coefficients[2:16]),
                coefficients = fit.qgcomp.c_0$fit$coefficients[2:16])
wgts2 <- tibble(exposure = c(names(fit.qgcomp.c_0$pos.weights), names(fit.qgcomp.c_0$neg.weights)),
                weights = c(fit.qgcomp.c_0$pos.weights, fit.qgcomp.c_0$neg.weights * -1))
temp2 <- left_join(temp2, wgts2, by = "exposure") %>%
  mutate(abs_weights = abs(weights))
gcomp_weights_results <- bind_rows(gcomp_weights_results, temp2)

#' Plot
plot_name <- paste0("QGCOMP_", "t1", "_0", "_zbmi_md.jpeg")
graphics.off()
jpeg(file = here::here("figs/qgcomp_plots/", plot_name))
plot(fit.qgcomp.b_0)
dev.off()

#' Save the model
model_name <- "QGComp_t1_0_md.rdata"
save(fit.qgcomp.c_0, fit.qgcomp.b_0,
     file = here::here("results/qgcomp_results", model_name))

#' --------------------------------------
#' Outcomes at 6 months
#' out_period = 6
#' --------------------------------------

df_6 <- md_data %>%
  filter(exp_period == "t1" & out_period == 6) %>%
  select(all_of(exps), all_of(covars), zbmi) %>%
  na.omit() %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Fit the quantile g computation models adjusted for all covariates
fit.qgcomp.c_6 <- qgcomp(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + dist_roads_m + 
                           ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                           tree_cover_500m + impervious_500m + lila + 
                           pct_hh_poverty + pct_less_hs_edu +
                           maternal_age + maternal_partner +
                           ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                           ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                           concep_spring + concep_summer + concep_fall + male,
                         data = df_6, expnms = exps, family = gaussian(),
                         bayes = T)
fit.qgcomp.b_6 <- qgcomp.boot(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + dist_roads_m + 
                                ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                                tree_cover_500m + impervious_500m + lila + 
                                pct_hh_poverty + pct_less_hs_edu +
                                maternal_age + maternal_partner +
                                ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                                ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                                concep_spring + concep_summer + concep_fall + male, 
                              data = df_6, expnms = exps, family = gaussian(),
                              bayes = T, seed = 123, B = 1000)

#' Summarize results
temp <- tibble(exposure_period = "t1",
               time_point = 6,
               outcome = "zbmi",
               group = "gis variables",
               model = c("qgcomp", "qgcomp.boot"),
               psi = c(fit.qgcomp.c_6$psi, fit.qgcomp.b_6$psi),
               psi_ll = c(fit.qgcomp.c_6$ci[1], fit.qgcomp.b_6$ci[1]),
               psi_ul = c(fit.qgcomp.c_6$ci[2], fit.qgcomp.b_6$ci[2]))
gcomp_psi_results <- bind_rows(gcomp_psi_results, temp)

temp2 <- tibble(exposure_period = "t1",
                time_point = 6,
                outcome = "zbmi",
                group = "gis variables",
                exposure = names(fit.qgcomp.c_6$fit$coefficients[2:16]),
                coefficients = fit.qgcomp.c_6$fit$coefficients[2:16])
wgts2 <- tibble(exposure = c(names(fit.qgcomp.c_6$pos.weights), names(fit.qgcomp.c_6$neg.weights)),
                weights = c(fit.qgcomp.c_6$pos.weights, fit.qgcomp.c_6$neg.weights * -1))
temp2 <- left_join(temp2, wgts2, by = "exposure") %>%
  mutate(abs_weights = abs(weights))

gcomp_weights_results <- bind_rows(gcomp_weights_results, temp2)

#' Plot
plot_name <- paste0("QGCOMP_", "t1", "_6", "_zbmi_md.jpeg")
graphics.off()
jpeg(file = here::here("figs/qgcomp_plots/", plot_name))
plot(fit.qgcomp.b_6)
dev.off()

#' Save the model
model_name <- "QGComp_t1_6_md.rdata"
save(fit.qgcomp.c_6, fit.qgcomp.b_6,
     file = here::here("results/qgcomp_results", model_name))

#' --------------------------------------
#' Outcomes at 12 months
#' out_period = 12
#' --------------------------------------

df_12 <- md_data %>%
  filter(exp_period == "t1" & out_period == 12) %>%
  select(all_of(exps), all_of(covars), zbmi) %>%
  na.omit() %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Fit the quantile g computation models adjusted for all covariates
fit.qgcomp.c_12 <- qgcomp(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + dist_roads_m + 
                           ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                           tree_cover_500m + impervious_500m + lila + 
                           pct_hh_poverty + pct_less_hs_edu +
                           maternal_age + maternal_partner +
                           ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                           ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                           concep_spring + concep_summer + concep_fall + male,
                         data = df_12, expnms = exps, family = gaussian(),
                         bayes = T)
fit.qgcomp.b_12 <- qgcomp.boot(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + dist_roads_m + 
                                ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                                tree_cover_500m + impervious_500m + lila + 
                                pct_hh_poverty + pct_less_hs_edu +
                                maternal_age + maternal_partner +
                                ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                                ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                                concep_spring + concep_summer + concep_fall + male, 
                              data = df_12, expnms = exps, family = gaussian(),
                              bayes = T, seed = 123, B = 1000)

#' Summarize results
temp <- tibble(exposure_period = "t1",
               time_point = 12,
               outcome = "zbmi",
               group = "gis variables",
               model = c("qgcomp", "qgcomp.boot"),
               psi = c(fit.qgcomp.c_12$psi, fit.qgcomp.b_12$psi),
               psi_ll = c(fit.qgcomp.c_12$ci[1], fit.qgcomp.b_12$ci[1]),
               psi_ul = c(fit.qgcomp.c_12$ci[2], fit.qgcomp.b_12$ci[2]))
gcomp_psi_results <- bind_rows(gcomp_psi_results, temp)

temp2 <- tibble(exposure_period = "t1",
                time_point = 12,
                outcome = "zbmi",
                group = "gis variables",
                exposure = names(fit.qgcomp.c_12$fit$coefficients[2:16]),
                coefficients = fit.qgcomp.c_12$fit$coefficients[2:16])
wgts2 <- tibble(exposure = c(names(fit.qgcomp.c_12$pos.weights), names(fit.qgcomp.c_12$neg.weights)),
                weights = c(fit.qgcomp.c_12$pos.weights, fit.qgcomp.c_12$neg.weights * -1))
temp2 <- left_join(temp2, wgts2, by = "exposure") %>%
  mutate(abs_weights = abs(weights))

gcomp_weights_results <- bind_rows(gcomp_weights_results, temp2)

#' Plot
plot_name <- paste0("QGCOMP_", "t1", "_12", "_zbmi_md.jpeg")
graphics.off()
jpeg(file = here::here("figs/qgcomp_plots/", plot_name))
plot(fit.qgcomp.b_12)
dev.off()

#' Save the model
model_name <- "QGComp_t1_12_md.rdata"
save(fit.qgcomp.c_12, fit.qgcomp.b_12,
     file = here::here("results/qgcomp_results", model_name))

#' --------------------------------------
#' Outcomes at 24 months
#' out_period = 24
#' --------------------------------------

df_24 <- md_data %>%
  filter(exp_period == "t1" & out_period == 24) %>%
  select(all_of(exps), all_of(covars), zbmi) %>%
  na.omit() %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Fit the quantile g computation models adjusted for all covariates
fit.qgcomp.c_24 <- qgcomp(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + dist_roads_m + 
                            ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                            tree_cover_500m + impervious_500m + lila + 
                            pct_hh_poverty + pct_less_hs_edu +
                            maternal_age + maternal_partner +
                            ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                            ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                            concep_spring + concep_summer + concep_fall + male,
                          data = df_24, expnms = exps, family = gaussian(),
                          bayes = T)
fit.qgcomp.b_24 <- qgcomp.boot(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + dist_roads_m + 
                                 ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                                 tree_cover_500m + impervious_500m + lila + 
                                 pct_hh_poverty + pct_less_hs_edu +
                                 maternal_age + maternal_partner +
                                 ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                                 ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                                 concep_spring + concep_summer + concep_fall + male, 
                               data = df_24, expnms = exps, family = gaussian(),
                               bayes = T, seed = 123, B = 1000)

#' Summarize results
temp <- tibble(exposure_period = "t1",
               time_point = 24,
               outcome = "zbmi",
               group = "gis variables",
               model = c("qgcomp", "qgcomp.boot"),
               psi = c(fit.qgcomp.c_24$psi, fit.qgcomp.b_24$psi),
               psi_ll = c(fit.qgcomp.c_24$ci[1], fit.qgcomp.b_24$ci[1]),
               psi_ul = c(fit.qgcomp.c_24$ci[2], fit.qgcomp.b_24$ci[2]))
gcomp_psi_results <- bind_rows(gcomp_psi_results, temp)

temp2 <- tibble(exposure_period = "t1",
                time_point = 24,
                outcome = "zbmi",
                group = "gis variables",
                exposure = names(fit.qgcomp.c_24$fit$coefficients[2:16]),
                coefficients = fit.qgcomp.c_24$fit$coefficients[2:16])
wgts2 <- tibble(exposure = c(names(fit.qgcomp.c_24$pos.weights), names(fit.qgcomp.c_24$neg.weights)),
                weights = c(fit.qgcomp.c_24$pos.weights, fit.qgcomp.c_24$neg.weights * -1))
temp2 <- left_join(temp2, wgts2, by = "exposure") %>%
  mutate(abs_weights = abs(weights))

gcomp_weights_results <- bind_rows(gcomp_weights_results, temp2)

#' Plot
plot_name <- paste0("QGCOMP_", "t1", "_24", "_zbmi_md.jpeg")
graphics.off()
jpeg(file = here::here("figs/qgcomp_plots/", plot_name))
plot(fit.qgcomp.b_24)
dev.off()

#' Save the model
model_name <- "QGComp_t1_24_md.rdata"
save(fit.qgcomp.c_24, fit.qgcomp.b_24,
     file = here::here("results/qgcomp_results", model_name))

#' -----------------------------------------------------------------------------
#' Fit the g-computation models using the spatiotemporal based traffic measures
#' -----------------------------------------------------------------------------

scale_exp <- c("pm", "no2", "o3", "rhavg", "tavg_c", "nox_pred",
               "ndvi_500m", "dist_parks_m", "park_area_500m", "tree_cover_500m",
               "impervious_500m", "pct_hh_poverty", "pct_less_hs_edu")

#' Exposures
exps <- c("pm", "no2", "o3", "rhavg", "tavg_c", "nox_pred",
          "ndvi_500m", "dist_parks_m", "park_area_500m", "park_count_500m",
          "tree_cover_500m", "impervious_500m", "lila", 
          "pct_hh_poverty", "pct_less_hs_edu")

#'Covariates
covars <- c("maternal_age", "maternal_partner", 
            "ed_no_hs", "ed_hs", "ed_aa", "ed_4yr", "smokesh", "gestsmoking", 
            "ppbmi_underweight", "ppbmi_overweight", "ppbmi_obese_i", "ppbmi_obese_ii", "ppbmi_obese_iii", 
            "concep_spring", "concep_summer", "concep_fall", "male")

#' --------------------------------------
#' Outcomes at delivery
#' out_period = 0
#' --------------------------------------

df_0 <- md_data %>%
  filter(exp_period == "t1" & out_period == 0) %>%
  select(all_of(exps), all_of(covars), zbmi) %>%
  na.omit() %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Fit the quantile g computation models adjusted for all covariates
fit.qgcomp.c_0 <- qgcomp(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + nox_pred + 
                           ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                           tree_cover_500m + impervious_500m + lila + 
                           pct_hh_poverty + pct_less_hs_edu +
                           maternal_age + maternal_partner +
                           ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                           ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                           concep_spring + concep_summer + concep_fall + male,
                         data = df_0, expnms = exps, family = gaussian(),
                         bayes = T)
fit.qgcomp.b_0 <- qgcomp.boot(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + nox_pred + 
                                ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                                tree_cover_500m + impervious_500m + lila + 
                                pct_hh_poverty + pct_less_hs_edu +
                                maternal_age + maternal_partner +
                                ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                                ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                                concep_spring + concep_summer + concep_fall + male, 
                              data = df_0, expnms = exps, family = gaussian(),
                              bayes = T, seed = 123, B = 1000)

#' Summarize results
temp <- tibble(exposure_period = "t1",
               time_point = 0,
               outcome = "zbmi",
               group = "st variables",
               model = c("qgcomp", "qgcomp.boot"),
               psi = c(fit.qgcomp.c_0$psi, fit.qgcomp.b_0$psi),
               psi_ll = c(fit.qgcomp.c_0$ci[1], fit.qgcomp.b_0$ci[1]),
               psi_ul = c(fit.qgcomp.c_0$ci[2], fit.qgcomp.b_0$ci[2]))
gcomp_psi_results <- bind_rows(gcomp_psi_results, temp)

temp2 <- tibble(exposure_period = "t1",
                time_point = 0,
                outcome = "zbmi",
                group = "st variables",
                exposure = names(fit.qgcomp.c_0$fit$coefficients[2:16]),
                coefficients = fit.qgcomp.c_0$fit$coefficients[2:16])
wgts2 <- tibble(exposure = c(names(fit.qgcomp.c_0$pos.weights), names(fit.qgcomp.c_0$neg.weights)),
                weights = c(fit.qgcomp.c_0$pos.weights, fit.qgcomp.c_0$neg.weights * -1))
temp2 <- left_join(temp2, wgts2, by = "exposure") %>%
  mutate(abs_weights = abs(weights))

gcomp_weights_results <- bind_rows(gcomp_weights_results, temp2)

#' Plot
plot_name <- paste0("QGCOMP_", "t1", "_0", "_zbmi_st_md.jpeg")
graphics.off()
jpeg(file = here::here("figs/qgcomp_plots/", plot_name))
plot(fit.qgcomp.b_0)
dev.off()

#' Save the model
model_name <- "QGComp_t1_0_st_md.rdata"
save(fit.qgcomp.c_0, fit.qgcomp.b_0,
     file = here::here("results/qgcomp_results", model_name))

#' --------------------------------------
#' Outcomes at 6 months
#' out_period = 6
#' --------------------------------------

df_6 <- md_data %>%
  filter(exp_period == "t1" & out_period == 6) %>%
  select(all_of(exps), all_of(covars), zbmi) %>%
  na.omit() %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Fit the quantile g computation models adjusted for all covariates
fit.qgcomp.c_6 <- qgcomp(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + nox_pred + 
                           ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                           tree_cover_500m + impervious_500m + lila + 
                           pct_hh_poverty + pct_less_hs_edu +
                           maternal_age + maternal_partner +
                           ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                           ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                           concep_spring + concep_summer + concep_fall + male,
                         data = df_6, expnms = exps, family = gaussian(),
                         bayes = T)
fit.qgcomp.b_6 <- qgcomp.boot(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + nox_pred + 
                                ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                                tree_cover_500m + impervious_500m + lila + 
                                pct_hh_poverty + pct_less_hs_edu +
                                maternal_age + maternal_partner +
                                ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                                ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                                concep_spring + concep_summer + concep_fall + male, 
                              data = df_6, expnms = exps, family = gaussian(),
                              bayes = T, seed = 123, B = 1000)

#' Summarize results
temp <- tibble(exposure_period = "t1",
               time_point = 6,
               outcome = "zbmi",
               group = "st variables",
               model = c("qgcomp", "qgcomp.boot"),
               psi = c(fit.qgcomp.c_6$psi, fit.qgcomp.b_6$psi),
               psi_ll = c(fit.qgcomp.c_6$ci[1], fit.qgcomp.b_6$ci[1]),
               psi_ul = c(fit.qgcomp.c_6$ci[2], fit.qgcomp.b_6$ci[2]))
gcomp_psi_results <- bind_rows(gcomp_psi_results, temp)

temp2 <- tibble(exposure_period = "t1",
                time_point = 6,
                outcome = "zbmi",
                group = "st variables",
                exposure = names(fit.qgcomp.c_6$fit$coefficients[2:16]),
                coefficients = fit.qgcomp.c_6$fit$coefficients[2:16])
wgts2 <- tibble(exposure = c(names(fit.qgcomp.c_6$pos.weights), names(fit.qgcomp.c_6$neg.weights)),
                weights = c(fit.qgcomp.c_6$pos.weights, fit.qgcomp.c_6$neg.weights * -1))
temp2 <- left_join(temp2, wgts2, by = "exposure") %>%
  mutate(abs_weights = abs(weights))

gcomp_weights_results <- bind_rows(gcomp_weights_results, temp2)

#' Plot
plot_name <- paste0("QGCOMP_", "t1", "_6", "_zbmi_st_md.jpeg")
graphics.off()
jpeg(file = here::here("figs/qgcomp_plots/", plot_name))
plot(fit.qgcomp.b_6)
dev.off()

#' Save the model
model_name <- "QGComp_t1_6_st_md.rdata"
save(fit.qgcomp.c_6, fit.qgcomp.b_6,
     file = here::here("results/qgcomp_results", model_name))

#' --------------------------------------
#' Outcomes at 12 months
#' out_period = 12
#' --------------------------------------

df_12 <- md_data %>%
  filter(exp_period == "t1" & out_period == 12) %>%
  select(all_of(exps), all_of(covars), zbmi) %>%
  na.omit() %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Fit the quantile g computation models adjusted for all covariates
fit.qgcomp.c_12 <- qgcomp(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + nox_pred + 
                            ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                            tree_cover_500m + impervious_500m + lila + 
                            pct_hh_poverty + pct_less_hs_edu +
                            maternal_age + maternal_partner +
                            ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                            ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                            concep_spring + concep_summer + concep_fall + male,
                          data = df_12, expnms = exps, family = gaussian(),
                          bayes = T)
fit.qgcomp.b_12 <- qgcomp.boot(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + nox_pred + 
                                 ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                                 tree_cover_500m + impervious_500m + lila + 
                                 pct_hh_poverty + pct_less_hs_edu +
                                 maternal_age + maternal_partner +
                                 ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                                 ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                                 concep_spring + concep_summer + concep_fall + male, 
                               data = df_12, expnms = exps, family = gaussian(),
                               bayes = T, seed = 123, B = 1000)

#' Summarize results
temp <- tibble(exposure_period = "t1",
               time_point = 12,
               outcome = "zbmi",
               group = "st variables",
               model = c("qgcomp", "qgcomp.boot"),
               psi = c(fit.qgcomp.c_12$psi, fit.qgcomp.b_12$psi),
               psi_ll = c(fit.qgcomp.c_12$ci[1], fit.qgcomp.b_12$ci[1]),
               psi_ul = c(fit.qgcomp.c_12$ci[2], fit.qgcomp.b_12$ci[2]))
gcomp_psi_results <- bind_rows(gcomp_psi_results, temp)

temp2 <- tibble(exposure_period = "t1",
                time_point = 12,
                outcome = "zbmi",
                group = "st variables",
                exposure = names(fit.qgcomp.c_12$fit$coefficients[2:16]),
                coefficients = fit.qgcomp.c_12$fit$coefficients[2:16])
wgts2 <- tibble(exposure = c(names(fit.qgcomp.c_12$pos.weights), names(fit.qgcomp.c_12$neg.weights)),
                weights = c(fit.qgcomp.c_12$pos.weights, fit.qgcomp.c_12$neg.weights * -1))
temp2 <- left_join(temp2, wgts2, by = "exposure") %>%
  mutate(abs_weights = abs(weights))

gcomp_weights_results <- bind_rows(gcomp_weights_results, temp2)

#' Plot
plot_name <- paste0("QGCOMP_", "t1", "_12", "_zbmi_st_md.jpeg")
graphics.off()
jpeg(file = here::here("figs/qgcomp_plots/", plot_name))
plot(fit.qgcomp.b_12)
dev.off()

#' Save the model
model_name <- "QGComp_t1_12_st_md.rdata"
save(fit.qgcomp.c_12, fit.qgcomp.b_12,
     file = here::here("results/qgcomp_results", model_name))

#' --------------------------------------
#' Outcomes at 24 months
#' out_period = 24
#' --------------------------------------

df_24 <- md_data %>%
  filter(exp_period == "t1" & out_period == 24) %>%
  select(all_of(exps), all_of(covars), zbmi) %>%
  na.omit() %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Fit the quantile g computation models adjusted for all covariates
fit.qgcomp.c_24 <- qgcomp(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + nox_pred + 
                            ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                            tree_cover_500m + impervious_500m + lila + 
                            pct_hh_poverty + pct_less_hs_edu +
                            maternal_age + maternal_partner +
                            ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                            ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                            concep_spring + concep_summer + concep_fall + male,
                          data = df_24, expnms = exps, family = gaussian(),
                          bayes = T)
fit.qgcomp.b_24 <- qgcomp.boot(zbmi ~ pm + no2 + o3 + rhavg + tavg_c + nox_pred + 
                                 ndvi_500m + dist_parks_m + park_area_500m + park_count_500m +
                                 tree_cover_500m + impervious_500m + lila + 
                                 pct_hh_poverty + pct_less_hs_edu +
                                 maternal_age + maternal_partner +
                                 ed_no_hs + ed_hs + ed_aa + ed_4yr + smokesh + gestsmoking +
                                 ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + 
                                 concep_spring + concep_summer + concep_fall + male, 
                               data = df_24, expnms = exps, family = gaussian(),
                               bayes = T, seed = 123, B = 1000)

#' Summarize results
temp <- tibble(exposure_period = "t1",
               time_point = 24,
               outcome = "zbmi",
               group = "st variables",
               model = c("qgcomp", "qgcomp.boot"),
               psi = c(fit.qgcomp.c_24$psi, fit.qgcomp.b_24$psi),
               psi_ll = c(fit.qgcomp.c_24$ci[1], fit.qgcomp.b_24$ci[1]),
               psi_ul = c(fit.qgcomp.c_24$ci[2], fit.qgcomp.b_24$ci[2]))
gcomp_psi_results <- bind_rows(gcomp_psi_results, temp)

temp2 <- tibble(exposure_period = "t1",
                time_point = 24,
                outcome = "zbmi",
                group = "st variables",
                exposure = names(fit.qgcomp.c_24$fit$coefficients[2:16]),
                coefficients = fit.qgcomp.c_24$fit$coefficients[2:16])
wgts2 <- tibble(exposure = c(names(fit.qgcomp.c_24$pos.weights), names(fit.qgcomp.c_24$neg.weights)),
                weights = c(fit.qgcomp.c_24$pos.weights, fit.qgcomp.c_24$neg.weights * -1))
temp2 <- left_join(temp2, wgts2, by = "exposure") %>%
  mutate(abs_weights = abs(weights))

gcomp_weights_results <- bind_rows(gcomp_weights_results, temp2)

#' Plot
plot_name <- paste0("QGCOMP_", "t1", "_24", "_zbmi_st_md.jpeg")
graphics.off()
jpeg(file = here::here("figs/qgcomp_plots/", plot_name))
plot(fit.qgcomp.b_24)
dev.off()

#' Save the model
model_name <- "QGComp_t1_24_st_md.rdata"
save(fit.qgcomp.c_24, fit.qgcomp.b_24,
     file = here::here("results/qgcomp_results", model_name))

#' -----------------------------------------------------------------------------
#' Write out the final results
#' -----------------------------------------------------------------------------

write_csv(gcomp_psi_results, file = here::here("results/qgcomp_results/psi_est", "qgcomp_t1_psi_estimates_md.csv"))
write_csv(gcomp_weights_results, file = here::here("results/qgcomp_results/weights", "qgcomp_t1_weights_md.csv"))