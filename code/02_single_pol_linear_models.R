#' -----------------------------------------------------------------------------
#' Date created: November 2, 2022
#' Date updated: January 5, 2023
#' Author: Sheena Martenies
#' Contact: smarte4@illinois.edu
#' 
#' Description: Linear models between exposures and continuous BMI by cohort
#' ----------------------------------------------------------------------------

#' Load required libraries
library(ggplot2)
library(viridis)
library(ggthemes)
library(tidyverse)
library(lubridate)
library(writexl)
library(readxl)
library(lme4)

#' Some coordinate reference systems
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm_13N <- "+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
ll_nad83 <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ll_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#' [Box Health] directory where the data MUST live
#' These are PHI and must be protected as such!
#' Do not save any object containing these data in any other directory!!
box_path <- "/Users/sheenamartenies/Library/CloudStorage/Box-Box/[Box Health - External] Echo_oif4d"

scale_this <- function(x) as.vector(scale(x))

tab_results <- function(lm, exp_period, outcome, out_period, group) {
  temp <- data.frame(exposure_period = exp_period,
                     time_point = out_period,
                     outcome = outcome,
                     exposure = rownames(summary(lm)$coefficients)[2],
                     group = group,
                     beta = summary(lm)$coefficients[2,1],
                     beta.se = summary(lm)$coefficients[2,2],
                     p.value = summary(lm)$coefficients[2,4])
  return(temp)
}

#' -----------------------------------------------------------------------------
#' Read in HS and MADRES data sets
#' Note: these are sensitive data and can ONLY live in the [Box Health] folder
#' -----------------------------------------------------------------------------

hs_data <- read_csv(paste0(box_path, "/clean_data/hs_analytic_data_clean.csv")) %>%
  mutate(cohort = "hs") 

md_data <- read_csv(paste0(box_path, "/clean_data/md_analytic_data_clean.csv")) %>%
  mutate(cohort = "md") 

c_vars <- "+ maternal_age + maternal_partner + ed_no_hs + ed_hs + ed_aa + ed_4yr + ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + smokesh + gestsmoking + concep_spring + concep_summer + concep_fall + male"

#' -----------------------------------------------------------------------------
#' Linear regression models for each exposure and outcome pair: Healthy Start
#' -----------------------------------------------------------------------------

#' Linear regression models
lm_results <- data.frame()

#' Continuous exposures
scale_exp <- c("pm", "no2", "o3", "rhavg", "tavg_c", "bc_pred", 
               "dist_roads_m", "ndvi_500m", "dist_parks_m", "park_area_500m", 
               "tree_cover_500m","impervious_500m", "pct_hh_poverty", "pct_less_hs_edu")

#' ----------------------------------------
#' Model controls
#' 
#' Exposures: pm, o3, no2, rhavg, tavg_c, bc_pred, dist_roads_m,
#' ndvi_500m, dist_parks_m, park_area_500m, 
#' park_count_500m, tree_cover_500m, impervious_500m,
#' lila, pct_hh_poverty, pct_less_hs_edu
#' 
#' Outcomes: zbmi
#' 
#' Exposure periods: pregnancy, t1, t2, t3
#' 
#' Outcome periods: 0, 6, 12, 24
#' ----------------------------------------

ex <- "pm"
out <- "zbmi"
ex_p <- "pregnancy"
out_p <- 0

#' ----------------------------------------
#' Regression code
#' ----------------------------------------

loop_df <- filter(hs_data, exp_period == ex_p & out_period == out_p) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Scatter plot of relationship
ggplot(data = loop_df, aes(x = .data[[ex]], y = .data[[out]])) + 
  geom_point() +
  geom_smooth(method='lm', se = T)
scatter_name <- paste0("scatter_", ex_p, "_", ex, "_", 
                       out_p, "_", out, "_hs.jpeg")
ggsave(filename = here::here("figs/lm_plots", scatter_name),
       device = "jpeg")

#' full cohort
lm_form <- paste(out, "~", ex, c_vars)
lm <- lm(as.formula(lm_form), data = loop_df)

#' append results to table 
lm_results <- bind_rows(lm_results, tab_results(lm, ex_p, out, out_p, "full cohort"))

#' diagnostic plots for linear regression
graphics.off()
jpeg(file = here::here("figs/lm_plots", paste0("diagnostic_", ex_p, "_", ex, "_", 
                                               out_p, "_", out, "_hs.jpeg")))
par(mfrow=c(2,2))
plot(lm)
dev.off()
par(mfrow=c(1,1))

#' Save results
lm_name <- paste0("LM_", ex_p, "_", ex, "_", out_p, "_", out, "_hs.rdata")
save(lm, file = here::here("results/lm_results", lm_name))

#' Write out summary files
lm_results$cohort <- "Healthy Start"
write_csv(lm_results, file = here::here("results", "tables", "linear_model_results_hs.csv"))

#' -----------------------------------------------------------------------------
#' Linear regression models for each exposure and outcome pair: MADRES
#' -----------------------------------------------------------------------------

#' Linear regression models
lm_results <- data.frame()

#' Continuous exposures
scale_exp <- c("pm", "no2", "o3", "rhavg", "tavg_c", "nox_pred", 
               "dist_roads_m", "ndvi_500m", "dist_parks_m", "park_area_500m", 
               "tree_cover_500m","impervious_500m", "pct_hh_poverty", "pct_less_hs_edu")

#' ----------------------------------------
#' Model controls
#' 
#' Exposures: pm, o3, no2, rhavg, tavg_c, nox_pred, dist_roads_m,
#' ndvi_500m, dist_parks_m, park_area_500m, 
#' park_count_500m, tree_cover_500m, impervious_500m,
#' lila, pct_hh_poverty, pct_less_hs_edu
#' 
#' Outcomes: zbmi
#' 
#' Exposure periods: pregnancy, t1, t2, t3
#' 
#' Outcome periods: 0, 6, 12, 24
#' ----------------------------------------

ex <- "pm"
out <- "zbmi"
ex_p <- "pregnancy"
out_p <- 0

#' ----------------------------------------
#' Regression code
#' ----------------------------------------

loop_df <- filter(md_data, exp_period == ex_p & out_period == out_p) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Scatter plot of relationship
ggplot(data = loop_df, aes(x = .data[[ex]], y = .data[[out]])) + 
  geom_point() +
  geom_smooth(method='lm', se = T)
scatter_name <- paste0("scatter_", ex_p, "_", ex, "_", 
                       out_p, "_", out, "_md.jpeg")
ggsave(filename = here::here("figs/lm_plots", scatter_name),
       device = "jpeg")

#' full cohort
lm_form <- paste(out, "~", ex, c_vars)
lm <- lm(as.formula(lm_form), data = loop_df)

#' append results to table 
lm_results <- bind_rows(lm_results, tab_results(lm, ex_p, out, out_p, "full cohort"))

#' diagnostic plots for linear regression
graphics.off()
jpeg(file = here::here("figs/lm_plots", paste0("diagnostic_", ex_p, "_", ex, "_", 
                                               out_p, "_", out, "_md.jpeg")))
par(mfrow=c(2,2))
plot(lm)
dev.off()
par(mfrow=c(1,1))

#' Save results
lm_name <- paste0("LM_", ex_p, "_", ex, "_", out_p, "_", out, "_md.rdata")
save(lm, file = here::here("results/lm_results", lm_name))

#' Write out summary files
lm_results$cohort <- "MADRES"
write_csv(lm_results, file = here::here("results", "tables", "linear_model_results_md.csv"))