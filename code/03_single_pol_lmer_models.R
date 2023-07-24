#' -----------------------------------------------------------------------------
#' Date created: February 28, 2023
#' Date updated: June 15, 2023
#' Author: Sheena Martenies
#' Contact: smarte4@illinois.edu
#' 
#' Description: Linear mixed models between exposures and continuous BMI 
#' Includes a random intercept for cohort
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

tab_results <- function(model, exp_period, outcome, out_period, group) {
  
  ci <- confint(model, method = "boot")
  
  temp <- data.frame(exposure_period = exp_period,
                     time_point = out_period,
                     outcome = outcome,
                     exposure = rownames(summary(model)$coefficients)[2],
                     group = group,
                     beta = summary(model)$coefficients[2,1],
                     beta_lcl = ci[4,1],
                     beta_ucl = ci[4,2])
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

comb_data <- bind_rows(hs_data, md_data) %>%
  select(-bc_pred, -nox_pred) %>%
  na.omit()

c_vars <- "+ maternal_age + maternal_partner + ed_no_hs + ed_hs + ed_aa + ed_4yr + ppbmi_underweight + ppbmi_overweight + ppbmi_obese_i + ppbmi_obese_ii + ppbmi_obese_iii + smokesh + gestsmoking + concep_spring + concep_summer + concep_fall + male"

#' -----------------------------------------------------------------------------
#' Linear mixed effect models for each exposure and outcome pair
#' Random effect for cohort
#' -----------------------------------------------------------------------------

#' Linear regression models
lmer_results <- data.frame()

#' Continuous exposures
scale_exp <- c("pm", "no2", "o3", "rhavg", "tavg_c", 
               "dist_roads_m", "ndvi_500m", "dist_parks_m", "park_area_500m", 
               "tree_cover_500m","impervious_500m", "pct_hh_poverty", "pct_less_hs_edu")

#' ----------------------------------------
#' Model controls
#' 
#' Exposures: pm, o3, no2, rhavg, tavg_c, dist_roads_m,
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

loop_df <- filter(comb_data, out_period == out_p & exp_period == ex_p) %>%
  mutate_at(.vars = vars(all_of(scale_exp)), .funs = scale_this)

#' Scatter plot of relationship
ggplot(data = loop_df, aes(x = .data[[ex]], y = .data[[out]],
                           color = cohort)) + 
  geom_point() +
  geom_smooth(method='lm', se = T)
ggsave(filename = here::here("figs/lmer_plots", 
                             paste0("scatter_", ex_p, "_", ex, "_", out_p, "_",
                                    out, "_combined.jpeg")),
       device = "jpeg")

#' full cohort
lmer_form <- paste(out, "~", ex, c_vars, "+ (1 | cohort)") #random intercept
lmer <- lmer(as.formula(lmer_form), data = loop_df)

#' append results to table 
lmer_results <- bind_rows(lmer_results, tab_results(lmer, ex_p, out, out_p, "full cohort"))

#' diagnostic plots for regression
graphics.off()
jpeg(file = here::here("figs/lmer_plots/", paste0("diagnostic_", ex_p, "_", ex, 
                                                  "_", out_p, "_", out, "_combined.jpeg")))
plot(lmer)
dev.off()

#' Save results
lmer_name <- paste0("lmer_", ex_p, "_", ex, "_", out_p, "_", out, "_combined.rdata")
save(lmer, file = here::here("results/lmer_results", lmer_name))

#' Write out summary file
write_csv(lmer_results, file = here::here("results", "tables", "lmer_model_results_combined.csv"))
