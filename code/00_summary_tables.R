#' -----------------------------------------------------------------------------
#' Date created: August 1, 2023
#' Author: Sheena Martenies
#' Contact: smarte4@illinois.edu
#' 
#' Description: Summary tables for the manuscript and supplemental materials
#' -----------------------------------------------------------------------------

#' Load required libraries
library(ggplot2)
library(viridis)
library(ggthemes)
library(tidyverse)
library(lubridate)
library(writexl)

#' Some coordinate reference systems
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm_13N <- "+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
ll_nad83 <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ll_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#' [Box Health] directory where the data MUST live
#' These are PHI and must be protected as such!
#' Do not save any object containing these data in any other director!!
box_path <- "/Users/sheenamartenies/Library/CloudStorage/Box-Box/[Box Health - External] Echo_oif4d"

#' Summary functions
mean_sd <- function(x) {
  m <- mean(x, na.rm = T)
  s <- sd(x, na.rm = T)
  miss <- sum(is.na(x))
  m_sd <- paste0(format(round(m, digits = 1), nsmall = 1), " (", format(round(s, digits = 1), nsmall = 1), ")")
  return(m_sd)
}
n_pct <- function(x, test) {
  n <- sum(x == test, na.rm = T)
  pct <- (sum(x == test, na.rm = T) / sum(!is.na(x))) * 100
  miss <- sum(is.na(x))
  n_pct <- paste0(format(round(n, digits = 1), nsmall = 0), " (", format(round(pct, digits = 1), nsmall = 1), ")")
}
beta_95ci <- function(beta, se) {
  ucl <- beta + 1.96*se
  lcl <- beta - 1.96*se
  beta_95ci <- paste0(format(round(beta, digits = 2), nsmall = 2), 
                      " (", format(round(lcl, digits = 2), nsmall = 2), 
                      ", ", format(round(ucl, digits = 2), nsmall = 2), ")")
  return(beta_95ci)
}
beta_90ci <- function(beta, se) {
  ucl <- beta + 1.64*se
  lcl <- beta - 1.64*se
  beta_90ci <- paste0(format(round(beta, digits = 2), nsmall = 2), 
                      " (", format(round(lcl, digits = 2), nsmall = 2), 
                      ", ", format(round(ucl, digits = 2), nsmall = 2), ")")
  return(beta_90ci)
}
beta_se <- function(beta, se) {
  beta_se <- paste0(format(round(beta, digits = 2), nsmall = 2), 
                      " (", format(round(se, digits = 2), nsmall = 2), ")")
  return(beta_se)
}


#' -----------------------------------------------------------------------------
#' Read in HS and MADRES data
#' Note: these are sensitive data and can ONLY live in the [Box Health] folder
#' -----------------------------------------------------------------------------

hs_data <- read_csv(paste0(box_path, "/clean_data/hs_analytic_data_clean.csv")) %>%
  mutate(cohort = "hs")
md_data <- read_csv(paste0(box_path, "/clean_data/md_analytic_data_clean.csv")) %>%
  mutate(cohort = "md")

#' -----------------------------------------------------------------------------
#' Table 1: Demographic and outcomes summary for each study and outcome time period
#' -----------------------------------------------------------------------------

hs_summary <- hs_data %>%
  filter(exp_period == "pregnancy") %>%
  na.omit() %>%
  group_by(out_period) %>%
  summarize(n = n(),
            maternal_age = mean_sd(maternal_age),
            maternal_partner = n_pct(maternal_partner, 1),
            hispanic_re = n_pct(hispanic_re, 1),
            nhblack_re = n_pct(nhblack_re, 1),
            nhwhite_re = n_pct(nhwhite_re, 1),
            other_re = n_pct(other_re, 1),
            ed_no_hs = n_pct(ed_no_hs, 1),
            ed_hs = n_pct(ed_hs, 1),
            ed_aa = n_pct(ed_aa, 1),
            ed_4yr = n_pct(ed_4yr, 1),
            ed_grad = n_pct(ed_grad, 1),
            ppbmi_underweight = n_pct(ppbmi_underweight, 1),
            ppbmi_normal = n_pct(ppbmi_normal, 1),
            ppbmi_overweight = n_pct(ppbmi_overweight, 1),
            ppbmi_obese_i = n_pct(ppbmi_obese_i, 1),
            ppbmi_obese_ii = n_pct(ppbmi_obese_ii, 1),
            ppbmi_obese_iii = n_pct(ppbmi_obese_iii, 1),
            gestsmoking = n_pct(gestsmoking, 1),
            smokesh = n_pct(smokesh, 1),
            male_infant = n_pct(male, 1),
            winter_conception = n_pct(concep_winter, 1),
            spring_conception = n_pct(concep_spring, 1),
            summer_conception = n_pct(concep_summer, 1),
            fall_conception = n_pct(concep_fall, 1),
            bmi_z = mean_sd(zbmi))
hs_summary_t <- t(hs_summary[,-1])
colnames(hs_summary_t) <- t(hs_summary[,1])
hs_summary_t <- as.data.frame(hs_summary_t)
hs_summary_t$var <- rownames(hs_summary_t)
hs_summary_t <- hs_summary_t[,c(6,1:5)]
colnames(hs_summary_t)[-1] <- paste("hs", colnames(hs_summary_t)[-1], sep = "_") 

md_summary <- md_data %>%
  filter(exp_period == "pregnancy") %>%
  na.omit() %>%
  group_by(out_period) %>%
  summarize(n = n(),
            maternal_age = mean_sd(maternal_age),
            maternal_partner = n_pct(maternal_partner, 1),
            hispanic_re = n_pct(hispanic_re, 1),
            nhblack_re = n_pct(nhblack_re, 1),
            nhwhite_re = n_pct(nhwhite_re, 1),
            other_re = n_pct(other_re, 1),
            ed_no_hs = n_pct(ed_no_hs, 1),
            ed_hs = n_pct(ed_hs, 1),
            ed_aa = n_pct(ed_aa, 1),
            ed_4yr = n_pct(ed_4yr, 1),
            ed_grad = n_pct(ed_grad, 1),
            ppbmi_underweight = n_pct(ppbmi_underweight, 1),
            ppbmi_normal = n_pct(ppbmi_normal, 1),
            ppbmi_overweight = n_pct(ppbmi_overweight, 1),
            ppbmi_obese_i = n_pct(ppbmi_obese_i, 1),
            ppbmi_obese_ii = n_pct(ppbmi_obese_ii, 1),
            ppbmi_obese_iii = n_pct(ppbmi_obese_iii, 1),
            gestsmoking = n_pct(gestsmoking, 1),
            smokesh = n_pct(smokesh, 1),
            male_infant = n_pct(male, 1),
            winter_conception = n_pct(concep_winter, 1),
            spring_conception = n_pct(concep_spring, 1),
            summer_conception = n_pct(concep_summer, 1),
            fall_conception = n_pct(concep_fall, 1),
            bmi_z = mean_sd(zbmi))
md_summary_t <- t(md_summary[,-1])
colnames(md_summary_t) <- t(md_summary[,1])
md_summary_t <- as.data.frame(md_summary_t)
md_summary_t$var <- rownames(md_summary_t)
md_summary_t <- md_summary_t[,c(5,1:4)]
colnames(md_summary_t)[-1] <- paste("md", colnames(md_summary_t)[-1], sep = "_") 

sum_table <- bind_cols(hs_summary_t, md_summary_t[,c(-1)]) %>%
  select(var, hs_0, md_0, hs_6, md_6, hs_12, md_12, hs_24, md_24)
write_csv(sum_table, here::here("manuscripts", "ec0548a_manuscript", 
                                "table_1_demographics.csv"))

#' -----------------------------------------------------------------------------
#' Table 2: Exposure summary for each study and exposure time period
#' Using data for participants with data at delivery
#' 
#' Table S1: Exposure summary by trimester
#' -----------------------------------------------------------------------------

#' -------------------------------------
#' Table 2
#' -------------------------------------

hs_exp_summary <- hs_data %>%
  filter(out_period == 0) %>%
  na.omit() %>%
  group_by(exp_period) %>%
  summarize(n = n(),
            pm2.5 = mean_sd(pm),
            no2 = mean_sd(no2),
            o3 = mean_sd(o3),
            rh = mean_sd(rhavg),
            temp_C = mean_sd(tavg_c),
            bc = mean_sd(bc_pred),
            nox = NA,
            dist_roads_m = mean_sd(dist_roads_m),
            ndvi_500m = mean_sd(ndvi_500m),
            dist_parks_m = mean_sd(dist_parks_m),
            park_area_500m = mean_sd(park_area_500m),
            park_count_500m = mean_sd(park_count_500m),
            tree_cover_500m = mean_sd(tree_cover_500m),
            impervious_500m = mean_sd(impervious_500m),
            food_desert = n_pct(lila, 1),
            pct_hh_poverty = mean_sd(pct_hh_poverty),
            pct_less_hs_edu = mean_sd(pct_less_hs_edu))
hs_exp_summary_t <- t(hs_exp_summary[,-1])
colnames(hs_exp_summary_t) <- t(hs_exp_summary[,1])
hs_exp_summary_t <- as.data.frame(hs_exp_summary_t)
hs_exp_summary_t$var <- rownames(hs_exp_summary_t)
hs_exp_summary_t <- hs_exp_summary_t[,c(5,1:4)]
colnames(hs_exp_summary_t)[-1] <- paste("hs", colnames(hs_exp_summary_t)[-1], sep = "_") 

md_exp_summary <- md_data %>%
  filter(out_period == 0) %>%
  na.omit() %>%
  group_by(exp_period) %>%
  summarize(n = n(),
            pm2.5 = mean_sd(pm),
            no2 = mean_sd(no2),
            o3 = mean_sd(o3),
            rh = mean_sd(rhavg),
            temp_C = mean_sd(tavg_c),
            bc = NA,
            nox = mean_sd(nox_pred),
            dist_roads_m = mean_sd(dist_roads_m),
            ndvi_500m = mean_sd(ndvi_500m),
            dist_parks_m = mean_sd(dist_parks_m),
            park_area_500m = mean_sd(park_area_500m),
            park_count_500m = mean_sd(park_count_500m),
            tree_cover_500m = mean_sd(tree_cover_500m),
            impervious_500m = mean_sd(impervious_500m),
            food_desert = n_pct(lila, 1),
            pct_hh_poverty = mean_sd(pct_hh_poverty),
            pct_less_hs_edu = mean_sd(pct_less_hs_edu))
md_exp_summary_t <- t(md_exp_summary[,-1])
colnames(md_exp_summary_t) <- t(md_exp_summary[,1])
md_exp_summary_t <- as.data.frame(md_exp_summary_t)
md_exp_summary_t$var <- rownames(md_exp_summary_t)
md_exp_summary_t <- md_exp_summary_t[,c(5,1:4)]
colnames(md_exp_summary_t)[-1] <- paste("md", colnames(md_exp_summary_t)[-1], sep = "_") 

exp_sum <- bind_cols(hs_exp_summary_t, md_exp_summary_t[,-1]) %>%
  select(var, hs_pregnancy, md_pregnancy)
write_csv(exp_sum, here::here("manuscripts", "ec0548a_manuscript", 
                              "table_2_exposures_pregnancy.csv"))

#' -------------------------------------
#' Table S1
#' -------------------------------------

hs_exp_summary <- hs_data %>%
  filter(out_period == 0) %>%
  filter(exp_period != "pregnancy") %>%
  na.omit() %>%
  group_by(exp_period) %>%
  summarize(n = n(),
            pm2.5 = mean_sd(pm),
            no2 = mean_sd(no2),
            o3 = mean_sd(o3),
            rh = mean_sd(rhavg),
            temp_C = mean_sd(tavg_c),
            bc = mean_sd(bc_pred),
            dist_roads_m = mean_sd(dist_roads_m),
            ndvi_500m = mean_sd(ndvi_500m),
            dist_parks_m = mean_sd(dist_parks_m),
            park_area_500m = mean_sd(park_area_500m),
            park_count_500m = mean_sd(park_count_500m),
            tree_cover_500m = mean_sd(tree_cover_500m),
            impervious_500m = mean_sd(impervious_500m),
            food_desert = n_pct(lila, 1),
            pct_hh_poverty = mean_sd(pct_hh_poverty),
            pct_less_hs_edu = mean_sd(pct_less_hs_edu))
hs_exp_summary_t <- t(hs_exp_summary[,-1])
colnames(hs_exp_summary_t) <- t(hs_exp_summary[,1])
hs_exp_summary_t <- as.data.frame(hs_exp_summary_t)
hs_exp_summary_t$var <- rownames(hs_exp_summary_t)
hs_exp_summary_t <- hs_exp_summary_t[,c(4,1:3)]
hs_exp_summary_t$cohort <- "hs"

md_exp_summary <- md_data %>%
  filter(out_period == 0) %>%
  filter(exp_period != "pregnancy") %>%
  na.omit() %>%
  group_by(exp_period) %>%
  summarize(n = n(),
            pm2.5 = mean_sd(pm),
            no2 = mean_sd(no2),
            o3 = mean_sd(o3),
            rh = mean_sd(rhavg),
            temp_C = mean_sd(tavg_c),
            nox = mean_sd(nox_pred),
            dist_roads_m = mean_sd(dist_roads_m),
            ndvi_500m = mean_sd(ndvi_500m),
            dist_parks_m = mean_sd(dist_parks_m),
            park_area_500m = mean_sd(park_area_500m),
            park_count_500m = mean_sd(park_count_500m),
            tree_cover_500m = mean_sd(tree_cover_500m),
            impervious_500m = mean_sd(impervious_500m),
            food_desert = n_pct(lila, 1),
            pct_hh_poverty = mean_sd(pct_hh_poverty),
            pct_less_hs_edu = mean_sd(pct_less_hs_edu))
md_exp_summary_t <- t(md_exp_summary[,-1])
colnames(md_exp_summary_t) <- t(md_exp_summary[,1])
md_exp_summary_t <- as.data.frame(md_exp_summary_t)
md_exp_summary_t$var <- rownames(md_exp_summary_t)
md_exp_summary_t <- md_exp_summary_t[,c(4,1:3)]
md_exp_summary_t$cohort <- "md"

exp_sum <- bind_rows(hs_exp_summary_t, md_exp_summary_t) 
write_csv(exp_sum, here::here("manuscripts", "ec0548a_manuscript", 
                              "table_s1_exposures_trimesters.csv"))

#' -----------------------------------------------------------------------------
#' Table 3: Linear regression results: pregnancy-wide exposures
#' Table S2 and S3: Linear regression results by trimester
#' -----------------------------------------------------------------------------

hs_lm <- read_csv(here::here("results", "tables", "linear_model_results_hs.csv"))
md_lm <- read_csv(here::here("results", "tables", "linear_model_results_md.csv")) 

#' -------------------------------------
#' Table 3
#' -------------------------------------

df_lm_p <- bind_rows(hs_lm, md_lm) %>%
  filter(exposure_period == "pregnancy") %>%
  mutate(beta_95ci = beta_95ci(beta, beta.se),
         p_value = format(round(p.value, digits = 2), nsmall = 2)) %>%
  select(exposure, time_point, cohort, beta_95ci, p_value) %>%
  pivot_wider(names_from = time_point, values_from = c(beta_95ci, p_value)) %>%
  select(exposure, cohort, 
         beta_95ci_0,
         beta_95ci_6, 
         beta_95ci_12, 
         beta_95ci_24) %>%
  arrange(match(cohort, c("Healthy Start", "MADRES")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c", 
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m", 
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(df_lm_p, here::here("manuscripts", "ec0548a_manuscript", 
                              "table_3_linear_reg_summary.csv"))

#' -------------------------------------
#' Table S2
#' -------------------------------------

hs_lm_t <- bind_rows(hs_lm) %>%
  filter(exposure_period != "pregnancy") %>%
  mutate(beta_95ci = beta_95ci(beta, beta.se),
         p_value = format(round(p.value, digits = 2), nsmall = 2)) %>%
  select(exposure, time_point, exposure_period, beta_95ci, p_value) %>%
  pivot_wider(names_from = time_point, values_from = c(beta_95ci, p_value)) %>%
  select(exposure, exposure_period, 
         beta_95ci_0, 
         beta_95ci_6, 
         beta_95ci_12, 
         beta_95ci_24) %>%
  arrange(match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c", 
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m", 
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")),
          match(exposure_period, c("t1", "t2", "t3")))
write_csv(hs_lm_t, here::here("manuscripts", "ec0548a_manuscript", 
                              "table_s2_linear_reg_summary_hs_trimesters.csv"))

#' -------------------------------------
#' Table S3
#' -------------------------------------

md_lm_t <- bind_rows(md_lm) %>%
  filter(exposure_period != "pregnancy") %>%
  mutate(beta_95ci = beta_95ci(beta, beta.se),
         p_value = format(round(p.value, digits = 2), nsmall = 2)) %>%
  select(exposure, time_point, exposure_period, beta_95ci, p_value) %>%
  pivot_wider(names_from = time_point, values_from = c(beta_95ci, p_value)) %>%
  select(exposure, exposure_period, 
         beta_95ci_0, 
         beta_95ci_6, 
         beta_95ci_12, 
         beta_95ci_24) %>%
  arrange(match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c", 
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m", 
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")),
          match(exposure_period, c("t1", "t2", "t3")))
write_csv(md_lm_t, here::here("manuscripts", "ec0548a_manuscript", 
                              "table_s3_linear_reg_summary_md_trimesters.csv"))

#' -----------------------------------------------------------------------------
#' Table 4: g-computation results (marginal effect estimates)
#' Table S4-S7: g-computation coefficients and weights for Healthy Start
#' Table S8-S11: g-computation coefficients and weights for Healthy Start MADRES
#' -----------------------------------------------------------------------------

#' -------------------------------------
#' Table 4
#' -------------------------------------

#' Summary of overall effect estimates
preg_hs <- read_csv(here::here("results/qgcomp_results", "psi_est", 
                               "qgcomp_pregnancy_psi_estimates_hs.csv")) %>%
  mutate(cohort = "hs")
t1_hs <- read_csv(here::here("results/qgcomp_results", "psi_est", 
                               "qgcomp_t1_psi_estimates_hs.csv")) %>%
  mutate(cohort = "hs")
t2_hs <- read_csv(here::here("results/qgcomp_results", "psi_est", 
                             "qgcomp_t2_psi_estimates_hs.csv")) %>%
  mutate(cohort = "hs")
t3_hs <- read_csv(here::here("results/qgcomp_results", "psi_est", 
                             "qgcomp_t3_psi_estimates_hs.csv")) %>%
  mutate(cohort = "hs")

preg_md <- read_csv(here::here("results/qgcomp_results", "psi_est", 
                               "qgcomp_pregnancy_psi_estimates_md.csv")) %>%
  mutate(cohort = "md")
t1_md <- read_csv(here::here("results/qgcomp_results", "psi_est", 
                             "qgcomp_t1_psi_estimates_md.csv")) %>%
  mutate(cohort = "md")
t2_md <- read_csv(here::here("results/qgcomp_results", "psi_est", 
                             "qgcomp_t2_psi_estimates_md.csv")) %>%
  mutate(cohort = "md")
t3_md <- read_csv(here::here("results/qgcomp_results", "psi_est", 
                             "qgcomp_t3_psi_estimates_md.csv")) %>%
  mutate(cohort = "md")

qgcomp <- bind_rows(preg_hs, t1_hs, t2_hs, t3_hs,
                    preg_md, t1_md, t2_md, t3_md) %>%
  filter(model == "qgcomp.boot") %>%
  mutate(psi_95ci = paste0(format(round(psi, digits = 2), nsmall = 2), 
                           " (", format(round(psi_ll, digits = 2), nsmall = 2), 
                           ", ", format(round(psi_ul, digits = 2), nsmall = 2), ")")) %>%
  select(time_point, exposure_period, group, cohort, psi_95ci) %>%
  pivot_wider(names_from = time_point, values_from = psi_95ci) %>%
  arrange(match(cohort, c("hs", "md")),
          match(group, c("gis variables", "st variables")),
          match(exposure_period, c("pregnancy", "t1", "t2", "t3")))
write_csv(qgcomp, here::here("manuscripts", "ec0548a_manuscript", 
                             "table_4_gcomp_psi_summary.csv"))

#' -------------------------------------
#' Table S4: HS exposure weights: pregnancy
#' -------------------------------------

#' Summary of exposure coefficients and weights
qgcomp_hs_preg_wts <- read_csv(here::here("results/qgcomp_results", "weights", 
                                      "qgcomp_pregnancy_weights_hs.csv")) %>%
  mutate(cohort = "hs") %>%
  mutate(coef = format(round(coefficients, digits = 2), nsmall = 2),
         wt = format(round(weights, digits = 2), nsmall = 2),
         coef_wt = paste0(coef, " (", wt, ")")) %>%
  select(exposure, group, time_point, coef_wt) %>%
  pivot_wider(names_from = time_point, values_from = coef_wt) %>%
  arrange(match(group, c("gis variables", "st variables")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(qgcomp_hs_preg_wts, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s4_gcomp_weights_pregnancy_hs.csv"))

#' -------------------------------------
#' Table S5: HS exposure weights: T1
#' -------------------------------------

#' Summary of exposure coefficients and weights
qgcomp_hs_t1_wts <- read_csv(here::here("results/qgcomp_results", "weights", 
                                          "qgcomp_t1_weights_hs.csv")) %>%
  mutate(cohort = "hs") %>%
  mutate(coef = format(round(coefficients, digits = 2), nsmall = 2),
         wt = format(round(weights, digits = 2), nsmall = 2),
         coef_wt = paste0(coef, " (", wt, ")")) %>%
  select(exposure, group, time_point, coef_wt) %>%
  pivot_wider(names_from = time_point, values_from = coef_wt) %>%
  arrange(match(group, c("gis variables", "st variables")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(qgcomp_hs_t1_wts, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s5_gcomp_weights_t1_hs.csv"))

#' -------------------------------------
#' Table S6: HS exposure weights: T2
#' -------------------------------------

#' Summary of exposure coefficients and weights
qgcomp_hs_t2_wts <- read_csv(here::here("results/qgcomp_results", "weights", 
                                        "qgcomp_t2_weights_hs.csv")) %>%
  mutate(cohort = "hs") %>%
  mutate(coef = format(round(coefficients, digits = 2), nsmall = 2),
         wt = format(round(weights, digits = 2), nsmall = 2),
         coef_wt = paste0(coef, " (", wt, ")")) %>%
  select(exposure, group, time_point, coef_wt) %>%
  pivot_wider(names_from = time_point, values_from = coef_wt) %>%
  arrange(match(group, c("gis variables", "st variables")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(qgcomp_hs_t2_wts, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s6_gcomp_weights_t2_hs.csv"))

#' -------------------------------------
#' Table S7: HS exposure weights: T3
#' -------------------------------------

#' Summary of exposure coefficients and weights
qgcomp_hs_t3_wts <- read_csv(here::here("results/qgcomp_results", "weights", 
                                        "qgcomp_t3_weights_hs.csv")) %>%
  mutate(cohort = "hs") %>%
  mutate(coef = format(round(coefficients, digits = 2), nsmall = 2),
         wt = format(round(weights, digits = 2), nsmall = 2),
         coef_wt = paste0(coef, " (", wt, ")")) %>%
  select(exposure, group, time_point, coef_wt) %>%
  pivot_wider(names_from = time_point, values_from = coef_wt) %>%
  arrange(match(group, c("gis variables", "st variables")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(qgcomp_hs_t3_wts, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s7_gcomp_weights_t3_hs.csv"))

#' -------------------------------------
#' Table S8: HS exposure weights: pregnancy
#' -------------------------------------

#' Summary of exposure coefficients and weights
qgcomp_md_preg_wts <- read_csv(here::here("results/qgcomp_results", "weights", 
                                          "qgcomp_pregnancy_weights_md.csv")) %>%
  mutate(cohort = "md") %>%
  mutate(coef = format(round(coefficients, digits = 2), nsmall = 2),
         wt = format(round(weights, digits = 2), nsmall = 2),
         coef_wt = paste0(coef, " (", wt, ")")) %>%
  select(exposure, group, time_point, coef_wt) %>%
  pivot_wider(names_from = time_point, values_from = coef_wt) %>%
  arrange(match(group, c("gis variables", "st variables")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(qgcomp_md_preg_wts, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s8_gcomp_weights_pregnancy_md.csv"))

#' -------------------------------------
#' Table S9: md exposure weights: T1
#' -------------------------------------

#' Summary of exposure coefficients and weights
qgcomp_md_t1_wts <- read_csv(here::here("results/qgcomp_results", "weights", 
                                        "qgcomp_t1_weights_md.csv")) %>%
  mutate(cohort = "md") %>%
  mutate(coef = format(round(coefficients, digits = 2), nsmall = 2),
         wt = format(round(weights, digits = 2), nsmall = 2),
         coef_wt = paste0(coef, " (", wt, ")")) %>%
  select(exposure, group, time_point, coef_wt) %>%
  pivot_wider(names_from = time_point, values_from = coef_wt) %>%
  arrange(match(group, c("gis variables", "st variables")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(qgcomp_md_t1_wts, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s9_gcomp_weights_t1_md.csv"))

#' -------------------------------------
#' Table S10: md exposure weights: T2
#' -------------------------------------

#' Summary of exposure coefficients and weights
qgcomp_md_t2_wts <- read_csv(here::here("results/qgcomp_results", "weights", 
                                        "qgcomp_t2_weights_md.csv")) %>%
  mutate(cohort = "md") %>%
  mutate(coef = format(round(coefficients, digits = 2), nsmall = 2),
         wt = format(round(weights, digits = 2), nsmall = 2),
         coef_wt = paste0(coef, " (", wt, ")")) %>%
  select(exposure, group, time_point, coef_wt) %>%
  pivot_wider(names_from = time_point, values_from = coef_wt) %>%
  arrange(match(group, c("gis variables", "st variables")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(qgcomp_md_t2_wts, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s10_gcomp_weights_t2_md.csv"))

#' -------------------------------------
#' Table S11: md exposure weights: T3
#' -------------------------------------

#' Summary of exposure coefficients and weights
qgcomp_md_t3_wts <- read_csv(here::here("results/qgcomp_results", "weights", 
                                        "qgcomp_t3_weights_md.csv")) %>%
  mutate(cohort = "md") %>%
  mutate(coef = format(round(coefficients, digits = 2), nsmall = 2),
         wt = format(round(weights, digits = 2), nsmall = 2),
         coef_wt = paste0(coef, " (", wt, ")")) %>%
  select(exposure, group, time_point, coef_wt) %>%
  pivot_wider(names_from = time_point, values_from = coef_wt) %>%
  arrange(match(group, c("gis variables", "st variables")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(qgcomp_md_t3_wts, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s11_gcomp_weights_t3_md.csv"))

#' -----------------------------------------------------------------------------
#' Table S12-S14: BKMR posterior means and SD, PIPs for Healthy Start
#' Table S15-S17: BKMR posterior means and SD, PIPs for MADRES
#' -----------------------------------------------------------------------------

#' -------------------------------------
#' Table S12: Healthy Start no var selection
#' -------------------------------------

time_points <- c(0, 6, 12, 24)
hs_bkmr_results <- tibble() 

for(i in 1:length(time_points)) {
  # GIS
  csv_name1a <- paste0("BKMR_singvar_pregnancy_", time_points[i], "_zbmi_hs.csv")
  df_risks1 <- read_csv(here::here("results", "bkmr_results", "risks", csv_name1a)) %>%
    filter(q.fixed == 0.5)
  
  csv_name1b <- paste0("BKMR_PIPs_pregnancy_", time_points[i], "_zbmi_hs.csv")
  df_pips1 <- read_csv(here::here("results", "bkmr_results", "pips", csv_name1b))
  
  df1 <- left_join(df_risks1, select(df_pips1, variable, PIP), by = "variable") %>%
    mutate(pm_psd_pip = paste0(beta_se(est, sd), "; ",
                               format(round(PIP, digits = 2), nsmall = 2))) %>%
    mutate(group = "gis variables",
           cohort = "hs",
           time_point = time_points[i])
  
  #ST
  csv_name2a <- paste0("BKMR_singvar_pregnancy_", time_points[i], "_zbmi_st_hs.csv")
  df_risks2 <- read_csv(here::here("results", "bkmr_results", "risks", csv_name2a)) %>%
    filter(q.fixed == 0.5)
  
  csv_name2b <- paste0("BKMR_PIPs_pregnancy_", time_points[i], "_zbmi_st_hs.csv")
  df_pips2 <- read_csv(here::here("results", "bkmr_results", "pips", csv_name2b))
  
  df2 <- left_join(df_risks2, select(df_pips2, variable, PIP), by = "variable") %>%
    mutate(pm_psd_pip = paste0(beta_se(est, sd), "; ",
                               format(round(PIP, digits = 2), nsmall = 2))) %>%
    mutate(group = "st variables",
           cohort = "hs",
           time_point = time_points[i])
  
  hs_bkmr_results <- bind_rows(hs_bkmr_results, df1, df2)
}

bkmr_hs_tab <- hs_bkmr_results %>%
  mutate(group = ifelse(group == "gis variables", "GIS", "ST")) %>%
  select(variable, pm_psd_pip, group, time_point) %>%
  pivot_wider(names_from = time_point, values_from = pm_psd_pip) %>%
  arrange(match(group, c("GIS", "ST")),
          match(variable, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(bkmr_hs_tab, here::here("manuscripts", "ec0548a_manuscript", 
                                  "table_s12_bkmr_pips_pregnancy_hs.csv"))

#' -------------------------------------
#' Trimester variable selection
#' -------------------------------------

time_points <- c(0, 6, 12, 24)
hs_bkmr_results_vs <- tibble() 

for(i in 1:length(time_points)) {
  # GIS
  csv_name1a <- paste0("BKMR_singvar_trimesters_", time_points[i], "_zbmi_hs.csv")
  df_risks1 <- read_csv(here::here("results", "bkmr_results", "risks", csv_name1a)) %>%
    mutate(group = "gis variables", 
           time_point = time_points[i])
  
  #ST
  csv_name2a <- paste0("BKMR_singvar_trimesters_", time_points[i], "_zbmi_st_hs.csv")
  df_risks2 <- read_csv(here::here("results", "bkmr_results", "risks", csv_name2a)) %>%
    mutate(group = "st variables", 
           time_point = time_points[i])
  
  hs_bkmr_results_vs <- bind_rows(hs_bkmr_results_vs, df_risks1, df_risks2)
}

hs_bkmr_results_vs2 <- hs_bkmr_results_vs %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(variable = trimws(gsub(paste(c("_t1", "_t2", "_t3"), collapse = "|"), 
                                "\\1", variable))) %>%
  filter(q.fixed == 0.5) %>%
  mutate(pm_psd = beta_se(est, sd)) %>%
  select(variable, group, trimester, time_point, pm_psd)

#' -------------------------------------
#' Table S13: Healthy Start var selection- GIS
#' Table S14: Healthy Start var selection- ST
#' -------------------------------------

#PIPS table: trimesters 
bkmr_tri_0_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                     "BKMR_PIPs_trimesters_0_zbmi_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_tri_6_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                     "BKMR_PIPs_trimesters_6_zbmi_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_tri_12_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_trimesters_12_zbmi_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_tri_24_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_trimesters_24_zbmi_hs.csv")) %>%
  mutate(cohort = "hs")

bkmr_tri_0_st_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                        "BKMR_PIPs_trimesters_0_zbmi_st_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_tri_6_st_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                        "BKMR_PIPs_trimesters_6_zbmi_st_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_tri_12_st_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                         "BKMR_PIPs_trimesters_12_zbmi_st_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_tri_24_st_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                         "BKMR_PIPs_trimesters_24_zbmi_st_hs.csv")) %>%
  mutate(cohort = "hs")

bkmr_pips_tri_hs <- bind_rows(bkmr_tri_0_hs, bkmr_tri_6_hs, bkmr_tri_12_hs, bkmr_tri_24_hs,
                              bkmr_tri_0_st_hs, bkmr_tri_6_st_hs, bkmr_tri_12_st_hs, bkmr_tri_24_st_hs)

hs_bkmr_results_gp <- select(bkmr_pips_tri_hs, -condPIP) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  select(-variable) %>%
  mutate(groupPIP = format(round(groupPIP, digits = 2), nsmall = 2)) %>%
  distinct() %>%
  pivot_wider(names_from = out_period, values_from = groupPIP)
write_csv(hs_bkmr_results_gp, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s13_s14_bkmr_gpips_hs.csv"))

hs_bkmr_results_cp <- select(bkmr_pips_tri_hs, -groupPIP) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(variable = trimws(gsub(paste(c("_t1", "_t2", "_t3"), collapse = "|"), 
                                "\\1", variable))) %>%
  mutate(condPIP = format(round(condPIP, digits = 2), nsmall = 2)) %>%
  select(variable, group, time_point = out_period, trimester, condPIP)

hs_bkmr_var_sel_gis <- hs_bkmr_results_vs2 %>%
  left_join(hs_bkmr_results_cp, by = c("variable", "trimester", "group", "time_point")) %>%
  mutate(pm_psd_pip = paste0(pm_psd, "; ", condPIP)) %>%
  filter(group == "gis variables") %>%
  select(variable, trimester, time_point, pm_psd_pip) %>%
  pivot_wider(names_from = time_point, values_from = pm_psd_pip) %>%
  arrange(match(trimester, c("t1", "t2", "t3")),
          match(variable, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(hs_bkmr_var_sel_gis, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s13_bkmr_means_pips_gis_hs.csv"))

hs_bkmr_var_sel_st <- hs_bkmr_results_vs2 %>%
  left_join(hs_bkmr_results_cp, by = c("variable", "trimester", "group", "time_point")) %>%
  mutate(pm_psd_pip = paste0(pm_psd, "; ", condPIP)) %>%
  filter(group == "st variables") %>%
  select(variable, trimester, time_point, pm_psd_pip) %>%
  pivot_wider(names_from = time_point, values_from = pm_psd_pip) %>%
  arrange(match(trimester, c("t1", "t2", "t3")),
          match(variable, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(hs_bkmr_var_sel_st, here::here("manuscripts", "ec0548a_manuscript", 
                                          "table_s14_bkmr_means_pips_st_hs.csv"))

#' -------------------------------------
#' Table S15: MADRES no var selection
#' -------------------------------------

time_points <- c(0, 6, 12)
md_bkmr_results <- tibble() 

for(i in 1:length(time_points)) {
  # GIS
  csv_name1a <- paste0("BKMR_singvar_pregnancy_", time_points[i], "_zbmi_md.csv")
  df_risks1 <- read_csv(here::here("results", "bkmr_results", "risks", csv_name1a)) %>%
    filter(q.fixed == 0.5)
  
  csv_name1b <- paste0("BKMR_PIPs_pregnancy_", time_points[i], "_zbmi_md.csv")
  df_pips1 <- read_csv(here::here("results", "bkmr_results", "pips", csv_name1b))
  
  df1 <- left_join(df_risks1, select(df_pips1, variable, PIP), by = "variable") %>%
    mutate(pm_psd_pip = paste0(beta_se(est, sd), "; ",
                               format(round(PIP, digits = 2), nsmall = 2))) %>%
    mutate(group = "gis variables",
           cohort = "md",
           time_point = time_points[i])
  
  #ST
  csv_name2a <- paste0("BKMR_singvar_pregnancy_", time_points[i], "_zbmi_st_md.csv")
  df_risks2 <- read_csv(here::here("results", "bkmr_results", "risks", csv_name2a)) %>%
    filter(q.fixed == 0.5)
  
  csv_name2b <- paste0("BKMR_PIPs_pregnancy_", time_points[i], "_zbmi_st_md.csv")
  df_pips2 <- read_csv(here::here("results", "bkmr_results", "pips", csv_name2b))
  
  df2 <- left_join(df_risks2, select(df_pips2, variable, PIP), by = "variable") %>%
    mutate(pm_psd_pip = paste0(beta_se(est, sd), "; ",
                               format(round(PIP, digits = 2), nsmall = 2))) %>%
    mutate(group = "st variables",
           cohort = "md",
           time_point = time_points[i])
  
  md_bkmr_results <- bind_rows(md_bkmr_results, df1, df2)
}

bkmr_md_tab <- md_bkmr_results %>%
  mutate(group = ifelse(group == "gis variables", "GIS", "ST")) %>%
  select(variable, pm_psd_pip, group, time_point) %>%
  pivot_wider(names_from = time_point, values_from = pm_psd_pip) %>%
  arrange(match(group, c("GIS", "ST")),
          match(variable, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(bkmr_md_tab, here::here("manuscripts", "ec0548a_manuscript", 
                                  "table_s15_bkmr_pips_pregnancy_md.csv"))

#' -------------------------------------
#' Trimester variable selection
#' -------------------------------------

time_points <- c(0, 6, 12)
md_bkmr_results_vs <- tibble() 

for(i in 1:length(time_points)) {
  # GIS
  csv_name1a <- paste0("BKMR_singvar_trimesters_", time_points[i], "_zbmi_md.csv")
  df_risks1 <- read_csv(here::here("results", "bkmr_results", "risks", csv_name1a)) %>%
    mutate(group = "gis variables", 
           time_point = time_points[i])
  
  #ST
  csv_name2a <- paste0("BKMR_singvar_trimesters_", time_points[i], "_zbmi_st_md.csv")
  df_risks2 <- read_csv(here::here("results", "bkmr_results", "risks", csv_name2a)) %>%
    mutate(group = "st variables", 
           time_point = time_points[i])
  
  md_bkmr_results_vs <- bind_rows(md_bkmr_results_vs, df_risks1, df_risks2)
}

md_bkmr_results_vs2 <- md_bkmr_results_vs %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(variable = trimws(gsub(paste(c("_t1", "_t2", "_t3"), collapse = "|"), 
                                "\\1", variable))) %>%
  filter(q.fixed == 0.5) %>%
  mutate(pm_psd = beta_se(est, sd)) %>%
  select(variable, group, trimester, time_point, pm_psd)

#' -------------------------------------
#' Table S16: MADRES var selection- GIS
#' Table S17: MADRES Start var selection- ST
#' -------------------------------------

#PIPS table: trimesters 
bkmr_tri_0_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                     "BKMR_PIPs_trimesters_0_zbmi_md.csv")) %>%
  mutate(cohort = "md")
bkmr_tri_6_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                     "BKMR_PIPs_trimesters_6_zbmi_md.csv")) %>%
  mutate(cohort = "md")
bkmr_tri_12_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_trimesters_12_zbmi_md.csv")) %>%
  mutate(cohort = "md")

bkmr_tri_0_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                        "BKMR_PIPs_trimesters_0_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")
bkmr_tri_6_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                        "BKMR_PIPs_trimesters_6_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")
bkmr_tri_12_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                         "BKMR_PIPs_trimesters_12_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")

bkmr_pips_tri_md <- bind_rows(bkmr_tri_0_md, bkmr_tri_6_md, bkmr_tri_12_md, 
                              bkmr_tri_0_st_md, bkmr_tri_6_st_md, bkmr_tri_12_st_md)

md_bkmr_results_gp <- select(bkmr_pips_tri_md, -condPIP) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  select(-variable) %>%
  mutate(groupPIP = format(round(groupPIP, digits = 2), nsmall = 2)) %>%
  distinct() %>%
  pivot_wider(names_from = out_period, values_from = groupPIP)
write_csv(md_bkmr_results_gp, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s16_s17_bkmr_gpips_md.csv"))

md_bkmr_results_cp <- select(bkmr_pips_tri_md, -groupPIP) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(variable = trimws(gsub(paste(c("_t1", "_t2", "_t3"), collapse = "|"), 
                                "\\1", variable))) %>%
  mutate(condPIP = format(round(condPIP, digits = 2), nsmall = 2)) %>%
  select(variable, group, time_point = out_period, trimester, condPIP)

md_bkmr_var_sel_gis <- md_bkmr_results_vs2 %>%
  left_join(md_bkmr_results_cp, by = c("variable", "trimester", "group", "time_point")) %>%
  mutate(pm_psd_pip = paste0(pm_psd, "; ", condPIP)) %>%
  filter(group == "gis variables") %>%
  select(variable, trimester, time_point, pm_psd_pip) %>%
  pivot_wider(names_from = time_point, values_from = pm_psd_pip) %>%
  arrange(match(trimester, c("t1", "t2", "t3")),
          match(variable, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(md_bkmr_var_sel_gis, here::here("manuscripts", "ec0548a_manuscript", 
                                          "table_s16_bkmr_means_pips_gis_md.csv"))

md_bkmr_var_sel_st <- md_bkmr_results_vs2 %>%
  left_join(md_bkmr_results_cp, by = c("variable", "trimester", "group", "time_point")) %>%
  mutate(pm_psd_pip = paste0(pm_psd, "; ", condPIP)) %>%
  filter(group == "st variables") %>%
  select(variable, trimester, time_point, pm_psd_pip) %>%
  pivot_wider(names_from = time_point, values_from = pm_psd_pip) %>%
  arrange(match(trimester, c("t1", "t2", "t3")),
          match(variable, c("pm", "o3", "no2", "rhavg", "tavg_c",
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m",
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(md_bkmr_var_sel_st, here::here("manuscripts", "ec0548a_manuscript", 
                                         "table_s17_bkmr_means_pips_st_md.csv"))


