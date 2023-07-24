
beta_95ci <- function(beta, se) {
  ucl <- beta + 1.96*se
  lcl <- beta - 1.96*se
  beta_95ci <- paste0(format(round(beta, digits = 2), nsmall = 2), 
                      " (", format(round(lcl, digits = 2), nsmall = 2), 
                      ", ", format(round(ucl, digits = 2), nsmall = 2), ")")
  return(beta_95ci)
}

or_95ci <- function(beta, se) {
  or <- exp(beta)
  ucl <- exp(beta + 1.96*se)
  lcl <- exp(beta - 1.96*se)
  or_95ci <- paste0(format(round(or, digits = 2), nsmall = 2, scientific = F), 
                      " (", format(round(lcl, digits = 2), nsmall = 2, scientific = F), 
                      ", ", format(round(ucl, digits = 2), nsmall = 2, scientific = F), ")")
  return(or_95ci)
}

library(tidyverse)
library(writexl)

#' -----------------------------------------------------------------------------
#' Linear regression results: pregnancy-wide exposures
#' -----------------------------------------------------------------------------

hs <- read_csv(here::here("results", "tables", "linear_model_results_hs.csv"))
md <- read_csv(here::here("results", "tables", "linear_model_results_md.csv")) 
comb <- read_csv(here::here("results", "tables", "lmer_model_results_combined.csv")) %>%
  mutate(beta_95ci = paste0(format(round(beta, digits = 2), nsmall = 2), 
                            " (", format(round(beta_lcl, digits = 2), nsmall = 2), 
                            ", ", format(round(beta_ucl, digits = 2), nsmall = 2), ")")) %>%
  mutate(cohort = "Combined") %>%
  filter(group == "full cohort") %>%
  filter(exposure_period == "pregnancy") %>%
  select(time_point, exposure_period, exposure, cohort, group, beta_95ci)

#' Main model results
df <- bind_rows(hs, md) %>%
  filter(group == "full cohort") %>%
  filter(exposure_period == "pregnancy") %>%
  mutate(beta_95ci = beta_95ci(beta, beta.se)) %>%
  select(time_point, exposure_period, exposure, cohort, group, beta_95ci) %>%
  pivot_wider(names_from = time_point, values_from = beta_95ci) %>%
  arrange(match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c", 
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m", 
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")),
          match(cohort, c("Healthy Start", "MADRES"))) %>%
  select(exposure, exposure_period, cohort, group, "0":"24")
write_csv(df, here::here("manuscripts", "ec0548a_manuscript", 
                         "linear_regression_results_by_cohort_pregnancy.csv"))

df2 <- comb %>%
  pivot_wider(names_from = time_point, values_from = beta_95ci) %>%
  arrange(match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c", 
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m", 
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(df2, here::here("manuscripts", "ec0548a_manuscript", 
                          "linear_regression_results_combined_cohort_pregnancy.csv"))

#' -----------------------------------------------------------------------------
#' Linear regression results: trimester specific exposures
#' -----------------------------------------------------------------------------

hs <- read_csv(here::here("results", "tables", "linear_model_results_hs.csv"))
md <- read_csv(here::here("results", "tables", "linear_model_results_md.csv")) 
comb <- read_csv(here::here("results", "tables", "lmer_model_results_combined.csv")) %>%
  mutate(beta_95ci = paste0(format(round(beta, digits = 2), nsmall = 2), 
                            " (", format(round(beta_lcl, digits = 2), nsmall = 2), 
                            ", ", format(round(beta_ucl, digits = 2), nsmall = 2), ")")) %>%
  mutate(cohort = "Combined") %>%
  filter(group == "full cohort") %>%
  filter(exposure_period != "pregnancy") %>%
  select(time_point, exposure_period, exposure, cohort, group, beta_95ci)

#' Main model results
df <- bind_rows(hs, md) %>%
  filter(group == "full cohort") %>%
  filter(exposure_period != "pregnancy") %>%
  mutate(beta_95ci = beta_95ci(beta, beta.se)) %>%
  select(time_point, exposure_period, exposure, cohort, group, beta_95ci) %>%
  pivot_wider(names_from = time_point, values_from = beta_95ci) %>%
  arrange(match(cohort, c("Healthy Start", "MADRES")),
          match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c", 
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m", 
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")),
          match(exposure_period, c("t1", "t2", "t3"))) %>%
  select(exposure, exposure_period, cohort, group, "0":"24")
write_csv(df, here::here("manuscripts", "ec0548a_manuscript", 
                         "linear_regression_results_by_cohort_trimesters.csv"))

df2 <- comb %>%
  pivot_wider(names_from = time_point, values_from = beta_95ci) %>%
  arrange(match(exposure, c("pm", "o3", "no2", "rhavg", "tavg_c", 
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m", 
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))
write_csv(df2, here::here("manuscripts", "ec0548a_manuscript", 
                          "linear_regression_results_combined_cohort_trimesters.csv"))

#' -----------------------------------------------------------------------------
#' Quantile-based g-computation results
#' -----------------------------------------------------------------------------

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
  pivot_wider(names_from = time_point, values_from = psi_95ci)
write_csv(qgcomp, here::here("manuscripts", "ec0548a_manuscript", 
                             "qgcomp_results_by_cohort_pregnancy_trimesters.csv"))

#' -----------------------------------------------------------------------------
#' BKMR results
#' -----------------------------------------------------------------------------

#PIPS table: pregnancy wide
bkmr_preg_0_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                               "BKMR_PIPs_pregnancy_0_zbmi_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_preg_6_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_6_zbmi_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_preg_12_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_12_zbmi_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_preg_24_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_24_zbmi_hs.csv")) %>%
  mutate(cohort = "hs")

bkmr_preg_0_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_0_zbmi_md.csv")) %>%
  mutate(cohort = "md")
bkmr_preg_6_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_6_zbmi_md.csv")) %>%
  mutate(cohort = "md")
bkmr_preg_12_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                       "BKMR_PIPs_pregnancy_12_zbmi_md.csv")) %>%
  mutate(cohort = "md")

bkmr_preg_0_st_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_0_zbmi_st_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_preg_6_st_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_6_zbmi_st_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_preg_12_st_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                       "BKMR_PIPs_pregnancy_12_zbmi_st_hs.csv")) %>%
  mutate(cohort = "hs")
bkmr_preg_24_st_hs <- read_csv(here::here("results/bkmr_results", "pips", 
                                       "BKMR_PIPs_pregnancy_24_zbmi_st_hs.csv")) %>%
  mutate(cohort = "hs")

bkmr_preg_0_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_0_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")
bkmr_preg_6_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_pregnancy_6_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")
bkmr_preg_12_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                       "BKMR_PIPs_pregnancy_12_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")

bkmr_pips_preg <- bind_rows(bkmr_preg_0_hs, bkmr_preg_6_hs, bkmr_preg_12_hs, bkmr_preg_24_hs,
                               bkmr_preg_0_md, bkmr_preg_6_md, bkmr_preg_12_md,
                               bkmr_preg_0_st_hs, bkmr_preg_6_st_hs, bkmr_preg_12_st_hs, bkmr_preg_24_st_hs,
                               bkmr_preg_0_st_md, bkmr_preg_6_st_md, bkmr_preg_12_st_md) %>%  
  mutate(PIP = format(round(PIP, digits = 2), nsmall = 2)) %>%
  pivot_wider(names_from = c(group, out_period), values_from = PIP) %>%
  select(variable:cohort, 'gis variables_0', 'st variables_0',
         'gis variables_6', 'st variables_6', 'gis variables_12', 'st variables_12',
         'gis variables_24', 'st variables_24') %>% 
  arrange(match(cohort, c("hs", "md")),
          match(variable, c("pm", "o3", "no2", "rhavg", "tavg_c", 
                            "bc_pred", "nox_pred", "dist_roads_m",
                            "ndvi_500m", "dist_parks_m", "park_area_500m", 
                            "park_count_500m", "tree_cover_500m", "impervious_500m",
                            "lila", "pct_hh_poverty", "pct_less_hs_edu")))

write_csv(bkmr_pips_preg, here::here("manuscripts", "ec0548a_manuscript", 
                                     "bkmr_pips_by_cohort_pregnancy.csv"))

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

bkmr_tri_0_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_trimesters_0_zbmi_md.csv")) %>%
  mutate(cohort = "md")
bkmr_tri_6_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                      "BKMR_PIPs_trimesters_6_zbmi_md.csv")) %>%
  mutate(cohort = "md")
bkmr_tri_12_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                       "BKMR_PIPs_trimesters_12_zbmi_md.csv")) %>%
  mutate(cohort = "md")

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

bkmr_tri_0_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                        "BKMR_PIPs_trimesters_0_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")
bkmr_tri_6_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                        "BKMR_PIPs_trimesters_6_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")
bkmr_tri_12_st_md <- read_csv(here::here("results/bkmr_results", "pips", 
                                         "BKMR_PIPs_trimesters_12_zbmi_st_md.csv")) %>%
  mutate(cohort = "md")

bkmr_pips_tri <- bind_rows(bkmr_tri_0_hs, bkmr_tri_6_hs, bkmr_tri_12_hs, bkmr_tri_24_hs,
                               bkmr_tri_0_md, bkmr_tri_6_md, bkmr_tri_12_md,
                               bkmr_tri_0_st_hs, bkmr_tri_6_st_hs, bkmr_tri_12_st_hs, bkmr_tri_24_st_hs,
                               bkmr_tri_0_st_md, bkmr_tri_6_st_md, bkmr_tri_12_st_md)

bkmr_group_pips_tri <- select(bkmr_pips_tri, -condPIP) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  select(-variable) %>%
  mutate(groupPIP = format(round(groupPIP, digits = 2), nsmall = 2)) %>%
  distinct() %>%
  pivot_wider(names_from = c(group, out_period), values_from = groupPIP) %>%
  select(exp_period:trimester, 'gis variables_0', 'st variables_0',
         'gis variables_6', 'st variables_6', 'gis variables_12', 'st variables_12',
         'gis variables_24', 'st variables_24') %>%
  arrange(match(cohort, c("hs", "md")))
write_csv(bkmr_group_pips_tri, here::here("manuscripts", "ec0548a_manuscript", 
                                          "bkmr_group_pips_by_cohort_trimesters.csv"))

vs1 <- c("pm", "o3", "no2", "rhavg", "tavg_c", 
         "bc_pred", "nox_pred", "dist_roads_m",
         "ndvi_500m", "dist_parks_m", "park_area_500m", 
         "park_count_500m", "tree_cover_500m", "impervious_500m",
         "lila", "pct_hh_poverty", "pct_less_hs_edu")
vs2 <- c("t1", "t2", "t3")
var_sort <- apply(expand.grid(vs1, vs2), 1, paste, collapse="_")

bkmr_cond_pips_tri <- select(bkmr_pips_tri, -groupPIP) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(condPIP = format(round(condPIP, digits = 2), nsmall = 2)) %>%
  pivot_wider(names_from = c(group, out_period), values_from = condPIP) %>%
  select(variable:trimester, 'gis variables_0', 'st variables_0',
         'gis variables_6', 'st variables_6', 'gis variables_12', 'st variables_12',
         'gis variables_24', 'st variables_24') %>%
  arrange(match(trimester, c("t1", "t2", "t3")),
          match(cohort, c("hs", "md")),
          match(variable, var_sort))
write_csv(bkmr_cond_pips_tri, here::here("manuscripts", "ec0548a_manuscript", 
                                          "bkmr_cond_pips_by_cohort_trimesters.csv"))
