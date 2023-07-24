#' -----------------------------------------------------------------------------
#' Date created: December 6, 2022
#' Author: Sheena Martenies
#' Contact: smarte4@illinois.edu
#' 
#' Description: summarize demographic, outcome, and exposure data for each 
#' cohort and period
#' ----------------------------------------------------------------------------

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

#' -----------------------------------------------------------------------------
#' Read in HS and MADRES data
#' Note: these are sensitive data and can ONLY live in the [Box Health] folder
#' -----------------------------------------------------------------------------

hs_data <- read_csv(paste0(box_path, "/clean_data/hs_analytic_data_clean.csv")) %>%
  mutate(cohort = "hs")
md_data <- read_csv(paste0(box_path, "/clean_data/md_analytic_data_clean.csv")) %>%
  mutate(cohort = "md")

comb_data <- bind_rows(hs_data, md_data) %>%
  mutate(cohort = as.factor(cohort))

#' -----------------------------------------------------------------------------
#' Sample sizes at each age time point
#' -----------------------------------------------------------------------------

hs_df <- filter(comb_data, cohort == "hs") %>%
  select(-nox_pred)

hs_birth <- filter(hs_df, exp_period == "pregnancy" & out_period == 0) %>%
  na.omit()
hs_age_4_6_mo <- filter(hs_df, exp_period == "pregnancy" & out_period == 6) %>%
  na.omit()
hs_age_10_12_mo <- filter(hs_df, exp_period == "pregnancy" & out_period == 12) %>%
  na.omit()
hs_age_22_24_mo <- filter(hs_df, exp_period == "pregnancy" & out_period == 24) %>%
  na.omit()
hs_age_4_6_yr <- filter(hs_df, exp_period == "pregnancy" & out_period == 72) %>%
  na.omit()

md_df <- filter(comb_data, cohort == "md") %>%
  select(-bc_pred)

md_birth <- filter(md_df, exp_period == "pregnancy" & out_period == 0) %>%
  na.omit()
md_age_4_6_mo <- filter(md_df, exp_period == "pregnancy" & out_period == 6) %>%
  na.omit()
md_age_10_12_mo <- filter(md_df, exp_period == "pregnancy" & out_period == 12) %>%
  na.omit()
md_age_22_24_mo <- filter(md_df, exp_period == "pregnancy" & out_period == 24) %>%
  na.omit()

#' -----------------------------------------------------------------------------
#' Demographic and outcomes summary for each study and outcome time period
#' -----------------------------------------------------------------------------

hs_summary <- comb_data %>%
  filter(cohort == "hs") %>%
  select(-nox_pred) %>%
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
write_csv(hs_summary_t, here::here("results", "tables", "hs_demographic_summary.csv"))

md_summary <- comb_data %>%
  filter(cohort == "md") %>%
  select(-bc_pred) %>%
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
write_csv(md_summary_t, here::here("results", "tables", "md_demographic_summary.csv"))

hs_ids <- unique(hs_birth$pid)
md_ids <- unique(md_birth$pid)

comb_summary <- comb_data %>%
  filter(pid %in% c(hs_ids, md_ids)) %>%
  select(-bc_pred, -nox_pred) %>%
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
comb_summary_t <- t(comb_summary[,-1])
colnames(comb_summary_t) <- t(comb_summary[,1])
comb_summary_t <- as.data.frame(comb_summary_t)
comb_summary_t$var <- rownames(comb_summary_t)
comb_summary_t <- comb_summary_t[,c(6,1:5)]
write_csv(comb_summary_t, here::here("results", "tables", "comb_demographic_summary.csv"))


#' -----------------------------------------------------------------------------
#' Exposure summary for each study and exposure time period
#' Using data for participants with data at delivery
#' -----------------------------------------------------------------------------

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
write_csv(hs_exp_summary_t, here::here("results", "tables", "hs_exposure_summary.csv"))

hs_bc <- filter(hs_data, out_period == 0) %>%
  filter(!is.na(exp_period)) %>%
  select(pid, exp_period, bc_pred) %>%
  pivot_wider(id_cols = pid, names_from = exp_period, values_from = bc_pred)
cor(hs_bc[,c(2:5)])

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
write_csv(md_exp_summary_t, here::here("results", "tables", "md_exposure_summary.csv"))

md_nox <- filter(md_data, out_period == 0) %>%
  filter(!is.na(exp_period)) %>%
  select(pid, exp_period, nox_pred) %>%
  pivot_wider(id_cols = pid, names_from = exp_period, values_from = nox_pred) %>%
  na.omit()
cor(md_nox[,c(2:5)])

comb_exp_summary <- comb_data %>%
  filter(pid %in% c(hs_ids, md_ids)) %>%
  select(-bc_pred, -nox_pred) %>%
  filter(out_period == 0) %>%
  na.omit() %>%
  group_by(exp_period) %>%
  summarize(n = n(),
            pm2.5 = mean_sd(pm),
            no2 = mean_sd(no2),
            o3 = mean_sd(o3),
            rh = mean_sd(rhavg),
            temp_C = mean_sd(tavg_c),
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
comb_exp_summary_t <- t(comb_exp_summary[,-1])
colnames(comb_exp_summary_t) <- t(comb_exp_summary[,1])
comb_exp_summary_t <- as.data.frame(comb_exp_summary_t)
comb_exp_summary_t$var <- rownames(comb_exp_summary_t)
comb_exp_summary_t <- comb_exp_summary_t[,c(5,1:4)]
write_csv(comb_exp_summary_t, here::here("results", "tables", "comb_exposure_summary.csv"))


#' -----------------------------------------------------------------------------
#' comparisons between HS and MADRES exposure variables
#' -----------------------------------------------------------------------------

hs_exp <- filter(comb_data, cohort == "hs") %>%
  select(-nox_pred) %>%
  filter(out_period == 0) %>%
  na.omit() %>%
  select(pm:tavg_c, dist_roads_m, ndvi_500m, 
         dist_parks_m, park_area_500m, park_count_500m,
         tree_cover_500m, impervious_500m,
         pct_hh_poverty, pct_less_hs_edu) %>%
  mutate(study = "hs")

md_exp <- filter(comb_data, cohort == "md") %>%
  select(-bc_pred) %>%
  filter(out_period == 0) %>%
  na.omit() %>%
  select(pm:tavg_c, dist_roads_m, ndvi_500m, 
         dist_parks_m, park_area_500m, park_count_500m,
         tree_cover_500m, impervious_500m,
         pct_hh_poverty, pct_less_hs_edu) %>%
  mutate(study = "md")

exp_df <- bind_rows(hs_exp, md_exp) %>%
  pivot_longer(cols = !study, names_to = "var", values_to = "value")

ggplot(exp_df) +
  geom_boxplot(aes(x = study, y = value, group = study, color = study)) +
  facet_wrap(~var, scales = "free") +
  scale_color_manual(name = "Study",
                     values = c("hs" = "red", "md" = "blue"),
                     labels = c("hs" = "Healthy Start", "md" = "MADRES")) +
  theme(legend.position = "bottom")
ggsave(filename = here::here("figs", "exposure_comparisons.jpeg"),
       device = "jpeg")











