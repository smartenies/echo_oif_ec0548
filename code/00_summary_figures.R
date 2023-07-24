library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(patchwork)
library(ggh4x)
library(viridis)

#' [Box Health] directory where the data MUST live
#' These are PHI and must be protected as such!
#' Do not save any object containing these data in any other directory!!
box_path <- "/Users/sheenamartenies/Library/CloudStorage/Box-Box/[Box Health - External] Echo_oif4d"

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
#' Figure 1: Boxplots of BMI z-scores by cohort
#' -----------------------------------------------------------------------------

hs_bmi <- filter(hs_data, exp_period == "pregnancy") %>%
  select(out_period, zbmi, zbmi_risk) %>%
  mutate(study = "hs")
md_bmi <- filter(md_data, exp_period == "pregnancy") %>%
  select(out_period, zbmi, zbmi_risk) %>%
  mutate(study = "md")

bmi_df <- bind_rows(hs_bmi, md_bmi) %>%
  filter(!is.na(out_period)) %>%
  filter(out_period != 72) %>%
  mutate(out_period = as.factor(out_period))
levels(bmi_df$out_period) <- c("Birth", "6 months", "12 months", "24 months")

ggplot(bmi_df) +
  geom_boxplot(aes(x = study, y = zbmi, group = study, color = study),
               show.legend = F) +
  facet_grid(. ~ out_period, scales = "free") +
  scale_x_discrete(labels = c("Healthy Start", "MADRES")) +
  scale_color_manual(name = "Cohort", 
                     values = c("hs" = "#440154FF", "md" = "#35B779FF"),
                     labels = c("hs" = "Healthy Start", "md" = "MADRES")) +
  xlab("") + ylab("BMI z-score")
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "figure1.jpeg"),
       device = "jpeg", units = "in", height = 3, width = 7)


#' -----------------------------------------------------------------------------
#' Figures 2 and 3: Correlations between variables for each cohort
#' -----------------------------------------------------------------------------

#' Exposures
exps_hs <- c("pm", "no2", "o3", "rhavg", "tavg_c", "dist_roads_m", "bc_pred",
             "ndvi_500m", "dist_parks_m", "park_area_500m", "park_count_500m",
             "tree_cover_500m", "impervious_500m", "lila", 
             "pct_hh_poverty", "pct_less_hs_edu")
exps_md <- c("pm", "no2", "o3", "rhavg", "tavg_c", "dist_roads_m", "nox_pred",
             "ndvi_500m", "dist_parks_m", "park_area_500m", "park_count_500m",
             "tree_cover_500m", "impervious_500m", "lila", 
             "pct_hh_poverty", "pct_less_hs_edu")

names_hs <- c("PM2.5", "NO2", "O3", "Relative humidity",
              "Temperature", "Distance to major roads", "Black carbon",
              "NDVI (500 m)",
              "Distance to parks", "Park area (500 m)", "Park count (500 m)",
              "% Tree cover (500 m)", "% Impervious (500 m)", 
              "LILA census tract", "% Households in poverty", 
              "% less than HS education")
names_md <- c("PM2.5", "NO2", "O3", "Relative humidity",
              "Temperature", "Distance to major roads", "Oxides of nitrogen",
              "NDVI (500 m)",
              "Distance to parks", "Park area (500 m)", "Park count (500 m)",
              "% Tree cover (500 m)", "% Impervious (500 m)", 
              "LILA census tract", "% Households in poverty", 
              "% less than HS education")

#'Covariates
covars <- c("maternal_age", "maternal_partner", 
            "ed_no_hs", "ed_hs", "ed_aa", "ed_4yr", "smokesh", "gestsmoking", 
            "ppbmi_underweight", "ppbmi_overweight", "ppbmi_obese_i", "ppbmi_obese_ii", "ppbmi_obese_iii", 
            "concep_spring", "concep_summer", "concep_fall", "male")

#' --------------------------------------
#' Pregnancy-wide exposures
#' --------------------------------------

df_0_hs <- hs_data %>%
  filter(exp_period == "pregnancy" & out_period == 0) %>%
  select(all_of(exps_hs), all_of(covars), zbmi) %>%
  na.omit()

#' Examine correlations between exposures (scaled)
ps_0_hs <- select(df_0_hs, all_of(exps_hs))
names(ps_0_hs) <- names_hs
corr_0_hs <- round(cor(ps_0_hs), 1)

plot_0_hs <- ggcorrplot(corr_0_hs, type = "upper", legend.title = "Correlation",
                        ggtheme = ggplot2::theme_gray,
                        outline.color = "white",
                        lab = T, lab_size = 2) +
  ggtitle("A) Healthy Start") +
  theme(axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = c(0.8, 0.3))


df_0_md <- md_data %>%
  filter(exp_period == "pregnancy" & out_period == 0) %>%
  select(all_of(exps_md), all_of(covars), zbmi) %>%
  na.omit()

#' Examine correlations between exposures (scaled)
ps_0_md <- select(df_0_md, all_of(exps_md))
names(ps_0_md) <- names_md
corr_0_md <- round(cor(ps_0_md), 1)

plot_0_md <- ggcorrplot(corr_0_md, type = "upper", legend.title = "Correlation",
                        ggtheme = ggplot2::theme_gray,
                        outline.color = "white",
                        lab = T, lab_size = 2) +
  ggtitle("B) MADRES") +
  theme(axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = c(0.8, 0.3))

plot_0_hs + plot_0_md
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "figure2.jpeg"),
       device = "jpeg", units = "in", height = 5, width = 10)

#' --------------------------------------
#' Trimester-specific exposures
#' --------------------------------------

#' T1
df_t1_hs <- hs_data %>%
  filter(exp_period == "t1" & out_period == 0) %>%
  select(all_of(exps_hs), all_of(covars), zbmi) %>%
  na.omit()

#' Examine correlations between exposures (scaled)
ps_t1_hs <- select(df_t1_hs, all_of(exps_hs))
names(ps_t1_hs) <- names_hs
corr_t1_hs <- round(cor(ps_t1_hs), 1)

plot_t1_hs <- ggcorrplot(corr_t1_hs, type = "upper", legend.title = "Correlation",
                         ggtheme = ggplot2::theme_gray,
                         outline.color = "white",
                         lab = T, lab_size = 2) +
  ggtitle("A) Healthy Start: First Trimester") +
  theme(axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = c(0.8, 0.3))


df_t1_md <- md_data %>%
  filter(exp_period == "t1" & out_period == 0) %>%
  select(all_of(exps_md), all_of(covars), zbmi) %>%
  na.omit()

#' Examine correlations between exposures (scaled)
ps_t1_md <- select(df_t1_md, all_of(exps_md))
names(ps_t1_md) <- names_md
corr_t1_md <- round(cor(ps_t1_md), 1)

plot_t1_md <- ggcorrplot(corr_t1_md, type = "upper", legend.title = "Correlation",
                         ggtheme = ggplot2::theme_gray,
                         outline.color = "white",
                         lab = T, lab_size = 2) +
  ggtitle("B) MADRES: First Trimester") +
  theme(axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = c(0.8, 0.3))

#' T2
df_t2_hs <- hs_data %>%
  filter(exp_period == "t2" & out_period == 0) %>%
  select(all_of(exps_hs), all_of(covars), zbmi) %>%
  na.omit()

#' Examine correlations between exposures (scaled)
ps_t2_hs <- select(df_t2_hs, all_of(exps_hs))
names(ps_t2_hs) <- names_hs
corr_t2_hs <- round(cor(ps_t2_hs), 1)

plot_t2_hs <- ggcorrplot(corr_t2_hs, type = "upper", legend.title = "Correlation",
                         ggtheme = ggplot2::theme_gray,
                         outline.color = "white",
                         lab = T, lab_size = 2) +
  ggtitle("C) Healthy Start: Second Trimester") +
  theme(axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = c(0.8, 0.3))


df_t2_md <- md_data %>%
  filter(exp_period == "t2" & out_period == 0) %>%
  select(all_of(exps_md), all_of(covars), zbmi) %>%
  na.omit()

#' Examine correlations between exposures (scaled)
ps_t2_md <- select(df_t2_md, all_of(exps_md))
names(ps_t2_md) <- names_md
corr_t2_md <- round(cor(ps_t2_md), 1)

plot_t2_md <- ggcorrplot(corr_t2_md, type = "upper", legend.title = "Correlation",
                         ggtheme = ggplot2::theme_gray,
                         outline.color = "white",
                         lab = T, lab_size = 2) +
  ggtitle("D) MADRES: Second Trimester") +
  theme(axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = c(0.8, 0.3))

#' T3
df_t3_hs <- hs_data %>%
  filter(exp_period == "t3" & out_period == 0) %>%
  select(all_of(exps_hs), all_of(covars), zbmi) %>%
  na.omit()

#' Examine correlations between exposures (scaled)
ps_t3_hs <- select(df_t3_hs, all_of(exps_hs))
names(ps_t3_hs) <- names_hs
corr_t3_hs <- round(cor(ps_t3_hs), 1)

plot_t3_hs <- ggcorrplot(corr_t3_hs, type = "upper", legend.title = "Correlation",
                         ggtheme = ggplot2::theme_gray,
                         outline.color = "white",
                         lab = T, lab_size = 2) +
  ggtitle("E) Healthy Start: Third Trimester") +
  theme(axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = c(0.8, 0.3))


df_t3_md <- md_data %>%
  filter(exp_period == "t3" & out_period == 0) %>%
  select(all_of(exps_md), all_of(covars), zbmi) %>%
  na.omit()

#' Examine correlations between exposures (scaled)
ps_t3_md <- select(df_t3_md, all_of(exps_md))
names(ps_t3_md) <- names_md
corr_t3_md <- round(cor(ps_t3_md), 1)

plot_t3_md <- ggcorrplot(corr_t3_md, type = "upper", legend.title = "Correlation",
                         ggtheme = ggplot2::theme_gray,
                         outline.color = "white",
                         lab = T, lab_size = 2) +
  ggtitle("F) MADRES: Third Trimester") +
  theme(axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = c(0.8, 0.3))


(plot_t1_hs + plot_t1_md) / (plot_t2_hs + plot_t2_md) / (plot_t3_hs + plot_t3_md)
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "figure3.jpeg"),
       device = "jpeg", units = "in", height = 14, width = 12)

#' -----------------------------------------------------------------------------
#' Figures 4, 5, 6, and 7: Exposure weights from quantile-based g-comp models
#' -----------------------------------------------------------------------------

exps_gis <- c("pm", "no2", "o3", "rhavg", "tavg_c", "dist_roads_m", 
             "ndvi_500m", "dist_parks_m", "park_area_500m", "park_count_500m",
             "tree_cover_500m", "impervious_500m", 
             "pct_hh_poverty", "pct_less_hs_edu")
exps_st <- c("pm", "no2", "o3", "rhavg", "tavg_c", 
             "bc_pred", "nox_pred",
             "ndvi_500m", "dist_parks_m", "park_area_500m", "park_count_500m",
             "tree_cover_500m", "impervious_500m", 
             "pct_hh_poverty", "pct_less_hs_edu")

names_gis <- c("PM2.5", "NO2", "O3", "Relative humidity",
              "Temperature", "Distance to major roads", 
              "NDVI (500 m)",
              "Distance to parks", "Park area (500 m)", "Park count (500 m)",
              "% Tree cover (500 m)", "% Impervious (500 m)", 
              "% Households in poverty", 
              "% less than HS education")
names_st <- c("PM2.5", "NO2", "O3", "Relative humidity",
              "Temperature", "Distance to major roads", 
              "Black carbon", "Oxides of nitrogen",
              "NDVI (500 m)",
              "Distance to parks", "Park area (500 m)", "Park count (500 m)",
              "% Tree cover (500 m)", "% Impervious (500 m)", 
              "% Households in poverty", 
              "% less than HS education")

#' ------------------------
#' Pregnancy-wide exposures, GIS variables
#' ------------------------

weights_hs <- read_csv(here::here("results", "qgcomp_results", "weights", 
                                  "qgcomp_pregnancy_weights_hs.csv")) %>%
  mutate(cohort = "hs")
weights_md <- read_csv(here::here("results", "qgcomp_results", "weights", 
                                  "qgcomp_pregnancy_weights_md.csv")) %>%
  mutate(cohort = "md")

weights <- bind_rows(weights_hs, weights_md)

weights_0_gis_pos <- filter(weights, group == "gis variables") %>%
  filter(weights > 0) %>%
  mutate(time_point = as.factor(time_point))
weights_0_gis_neg <- filter(weights, group == "gis variables") %>%
  filter(weights < 0) %>%
  mutate(time_point = as.factor(time_point))

levels(weights_0_gis_pos$time_point) <- c("Birth", "6 months", "12 months", "24 months")
levels(weights_0_gis_neg$time_point) <- c("Birth", "6 months", "12 months", "24 months")

ggplot() +
  geom_linerange(data = weights_0_gis_pos, aes(x = exposure, ymin = 0, ymax = weights, 
                                      group = cohort, color = cohort),
                 linewidth = 2,
                 position = position_dodge(width = 0.6)) +
  geom_linerange(data = weights_0_gis_neg, aes(x = exposure, ymin = 0, ymax = weights, 
                                      group = cohort, color = cohort),
                 linewidth = 2,
                 position = position_dodge(width = 0.6)) +
  coord_flip() +
  xlab("") +
  scale_color_manual(name = "Cohort",
                     values = c("hs" = "#440154FF", "md" = "#35B779FF"),
                     labels = c("hs" = "Healthy Start", "md" = "MADRES")) +
  scale_x_discrete(limits = rev(exps_gis), labels = rev(names_gis)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  facet_grid(. ~ time_point) +
  theme(legend.position = "bottom")
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "figure4.jpeg"),
       device = "jpeg", units = "in", height = 4, width = 8)

#' ------------------------
#' Trimester-specific exposures, GIS variables
#' ------------------------

weights_hs_t1 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t1_weights_hs.csv")) %>%
  mutate(cohort = "hs")
weights_hs_t2 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t2_weights_hs.csv")) %>%
  mutate(cohort = "hs")
weights_hs_t3 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t3_weights_hs.csv")) %>%
  mutate(cohort = "hs")

weights_md_t1 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t1_weights_md.csv")) %>%
  mutate(cohort = "md")
weights_md_t2 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t2_weights_md.csv")) %>%
  mutate(cohort = "md")
weights_md_t3 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t3_weights_md.csv")) %>%
  mutate(cohort = "md")

weights <- bind_rows(weights_hs_t1, weights_hs_t2, weights_hs_t3, 
                     weights_md_t1, weights_md_t2, weights_md_t3)

weights_0_gis_pos <- filter(weights, group == "gis variables") %>%
  filter(weights > 0) %>%
  mutate(time_point = as.factor(time_point),
         exposure_period = as.factor(exposure_period))
weights_0_gis_neg <- filter(weights, group == "gis variables") %>%
  filter(weights < 0) %>%
  mutate(time_point = as.factor(time_point),
         exposure_period = as.factor(exposure_period))

levels(weights_0_gis_pos$time_point) <- c("Birth", "6 months", "12 months", "24 months")
levels(weights_0_gis_neg$time_point) <- c("Birth", "6 months", "12 months", "24 months")
levels(weights_0_gis_pos$exposure_period) <- c("Trimester 1", "Trimester 2", "Trimester 3")
levels(weights_0_gis_neg$exposure_period) <- c("Trimester 1", "Trimester 2", "Trimester 3")

ggplot() +
  geom_linerange(data = weights_0_gis_pos, aes(x = exposure, ymin = 0, ymax = weights, 
                                               group = cohort, color = cohort),
                 linewidth = 2,
                 position = position_dodge(width = 0.6)) +
  geom_linerange(data = weights_0_gis_neg, aes(x = exposure, ymin = 0, ymax = weights, 
                                               group = cohort, color = cohort),
                 linewidth = 2,
                 position = position_dodge(width = 0.6)) +
  coord_flip() +
  xlab("") +
  scale_color_manual(name = "Cohort",
                     values = c("hs" = "#440154FF", "md" = "#35B779FF"),
                     labels = c("hs" = "Healthy Start", "md" = "MADRES")) +
  scale_x_discrete(limits = rev(exps_gis), labels = rev(names_gis)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  facet_grid(exposure_period ~ time_point) +
  theme(legend.position = "bottom")
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "figure5.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

#' ------------------------
#' Pregnancy-wide exposures, GIS variables
#' ------------------------

weights_hs <- read_csv(here::here("results", "qgcomp_results", "weights", 
                                  "qgcomp_pregnancy_weights_hs.csv")) %>%
  mutate(cohort = "hs")
weights_md <- read_csv(here::here("results", "qgcomp_results", "weights", 
                                  "qgcomp_pregnancy_weights_md.csv")) %>%
  mutate(cohort = "md")

weights <- bind_rows(weights_hs, weights_md)

weights_0_gis_pos <- filter(weights, group == "st variables") %>%
  filter(weights > 0) %>%
  mutate(time_point = as.factor(time_point))
weights_0_gis_neg <- filter(weights, group == "st variables") %>%
  filter(weights < 0) %>%
  mutate(time_point = as.factor(time_point))

levels(weights_0_gis_pos$time_point) <- c("Birth", "6 months", "12 months", "24 months")
levels(weights_0_gis_neg$time_point) <- c("Birth", "6 months", "12 months", "24 months")

ggplot() +
  geom_linerange(data = weights_0_gis_pos, aes(x = exposure, ymin = 0, ymax = weights, 
                                               group = cohort, color = cohort),
                 linewidth = 2,
                 position = position_dodge(width = 0.6)) +
  geom_linerange(data = weights_0_gis_neg, aes(x = exposure, ymin = 0, ymax = weights, 
                                               group = cohort, color = cohort),
                 linewidth = 2,
                 position = position_dodge(width = 0.6)) +
  coord_flip() +
  xlab("") +
  scale_color_manual(name = "Cohort",
                     values = c("hs" = "#440154FF", "md" = "#35B779FF"),
                     labels = c("hs" = "Healthy Start", "md" = "MADRES")) +
  scale_x_discrete(limits = rev(exps_st), labels = rev(names_st)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  facet_grid(. ~ time_point) +
  theme(legend.position = "bottom")
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "figure6.jpeg"),
       device = "jpeg", units = "in", height = 4, width = 8)

#' ------------------------
#' Trimester-specific exposures, GIS variables
#' ------------------------

weights_hs_t1 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t1_weights_hs.csv")) %>%
  mutate(cohort = "hs")
weights_hs_t2 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t2_weights_hs.csv")) %>%
  mutate(cohort = "hs")
weights_hs_t3 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t3_weights_hs.csv")) %>%
  mutate(cohort = "hs")

weights_md_t1 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t1_weights_md.csv")) %>%
  mutate(cohort = "md")
weights_md_t2 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t2_weights_md.csv")) %>%
  mutate(cohort = "md")
weights_md_t3 <- read_csv(here::here("results", "qgcomp_results", "weights",
                                     "qgcomp_t3_weights_md.csv")) %>%
  mutate(cohort = "md")

weights <- bind_rows(weights_hs_t1, weights_hs_t2, weights_hs_t3, 
                     weights_md_t1, weights_md_t2, weights_md_t3)

weights <- bind_rows(weights_hs_t1, weights_hs_t2, weights_hs_t3, 
                     weights_md_t1, weights_md_t2, weights_md_t3)

weights_0_gis_pos <- filter(weights, group == "st variables") %>%
  filter(weights > 0) %>%
  mutate(time_point = as.factor(time_point),
         exposure_period = as.factor(exposure_period))
weights_0_gis_neg <- filter(weights, group == "st variables") %>%
  filter(weights < 0) %>%
  mutate(time_point = as.factor(time_point),
         exposure_period = as.factor(exposure_period))

levels(weights_0_gis_pos$time_point) <- c("Birth", "6 months", "12 months", "24 months")
levels(weights_0_gis_neg$time_point) <- c("Birth", "6 months", "12 months", "24 months")
levels(weights_0_gis_pos$exposure_period) <- c("Trimester 1", "Trimester 2", "Trimester 3")
levels(weights_0_gis_neg$exposure_period) <- c("Trimester 1", "Trimester 2", "Trimester 3")

ggplot() +
  geom_linerange(data = weights_0_gis_pos, aes(x = exposure, ymin = 0, ymax = weights, 
                                               group = cohort, color = cohort),
                 linewidth = 2,
                 position = position_dodge(width = 0.6)) +
  geom_linerange(data = weights_0_gis_neg, aes(x = exposure, ymin = 0, ymax = weights, 
                                               group = cohort, color = cohort),
                 linewidth = 2,
                 position = position_dodge(width = 0.6)) +
  coord_flip() +
  xlab("") +
  scale_color_manual(name = "Cohort",
                     values = c("hs" = "#440154FF", "md" = "#35B779FF"),
                     labels = c("hs" = "Healthy Start", "md" = "MADRES")) +
  scale_x_discrete(limits = rev(exps_st), labels = rev(names_st)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  facet_grid(exposure_period ~ time_point) +
  theme(legend.position = "bottom")
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "figure7.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

#' -----------------------------------------------------------------------------
#' BKMR exposure response functions
#' A) Full mixture
#' B) By exposure
#' -----------------------------------------------------------------------------

library(bkmr)

exps <- c("pm", "o3", "no2", "rhavg", "tavg_c", 
          "bc_pred", "nox_pred", "dist_roads_m",
          "ndvi_500m", "dist_parks_m", "park_area_500m", 
          "park_count_500m", "tree_cover_500m", "impervious_500m",
          "lila", "pct_hh_poverty", "pct_less_hs_edu")
exp_labs <- c("PM2.5", "O3", "NO2", "Rel. Humidity", "Temperature",
              "Black Carbon", "NOx", "Dist. to Roads", "NDVI",
              "Dist. to Parks", "Park Area", "Park Count",
              "% Tree Cover", "% Impervious", "LILA Tracts",
              "% HH Poverty", "% Low Education")
names(exp_labs) <- c("pm", "o3", "no2", "rhavg", "tavg_c", 
                     "bc_pred", "nox_pred", "dist_roads_m",
                     "ndvi_500m", "dist_parks_m", "park_area_500m", 
                     "park_count_500m", "tree_cover_500m", "impervious_500m",
                     "lila", "pct_hh_poverty", "pct_less_hs_edu")

#' --------------------------------------
#' Pregnancy-wide Healthy Start: GIS variables
#' Pregnancy-wide Healthy Start: ST variables
#' --------------------------------------

#' Birth: GIS
mod_p_0 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_results_hs.rdata")
load(here::here("results/bkmr_results", mod_p_0))

pred_univariate_50_p_0 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_0 <- o_risks %>%
  mutate(out_period = "Birth",
         exp_set = "GIS TRAP")

s_risks_p_0 <- s_risks %>%
  mutate(out_period = "Birth",
         exp_set = "GIS TRAP")

er_p_0 <- ggplot(pred_univariate_50_p_0, 
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_0

#' Birth: ST
mod_p_st_0 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_results_st_hs.rdata")
load(here::here("results/bkmr_results", mod_p_st_0))

pred_univariate_50_p_st_0 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_st_0 <- o_risks %>%
  mutate(out_period = "Birth",
         exp_set = "ST TRAP")

s_risks_p_st_0 <- s_risks %>%
  mutate(out_period = "Birth",
         exp_set = "ST TRAP")

er_p_st_0 <- ggplot(pred_univariate_50_p_st_0, 
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_st_0

#' 6 Months: GIS
mod_p_6 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_results_hs.rdata")
load(here::here("results/bkmr_results", mod_p_6))

pred_univariate_50_p_6 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_6 <- o_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "GIS TRAP")

s_risks_p_6 <- s_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "GIS TRAP")

er_p_6 <- ggplot(pred_univariate_50_p_6, 
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_6

#' 6 Months: ST
mod_p_st_6 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_results_st_hs.rdata")
load(here::here("results/bkmr_results", mod_p_st_6))

pred_univariate_50_p_st_6 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_st_6 <- o_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "ST TRAP")

s_risks_p_st_6 <- s_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "ST TRAP")

er_p_st_6 <- ggplot(pred_univariate_50_p_st_6, 
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_st_6

#' 12 Months: GIS
mod_p_12 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_results_hs.rdata")
load(here::here("results/bkmr_results", mod_p_12))

pred_univariate_50_p_12 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_12 <- o_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "GIS TRAP")

s_risks_p_12 <- s_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "GIS TRAP")

er_p_12 <- ggplot(pred_univariate_50_p_12, 
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_12

#' 12 Months: ST
mod_p_st_12 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_results_st_hs.rdata")
load(here::here("results/bkmr_results", mod_p_st_12))

pred_univariate_50_p_st_12 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_st_12 <- o_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "ST TRAP")

s_risks_p_st_12 <- s_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "ST TRAP")

er_p_st_12 <- ggplot(pred_univariate_50_p_st_12, 
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_st_12

#' 24 Months: GIS
mod_p_24 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_results_hs.rdata")
load(here::here("results/bkmr_results", mod_p_24))

pred_univariate_50_p_24 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_24 <- o_risks %>%
  mutate(out_period = "24 Months",
         exp_set = "GIS TRAP")

s_risks_p_24 <- s_risks %>%
  mutate(out_period = "24 Months",
         exp_set = "GIS TRAP")

er_p_24 <- ggplot(pred_univariate_50_p_24, 
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_24

#' 24 Months: ST
mod_p_st_24 <- paste0("BKMR_", "pregnancy", "_", 24, "_zbmi_results_st_hs.rdata")
load(here::here("results/bkmr_results", mod_p_st_24))

pred_univariate_50_p_st_24 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_st_24 <- o_risks %>%
  mutate(out_period = "24 Months",
         exp_set = "ST TRAP")

s_risks_p_st_24 <- s_risks %>%
  mutate(out_period = "24 Months",
         exp_set = "ST TRAP")

er_p_st_24 <- ggplot(pred_univariate_50_p_st_24, 
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_st_24

#' Summary Figures: Healthy Start
library(patchwork)

#' Exposure Response
er_p_0 / er_p_st_0 +
  plot_annotation(title = "Healthy Start Cohort",
                  tag_levels = list(c("A: GIS TRAP", "B: ST TRAP")))
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_p_0_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

er_p_6 / er_p_st_6 +
  plot_annotation(title = "Healthy Start Cohort",
                  tag_levels = list(c("A: GIS TRAP", "B: ST TRAP")))
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_p_6_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

er_p_12 / er_p_st_12 +
  plot_annotation(title = "Healthy Start Cohort",
                  tag_levels = list(c("A: GIS TRAP", "B: ST TRAP")))
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_p_12_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

er_p_24 / er_p_st_24 +
  plot_annotation(title = "Healthy Start Cohort",
                  tag_levels = list(c("A: GIS TRAP", "B: ST TRAP")))
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_p_24_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

#' Overall Risk
o_risks <- bind_rows(o_risks_p_0, o_risks_p_6, o_risks_p_12, o_risks_p_24,
                           o_risks_p_st_0, o_risks_p_st_6, o_risks_p_st_12, o_risks_p_st_24) %>%
  mutate(out_period = factor(out_period, levels = c("Birth", "6 Months", "12 Months", "24 Months")),
         exp_set = factor(exp_set, levels = c("GIS TRAP", "ST TRAP")))

ro_p_hs <- ggplot(o_risks, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  ggtitle("Healthy Start Cohort") +
  facet_grid(out_period ~ exp_set) +
  xlab("Exposure quantile") +
  ylab("Change in BMI z-score") +
  geom_pointrange()
ro_p_hs
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_ro_p_hs.jpeg"),
       device = "jpeg", units = "in", height = 8, width = 10)

#' Single variable risks
s_risks <- bind_rows(s_risks_p_0, s_risks_p_6, s_risks_p_12, s_risks_p_24,
                           s_risks_p_st_0, s_risks_p_st_6, s_risks_p_st_12, s_risks_p_st_24) %>%
  mutate(out_period = factor(out_period, levels = c("Birth", "6 Months", "12 Months", "24 Months")),
         exp_set = factor(exp_set, levels = c("GIS TRAP", "ST TRAP")),
         variable = factor(variable, levels = exps, labels = exp_labs))

rs_p_hs <- ggplot(s_risks, aes(variable, est, ymin = est - 1.96*sd, 
                                        ymax = est + 1.96*sd, col = q.fixed)) +
  ggtitle("Healthy Start Cohort") +
  facet_grid(out_period ~ exp_set, scales = "free") +
  scale_color_manual(values = c("0.25" = viridis(n = 10)[1],
                                "0.5" = viridis(n = 10)[5],
                                "0.75" = viridis(n = 10)[9]),
                     name = "Quantile of other exposures") +
  geom_pointrange(position = position_dodge(width = 0.75)) +
  xlab("") + ylab("Change in BMI z-score") +
  scale_x_discrete(labels = exp_labs, drop = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom")
rs_p_hs
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_rs_p_hs.jpeg"),
       device = "jpeg", units = "in", height = 8, width = 10)

#' --------------------------------------
#' Trimester-specific Healthy Start: GIS variables
#' Trimester-specific Healthy Start: ST variables
#' --------------------------------------

t_names <- c("_t1", "_t2", "_t3")

vs1 <- c("pm", "o3", "no2", "rhavg", "tavg_c", 
         "bc_pred", "nox_pred", "dist_roads_m",
         "ndvi_500m", "dist_parks_m", "park_area_500m", 
         "park_count_500m", "tree_cover_500m", "impervious_500m",
         "lila", "pct_hh_poverty", "pct_less_hs_edu")
vs2 <- c("t1", "t2", "t3")
var_sort <- apply(expand.grid(vs1, vs2), 1, paste, collapse="_")

#' Birth: GIS
mod_t_0 <- paste0("BKMR_", "trimesters", "_", 0, "_zbmi_results_hs.rdata")
load(here::here("results/bkmr_results", mod_t_0))

pred_univariate_50_t_0 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_0 <- o_risks %>%
  mutate(out_period = "Birth",
         exp_set = "GIS TRAP") 

s_risks_t_0 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "Birth",
         exp_set = "GIS TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_0 <- ggplot(pred_univariate_50_t_0,
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_0
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_0_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' Birth: ST
mod_t_st_0 <- paste0("BKMR_", "trimesters", "_", 0, "_zbmi_results_st_hs.rdata")
load(here::here("results/bkmr_results", mod_t_st_0))

pred_univariate_50_t_st_0 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_st_0 <- o_risks %>%
  mutate(out_period = "Birth",
         exp_set = "ST TRAP") 

s_risks_t_st_0 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "Birth",
         exp_set = "ST TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_st_0 <- ggplot(pred_univariate_50_t_st_0,
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_st_0
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_st_0_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 6 Months: GIS
mod_t_6 <- paste0("BKMR_", "trimesters", "_", 6, "_zbmi_results_hs.rdata")
load(here::here("results/bkmr_results", mod_t_6))

pred_univariate_50_t_6 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_6 <- o_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "GIS TRAP") 

s_risks_t_6 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  
  mutate(out_period = "6 Months",
         exp_set = "GIS TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_6 <- ggplot(pred_univariate_50_t_6,
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_6
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_6_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 6 Months: ST
mod_t_st_6 <- paste0("BKMR_", "trimesters", "_", 6, "_zbmi_results_st_hs.rdata")
load(here::here("results/bkmr_results", mod_t_st_6))

pred_univariate_50_t_st_6 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_st_6 <- o_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "ST TRAP") 

s_risks_t_st_6 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "6 Months",
         exp_set = "ST TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_st_6 <- ggplot(pred_univariate_50_t_st_6,
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_st_6
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_st_6_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 12 Months: GIS
mod_t_12 <- paste0("BKMR_", "trimesters", "_", 12, "_zbmi_results_hs.rdata")
load(here::here("results/bkmr_results", mod_t_12))

pred_univariate_50_t_12 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_12 <- o_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "GIS TRAP") 

s_risks_t_12 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "12 Months",
         exp_set = "GIS TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_12 <- ggplot(pred_univariate_50_t_12, 
                  aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_12
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_12_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 12 Months: ST
mod_t_st_12 <- paste0("BKMR_", "trimesters", "_", 12, "_zbmi_results_st_hs.rdata")
load(here::here("results/bkmr_results", mod_t_st_12))

pred_univariate_50_t_st_12 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_st_12 <- o_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "ST TRAP") 

s_risks_t_st_12 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "12 Months",
         exp_set = "ST TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_st_12 <- ggplot(pred_univariate_50_t_st_12, 
                     aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_st_12
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_st_12_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 24 Months: GIS
mod_t_24 <- paste0("BKMR_", "trimesters", "_", 24, "_zbmi_results_hs.rdata")
load(here::here("results/bkmr_results", mod_t_24))

pred_univariate_50_t_24 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_24 <- o_risks %>%
  mutate(out_period = "24 Months",
         exp_set = "GIS TRAP") 

s_risks_t_24 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "24 Months",
         exp_set = "GIS TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_24 <- ggplot(pred_univariate_50_t_24,
                  aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_24
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_24_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 24 Months: ST
mod_t_st_24 <- paste0("BKMR_", "trimesters", "_", 24, "_zbmi_results_st_hs.rdata")
load(here::here("results/bkmr_results", mod_t_st_24))

pred_univariate_50_t_st_24 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_st_24 <- o_risks %>%
  mutate(out_period = "24 Months",
         exp_set = "ST TRAP") 

s_risks_t_st_24 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "24 Months",
         exp_set = "ST TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_st_24 <- ggplot(pred_univariate_50_t_st_24, 
                     aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_st_24
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_st_24_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' Summary Figures: Healthy Start

#' Overall Risk
o_risks <- bind_rows(o_risks_t_0, o_risks_t_6, o_risks_t_12, o_risks_t_24,
                           o_risks_t_st_0, o_risks_t_st_6, o_risks_t_st_12, o_risks_t_st_24) %>%
  mutate(out_period = factor(out_period, levels = c("Birth", "6 Months", "12 Months", "24 Months")),
         exp_set = factor(exp_set, levels = c("GIS TRAP", "ST TRAP")))

ro_t_hs <- ggplot(o_risks, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  ggtitle("Healthy Start Cohort") +
  facet_grid(out_period ~ exp_set) +
  xlab("Exposure quantile") +
  ylab("Change in BMI z-score") +
  geom_pointrange()
ro_t_hs
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_ro_t_hs.jpeg"),
       device = "jpeg", units = "in", height = 8, width = 10)

#' Single variable risks
s_risks <- bind_rows(s_risks_t_0, s_risks_t_6, s_risks_t_12, s_risks_t_24,
                           s_risks_t_st_0, s_risks_t_st_6, s_risks_t_st_12, s_risks_t_st_24) %>%
  mutate(out_period = factor(out_period, levels = c("Birth", "6 Months", "12 Months", "24 Months")),
         exp_set = factor(exp_set, levels = c("GIS TRAP", "ST TRAP")),
         #variable = factor(variable, levels = exps, labels = exp_labs),
         trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("T1", "T2", "T3")))

rs_t_hs <- ggplot(s_risks, 
                  aes(variable, est, ymin = est - 1.96*sd, 
                      ymax = est + 1.96*sd, col = q.fixed)) +
  geom_pointrange(position = position_dodge(width = 0.75)) +
  facet_nested(out_period + trimester ~ exp_set , scales = "free") +
  scale_color_manual(values = c("0.25" = viridis(n = 10)[1],
                                "0.5" = viridis(n = 10)[5],
                                "0.75" = viridis(n = 10)[9]),
                     name = "Quantile of other exposures") +
  xlab("") + ylab("Change in BMI z-score") +
  scale_x_discrete(labels = exp_labs, drop = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom")
rs_t_hs
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_rs_t_hs.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

#' --------------------------------------
#' Pregnancy-wide MADRES: GIS variables
#' Pregnancy-wide MADRES: ST variables
#' --------------------------------------

#' Birth: GIS
mod_p_0 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_results_md.rdata")
load(here::here("results/bkmr_results", mod_p_0))

pred_univariate_50_p_0 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_0 <- o_risks %>%
  mutate(out_period = "Birth",
         exp_set = "GIS TRAP")

s_risks_p_0 <- s_risks %>%
  mutate(out_period = "Birth",
         exp_set = "GIS TRAP")

er_p_0 <- ggplot(pred_univariate_50_p_0, 
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_0

#' Birth: ST
mod_p_st_0 <- paste0("BKMR_", "pregnancy", "_", 0, "_zbmi_results_st_md.rdata")
load(here::here("results/bkmr_results", mod_p_st_0))

pred_univariate_50_p_st_0 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_st_0 <- o_risks %>%
  mutate(out_period = "Birth",
         exp_set = "ST TRAP")

s_risks_p_st_0 <- s_risks %>%
  mutate(out_period = "Birth",
         exp_set = "ST TRAP")

er_p_st_0 <- ggplot(pred_univariate_50_p_st_0, 
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_st_0

#' 6 Months: GIS
mod_p_6 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_results_md.rdata")
load(here::here("results/bkmr_results", mod_p_6))

pred_univariate_50_p_6 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_6 <- o_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "GIS TRAP")

s_risks_p_6 <- s_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "GIS TRAP")

er_p_6 <- ggplot(pred_univariate_50_p_6, 
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_6

#' 6 Months: ST
mod_p_st_6 <- paste0("BKMR_", "pregnancy", "_", 6, "_zbmi_results_st_md.rdata")
load(here::here("results/bkmr_results", mod_p_st_6))

pred_univariate_50_p_st_6 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_st_6 <- o_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "ST TRAP")

s_risks_p_st_6 <- s_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "ST TRAP")

er_p_st_6 <- ggplot(pred_univariate_50_p_st_6, 
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_st_6

#' 12 Months: GIS
mod_p_12 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_results_md.rdata")
load(here::here("results/bkmr_results", mod_p_12))

pred_univariate_50_p_12 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_12 <- o_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "GIS TRAP")

s_risks_p_12 <- s_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "GIS TRAP")

er_p_12 <- ggplot(pred_univariate_50_p_12, 
                  aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_12

#' 12 Months: ST
mod_p_st_12 <- paste0("BKMR_", "pregnancy", "_", 12, "_zbmi_results_st_md.rdata")
load(here::here("results/bkmr_results", mod_p_st_12))

pred_univariate_50_p_st_12 <- pred_univariate_50 %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_p_st_12 <- o_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "ST TRAP")

s_risks_p_st_12 <- s_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "ST TRAP")

er_p_st_12 <- ggplot(pred_univariate_50_p_st_12, 
                     aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable, scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_p_st_12

#' Summary Figures: Healthy Start
library(patchwork)

#' Exposure Response
er_p_0 / er_p_st_0 +
  plot_annotation(title = "MADRES Cohort",
                  tag_levels = list(c("A: GIS TRAP", "B: ST TRAP")))
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_p_0_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

er_p_6 / er_p_st_6 +
  plot_annotation(title = "MADRES Cohort",
                  tag_levels = list(c("A: GIS TRAP", "B: ST TRAP")))
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_p_6_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

er_p_12 / er_p_st_12 +
  plot_annotation(title = "MADRES Cohort",
                  tag_levels = list(c("A: GIS TRAP", "B: ST TRAP")))
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_p_12_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 10)

#' Overall Risk
o_risks <- bind_rows(o_risks_p_0, o_risks_p_6, o_risks_p_12,
                           o_risks_p_st_0, o_risks_p_st_6, o_risks_p_st_12) %>%
  mutate(out_period = factor(out_period, levels = c("Birth", "6 Months", "12 Months")),
         exp_set = factor(exp_set, levels = c("GIS TRAP", "ST TRAP")))

ro_p_md <- ggplot(o_risks, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  ggtitle("MADRES Cohort") +
  facet_grid(out_period ~ exp_set) +
  xlab("Exposure quantile") +
  ylab("Change in BMI z-score") +
  geom_pointrange()
ro_p_md
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_ro_p_md.jpeg"),
       device = "jpeg", units = "in", height = 6, width = 10)

#' Single variable risks
s_risks <- bind_rows(s_risks_p_0, s_risks_p_6, s_risks_p_12,
                           s_risks_p_st_0, s_risks_p_st_6, s_risks_p_st_12) %>%
  mutate(out_period = factor(out_period, levels = c("Birth", "6 Months", "12 Months", "24 Months")),
         exp_set = factor(exp_set, levels = c("GIS TRAP", "ST TRAP")),
         variable = factor(variable, levels = exps, labels = exp_labs))

rs_p_md <- ggplot(s_risks, aes(variable, est, ymin = est - 1.96*sd, 
                                     ymax = est + 1.96*sd, col = q.fixed)) +
  ggtitle("MADRES Cohort") +
  facet_grid(out_period ~ exp_set, scales = "free") +
  scale_color_manual(values = c("0.25" = viridis(n = 10)[1],
                                "0.5" = viridis(n = 10)[5],
                                "0.75" = viridis(n = 10)[9]),
                     name = "Quantile of other exposures") +
  geom_pointrange(position = position_dodge(width = 0.75)) +
  xlab("") + ylab("Change in\nBMI z-score") +
  scale_x_discrete(labels = exp_labs, drop = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom")
rs_p_md
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_rs_p_md.jpeg"),
       device = "jpeg", units = "in", height = 6, width = 10)

#' --------------------------------------
#' Trimester-specific MADRES: GIS variables
#' Trimester-specific MADRES: ST variables
#' --------------------------------------

t_names <- c("_t1", "_t2", "_t3")

vs1 <- c("pm", "o3", "no2", "rhavg", "tavg_c", 
         "bc_pred", "nox_pred", "dist_roads_m",
         "ndvi_500m", "dist_parks_m", "park_area_500m", 
         "park_count_500m", "tree_cover_500m", "impervious_500m",
         "lila", "pct_hh_poverty", "pct_less_hs_edu")
vs2 <- c("t1", "t2", "t3")
var_sort <- apply(expand.grid(vs1, vs2), 1, paste, collapse="_")

#' Birth: GIS
mod_t_0 <- paste0("BKMR_", "trimesters", "_", 0, "_zbmi_results_md.rdata")
load(here::here("results/bkmr_results", mod_t_0))

pred_univariate_50_t_0 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_0 <- o_risks %>%
  mutate(out_period = "Birth",
         exp_set = "GIS TRAP") 

s_risks_t_0 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "Birth",
         exp_set = "GIS TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_0 <- ggplot(pred_univariate_50_t_0,
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_0
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_0_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' Birth: ST
mod_t_st_0 <- paste0("BKMR_", "trimesters", "_", 0, "_zbmi_results_st_md.rdata")
load(here::here("results/bkmr_results", mod_t_st_0))

pred_univariate_50_t_st_0 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_st_0 <- o_risks %>%
  mutate(out_period = "Birth",
         exp_set = "ST TRAP") 

s_risks_t_st_0 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "Birth",
         exp_set = "ST TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_st_0 <- ggplot(pred_univariate_50_t_st_0,
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_st_0
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_st_0_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 6 Months: GIS
mod_t_6 <- paste0("BKMR_", "trimesters", "_", 6, "_zbmi_results_md.rdata")
load(here::here("results/bkmr_results", mod_t_6))

pred_univariate_50_t_6 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_6 <- o_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "GIS TRAP") 

s_risks_t_6 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  
  mutate(out_period = "6 Months",
         exp_set = "GIS TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_6 <- ggplot(pred_univariate_50_t_6,
                 aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_6
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_6_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 6 Months: ST
mod_t_st_6 <- paste0("BKMR_", "trimesters", "_", 6, "_zbmi_results_st_md.rdata")
load(here::here("results/bkmr_results", mod_t_st_6))

pred_univariate_50_t_st_6 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_st_6 <- o_risks %>%
  mutate(out_period = "6 Months",
         exp_set = "ST TRAP") 

s_risks_t_st_6 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "6 Months",
         exp_set = "ST TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_st_6 <- ggplot(pred_univariate_50_t_st_6,
                    aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_st_6
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_st_6_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 12 Months: GIS
mod_t_12 <- paste0("BKMR_", "trimesters", "_", 12, "_zbmi_results_md.rdata")
load(here::here("results/bkmr_results", mod_t_12))

pred_univariate_50_t_12 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_12 <- o_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "GIS TRAP") 

s_risks_t_12 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "12 Months",
         exp_set = "GIS TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_12 <- ggplot(pred_univariate_50_t_12, 
                  aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_12
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_12_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' 12 Months: ST
mod_t_st_12 <- paste0("BKMR_", "trimesters", "_", 12, "_zbmi_results_st_md.rdata")
load(here::here("results/bkmr_results", mod_t_st_12))

pred_univariate_50_t_st_12 <- pred_univariate_50 %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("Trimester 1", "Trimester 2", "Trimester 3"))) %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

o_risks_t_st_12 <- o_risks %>%
  mutate(out_period = "12 Months",
         exp_set = "ST TRAP") 

s_risks_t_st_12 <- s_risks %>%
  mutate(variable = as.character(variable)) %>%
  mutate(trimester = sapply(strsplit(variable,"_"), tail, 1)) %>%
  mutate(out_period = "12 Months",
         exp_set = "ST TRAP") %>%
  mutate(variable = str_remove_all(variable, paste(t_names, collapse = "|"))) %>%
  mutate(variable = factor(variable, levels = exps, labels = exp_labs))

er_t_st_12 <- ggplot(pred_univariate_50_t_st_12, 
                     aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
  geom_smooth(stat = "identity") +
  facet_nested_wrap(vars(trimester, variable), scales = "fixed", ncol = 5) +
  ylab("h(z) holding other exposures at 50th percentile") +
  xlab("Scaled exposure") 
er_t_st_12
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_er_t_st_12_md.jpeg"),
       device = "jpeg", units = "in", height = 12, width = 8)

#' Summary Figures: MADRES

#' Overall Risk
o_risks <- bind_rows(o_risks_t_0, o_risks_t_6, o_risks_t_12, 
                           o_risks_t_st_0, o_risks_t_st_6, o_risks_t_st_12) %>%
  mutate(out_period = factor(out_period, levels = c("Birth", "6 Months", "12 Months")),
         exp_set = factor(exp_set, levels = c("GIS TRAP", "ST TRAP")))

ro_t_md <- ggplot(o_risks, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +
  ggtitle("MADRES Cohort") +
  facet_grid(out_period ~ exp_set) +
  xlab("Exposure quantile") +
  ylab("Change in BMI z-score") +
  geom_pointrange()
ro_t_md
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_ro_t_md.jpeg"),
       device = "jpeg", units = "in", height = 6, width = 10)

#' Single variable risks
s_risks <- bind_rows(s_risks_t_0, s_risks_t_6, s_risks_t_12,
                           s_risks_t_st_0, s_risks_t_st_6, s_risks_t_st_12) %>%
  mutate(out_period = factor(out_period, levels = c("Birth", "6 Months", "12 Months")),
         exp_set = factor(exp_set, levels = c("GIS TRAP", "ST TRAP")),
         #variable = factor(variable, levels = exps, labels = exp_labs),
         trimester = factor(trimester, levels = c("t1", "t2", "t3"),
                            labels = c("T1", "T2", "T3")))

rs_t_md <- ggplot(s_risks, 
                  aes(variable, est, ymin = est - 1.96*sd, 
                      ymax = est + 1.96*sd, col = q.fixed)) +
  geom_pointrange(position = position_dodge(width = 0.75)) +
  facet_nested(out_period + trimester ~ exp_set , scales = "free") +
  scale_color_manual(values = c("0.25" = viridis(n = 10)[1],
                                "0.5" = viridis(n = 10)[5],
                                "0.75" = viridis(n = 10)[9]),
                     name = "Quantile of other exposures") +
  xlab("") + ylab("Change in BMI z-score") +
  scale_x_discrete(labels = exp_labs, drop = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom")
rs_t_md
ggsave(filename = here::here("manuscripts", "ec0548a_manuscript", "bkmr_rs_t_md.jpeg"),
       device = "jpeg", units = "in", height = 9, width = 10)