
# EC0584: Examining the joint effects of prenatal exposures to traffic-related air pollutants, other environmental hazards, and social stressors on childhood obesity outcomes

EC0548 is an OIF-funded ECHO project investigating the effects of neighborhood-level environmental exposures on BMI in early childhood. The goals of EC0548 are to:

- Ascertain associations between individual neighborhood-level exposures and BMI z-scores in early childhood
- Explore the effects of the mixture of neighborhood-level exposures on BMI z-scores using machine learning methods
- Examine whether using study-specific models of traffic-related air pollution (TRAP) improves the models used to measure these associations relative to GIS-based measures of TRAP exposure.

This study brings together data from two ECHO cohorts: [Healthy Start](https://healthystartstudy.org/) (based in the Denver metropolitan area) and [MADRES](https://madres.usc.edu/) (based in Los Angeles). 

This repository contains code used to fit the regression models described in our manuscript. A ["manual of operations"](https://smartenies.github.io/echo_oif_ec0548/manual_of_operations) with information on how we derived our exposure variables has also been provided and details how we developed our exposure data sets. The operations manual includes links to the publicly available data sets used in this study. A [data dictionary]() defining each variable in the data set is also available.

Participant data from the Healthy Start and MADRES cohorts were shared with the primary investigator (Sheena Martenies) under data use agreements with their respective universities. Researchers interested in using these data should contact the PIs of each cohort: Dana Dabelea (Healthy Start, University of Colorado Anschutz Medical Campus) and Tracy Bastain or Carrie Breton (MADRES; University of Southern California).

We used several R packages in the course of this work. The work relies heavily on the `qgcomp` (Keil, 2022) and `bkmr` (Bobb, 2023) packages. Links to the vignettes for these packages are available in the scripts. 

The following scripts are available in this repository:

- 00_summary_figures.R: code to generate the figures used in the manuscript and supplemental materials
- 00_summary_tables.R: code to generate the results tables presented in the manuscript and supplemental materials
- 01_summary_statistics.R: code to generate tables 1 and 2 presented in the manuscript 
- 02_single_pol_linear_models.R: code to fit the linear regression models for each cohort individually
- 03_single_pol_lmer_models.R: code to fit the linear mixed effect regression models for the combined cohort
- 04a_qgcomp_models_healthy_start.R: code to fit the quantile-based g-computation models for the Healthy Start cohort using exposures averaged across pregnancy
- 04b_qgcomp_models_healthy_start_t1.R: code to fit the quantile-based g-computation models for the Healthy Start cohort using exposures averaged during the first trimester
- 04c_qgcomp_models_healthy_start_t2.R: code to fit the quantile-based g-computation models for the Healthy Start cohort using exposures averaged during the second trimester
- 04d_qgcomp_models_healthy_start_t3.R: code to fit the quantile-based g-computation models for the Healthy Start cohort using exposures averaged during the third trimester
- 05a_qgcomp_models_madres.R: code to fit the quantile-based g-computation models for the MADRES cohort using exposures averaged across pregnancy
- 05b_qgcomp_models_madres_t1.R: code to fit the quantile-based g-computation models for the MADRES cohort using exposures averaged during the first trimester
- 05c_qgcomp_models_madres_t2.R: code to fit the quantile-based g-computation models for the MADRES cohort using exposures averaged during the second trimester
- 05d_qgcomp_models_madres_t3.R: code to fit the quantile-based g-computation models for the MADRES cohort using exposures averaged during the third trimester
- 06a_bkmr_models_healthy_start.R: code to fit BKMR models for the Healthy Start cohort using pregnancy-wide exposures
- 06b_bkmr_models_healthy_start_trimesters.R: code to fit BKMR models for the Healthy Start cohort using hierarchical variable selection for exposures averaged by trimester
- 07a_bkmr_models_madres.R: code to fit BKMR models for the MADRES cohort using pregnancy-wide exposures
- 07b_bkmr_models_madres_trimesters.R: code to fit BKMR models for the MADRES cohort using hierarchical variable selection for exposures averaged by trimester

## Acknowledgements

This work was supported by the Office Of The Director, National Institutes Of Health of the National Institutes of Health under Award Number U2COD023375. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

## References

Bobb J (2023). _bkmr: Bayesian Kernel Machine Regression_. R package version 0.2.2,
  <https://github.com/jenfb/bkmr>.
  
Keil A (2022). _qgcomp: Quantile G-Computation_. R package version 2.10.1,
  <https://CRAN.R-project.org/package=qgcomp>.


