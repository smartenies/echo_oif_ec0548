# EC0548a Manual of Operations

Sheena Martenies
July 25, 2023

The following outlines the data manipulation steps used to assign exposures to environmental and social exposures at the neighborhood level for Healthy Start and MADRES participants. These exposures were obtained from secondary, publicly available data sets.

## A. Data Sources

The following table summarizes the publicly available data sources that were used to estimate residential and census tract exposures for the Healthy Start and MADRES participants. Differences by cohort are noted in the table below.

Data Set | Variables | Link to the Data Set
----------| ------ | ----
Criteria Pollutants | PM~2.5~, O~3~, NO~2~ | https://aqs.epa.gov/aqsweb/airdata
Meteorology^1^ | Temperature, Relative humidity | https://www.climatologylab.org/gridmet.html
Green Space | NDVI (MYD13Q1 v061)^2^ | https://lpdaac.usgs.gov/products/myd13q1v061/
Built Environment | Parks^3^ | https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2022&layergroup=Landmarks
Built Environment | Impervious Surface | https://www.mrlc.gov/data 
Built Environment | Tree Canopy | https://www.mrlc.gov/data 
Built Environment | Major Roads^4^ | https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2021&layergroup=Roads  
Food Access | LILA Tract^5^ | https://www.ers.usda.gov/data-products/food-access-research-atlas/download-the-data/
Neighborhood SES^6^ | Educational attainment^7^ | https://www.census.gov/data/developers/data-sets/acs-5year.html
Neighborhood SES^6^ | Households in poverty | https://www.census.gov/data/developers/data-sets/acs-5year.html

_^1^ For the Healthy Start cohort, these data were accessed via the climateR package for R: https://github.com/mikejohnson51/climateR_

_^2^ MODIS/Aqua Vegetation Indices 16-Day L3 Global 250 m SIN Grid_

_^3^ Parks were defined using MAF/TIGER Feature Class Codes (MTFCC) K2180 through K2190._

_^4^ For MADRES, estimates of distance to major roads were based on data a proprietary road network data set._

_^5^ Census tracts identified as low income and low access tracts (based on 0.5 miles for urban census tracts and 10 miles for rural tracts)._

_^6^ For Healthy Start, variables at the census tract level were accessed using the tidycensus package for R: https://walker-data.com/tidycensus/. Note that the use of tidycensus requires a free API key from the US Census Bureau._

_^7^ Defined as the percentage of the population ages 25 and older without a high school diploma or equivalent._

## B. Data Processing
The following describes the data cleaning and processing steps to obtain the exposure data sets for the Healthy Start and MADRES cohorts. 

Assigning exposures to cohort participants required the use of geocoded address information. For Healthy Start, addresses were geocoded using the Google geocoding API via the ggmap package in R (https://github.com/dkahle/ggmap). Note that the use of the Google geocoding API requires a key from Google. MADRES data were geocoded by analysts at the University. of Southern California. For additional details or to use MADRES data, please contact the study: https://madres.usc.edu/investigators/

For area-based variables, we selected a 500 m (0.3 mile) buffer to represent a reasonable walking distance during pregnancy. Although people who walk for transportation often walk farther than 0.25 miles regularly [1], physical activity research has found that pregnant people tend to walk less than recommended amounts [2,3]. Therefore, we opted to use the typical buffer distance to characterize the environment around the homes of participants during their pregnancies.

When residential history data were available, values were time-weighted based on the time spent at each location. Residential history data were available for MADRES participants; Healthy Start exposures were based on the address at the time of enrollment.

### B1. Criteria Pollutants and Meteorology
For each residential location, daily PM2.5, O3, and NO2 concentrations were estimated using inverse distance weighting (IDW) with a power of 2. PM2.5 and NO2 data were summarized as daily mean concentrations. O3 data were summarized as daily means and daily 8-hour max concentrations. Monitors within 50 km were included in the calculations. IDW estimates for the Healthy Start cohort were generated using the gstat package in R (https://github.com/r-spatial/gstat/)

Daily temperature (tmin, tmax) and relative humidity (rhmin, rhmax) were summarized at each residential location using gridded data (gridMET) from the Climatology Lab (www.climatologylab.org). For Healthy Start, data were obtained using the climateR package for R https://github.com/mikejohnson51/climateR. GridMET data were available at minimum and maximum values each day. We averaged these values to estimate the mean temperature and relative humidity for each location each day. 

Daily air pollutant and meteorology variables were summarized as pregnancy averages (average of all measurements from the estimated date of conception through the date of birth), as trimester specific averages, and as averages for the first year of life. We defined the estimated date of conception as the date of birth minus the gestational age (in days). Trimester 1 was defined as the date of conception to day 91 of gestation, trimester 2 was defined as day 92 to day 182 of gestation, and trimester 3 was defined as day 183 of gestation to the date of birth. Postnatal exposures were calculated as weekly averages based on the age of the child in weeks.

### B2. Traffic-related air pollution
We used two different approaches to estimate exposure to traffic-related air pollutants (TRAP). First, we used distance to major roadways as a proxy measure of TRAP exposure. For Healthy Start, we calculated the minimum distance to a major road using the US Census Bureau TIGER/Line shapefile for primary and secondary roads (https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2021&layergroup=Roads)[4]. For MADRES, we calculated the minimum distance to a class 1 or class 2 road using the Streetlytics model from Citilabs. 

Second, we estimated TRAP exposures using estimates from spatiotemporal prediction models developed for each cohort. For Healthy Start, we based TRAP exposures on a prediction model for black carbon (BC). Details on this model are available elsewhere (Martenies et al., 2021). Briefly, we collected weekly PM2.5 samples in using low-cost monitors and measured BC mass using transmissometry [6,7]. Our model was developed using the ‘SpatioTemporal’ package in R [8] and featured 13 spatial covariates, two spatiotemporal covariates, and two temporal basis functions were derived from regional NO2, PM2.5, and BC monitoring data. The prediction model generated weekly BC concentrations at Healthy Start residential locations. For MADRES, we based TRAP exposures on a prediction model for oxides of nitrogen (NOx). Details on this model are available elsewhere [9,10]. Briefly, the model uses a flexible, three-stage hierarchical framework with spatiotemporally referenced covariates and data from both long-term routine monitoring stations with high temporal resolution and short-term, sporadic measurement campaigns with high spatial resolution. Temporal basis functions are fit to capture seasonality and longer term temporal variation [11]. The second stage of the model utilizes ensemble learning (bootstrap aggregation) to increase the precision of estimates and output their standard deviation, while the third constrained optimization stage adjusts the parameter estimates based on real-life physical and chemical properties and trends. The prediction model generated weekly NOx concentrations at MADRES residential locations. Weekly TRAP exposures (BC or NOx) were averaged across the pregnancy and by trimester for each participant.

### B3. Green Space (NDVI)
Normalized difference vegetation index (NDVI) data (250 m resolution) were summarized for each residential location at different buffer distances. NDVI data are available every 16 days. We obtained all NDVI raster files for each year in our study period. Raster files were averaged to generate an annual NDVI estimate. Then, the average NDVI value within a given buffer were calculated for each residential location. 

### B4. Built Environment: Parks
We calculated three variables related to park space: distance to the closest park; total park area within a buffer, and number of parks within a buffer distance. Distance to the closest park (meters) was defined as the minimum distance between the residential location (point) and the boundary of any nearby park. Total park area (meters squared) was defined as the sum of park area that fell within the specified buffer distances. The number of parks was defined as the sum of park polygons that intersected with the buffer.

### B5. Built Environment: Impervious surface and tree canopy
We used two variables from the National Land Cover Database (NLCD): percent impervious surface and percent tree canopy. NLCD data are available at a 50 m resolution for the conterminous United States. Mean percent impervious surface and mean percent tree canopy in a 500 m buffer were calculated for Healthy Start and MADRES participants using NLCD data from 2011 and 2016, respectively. 

### B6. Neighborhood Socioeconomic Status
Three variables were used as measures of neighborhood socioeconomic status. First, we used the Food Access Research Atlas (FARA) to determine whether participants lived in a low income and low food access (LILA) census tract. Data from the 2015 and 20XX releases were used to assign values to Healthy Start and MADRES participants, respectively. We used the FARA variable that defined low income and low access using the half mile designation for urban areas and 10 miles for nonurban areas (variable name = LILATracts_halfAnd10).

Then, we derived two socioeconomic indicators for each census tract using data from the American Community Survey [12]: the percentage of the population aged 25 years or older without at least a high school diploma or equivalent (low educational attainment) and the percentage of households with past year income below the poverty level (households in poverty). We used the following variables from the ACS data set:

- B15003_001    Population 25 years and older
- B15003_002    No schooling
- B15003_003    Nursery school
- B15003_004    K
- B15003_005    1st grade
- B15003_006    2nd grade
- B15003_007    3rd grade
- B15003_008    4th grade
- B15003_009    5th grade
- B15003_010    6th grade
- B15003_011    7th grade
- B15003_012    8th grade
- B15003_013    9th grade
- B15003_014    10th grade
- B15003_015    11th grade
- B15003_016   12th grade, no diploma
- B17017_001   Total households
- B17017_002   Households with come in the past 12 months below poverty level

## References
1. Yang Y, Diez-Roux AV. Walking Distance by Trip Purpose and Population Subgroups. Am J Prev Med. 2012;43:11–9. 

2. Kim Y, Chung, Eunhee. Descriptive Epidemiology of Objectively Measured Walking Among US Pregnant Women: National Health and Nutrition Examination Survey, 2005–2006. Prev Chronic Dis [Internet]. 2015 [cited 2023 Jun 6];12. Available from: https://www.cdc.gov/pcd/issues/2015/15_0437.htm

3. Connolly CP, Conger SA, Montoye AHK, Marshall MR, Schlaff RA, Badon SE, et al. Walking for health during pregnancy: A literature review and considerations for future research. J Sport Health Sci. 2019;8:401–11. 

4. US Census Bureau. TIGER/Line Shapefiles [Internet]. Census.gov. 2022 [cited 2023 Jul 3]. Available from: https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html

5. Martenies SE, Keller JP, WeMott S, Kuiper G, Ross Z, Allshouse WB, et al. A Spatiotemporal Prediction Model for Black Carbon in the Denver Metropolitan Area, 2009–2020. Environ Sci Technol. 2021;55:3112–23. 

6. Ahmed T, Dutkiewicz VA, Shareef A, Tuncel G, Tuncel S, Husain L. Measurement of black carbon (BC) by an optical method and a thermal-optical method: Intercomparison for four sites. Atmospheric Environment. 2009;43:6305–11. 

7. Presler-Jur P, Doraiswamy P, Hammond O, Rice J. An evaluation of mass absorption cross-section for optical carbon analysis on Teflon filter media. Journal of the Air & Waste Management Association. 2017;67:1213–28. 

8. Lindstrom J, Szpiro A, Sampson PD, Bergen S, Oron AP. SpatioTemporal: Spatio-Temporal Model Estimation. R package version 1.1.9.1. [Internet]. 2019. Available from: https://CRAN.R-project.org/package=SpatioTemporal

9. Li L, Girguis M, Lurmann F, Wu J, Urman R, Rappaport E, et al. Cluster-based bagging of constrained mixed-effects models for high spatiotemporal resolution nitrogen oxides prediction over large regions. Environ Int. 2019;128:310–23. 

10. Li L, Lurmann F, Habre R, Urman R, Rappaport E, Ritz B, et al. Constrained Mixed-Effect Models with Ensemble Learning for Prediction of Nitrogen Oxides Concentrations at High Spatiotemporal Resolution. Environ Sci Technol. 2017;51:9920–9. 

11. Szpiro AA, Sampson PD, Sheppard L, Lumley T, Adar SD, Kaufman JD. Predicting intra-urban variation in air pollution concentrations with complex spatio-temporal dependencies. Environmetrics. 2010;21:606–31. 
12. US Census Bureau. 2010-2014 American Community Survey (ACS) 5-year Estimates [Internet]. 2014 [cited 2016 Oct 6]. Available from: https://www.census.gov/programs-surveys/acs/

