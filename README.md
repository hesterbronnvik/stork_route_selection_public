# stork_route_selection_public

This repository contains R scripts to reproduce the results and figures of "Experience reduces route selection for conspecifics by the collectively migrating white stork" by Bronnvik et al. 2024, _Current Biology_ (https://doi.org/10.1016/j.cub.2024.03.052).

# Abstract
Migration can be an energetically costly behavior with strong fitness consequences in terms of mortality and reproduction. The fact that costs and benefits may change with age and experience raises the question of whether experience changes the criteria that individuals use to select their migratory routes. Here we investigate the effect of age on route selection criteria in a collectively migrating soaring bird, the white stork (_Ciconia ciconia_). We perform step selection analysis on a longitudinal data set tracking 158 white storks over up to nine years to quantify how they select their routes on the basis of the social and atmospheric environments, and to examine how this selection changes with age. We find clear ontogenetic shifts in route selection criteria. Juveniles choose routes that have good atmospheric conditions and high conspecific densities. Yet, as they gain experience storks' selection on the availability of social information reduces---after their fifth migration experienced birds also choose routes with low conspecific densities. Thus, our results suggest that as individuals age, they gradually replace information gleaned from other individuals with information gained from experience, allowing them to shift their migration timing and increasing the time scale at which they select their routes. 

# Contents:
This repository consists of the R scripts:

01_data_location.R: accesses, cleans, and segments the data that we used for our analyses.

02_step_generation.R: provides the script necessary for generating alternative steps for step-selection analysis (adapted from Nourani et al. 2021 https://github.com/mahle68/global_seascape_public/blob/main/step_generation.R)

03_uplift_annotation.R: downloads ECMWF ERA-5 data and annotates the stork locations with the convective velocity scale.

04_utilization_distribution_estimation.R: builds auto-correlated kernel density estimates for each hour of the tracking data and annotates the stork locations with utilization distribution probability density.

05_step_selection_analysis.R: uses the data generated in the previous scripts to run the step selection analyses.

06_plotting_WS.R: plots exploratory and results figures from the paper and supplement.

how_to_estimate_uplift.pdf: a complete guide to estimating the convective velocity scale, w*, using R and ECMWF Reanalysis ERA5 data.

supplemental_wind_support_analysis.pdf: runs the step selection models including wind support as a predictor, which is non-significant.

All input data are available on Movebank.org
