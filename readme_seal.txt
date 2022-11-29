The pipeline for getting data to the point of a CHESS seal analysis with TMB is quite long and convoluted. Often, output from one script
is read into another. One will want to make sure that scripts point to the correct .RData files. Steps include

- In git/haulout/ringed, get_temp_grid.R produces spring_temp [average spring temperature] for 2005-2020 in each of the CHESS grid cells 
  (in Temps_cell_year.RDa).  This is required for ringed seal haulout analysis, and also for predicting haul-out for CHESS surveys

- In git/haulout/ringed, attach_snow_melt_grid.R attaches melt products (Markus et al. 2009, MDSDA) and remotely sensed snow depth data
  to the CHESS grid.  This is used in ringed seal haul-out analysis and in haul-out predictions for CHESS surveys (saved as snow_melt_grid_Feb2021.RDa)

- In git/haulout/ringed, conduct an analysis of ringed seal telemetry records using analyze_ringed_gam.R.  This produces 
  big_gam_ringed_ho_model.RData, which includes the fitted gam object needed to predict haul-out proportions for CHESS surveys

- In git/CHESS/, format_Chukchi_covs_all.R formats covariates specifically for the analysis grid, producing a list that includes covariates for each
  day of the survey.  Output is chess_grid_list_all_Feb2021.RData

- In git/CHESS/, format_effort_chukchi_large.R formats US aerial survey effort, producing CHESS_us_data_TMB_Feb2021.RData

- In git/CHESS/, format_russian_data_TMB_Feb2021.R formats Russian aerial survey effort and formats habitat covariates for the whole grid.  
  This produces CHESS_russia_dta_TMB_Feb2021.RData

- In git/CHESS, merge_data_TMB.R merges the US and Russian survey effort files into a format suitable for TMB analysis, producing CHESS_data_Feb2021.RData

- finally, run_CHESS_ho_pspline_Feb2021.R will conduct spatio-temporal abundance analysis from CHESS survey data