
initial: 2-8-2022 - Zach Brodeur

Overview: Contains data and code to fit a multivariate statistical model to HEFS ensemble forecast data of Lake Mendocino inflow for period of 10-01-1985 to 09-30-2010. Raw HEFS data are output at hourly forecast timesteps out to 14 days lead, while model framework condenses output to daily timescales to enable fitting between observations at daily scale and forecasts. Inputs are sequential .csv files produced by CNRFC and output is a m x e x n x l array, where m is the number of synthetic ensembles, e is the ensemble size (61 members in this case), n is the length of the observational data, and l is the number of leads.

Repository - 'hec-wat_ensemble':

	Sub directory - 'data': 
		
		Sub repository - 'hefs_lamc_act-meteo': Daily output from CNRFC HEFS hindcasts in .csv format, forecasts are hourly to 14 day lead
		'LAMC_local': .csv containing daily observed flow
		**'raw_data_process' script outputs various files to this sub-repository in R data structure (rds) format. Most importantly is the 		'lamc_hefs_ens_forc.rds' file which matches each daily observation with the forecasts at different lead times for model fitting
	
	Sub directory - 'fit':
		
		**'lamc_fit-model_cmean' script outputs fitted parameters and timeseries to this directory

	Sub directory - 'out':

		**'lamc_synthetic-gen_cmean' and 'syn_hefs_out' output synthetic ensembles and required metadata to this directory

	Sub directory - 'common':

	'GL_maineqs': Primary equations for generalized likelihood (GL) function used to fit distributions
	'GL_subeqs': Required sub-equations for GL functionality

	Sub directory - 'fitting_process':

	'raw_data_process': Converts daily .csv data files to R array for model fitting; outputs to 'data' sub-directory
	'lamc_fit-model_cmean': Fits a statistical model by month and lead time to the observations and forecast data by ensemble; output to 'fit' 	sub-directory
	'lamc_synthetic-gen_cmean': Generates synthetic ensembles of same dimension as HEFS input; output to 'out' sub-directory
	'syn_hefs_out': Reorganizes ensembles to original forward looking format and outputs to both .xlsx and .feather formats
		*Note: .xlsx file is a single workbook for each synthetic ensemble, .feather saves individual ensemble dataframes for each 		synthetic generation run in a single sub-directory

	Sub directory - 'output_process':
	'diagnostics.R': ggplot script to show synthetic ensemble skill against HEFS skill
	'fcst_reformatters.R': transforms forecast matrix into tidy array for outputs
	'ens_testing': A short script to plot results and verify normal operation

	Sub directory - 'wat_files'	
    This directory can be copied wholesale into the WAT watershed.
  'scripts' and 'scripting': 
    folders containing Jython script to run R process from WAT, WAT alternative file to launch script
  'syntheticFcsts' directory:
  	'wat_synthetic-gen_cmean.R': synthetics generation for ADOC, generalized from LAMC version.
  	'syn_hefs_out_tsensembles.R': writes forecasts to various file types including SQLite
  	'wat_launcher.r': launches process for WAT and manages inputs
  	'forecast_config.json': configuration parameters needed to run script in WAT.


	
		