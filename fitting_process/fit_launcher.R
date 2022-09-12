require(rjson)
require(stringr)

# scriptArgs = commandArgs(trailingOnly=TRUE)

# function to read JSON config file
parseConfigFile <- function(configFileName){
  jsonObj = fromJSON(file=configFileName)
  return(jsonObj)
}

# function to unpacks json values into namespace
setVars <- function(jsonConfig){
  for(i in 1:length(jsonConfig)){
    assign(names(jsonConfig[i]),jsonConfig[i],env=parent.frame(n=1))
  }
}

fitConfig = parseConfigFile(".\\wat_files\\synForecasts\\forecastConfig.json")$locations[[1]]

# TODO: for each location do this

# read these from fit config
location_name = fitConfig$name
# dimensions
ens_num = fitConfig$ens_num #no of ensemble members
leads = fitConfig$leads #total number of lead times
ar = fitConfig$ar #lags for VAR model
# HEFS params
hefs_folder = fitConfig$hefs_folder
skip_rows = fitConfig$skip_rows
# observed data
obs_csv_file = fitConfig$obs_csv_file
# start/end of fitting process, must exist for both
fit_start = fitConfig$start_date
fit_end = fitConfig$end_date
# convert start/end to dates
fit_start_date = as.Date(fit_start)
fit_end_date = as.Date(fit_end)

source("raw_data_process.R")
source("generalized_init-fit-model_v2.R")
source("generalized_fit-model_v2.R")