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

# read these from fit config, not just unpacking them so I know what to refactor later. ;)
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
fit_start = fitConfig$fit_start
fit_end = fitConfig$fit_end
# convert start/end to dates
fit_start_date = as.Date(fit_start, format="%m/%d/%Y %H:%M")
fit_end_date = as.Date(fit_end, format="%m/%d/%Y %H:%M")

#source("./fitting_process/raw_data_process.R")
#source("./fitting_process/generalized_init-fit-model_v2.R")
source("./fitting_process/generalized_fit-model_v2.R")
