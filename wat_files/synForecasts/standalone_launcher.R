# wat_launcher.R
# Evan Heisman, USACE IWR-HEC
# 
# Runs R from a WAT jython scripting plugin
#
.libPaths("C:\\programs\\r\\win-library\\4.1")
Sys.setenv(TEMP="C:\\temp\\fcst_gen")
# Functions to parses a configuration file from WAT's "RunRCmd" script and
# unpack into environment if needed
require(rjson)
require(stringr)

# just a check that we got here from WAT
#require(tcltk)
#msgBox <- tkmessageBox(title = "WAT Compute",
#                       message = "Hello, WAT World!", icon = "info", type = "ok")

scriptArgs = commandArgs(trailingOnly=TRUE)

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


# read metadata inputs from WAT
#eventConfig = parseConfigFile(scriptArgs[1])
# set synthetic flows file to generate forecasts for
syntheticFlowFile = scriptArgs[2]
print(syntheticFlowFile)
# read config file for script - this should specify a few folder names and other variables
# comment out defintions in in other scripts and add to this file
scriptConfig = parseConfigFile(".\\wat_files\\synForecasts\\forecastConfig.json")
eventConfig = scriptConfig$standalone_generator
# use this to unpack config variables that we turn off in wat_synthetics-...
# let's try turning this off
#setVars(scriptConfig$forecast_generator_config)

# seed needs to vary by integer - WAT's event random number isn't sufficient
set.seed(as.integer(eventConfig$seed))

# Set directory for output
outputDir = "fcst_gen_out" #str_replace_all(eventConfig$Outputs$`Run Directory`, fixed("\\Scripting"), "") # one level up, but being lazy

# Finally source WAT forecast generator
# Would like to do this all in the watershed at some point
setwd("C:\\projects\\Prado_WAT_FIRO_Dev\\hec-wat_ensemble")
# OR
#setwd(paste0(eventConfig$Outputs$`Watershed Directory`, scriptDir))

use_observed_flows = T # use obs dataset?  F for generating WAT events, T for historical data
source(".\\wat_files\\synForecasts\\synthetic-gen_v2_parallel.R")
#source(".\wat_files\synForecasts\wat_synthetics-gen_cmean.R")

# create plots
source(".\\common\\fcst_reformatters.R")
plotDir = "fcst_gen_out"
source(".\\output_process\\diagnostics.R")

# write out to sqlite file?
#source(".\\wat_files\\synForecasts\\syn-hefs_out_tsensembles.R")
