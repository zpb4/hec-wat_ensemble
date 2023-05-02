# wat_launcher.R
# Evan Heisman, USACE IWR-HEC
# 
# Runs R from a WAT jython scripting plugin
#

scriptArgs = commandArgs(trailingOnly=TRUE)
rDir = scriptArgs[1]


.libPaths(file.path(rDir, "library"))
Sys.setenv(TEMP="C:\\temp\\fcst_gen")

# Functions to parses a configuration file from WAT's "RunRCmd" script and
# unpack into environment if needed
require(rjson)
require(stringr)

# just a check that we got here from WAT
#require(tcltk)
#msgBox <- tkmessageBox(title = "WAT Compute",
#                       message = "Hello, WAT World!", icon = "info", type = "ok")

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
eventConfig = parseConfigFile(scriptArgs[2])
# set synthetic flows file to generate forecasts for
syntheticFlowFile = scriptArgs[3]
print(syntheticFlowFile)
# read config file for script - this should specify a few folder names and other variables
# comment out defintions in in other scripts and add to this file
scriptConfig = parseConfigFile(paste0(eventConfig$Outputs$`Watershed Directory`, "synForecasts\\forecastConfig.json"))
# use this to unpack config variables that we turn off in wat_synthetics-...
setVars(scriptConfig$forecast_generator_config)

# seed needs to vary by integer - WAT's event random number isn't sufficient
set.seed(as.integer(eventConfig$Indices$`Event Number` * eventConfig$Indices$`Lifecycle Number`))

# Set directory for output
outputDir = str_replace_all(eventConfig$Outputs$`Run Directory`, fixed("\\Scripting"), "") # one level up, but being lazy

# Finally source WAT forecast generator
# Would like to do this all in the watershed at some point
#setwd("C:\\projects\\Prado_WAT_FIRO_Dev\\hec-wat_ensemble")
# OR
setwd(paste0(eventConfig$Outputs$`Watershed Directory`, scriptDir))

use_observed_flows = F # use obs dataset?  F for generating WAT events, T for historical data
source("synthetic-gen_v2_parallel.R")
#source("wat_synthetics-gen_cmean.R")

# create plots?
source(".\\output_process\\fcst_reformatters.R")
#plotDir = paste(eventConfig$Outputs$`Watershed Directory`, "synForecasts", sep="\\")
#source("diagnostics.R")

# write out to sqlite file?
source("syn-hefs_out_tsensembles.R")

#remove variables and clean environment
#if(scriptConfig$useGC){
#  rm(list=ls());gc()
#}