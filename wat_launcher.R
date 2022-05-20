# wat_launcher.R
# Evan Heisman, USACE IWR-HEC
# 
# Runs R from a WAT jython scripting plugin
#
# Functions to parses a configuration file from WAT's "RunRCmd" script and
# unpack into environment if needed
require(rjson)

# just a check that we got here from WAT
#require(tcltk)
#msgBox <- tkmessageBox(title = "WAT Compute",
#                       message = "Hello, WAT World!", icon = "info", type = "ok")

scriptArgs = commandArgs(trailingOnly=TRUE)

# Don't do this except to check that the args are showing up
#s = paste(scriptArgs, sep="\n")
#msgBox <- tkmessageBox(title = "WAT Compute",
#                       message = s, icon = "info", type = "ok")


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
eventConfig = parseConfigFile(scriptArgs[1])
# set synthetic flows file to generate forecasts for
syntheticFlowFile = scriptArgs[2]
print(syntheticFlowFile)
# read config file for script - this should specify a few folder names and other variables
# comment out defintions in in other scripts and add to this file
scriptConfig = parseConfigFile(paste0(eventConfig$Outputs$`Watershed Directory`, "synForecasts\\forecastConfig.json"))
# use this to unpack config variables that we turn off in wat_synthetics-...
setVars(scriptConfig)

# seed needs to vary by integer - WAT's event random number isn't sufficient
set.seed(as.integer(eventConfig$Indices$`Event Number` * eventConfig$Indices$`Lifecycle Number`))

# Set directory for output
outputDir = paste0(eventConfig$Outputs$`Run Directory`, "..\\") # one level up, but being lazy

# Finally source WAT forecast generator
# Would like to do this all in the watershed at some point
#setwd("C:\\projects\\Prado_WAT_FIRO_Dev\\hec-wat_ensemble")
# OR
setwd(paste0(eventConfig$Outputs$`Watershed Directory`, scriptDir))
source("wat_synthetics-gen_cmean.R")

# write out to sqlite file?
source("syn-hefs_out_tsensembles.R")

#remove variables and clean environment
#if(scriptConfig$useGC){
#  rm(list=ls());gc()
#}