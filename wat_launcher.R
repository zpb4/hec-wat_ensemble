# wat_launcher.R
# Evan Heisman, USACE IWR-HEC
# 
# Runs R from a WAT jython scripting plugin
#
# Functions to parses a configuration file from WAT's "RunRCmd" script and
# unpack into environment if needed
require(rjson)
require(tcltk)

# just a check that we got here
#msgBox <- tkmessageBox(title = "WAT Compute",
#                       message = "Hello, WAT World!", icon = "info", type = "ok")

scriptArgs = commandArgs(trailingOnly=TRUE)

# Don't do this except to check that the args are showing up
#s = paste(scriptArgs, sep="\n")
#msgBox <- tkmessageBox(title = "WAT Compute",
#                       message = s, icon = "info", type = "ok")

parseConfigFile <- function(configFileName){
  jsonObj = fromJSON(file=configFileName)
  return(jsonObj)
}
config = parseConfigFile(scriptArgs[1])
syntheticFlowFile = scriptArgs[2]

# print(config)

# seed needs to vary by integer - WAT's event random number isn't sufficient
set.seed(config$Seeds$`Event Seed` * config$Indices$`Lifecycle Number`)

# Set directory for output
outputDir = paste0(config$Outputs$`Run Directory`, "..\\") # one level up, being lazy

# Finally source WAT forecast generator
# Would like to do this all in the watershed at some point
setwd("C:\\projects\\Prado_WAT_FIRO_Dev\\hec-wat_ensemble")
source("wat_synthetics-gen_cmean.R")