# wat_helper.R
# Evan Heisman, USACE IWR-HEC
# 
# Functions to parses a configuration file from WAT's "RunRCmd" script and
# unpack into environment if needed
require(argparser)
require(rjson)

readJsonConfig(file) <- function(configFileName){
  jsonObj = fromJSON(file=configFileName)
  return(jsonObj)
  
}

setVars <- function(jsonConfig){
  
  # for object in json's rvars list:
  #   assign(k,v,env=parent.frame(n=1))
  
  return(jsonObj)
}
