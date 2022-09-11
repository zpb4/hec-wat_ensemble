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


readObs <- function(fi, obsColumn=2){
  inf = read.csv(obsCSVFile)
  obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),1] # col 2 is obsv ADOC flow
  return(obs)
}


fitConfig = parseConfigFile(".\\wat_files\\synForecasts\\forecastConfig.json")$locations[[1]]