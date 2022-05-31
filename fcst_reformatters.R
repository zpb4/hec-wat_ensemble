# fcst_reformatters.R
# reformats forecast matrices into tidy formats

require(reshape2)
require(plyr)

# function to convert forecasts to useful formats
meltForecasts <- function(fcstMatrix){
  melt(fcstMatrix, varnames=c("ens_num", "day", "lead"), value.name="flow")
}

# compute percentiles
summarizeForecasts <- function(fcstDF, fcstLead=3, p=c(0, 0.05, 0.5, 0.95, 1), dropMembers=c()){
  #fcstDF$day = fcstDF$day+fcstLead
  keepMembers = unique(fcstDF$ens_num)
  keepMembers = keepMembers[which(!(keepMembers %in% dropMembers))]
  ddply(subset(fcstDF, lead==fcstLead & ens_num %in% keepMembers), .(day), 
        function(df){
          data.frame(percentiles=paste0("p", p), flow=quantile(df$flow, p))
        })
}

# tidy percentiles
flipForecastSummary <- function(sumFcst) {
  fcstDF = dcast(sumFcst, day ~ percentiles, value.var="flow")
  fcstDF$day = as.Date(ix_sim)
  fcstDF
}

