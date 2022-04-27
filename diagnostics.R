
# ensemble forecast generation diagnostics code
# this expects you to turn off environment cleaning in synthetic generation script

require(reshape2)
require(ggplot2)
require(plyr)

# function to convert forecasts to useful formats
meltForecasts <- function(fcstMatrix){
  melt(fcstMatrix, varnames=c("ens_num", "day", "lead"), value.name="flow")
}

# compute percentiles
summarizeForecasts <- function(fcstDF, fcstLead=3, p=c(0, 0.05, 0.5, 0.95, 1)){
  ddply(subset(fcstDF, lead==fcstLead), .(day), 
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

# This is a silly plot, but useful to validate that I have dimensions right - don't use it.
#plotEnsemble = 1
#ensFcst = melt(t(syn_hefs_flow[1,plotEnsemble,,]), varnames=c("lead", "day"), value.name="flow", )
#ensFcst$fcst.day = ensFcst$day + ensFcst$lead
#ggplot() + theme_bw() +
#  geom_line(data=ensFcst, aes(x=fcst.day, y=flow, group=day)) + 
#  geom_line(data=obsdf, aes(x=day, y=flow), color="blue")


# gather data
obsdf = data.frame(flow=new_obs, day=ix_sim)
# create tidy table of synthetics for ggploting
synFcsts = meltForecasts(syn_hefs_flow[1,,,])
# create same tidy table of HEFS
hefsFcsts = meltForecasts(hefs_mat[,which(ix2 %in% ixx_sim),])

# create summary plots
sumSynFcsts = summarizeForecasts(synFcsts)
sumHefsFcsts = summarizeForecasts(hefsFcsts)
  
rangeSynFcsts = flipForecastSummary(sumSynFcsts)
rangeHefsFcsts = flipForecastSummary(sumHefsFcsts)

fcstSkillPlot <- function(rangeFcsts, pltName){
  ggplot() + theme_bw() +
    geom_ribbon(data=rangeFcsts, aes(x=day, ymin=p0.05, ymax=p0.95), fill="lightgrey") + 
    geom_line(data=rangeFcsts, aes(x=day, y=p0.5), color="black") + scale_color_manual(labels="median") + 
    geom_line(data=rangeFcsts, aes(x=day, y=p1), color="grey", linetype='dashed') +
    geom_line(data=rangeFcsts, aes(x=day, y=p0), color="grey", linetype='dashed') +
    geom_line(data=obsdf, aes(x=day, y=flow), color="blue") + 
    scale_y_continuous(limits=c(-10,2e4)) +
    scale_x_date() + labs(title=pltName, y="flow [cfs]", x="date")
}

dateRange = paste0(ix_sim[1], " to ", tail(ix_sim,1))
print(fcstSkillPlot(rangeSynFcsts, paste0("Synthetic forecasts: ", dateRange)))
#fcstSkillPlot(rangeHefsFcsts, paste0("HEFS forecasts: ", dateRange))
