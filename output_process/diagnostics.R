
# ensemble forecast generation diagnostics code
# this expects you to turn off environment cleaning in synthetic generation script

require(ggplot2)


# This is a silly plot, but useful to validate that I have dimensions right - don't use it.
#plotEnsemble = 1
#ensFcst = melt(t(syn_hefs_flow[1,plotEnsemble,,]), varnames=c("lead", "day"), value.name="flow", )
#ensFcst$fcst.day = ensFcst$day + ensFcst$lead
#ggplot() + theme_bw() +
#  geom_line(data=ensFcst, aes(x=fcst.day, y=flow, group=day)) + 
#  geom_line(data=obsdf, aes(x=day, y=flow), color="blue")

source("./output_process/fcst_reformatters.R")

# gather data
obsdf = data.frame(flow=new_obs, day=ix_sim)
# create tidy table of synthetics for ggploting
synFcsts = meltForecasts(synflow_out) #syn_hefs_flow[1,,,])
# create same tidy table of HEFS
#hefsFcsts = meltForecasts(hefs_mat[,which(ix2 %in% ixx_sim),])

# create summary plots
sumSynFcsts = summarizeForecasts(synFcsts, fcstLead=3)
#sumHefsFcsts = summarizeForecasts(hefsFcsts, fcstLead=3)
  
rangeSynFcsts = flipForecastSummary(sumSynFcsts)
#rangeHefsFcsts = flipForecastSummary(sumHefsFcsts)

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
eventPlotName = paste(eventConfig$Outputs$`Watershed Directory`, "synForecasts", sprintf("plots-event_%03d.pdf", eventConfig$Indices$`Event Number`), sep="\\")
print(eventPlotName)
ggsave(eventPlotName, fcstSkillPlot(rangeSynFcsts, paste0("Synthetic forecasts: ", dateRange)))

#print(fcstSkillPlot(rangeHefsFcsts, paste0("HEFS forecasts: ", dateRange)))
