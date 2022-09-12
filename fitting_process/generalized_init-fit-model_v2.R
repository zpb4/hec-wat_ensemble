#Script to fit conditional mean, VAR, and GL model to forecast residuals and create
#Schaake shuffled empirical copula to sample from

#setwd('h:/firo_lamc/hec-wat_ensemble//')0

#ens_num <- 68 #no of ensemble members
#leads <- 14 #total number of lead times
#ar<-3 #lags for VAR model

#1) Read in observed and simulated inflows
inf<-read.csv('data/ADOC_inflow.csv')
colnames(inf) <- c("GMT", "ADOC")
obs<-inf[which(inf$GMT==fit_start):which(inf$GMT==fit_end),2] # col 2 is obsv ADOC flow


#set negative fnf observations to zero
obs[obs<0]<-0

#define matrices of observations that match forecast output for residual calculations
#e.g. repeat observations across 14 forecast leads, so this is a (n x 14) matrix
obs_mat<-cbind(matrix(rep(obs,leads),ncol=leads))

#date indices, fit is calculated starting 15 days after beginning of data since full forecast data for all leads required
ix<-seq(fit_start_date, fit_end_date,'day')
ix2<-as.POSIXlt(ix)
ix3<-seq(fit_start_date, fit_end_date,'day')
hefs_idx<-which(ix3==fit_start_date):which(ix3==fit_end_date)

#forecast matrices
fcst_hefs<-readRDS(paste0('data/', location_name, '_hefs_ens_forc.rds'))

#scale by 1000 to convert from kcfs to cfs
fcst_hefs<-fcst_hefs[,hefs_idx,]*1000

#for loop to do across all 61 ensembles, takes a while to fit all samples so 
#might want to comment this out to just play around with a single ensemble
print(Sys.time()) #start time 

obs_mat<-array(rep(obs),c(length(hefs_idx),ens_num))
rresids<-array(NA,c(ens_num,length(ix),leads))

#Fit a LOESS model to ensemble median as conditional mean estimate
loess_fit<-vector('list',leads)
loess_fit_sub<-vector('list',12)

for(i in 1:leads){
  loess_fit[[i]]<-loess_fit_sub
  for(j in 1:12){
    seas<-which(ix2$mon==(j-1))
    loess_fit[[i]][[j]]<-loess(apply(fcst_hefs[,seas,i],2,median)~obs_mat[seas,1],span=1,degree = 2,family='symmetric',control=loess.control(surface='direct'))
  }
}

saveRDS(loess_fit,paste0('fit/', location_name, '_loess_fit_v2.rds'))

#save an array of the LOESS estimated conditional mean
loess_mat<-array(NA,c(length(ix),leads))

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  obs_inf<-matrix(rep(obs[seas],leads),ncol=leads,byrow=F)
  loess_inf<-array(NA,dim(obs_inf))
  for(j in 1:leads){
    ob_pred<-obs_inf[,j]
    cmn<-predict(loess_fit[[j]][[i]],ob_pred)
    zero_idx<-which(cmn<0)
    cmn[zero_idx]<-ob_pred[zero_idx]
    loess_inf[,j]<-cmn
  }
  loess_mat[seas,]<-loess_inf
}

saveRDS(loess_mat,paste0('fit/', location_name, '_loess_mat_v2.rds'))

#store maximum forecast values in historical period for wild hair processing
max_val<-array(NA,c(12,leads))

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  max_val[i,]<-apply(apply(fcst_hefs,c(2,3),max)[seas,],2,function(x){max(x,na.rm=T)})
}

saveRDS(max_val,paste0('fit/', location_name, '_max_val_v2.rds'))

#calculate raw residuals
for(e in 1:ens_num){
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    rresid_mat<-matrix(ncol=leads,nrow=length(seas))
    sim_inf<-fcst_hefs[e,,]
    loess_inf<-loess_mat[seas,]
    rresid_mat <- loess_inf - sim_inf[seas,]
    rresids[e,seas,]<-rresid_mat
  }

}

saveRDS(rresids,paste0('fit/', location_name, '_rresids_v2.rds'))

#remove variables and clean environment
rm(list=ls());gc()

############################################END#########################################################
