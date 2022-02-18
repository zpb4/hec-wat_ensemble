#Script to process csv HEFS data from RFC
#lamc,ukac,hopc,hopcL data
setwd('h:/firo_lamc/ensemble-will/data/')
library(stringr)
library(abind)

#define folder, date/time index, and data dimensions
folder<-'hefs_lamc_act-meteo/' #enter lamc, ukac, hopc, hopcL
loc<-'LAMC1'
ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day') #same length
ixx3<-as.POSIXlt(ix3)
leads<-14
ens_num<-61

#
i<-1
fname<-paste(folder,ixx3[i]$year+1900,str_pad(ixx3[i]$mon+1,2,'left','0'),str_pad(ixx3[i]$mday,2,'left','0'),'12_',loc,'_SQIN_hourly.csv',sep='')
dat1 <- read.csv(fname)
dat1[2:337,2:62]<-NA

m<-length(ix3)
hefs_raw<-array(NA,c(m,336,61))

for(i in 1:m){
  fname<-paste(folder,ixx3[i]$year+1900,str_pad(ixx3[i]$mon+1,2,'left','0'),str_pad(ixx3[i]$mday,2,'left','0'),'12_',loc,'_SQIN_hourly.csv',sep='')
  dat <- tryCatch(read.csv(fname),error = function (e){print(paste(fname,'no data'));return(dat1)},warning = function (w){print(paste(fname,'no data'));return(dat1)})
  lst_dat<-dat[2:337,2:62]
  for(j in 1:ens_num){
    out<-levels(lst_dat[,j])[lst_dat[,j]]
    out2<-as.numeric(out)
    hefs_raw[i,,j]<-as.numeric(out2)
  }
}

saveRDS(hefs_raw,'lamc_hefs_raw.rds')

#------------------------------------------------------------------------------------------------
#process raw hourly hefs data into daily mean array

idx<-cbind(seq(1,336,24),seq(24,336,24))

hefs_dly_mean<-array(NA,c(m,leads,ens_num))

for(i in 1:leads){
  mn<-apply(hefs_raw[,idx[i,1]:idx[i,2],],c(1,3),function(x){mean(x,na.rm=T)})
  hefs_dly_mean[,i,]<-mn
}

saveRDS(hefs_dly_mean,'lamc_hefs_dly_mean.rds')

#----------------------------------------------------------------------------------------------------
#Rearrange HEFS to line up observations with forecast

hefs_ens_forc<-array(0,c(ens_num,m,leads))

for(i in 1:14){
  int<-array(0,c(61,m))
  int[,(i+1):m]<-t(hefs_dly_mean[,i,])[,-c((m-i+1):m)]
  hefs_ens_forc[,,i]<-int
}

#replace any missing data with linear interpolated data
na_fun<-function(x){for(i in 1:length(x)){
  if(is.na(x[i])==T){x[i]<-x[i-1]+((x[i+1]-x[i-1])/2)}}
  return(x)}

for(i in 1:ens_num){
  hefs_ens_forc[i,,]<-apply(hefs_ens_forc[i,,],2,na_fun)
}

saveRDS(hefs_ens_forc,'lamc_hefs_ens_forc.rds')


#calculate ensemble mean and save in array as needed
hefs_ens_mean_forc<-apply(hefs_ens_forc,c(2,3),mean)

saveRDS(hefs_ens_mean_forc,'lamc_hefs_ens_mean_forc.rds')

#remove variables and clean environment
rm(list=ls());gc()


##########################################END#################################################################
