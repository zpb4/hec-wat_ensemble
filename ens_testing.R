##Testing Script for ensemble output

setwd('h:/firo_lamc/ensemble-will/')

syn_lam <- readRDS("out/ADOC1_syn_lamc_out.rds")

syn_lam_out <- syn_lam[1,,,]

ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

roll_sum<-function(x){out<-c();out[1]<-x[1];
for(i in 1:(length(x)-1)){out[i+1]<-out[i]+x[i+1]};return(out)}
u90<-function(x){y<-sort(x);return(y[55])}
l10<-function(x){y<-sort(x);return(y[6])}

ylm1<-c(0,5)
ylm2<-c(50,50,50,40,20,5,5,5,5,5,5,50)

cfs_taf<-2.29568411*10**-5 * 86400 / 1000

ens_num<-61

#--------------------------------------------------------------------------------------------------------
#Combined Plots
yrs<-1986:2009
dy<-'01'
mth<-c('01','02','03','04','05','06','07','08','09','10','11','12')

par(mfrow=c(4,3))
for(k in 1:12){
  for(j in 1:length(yrs)){
    idx<-which(ix==paste(yrs[j],mth[k],dy,sep='-'))
    plot(0:14,roll_sum(syn_lam_out[1,idx:(idx+14),1]) * cfs_taf,type='l',lwd=3,xlab='Days',
         ylab='Cumulative Flow (TAF)',main=ix[idx],ylim=c(0,ylm2[k]))
    for(i in 1:ens_num){
      lines(1:15,roll_sum(syn_lam_out[i,idx,]) * cfs_taf,col='light gray')
    }
    mn<-apply(syn_lam_out[,idx,],1,roll_sum)
    lines(1:15,apply(mn,1,mean) * cfs_taf,col='green',lwd=3,lty=2)
    lines(1:15,apply(mn,1,median) * cfs_taf,col='orange',lwd=3,lty=3)
    lines(1:15,roll_sum(syn_lam_out[1,idx:(idx+14),1]) * cfs_taf,lwd=3)
    lines(1:15,apply(mn,1,u90) * cfs_taf,lwd=3,lty=4,col='dark gray')
    lines(1:15,apply(mn,1,l10) * cfs_taf,lwd=3,lty=4,col='dark gray')
  }
}

#------------------------------------------------------------------------------------------------------------

#4 panel plot
idx<-which(ix=='1986-02-18')

ylm1<-c(0,30)
ylm2<-c(0,110)

par(mfrow=c(1,4))
for(i in c(10,5,3,1)){
  plot(1:15,syn_lam_out[1,(idx-i):(idx+14-i),1] * cfs_taf,type='l',lwd=3,ylim=ylm1,xlab='Days',
       ylab='Flow (TAF/d)',main=paste(ix[idx-i],i,'day lead'))
  for(j in 1:ens_num){
    lines(1:15,syn_lam_out[j,idx-i,] * cfs_taf,col='light gray')
  }
  lines(1:15,apply(syn_lam_out[,idx-i,],2,mean) * cfs_taf,col='green',lwd=3,lty=2)
  lines(1:15,apply(syn_lam_out[,idx-i,],2,median) * cfs_taf,col='orange',lwd=3,lty=3)
  lines(1:15,syn_lam_out[1,(idx-i):(idx+14-i),1] * cfs_taf,lwd=3)
  legend('topright',c('Obs','syn-HEFS','syn-HEFS mean','syn-HEFS median'),
         col=c('black','light gray','green','orange'),lwd=c(1,1,3,3),lty=c(1,1,2,3))
}

for(i in c(10,5,3,1)){
  plot(1:15,roll_sum(syn_lam_out[1,(idx-i):(idx+14-i),1]) * cfs_taf,type='l',lwd=3,xlab='Days',
       ylab='Cumulative Flow (TAF)',main=paste(ix[idx-i],i,'day lead'),ylim=ylm2)
  for(j in 1:ens_num){
    lines(1:15,roll_sum(syn_lam_out[j,idx-i,]) * cfs_taf,col='light gray')
  }
  mn<-apply(syn_lam_out[,idx-i,],1,roll_sum)
  lines(1:15,apply(mn,1,mean) * cfs_taf,col='green',lwd=3,lty=2)
  lines(1:15,apply(mn,1,median) * cfs_taf,col='orange',lwd=3,lty=3)
  lines(1:15,roll_sum(syn_lam_out[1,(idx-i):(idx+14-i),1]) * cfs_taf,lwd=3)
  lines(1:15,apply(mn,1,u90) * cfs_taf,lwd=3,lty=4,col='dark gray')
  lines(1:15,apply(mn,1,l10) * cfs_taf,lwd=3,lty=4,col='dark gray')
  legend('topleft',c('Obs','syn-HEFS','syn-HEFS mean','syn-HEFS median','90/10 pcnt'),
         col=c('black','light gray','green','orange','dark gray'),lwd=c(1,1,3,3,3),lty=c(1,1,2,3,4))
}


#---------------------------------------------------------------------------------------------