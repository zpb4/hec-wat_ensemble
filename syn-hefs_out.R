#Script to convert files to 'normal' forward looking ensemble format and then store as 'feather' file

setwd('h:/firo_lamc/ensemble-will/')
library(feather)
library(openxlsx)

syn_forc<-readRDS('out/syn_hefs_flow_cm.rds')
inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),5]
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

ens_num<-61
n<-10
leads<-14

syn_lam_out<-array(NA,c(n,ens_num,dim(syn_forc)[3],(dim(syn_forc)[4]+1)))

for(m in 1:n){
  syn_lam<-syn_forc[m,,,]
  
  syn_lam_out[m,,,1]<-matrix(rep(obs,ens_num),nrow=ens_num,byrow=T)

  for(i in 1:leads){
    syn_lam_out[m,,,(i+1)]<-syn_lam[,c((i+1):dim(syn_lam)[2],rep(dim(syn_lam)[2],(i))),i]
  }
}

saveRDS(syn_lam_out,'out/syn_lamc_out.rds')

#output to single excel file
for(m in 1:n){
  wb<-createWorkbook()
  
  for(e in 1:ens_num){
    syn_lam_xl<-as.data.frame(syn_lam_out[m,e,,])
    colnames(syn_lam_xl)<-paste(0:14,'d')
    rownames(syn_lam_xl)<-as.Date(ix)
    addWorksheet(wb, sheetName= paste('Ens',e))
    writeData(wb,paste('Ens',e),syn_lam_xl,rowNames = T)
  }
  saveWorkbook(wb, file = paste('out/excel/syn-ensemble_',m,'.xlsx',sep=''), overwrite = TRUE)
}

#output to separate feather files
for(m in 1:n){
  if(dir.exists(paste('out/feather/syn-ens_',m,sep=''))==T){
    unlink(paste('out/feather/syn-ens_',m,sep=''),recursive = T)}
  dir.create(paste('out/feather/syn-ens_',m,sep=''))
  for(e in 1:ens_num){
    syn_lam_xl<-as.data.frame(syn_lam_out[m,e,,])
    colnames(syn_lam_xl)<-paste(0:14,'d')
    rownames(syn_lam_xl)<-as.Date(ix)
    write_feather(syn_lam_xl, paste('out/feather/syn-ens_',m,'/ens-mbr_',e,'.feather',sep=''))
  }
}


rm(list=ls());gc()

#############################################END############################################

