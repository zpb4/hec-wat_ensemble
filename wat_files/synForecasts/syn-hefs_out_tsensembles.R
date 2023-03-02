#Script to convert files to 'normal' forward looking ensemble format and then store as 'feather' file

#setwd('h:/firo_lamc/ensemble-will/')


#syn_forc<-readRDS('out/syn_hefs_flow_cm.rds')
#inf<-read.csv('data/LAMC_local.csv')
#obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),5]
#ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
#ix2<-as.POSIXlt(ix)

# dims required
#ens_num<-68
n<-1
#leads<-14
syn_forc = synflow_out #syn_hefs_flow


# returns syn_lam_out array used by Excel and feather outputs
fcstToSynLamOut <- function(){
  # returns array of dimensions:
  # n ensemble sets
  # ens_num ensemble members
  # # of timesteps (matches ix or new_obs)
  # 15?  1+leads?
  syn_lam_out<-array(NA,c(n,ens_num,dim(syn_forc)[3],(dim(syn_forc)[4]+1)))
  
  for(m in 1:n){
    syn_lam<-syn_forc[m,,,]
    
    # fully populate with observed data
    syn_lam_out[m,,,1]<-matrix(rep(new_obs,ens_num),nrow=ens_num,byrow=T)

    # replace all but one with ensembles, leaving obs row    
    for(i in 1:leads){
      syn_lam_out[m,,,(i+1)]<-syn_lam[,c((i+1):dim(syn_lam)[2],rep(dim(syn_lam)[2],(i))),i]
    }
  }
  return(syn_lam_out)
}

fcstOutToRDS <- function(){
  syn_lam_out = fcstToSynLamOut()
  saveRDS(syn_lam_out,'out/syn_lamc_out.rds')
}


# write out to Excel file
fcstOutToExcel <- function(){
  require(openxlsx)
  syn_lam_out = fcstToSynLamOut()
  
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
}

# write out to a feather file
fcstOutToFeather <- function(ensembleFilename){
  require(feather)
  syn_lam_out = fcstToSynLamOut()
  
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
  
}

# attempt to do this with rJava!
fcstOutToEnsembleFile <- function(ensembleFilename){
  
  # getting the correct Java configuration (JRE 11)
  Sys.setenv(JAVA_HOME=scriptConfig$java_config$java_home)
  # Security settings won't let JRE connect if rJava is in user home
  require(rJava, lib.loc=scriptConfig$java_config$rjava_libloc)
  # add TSEnsembles library and dependencies to class path
  # jvmArgs should include temp directory with "-Djava.io.tmpdir=C:\\Temp"
  # pick something better - https://devblogs.microsoft.com/oldnewthing/20121031-00/?p=6203
  .jinit(classpath=scriptConfig$java_config$java_libraries, parameters=scriptConfig$java_config$jvmArgs)

  # create fresh database file, remove if exists already
  if(file.exists(ensembleFilename)){
    file.remove(ensembleFilename)
  }
  db = .jnew("hec/SqliteDatabase", ensembleFilename, J("hec.SqliteDatabase")$CREATION_MODE$CREATE_NEW)
  # create an ensemble timeseries in Java
  recordID = .jnew("hec/RecordIdentifier", "ADOC", "FLOW")
  ensembleTS = .jnew("hec.ensemble.EnsembleTimeSeries", recordID, "cfs", "inst", "synthetic")
  # populate ensemble TS with ensembles
  # get a UTC ID
  int0 = as.integer(0)
  utc_zid = .jcall("java/time/ZoneId", returnSig="Ljava/time/ZoneId;", "of", "Z")
  dur = .jcall("java/time/Duration", returnSig="Ljava/time/Duration;", "ofDays", .jlong(1))
  print(paste0("len(ixx_sim) is ", length(ixx_sim)))
  print(paste0("dim(syn_forc) is ", paste0(dim(syn_forc), collapse="x")))
  for(i in 1:length(ixx_sim)){ # simulation timesteps
    simts = ixx_sim[i]
    issueDate = .jcall("java/time/ZonedDateTime", returnSig="Ljava/time/ZonedDateTime;", "of", 
                       as.integer(year(simts)), as.integer(month(simts)), as.integer(day(simts)),
                       int0, int0, int0, int0, utc_zid)
    startDate = issueDate
    # get the slice for this issue date
    # from Ensemble.java:
    # - row represents ensemble members
    # - columns are time steps
    # TODO: do I need to transpose this, code works in either way?
    ens_array = syn_forc[,which(simts == ixx_sim),]
    # convert to java float[][]
    ens_array = .jarray(.jfloat(ens_array), dispatch=TRUE)
    # add to ensemble TS
    .jcall(ensembleTS, returnSig="V", "addEnsemble", issueDate, ens_array, startDate, dur, "cfs")
  }
  # write to file
  db$write(ensembleTS)
  db$close()
}

fcstOutToEnsembleFile(paste0(outputDir, scriptConfig$forecast_generator_config
                             $sqlFilename))
#rm(list=ls());gc()

#############################################END############################################

