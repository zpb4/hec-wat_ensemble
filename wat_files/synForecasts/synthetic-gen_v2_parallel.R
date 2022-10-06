#7) Synthetic forecast generation
library(fGarch)
library(doParallel) #required package for parallel ops
library(abind) #required for parallel combining ops
library(stringr)
require(rjson)

# setwd('z:/hec-wat_ensemble/')

#parallelization set-up
parallel::detectCores()

n.cores <- parallel::detectCores()-2
my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

use_observed_flows = T # use obs dataset?
# if false, use this file
#syntheticFlowFile = "C:\\Projects\\Prado_WAT_FIRO_Dev\\Watersheds\\FIRO_Prado_Dev\\runs\\WCM_Ops\\RTestFRA\\realization 1\\lifecycle 1\\event 7\\obsTimeseries.csv"

#outputDir = "out\\" # local output
# function to read JSON config file
parseConfigFile <- function(configFileName){
  jsonObj = fromJSON(file=configFileName)
  return(jsonObj)
}
scriptConfig = parseConfigFile("./wat_Files/synForecasts/forecastConfig.json")
fcstConfig = scriptConfig$forecast_generator_config

#1. User defined inputs
#date vector for fitted data
ix<-seq(as.Date(fcstConfig$fit_start, format="%m/%d/%Y %H:%M"),
        as.Date(fcstConfig$fit_end, format="%m/%d/%Y %H:%M"),'day')
ix2<-as.POSIXlt(ix)

#1a. 'Standard' parameters that you likely won't want to change
leads <- fcstConfig$leads #daily leads, should stay as 14 for HEFS
ens_num <- fcstConfig$ens_num #no. of ensembles, model is currently fit to 61 members, so likely don't need to change
ar <- fcstConfig$ar #no. of lags in vector auto-regressive model; also don't recommend changing

#1b. Primary user defined parameters to change as desired
env_scale <- fcstConfig$env_scale #scaling of forecast envelope (wild hair removal)
no_samps <- fcstConfig$no_samps #number of ensemble sets (samples) desired

# directories and file names
locationName = fcstConfig$locationName
forecastFitDir = fcstConfig$forecastFitDir
#outputDir = "out\\"
scriptDataDir = fcstConfig$scriptDataDir

#2. Read in raw observed data
inf<-read.csv(fcstConfig$obs_csv_file)
colnames(inf)[1] <- "GMT" # fix header
obs<-inf[which(inf$GMT==fcstConfig$fit_start):which(inf$GMT==fcstConfig$fit_end),fcstConfig$obs_column]
obs[obs<0]<-0

#2a. Create matrix of observations matching forecast
obs_mat<-matrix(rep(obs,leads),ncol=leads)

#3. Define observed data matrix to create synthetic samples
if(use_observed_flows){
  
  #Define s start in year, month, and day; minimum 1948-10-01
  st_yr <- 1991 #4 digit year
  st_mo <- 10 #specify with leading zero for single digits, e.g. '01' instead of '1'
  st_dy <- 01 #specify with leading zero for single digits, e.g. '01' instead of '1'
  
  #Define simulation end in year, month, and day; maximum 2010-09-30
  end_yr <- 1999 #4 digit year
  end_mo <- 03 #specify with leading zero for single digits, e.g. '01' instead of '1'
  end_dy <- 15 #specify with leading zero for single digits, e.g. '01' instead of '1'
  
  st_date<-paste(str_remove(st_mo,'^0'),str_remove(st_dy,'^0'),st_yr,sep='/')
  end_date<-paste(str_remove(end_mo,'^0'),str_remove(end_dy,'^0'),end_yr,sep='/')
  
  new_obs<-inf[which(inf$GMT==paste(st_date,' 12:00',sep='')):which(inf$GMT==paste(end_date,' 12:00',sep='')),2]
  
} else { # WAT sampled synthetics
  new_obs_df = read.csv(syntheticFlowFile)
  new_obs = new_obs_df$Prado
  # overwrite timestamps used
  st_timestamp = as.POSIXlt(new_obs_df$GMT[1], format="%m/%d/%Y")
  st_date = str_replace(as.character(st_timestamp, format="%m/%d/%Y"), fixed(" 24:00"), "")
  st_yr = year(st_timestamp) 
  st_mo = month(st_timestamp)
  st_dy = day(st_timestamp)
  end_timestamp = as.POSIXlt(new_obs_df$GMT[nrow(new_obs_df)], format="%m/%d/%Y")
  end_date = str_replace(as.character(end_timestamp, format="%m/%d/%Y"), fixed(" 24:00"), "")
  end_yr = year(end_timestamp)
  end_mo = month(end_timestamp)
  end_dy = day(end_timestamp)
}

new_obs[new_obs<0]<-0
new_obs_mat<-matrix(rep(new_obs,leads),ncol=leads)

# Error checks on dates
err_code<-c()
if(end_yr<st_yr){stop(print('end year must be greater than or equal to start year'))} 
if(end_yr==st_yr & end_mo<st_mo){stop(print('end month must be greater than start month if same year'))}
if(end_yr==st_yr & end_mo==st_mo & end_dy<st_dy){stop(print('end day must be greater than start day if same year and month'))}

#change 1: define a 'sampling month sequence' that is flexible for shorter scenario generation
samp_mo_seq<-1:12 #baseline assumption

#change 1: if same year and not full year, ensure sample month sequence only includes desired months
if(st_yr==end_yr & length(st_mo:end_mo)<12){
  samp_mo_seq<-st_mo:end_mo
}

#change 1: if crossing into another year but no further, then define the appropriate crossover sequence of months
if((end_yr-st_yr)==1){
  samp_mo_seq<-c(st_mo:12,1:end_mo)
}


#5) Synthetic Generation
#date vector for simulated data
st<-paste(st_yr,st_mo,st_dy,sep='-')
ed<-paste(end_yr,end_mo,end_dy,sep='-')

ix_sim<-seq(as.Date(st),as.Date(ed),'day')
ixx_sim<-as.POSIXlt(ix_sim)

#this section estimates conditional mean from fitted LOESS model for syn generation
loess_fit<-readRDS(paste0(forecastFitDir, locationName, '_loess_fit_v2.rds'))

ls_sim<-array(NA,c(length(ix_sim),leads))

#each lead time has a separate LOESS model to estimate conditional mean
for(i in samp_mo_seq){
  seas_sim<-which(ixx_sim$mon==(i-1))
  nobs_inf<-matrix(rep(new_obs[seas_sim],leads),ncol=leads,byrow=F)
  loess_inf<-array(NA,dim(nobs_inf))
  for(j in 1:leads){
    ls<-predict(loess_fit[[j]][[i]],nobs_inf[,j])
    zero_idx<-which(ls<0)
    ls[zero_idx]<-nobs_inf[zero_idx,j]
    loess_inf[,j]<-ls
  }
ls_sim[seas_sim,]<-loess_inf
}

##NOTE--this section is a little wonky and hard to follow. It sets up an index to sample in a continuous
#fashion across the simulated timespan to maintain VAR model continuity...

yr_idx<-ixx_sim$year

#change 1: only include sampled months for shorter periods
yr_idx_lst<-vector('list',length(samp_mo_seq))

#change 1: only include sampled months for shorter periods
for(i in samp_mo_seq){  
  seas3<-which(ixx_sim$mon==(i-1))
  yr_idx_lst[[i]]<-yr_idx[seas3]
}

styr_idx<-st_yr-1900
edyr_idx<-end_yr-1900

yr_seq<-c(rep(styr_idx,length=(13-st_mo)),rep((styr_idx+1):(edyr_idx-1),each=12),rep((edyr_idx),length=end_mo))
mo_seq<-c(st_mo:12,rep(1:12,length((styr_idx+1):(edyr_idx-1))),1:end_mo)

#change 1: for generation of short periods in the same years
if(st_yr==end_yr & length(st_mo:end_mo)<12){
  yr_seq<-rep(styr_idx,length(samp_mo_seq))
  mo_seq<-samp_mo_seq
}

#change 1: same as above, but for crossover year situations
if((end_yr-st_yr)==1){
  yr_seq<-c(rep(styr_idx,length=(13-st_mo)),rep(edyr_idx,length=end_mo))
  mo_seq<-samp_mo_seq
}

#5b. knn set up
seas<-which(ix2$mon==0) #find length of january monthly subset of fitted data

knn<-round(sqrt(length(seas)))
tot<-sum(rep(1,knn) / 1:knn)
wts<-rep(1,knn) / 1:knn / rep(tot,knn)

#array of maximum historical forecast values to assist in wild-hair processing
max_val<-readRDS(paste0(forecastFitDir, locationName, '_max_val_v2.rds'))

#functions needed in generation
lin_mod<-function(intcpt,slope,x){
  y = intcpt + slope * x
  return(y)
}

abnd<-function(x,y){abind(x,y,along=4)}

#Synthetic Generation

#parallel for loop 
synflow_out <- foreach(m = 1:no_samps,.combine='abnd',.packages='fGarch') %dopar% {

syn_hefs_flow<-array(NA,c(leads,length(ixx_sim),ens_num))
  
knn_lst<-vector('list',length(samp_mo_seq))
  
for(i in samp_mo_seq){
  knn_vec<-c()
  seas<-which(ix2$mon==(i-1)) #define fitted data monthly index
  seas_sim<-which(ixx_sim$mon==(i-1)) #define synthetic timespan monthly index
  n_obs<-new_obs[seas_sim]
  ob<-obs[seas]
  id0 <- which(ob == 0)
  for(j in 1:length(n_obs)){
    if(n_obs[j] == 0 & length(id0)>=1) {
      s<-sample(id0,1); knn_vec[j]<-s} 
    else {
      nob_mag<-n_obs[j]
      y<-abs(nob_mag - ob) #find NEP closest by Euclidean distance
      x<-sort(y)
      x<-x[1:knn] #top 'knn' values are selected
      s<-sample(x,1,prob=wts)
      id<-which(y==s)
      if(length(id)>1) {
        id <- sample(id,1)} #resample the sample for any duplicated values
      knn_vec[j]<-id}
  }
  knn_lst[[i]]<-knn_vec
}


for(e in 1:ens_num){
  # TODO: move this to read fewer times  
  param_lst<-readRDS(paste0(forecastFitDir, locationName, '_param_lst_v2.rds'))
  #all parameters now compiled in a R list format
  norm_fit<-param_lst[[e]][[1]]
  sd_arr<-param_lst[[e]][[2]]
  var_coefs<-param_lst[[e]][[3]]
  gl_par_arr<-param_lst[[e]][[4]]
  # TODO: move this to read fewer times
  ats<-readRDS(paste0(forecastFitDir, locationName, '_ats_v2.rds'))
  
  #1) Create array of n Schaake Shuffled sequences for KNN
  syn_cop<-vector('list',length(samp_mo_seq))

  for(i in samp_mo_seq){
    seas<-which(ix2$mon==(i-1))
    mat<-ats[e,seas,]
    ecop<-apply(mat,2,function(x){rank(x,ties.method = 'random')})
    syn_ecop_mat<-array(NA,c(dim(mat)[1],(leads)))
    for(j in 1:leads){
      syn_at<-rsged(dim(mat)[1],mean=gl_par_arr[i,j,1],sd=gl_par_arr[i,j,2],nu=gl_par_arr[i,j,3],xi=gl_par_arr[i,j,4])
      r_syn_at<-rank(syn_at,ties.method = 'random')
      for(k in 1:length(r_syn_at)){
      syn_ecop_mat[k,j]<-syn_at[which(r_syn_at==ecop[k,j])]
      }
    }
    syn_cop[[i]]<-syn_ecop_mat
  }

#KNN process to generate array of 'ens_num' samples
  syn_cop_knn<-vector('list',length(samp_mo_seq))
 
  for(i in samp_mo_seq){
    synecop_knn<-syn_cop[[i]]
    syncop<-array(NA,c(length(n_obs),leads))
    id<-knn_lst[[i]]
    syncop<-synecop_knn[id,]
    syn_cop_knn[[i]]<-syncop
  }
    
  app_mat<-syn_cop_knn[[mo_seq[1]]][1:3,]
    
  for(i in 1:length(yr_seq)){
    seas_sim<-which(ixx_sim$mon==(mo_seq[i]-1))
    mat2<-syn_cop_knn[[mo_seq[i]]][which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
    coeff<-var_coefs[mo_seq[i],,]
    ls<-ls_sim[seas_sim,]
    ls_mat<-ls[which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
    seas_input<-which(ixx_sim$mon==(mo_seq[i]-1) & ixx_sim$year==yr_seq[i])
    ob<-ls_mat
    
    oos<-new_obs_mat[seas_sim,]
    oos_mat<-oos[which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
    syn_resid_mat<-matrix(0,ncol=leads,nrow=(dim(mat2)[1]+ar))
    syn_var_mat<-matrix(0,ncol=leads,nrow=(dim(mat2)[1]+ar))
    syn_var_mat[1:3,]<-app_mat
    for(j in 1:leads){
      for(k in (ar+1):(dim(mat2)[1]+ar)){
        bvar_res<-t(matrix(c(1,syn_var_mat[(k-1),],syn_var_mat[(k-2),],syn_var_mat[(k-3),]))) %*% matrix(coeff[j,]) + mat2[(k-ar),j]
        #aobs<-act_obs[(k-ar),j]
        o<-ob[(k-ar),j]
        ao<-oos_mat[(k-ar),j]
        fmax<-max_val[mo_seq[i],j]
        if(o > fmax | ao > fmax){fmax<-max(o,ao)}
        
        syn_var_mat[k,j]<-bvar_res[1,1]
        norm_val<-lin_mod(norm_fit[[mo_seq[i]]][[j]][1],norm_fit[[mo_seq[i]]][[j]][2],o)
        #ensure no normalization values less than zero (bad prediction from linear model)
        if(norm_val<=0){norm_val<-sd_arr[mo_seq[i],j]}
        
        res<-bvar_res[1,1] * norm_val 
        
        #INTERNAL WILD HAIR PROCESSING
        #if residual will exceed reasonable limits (1.5X the max observed forecast)
        #try first to use a the more stable normalization constant
        if((o - res)>(env_scale*fmax)){
          res<-bvar_res[1,1]*sd_arr[mo_seq[i],j]
        }
        
        #then try removing the VAR part to prevent runaway autocorrelation
        if((o - res)>(env_scale*fmax)){
          res<-mat2[(k-ar),j]*sd_arr[mo_seq[i],j]
          syn_var_mat[k,j]<-res}
        
        #finally, try resampling from the distribution and pick from samples that don't exceed reasonable limits
        if((o - res)>(env_scale*fmax)){
          new_at<-rsged(100,mean=gl_par_arr[mo_seq[i],j,1],sd=gl_par_arr[mo_seq[i],j,2],nu= gl_par_arr[mo_seq[i],j,3],xi=gl_par_arr[mo_seq[i],j,4])
          new_res<-new_at*sd_arr[mo_seq[i],j]
          samp_idx<-which((o-new_res)<(env_scale*fmax))
          res<-sample(new_res[samp_idx],1)
          syn_var_mat[k,j]<-res}
        
        
        syn_resid_mat[k,j]<-res
        
        #MINIMIZE 0 value forecasts
        #if residual from above will produce a negative value
        if(res > o){
          #generate a bunch of new samples from fitted distribution
          syn_at<-rsged(1000,mean=gl_par_arr[mo_seq[i],j,1],sd=gl_par_arr[mo_seq[i],j,2],nu= gl_par_arr[mo_seq[i],j,3],xi=gl_par_arr[mo_seq[i],j,4])
          #apply VAR model to all residuals
          var_res1<-t(matrix(c(1,syn_var_mat[(k-1),],syn_var_mat[(k-2),],syn_var_mat[(k-3),]))) %*% matrix(coeff[j,])
          #secondary set of residuals if can't get a suitable sample
          var_res2<-syn_at
          var_res<-(var_res1[1,1] + var_res2) * norm_val
          #determine positive samples which will produce a forecast flow value between the observation(cmean) and 0 flow
          var_idx<-which(var_res>=0 & var_res<=o)
          at_idx<-which(var_res2>=0 & var_res2<=o)
          #select first any residuals generated by full VAR model that will work
          if(length(at_idx)>0){syn_resid_mat[k,j]<-var_res2[sample(at_idx,1)]}
          #if no good samples there, try the raw a_t samples
          if(length(var_idx)>0){syn_resid_mat[k,j]<-var_res[sample(var_idx,1)]}
          #accept 0 as a last resort
          else{syn_resid_mat[k,j]<-o}
          #else{syn_resid_mat[k,j]<-(-min(abs(var_res2)))}
        }
        if((o - syn_resid_mat[k,j])>(env_scale*fmax)){print(paste(o - syn_resid_mat[k,j],env_scale*fmax,e,i,j,k,sep=','))}
      }
      syn_hefs_flow[j,seas_input,e]<-ob[,j] - syn_resid_mat[(ar+1):k,j]
    }
    app_mat<-syn_var_mat[(k-2):k,]
  }
}

print(paste(m,Sys.time()))
return(syn_hefs_flow)
}

print(paste('end',Sys.time()))

saveRDS(synflow_out, paste('out/syn_hefs_flow_',st_mo,'-',st_dy,'-',st_yr,'_',end_mo,'-',end_dy,'-',end_yr,'_v2_parallel.rds',sep=''))

rm(list=ls());gc()

###################################END######################################

