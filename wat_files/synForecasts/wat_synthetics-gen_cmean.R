#Script to generate synthetic forecasts with parameters from '...fit' script

library(fGarch)
library(BigVAR)
library(stringr)
library(lubridate)

source('common/GL_maineqs.R')
#source('wat_helper.R')

#date vector for fitted data
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

#date vector for hefs sub-indexing
ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day')
hefs_idx<-which(ix3=='1985-10-15'):which(ix3=='2010-09-30')

#1. User defined inputs

#1a. 'Standard' parameters that you likely won't want to change
leads <- 14 #daily leads, should stay as 14 for HEFS
ens_num <- 66 #no. of ensembles, model is currently fit to 61 members, so likely don't need to change
ar <- 3 #no. of lags in vector auto-regressive model; also don't recommend changing
loc <- "ADOC1"
#1b. Primary user defined parameters to change as desired
n <- 1 #no. of ensemble sets desired

use_observed_flows = T # use obs dataset?
# if false, use this file
#syntheticFlowFile = "C:\\Projects\\Prado_WAT_FIRO_Dev\\Watersheds\\FIRO_Prado_Dev\\runs\\WCM_Ops\\RTestFRA\\realization 1\\lifecycle 1\\event 7\\obsTimeseries.csv"
outputDir = "out\\" # local output


#2. Read in raw observed data
inf<-read.csv('data/adoc_inflow.csv')
colnames(inf) <- c("GMT", "ADOC", "OCWD")
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),2] # 2 is ADOC
obs[obs<0]<-0

#2a. Create matrix of observations matching forecast
obs_mat<-cbind(matrix(rep(obs,leads),ncol=leads))

#3. Define observed data matrix to create synthetic samples
if(use_observed_flows){
  forecastFitDir = "fit\\"
  outputDir = "out\\"
  scriptDataDir = "data\\"
  
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

# clean obs data
new_obs[new_obs<0]<-0
# create matrix by leads
new_obs_mat<-cbind(matrix(rep(new_obs,leads),ncol=leads))

## validation of start/end dates
#change 1: code alerts to erroneous entries but will not stop script
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

#change 1: if crossing into another year but no further, then define the appropriate crossover sequency of months
if((end_yr-st_yr)==1){
  samp_mo_seq<-c(st_mo:12,1:end_mo)
}



#-------------------------------------------------------------------------------------------------
#4) Create array of n (# of desired sample runs) Schaake Shuffled sequences for KNN
ats<-readRDS(paste0(forecastFitDir, loc, '_ats_cm.rds'))
gl_par_arr<-readRDS(paste0(forecastFitDir, loc, '_gl_par_arr_cm.rds'))
syn_ecop<-array(NA,c(n,ens_num,dim(obs_mat)))

print(paste(0,Sys.time()))

for(m in 1:n){
  for(e in 1:ens_num){
    for(i in 1:12){
      seas<-which(ix2$mon==(i-1))
      mat<-ats[e,seas,]
      ecop<-apply(mat,2,function(x){rank(x,ties.method = 'random')})
      syn_ecop_mat<-array(NA,c(dim(mat)[1],leads))
      #generate new a_t from fitted dist'n and reorder via Schaake shuffle
      for(j in 1:leads){
        syn_at<-rsged(dim(mat)[1],mean=0,sd=1,nu= (2 / (1 + gl_par_arr[e,i,j,3])),xi=gl_par_arr[e,i,j,4])
        r_syn_at<-rank(syn_at,ties.method = 'random')
        for(k in 1:length(r_syn_at)){
          syn_ecop_mat[k,j]<-syn_at[which(r_syn_at==ecop[k,j])]
        }
      }
      syn_ecop[m,e,seas,]<-syn_ecop_mat #resample a_t matrices indexed by sample ID 'm'
    }
  }
  print(paste(m,Sys.time()))
}

saveRDS(syn_ecop, paste0(outputDir, loc, '_syn_ecop_cm.rds'))

#-------------------------------------------------------------------------------------------------------
#5) Synthetic Generation
#date vector for simulated data
st<-paste(st_yr,st_mo,st_dy,sep='-')
ed<-paste(end_yr,end_mo,end_dy,sep='-')

ix_sim<-seq(as.Date(st),as.Date(ed),'day')
ixx_sim<-as.POSIXlt(ix_sim)

#5a. Create conditional mean matrix for simulated data

#conditional mean function for out of sample data
cond_mean<-function(Qfit,Qsim,f){
  cmean <- mean(f) + (cov(Qfit,f)/(sd(Qfit))^2) * (Qsim - mean(Qfit))
  return(cmean)
}

#forecast matrices for fitting
hefs_mat<-readRDS(paste0(scriptDataDir,loc,'_hefs_ens_forc.rds'))
#scale by 1000 to convert from kcfs to cfs
hefs_mat<-hefs_mat[,hefs_idx,]*1000

cmean_sim<-array(NA,c(ens_num,length(ix_sim),leads))

for(e in 1:ens_num){
  
  #1a. Raw Resid Matrix
  #calculate raw residuals based on conditional mean (cmean) estimation (by month)
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    seas_sim<-which(ixx_sim$mon==(i-1))
    sim_inf<-hefs_mat[e,seas,]
    obs_inf<-obs_mat[seas,]
    obs_sim<-new_obs_mat[seas_sim,]
    cmean_inf<-array(NA,dim(obs_sim))
    for(j in 1:leads){
      cmean_inf[,j]<-cond_mean(obs_inf[,j],obs_sim[,j],sim_inf[,j])
      cmean_inf[cmean_inf<0]<-0
    }
    cmean_sim[e,seas_sim,]<-cmean_inf
  }
}

#5a. set up index for sequential sampling and VAR continuity across months

##NOTE--this section is a little wonky and hard to follow. It sets up an index to sample in a continuous
#fashion across the simulated timespan to maintain VAR model continuity...
#I could explain in more detail in person, but recommend to just accept as is for now

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

knn<-round(sqrt(length(seas))) #knn set to square root of monthly # of samples (from fitted data)
tot<-sum(rep(1,knn) / 1:knn)
wts<-rep(1,knn) / 1:knn / rep(tot,knn) #weights for kernel weighted sampling

#5c. load required fit data
syn_ecop<-readRDS(paste0(outputDir, loc, '_syn_ecop_cm.rds'))
var_coefs<-readRDS( paste0(forecastFitDir, loc, '_var_coefs_cm.rds'))
gl_par_arr<-readRDS( paste0(forecastFitDir, loc, '_gl_par_arr_cm.rds'))
cmean<-readRDS( paste0(forecastFitDir, loc, '_cmean.rds'))

#5d. define matrices to store synthetic forecast residuals and forecasts themselves
syn_hefs_resid<-array(NA,c(n,ens_num,length(ixx_sim),leads))
syn_hefs_flow<-array(NA,c(n,ens_num,length(ixx_sim),leads))

print(Sys.time()) #start time


#Script generates 'n' ensembles of size 'ens_num'
for(m in 1:n){
  
  knn_lst<-vector('list',length(samp_mo_seq))
  
  #change 1: only include sampled months for shorter periods
  for(i in samp_mo_seq){ ###
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
        ob_val<-n_obs[j]
        y<-sqrt((ob_val - ob)^2) #find NEP closest by Euclidean distance
        x<-sort(y)
        x<-x[1:knn] #sort top values
        s<-sample(x,1,prob=wts) #distance weighted kNN sampling
        id<-which(y==s)
        if(length(id)>1) {
          id <- sample(id,1)} #resample the sample for any duplicated values
        knn_vec[j]<-id}
    }
    knn_lst[[i]]<-knn_vec
  }
  
  #loop across all 'ens_num' ensemble members
  for(e in 1:ens_num){
  
    #KNN process to generate array of 'ens_num' samples
    syn_cop_knn<-vector('list',length(samp_mo_seq))
    syn_ecop_knn<-syn_ecop[m,e,,]
    
    #change 1: only include sampled months for shorter periods
    for(i in samp_mo_seq){ ###
      seas<-which(ix2$mon==(i-1))
      seas_sim<-which(ixx_sim$mon==(i-1))
      synecop_knn<-syn_ecop_knn[seas,]
      syncop<-array(NA,c(length(seas_sim),leads))
      id<-knn_lst[[i]]
      syncop<-synecop_knn[id,]
      syn_cop_knn[[i]]<-syncop
    }
    
    # Minimum of 3-months of spin-up for the model? -EAH
    #initial set-up of 'app_mat' which is matrix of (t-3:t-1) previous month errors to append to current month
    #matrix to maintain VAR continuity
    app_mat<-syn_cop_knn[[mo_seq[1]]][1:3,]
    
    #proceed sequentially through 'yr_seq' index to maintain VAR continuity
    for(i in 1:length(yr_seq)){
      seas_sim<-which(ixx_sim$mon==(mo_seq[i]-1))
      mat<-syn_cop_knn[[mo_seq[i]]][which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
      cmn<-cmean[e,seas_sim,]
      cmean_mat<-cmn[which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
      coeff<-var_coefs[e,mo_seq[i],,]
      seas_input<-which(ixx_sim$mon==(mo_seq[i]-1) & ixx_sim$year==yr_seq[i])
      ob<-cmean_mat
      #create empty residual matrix with 3 extra spaces for
      syn_resid_mat<-matrix(0,ncol=leads,nrow=(dim(mat)[1]+ar))
      syn_resid_mat[1:3,]<-app_mat #add previous months errors to empty matrix
      #generate new errors using fitted VAR coefficients, sampled a_t, and heteroscedastic scaling
      for(j in 1:leads){
        for(k in (ar+1):(dim(mat)[1]+ar)){
          #VAR linear model error estimation
          syn_resid_mat[k,j]<-t(matrix(c(syn_resid_mat[(k-1),],syn_resid_mat[(k-2),],syn_resid_mat[(k-3),]))) %*% matrix(coeff[j,]) + 
            #plus heteroscedastic scaled random error term
            sigma_t(gl_par_arr[e,mo_seq[i],j,1],gl_par_arr[e,mo_seq[i],j,2],ob[(k-ar),j])*mat[(k-ar),j]
          #ensure residuals produce >= 0 output
          if(syn_resid_mat[k,j]>ob[(k-ar),j]) syn_resid_mat[k,j]<-ob[(k-ar),j]
        }
        syn_hefs_resid[m,e,seas_input,j]<-syn_resid_mat[(ar+1):k,j]
        #synthetic forecasts are observation (in this case conditional mean) - errors
        syn_hefs_flow[m,e,seas_input,j]<-ob[,j] - syn_resid_mat[(ar+1):k,j]
      }
      app_mat<-syn_resid_mat[(k-2):k,]
    }
  }
  print(paste(m,Sys.time())) #keep track of progress, verbose
  #saveRDS(syn_hefs_flow,  paste0(outputDir, loc, '_syn_hefs_flow_cm.rds'))
  #saveRDS(syn_hefs_resid, paste0(outputDir, loc, '_syn_hefs_resid_cm.rds') #commented out, you probably don't really need the forecast residuals for anything
}

if(anyNA(syn_hefs_flow)==T){stop('Bad Forecast Output')}

# create plots?
#source("diagnostics.R")
# write out to sqlite file?
#source("syn-hefs_out_tsensembles.R")

#remove variables and clean environment
#if(useGC){
#  rm(list=ls());gc()
#}

###################################END######################################


