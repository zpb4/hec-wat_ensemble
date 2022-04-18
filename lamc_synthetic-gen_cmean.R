#Script to generate synthetic forecasts with parameters from '...fit' script

setwd('h:/firo_lamc/hec-wat_ensemble/')
library(fGarch)
library(BigVAR)
library(stringr)

source('GL_maineqs.R')

#date vector for fitted data
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)
ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day')
hefs_idx<-which(ix3=='1985-10-15'):which(ix3=='2010-09-30')

#1. User defined inputs

#1a. 'Standard' parameters that you likely won't want to change
leads <- 14 #daily leads, should stay as 14 for HEFS
ens_num <- 61 #no. of ensembles, model is currently fit to 61 members, so likely don't need to change
ar <- 3 #no. of lags in vector auto-regressive model; also don't recommend changing

#1b. Primary user defined parameters to change as desired
n <- 1 #no. of ensemble sets desired
#Define simulation start in year, month, and day; minimum 1948-10-01
 st_yr <- 1950 #4 digit year
 st_mo <- 01 #specify with leading zero for single digits, e.g. '01' instead of '1'
 st_dy <- 01 #specify with leading zero for single digits, e.g. '01' instead of '1'
 
#Define simulation end in year, month, and day; maximum 2010-09-30
 end_yr <- 1952 #4 digit year
 end_mo <- 07 #specify with leading zero for single digits, e.g. '01' instead of '1'
 end_dy <- 09 #specify with leading zero for single digits, e.g. '01' instead of '1'


#2. Read in raw observed data
inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),5]
obs[obs<0]<-0

#2a. Create matrix of observations matching forecast
obs_mat<-cbind(matrix(rep(obs,leads),ncol=leads))

#3. Define observed data matrix to create synthetic samples
st_date<-paste(str_remove(st_mo,'^0'),str_remove(st_dy,'^0'),st_yr,sep='/')
end_date<-paste(str_remove(end_mo,'^0'),str_remove(end_dy,'^0'),end_yr,sep='/')
new_obs<-inf[which(inf$GMT==paste(st_date,' 12:00',sep='')):which(inf$GMT==paste(end_date,' 12:00',sep='')),5]
new_obs[new_obs<0]<-0
new_obs_mat<-cbind(matrix(rep(new_obs,leads),ncol=leads))


#-------------------------------------------------------------------------------------------------
#4) Create array of n (# of desired sample runs) Schaake Shuffled sequences for KNN
ats<-readRDS('fit/ats_cm.rds')
gl_par_arr<-readRDS('fit/gl_par_arr_cm.rds')
syn_ecop<-array(NA,c(n,ens_num,dim(obs_mat)))

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

saveRDS(syn_ecop,'out/syn_ecop_cm.rds')

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
lamc_hefs<-readRDS('data/lamc_hefs_ens_forc.rds')
#scale by 1000 to convert from kcfs to cfs
lamc_hefs<-lamc_hefs[,hefs_idx,]*1000

cmean_sim<-array(NA,c(ens_num,length(ix_sim),leads))

for(e in 1:ens_num){
  
  #1a. Raw Resid Matrix
  #calculate raw residuals based on conditional mean (cmean) estimation (by month)
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    seas_sim<-which(ixx_sim$mon==(i-1))
    sim_inf<-lamc_hefs[e,seas,]
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
yr_idx_lst<-vector('list',12)

for(i in 1:12){
  seas3<-which(ixx_sim$mon==(i-1))
  yr_idx_lst[[i]]<-yr_idx[seas3]
}

styr_idx<-st_yr-1900
edyr_idx<-end_yr-1900

yr_seq<-c(rep(styr_idx,length=(13-st_mo)),rep((styr_idx+1):(edyr_idx-1),each=12),rep((edyr_idx),length=end_mo))
mo_seq<-c(st_mo:12,rep(1:12,length((styr_idx+1):(edyr_idx-1))),1:end_mo)


#5b. knn set up
seas<-which(ix2$mon==0) #find length of january monthly subset of fitted data

knn<-round(sqrt(length(seas))) #knn set to square root of monthly # of samples (from fitted data)
tot<-sum(rep(1,knn) / 1:knn)
wts<-rep(1,knn) / 1:knn / rep(tot,knn) #weights for kernel weighted sampling

#5c. load required fit data
syn_ecop<-readRDS('out/syn_ecop_cm.rds')
var_coefs<-readRDS('fit/var_coefs_cm.rds')
gl_par_arr<-readRDS('fit/gl_par_arr_cm.rds')
cmean<-readRDS('fit/cmean.rds')

#5d. define matrices to store synthetic forecast residuals and forecasts themselves
syn_hefs_resid<-array(NA,c(n,ens_num,length(ixx_sim),leads))
syn_hefs_flow<-array(NA,c(n,ens_num,length(ixx_sim),leads))

print(Sys.time()) #start time


#Script generates 'n' ensembles of size 'ens_num'
for(m in 1:n){
  
  knn_lst<-vector('list',12)
  
  for(i in 1:12){
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
    syn_cop_knn<-vector('list',12)
    syn_ecop_knn<-syn_ecop[m,e,,]
    
    for(i in 1:12){
      seas<-which(ix2$mon==(i-1))
      seas_sim<-which(ixx_sim$mon==(i-1))
      synecop_knn<-syn_ecop_knn[seas,]
      syncop<-array(NA,c(length(seas_sim),leads))
      id<-knn_lst[[i]]
      syncop<-synecop_knn[id,]
      syn_cop_knn[[i]]<-syncop
    }
    
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
  saveRDS(syn_hefs_flow, 'out/syn_hefs_flow_cm.rds')
  #saveRDS(syn_hefs_resid, 'out/syn_hefs_resid_cm.rds') #commented out, you probably don't really need the forecast residuals for anything
}

#remove variables and clean environment
#rm(list=ls());gc()

###################################END######################################