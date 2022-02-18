#Script to generate synthetic forecasts with parameters from '...fit' script

setwd('h:/firo_lamc/ensemble-will/')
library(fGarch)
library(BigVAR)
source('GL_maineqs.R')

leads <- 14
ens_num <- 61
ar <- 3
n <- 10

#date vectors
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

#read in raw observed data
inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),5]
obs[obs<0]<-0

#create matrix of observations matching forecast
obs_mat<-cbind(matrix(rep(obs,leads),ncol=leads))

#Generate new out of sample data with KNN method
new_obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),5]
new_obs[new_obs<0]<-0
new_obs_mat<-cbind(matrix(rep(new_obs,leads),ncol=leads))

#1) Create array of n (# of desired sample runs) Schaake Shuffled sequences for KNN
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
##Synthetic Generation
#set up index for sequential sampling and VAR continuity across months
yr_idx<-ix2$year
yr_idx_lst<-vector('list',12)

for(i in 1:12){
  seas3<-which(ix2$mon==(i-1))
  yr_idx_lst[[i]]<-yr_idx[seas3]
}

yr_seq<-c(rep(85,3),rep(86:109,each=12),rep(110,9))
mo_seq<-c(10:12,rep(1:12,length(86:109)),1:9)


#knn set up
knn<-round(sqrt(length(seas3))) #knn set to square root of monthly # of samples
tot<-sum(rep(1,knn) / 1:knn)
wts<-rep(1,knn) / 1:knn / rep(tot,knn) #weights for kernel weighted sampling

syn_ecop<-readRDS('out/syn_ecop_cm.rds')
var_coefs<-readRDS('fit/var_coefs_cm.rds')
gl_par_arr<-readRDS('fit/gl_par_arr_cm.rds')
cmean<-readRDS('fit/cmean.rds')

#define matrices to store synthetic forecast residuals and forecasts themselves
syn_hefs_resid<-array(NA,c(n,ens_num,length(ix2),leads))
syn_hefs_flow<-array(NA,c(n,ens_num,length(ix2),leads))

print(Sys.time()) #start time


#Script generates 'n' ensembles of size 'ens_num'
for(m in 1:n){
  
  knn_lst<-vector('list',12)
  
  for(i in 1:12){
    knn_vec<-c()
    seas<-which(ix2$mon==(i-1))
    n_obs<-new_obs[seas]
    ob<-obs[seas]
    id0 <- which(ob == 0)
    for(j in 1:length(n_obs)){
      if(n_obs[j] == 0 & length(id0)>=1) {s<-sample(id0,1); knn_vec[j]<-s} 
      else {ob_val<-n_obs[j]
      y<-sqrt((ob_val - ob)^2) #find NEP closest by Euclidean distance
      x<-sort(y)
      x<-x[1:knn] #sort top values
      s<-sample(x,1,prob=wts)
      id<-which(y==s)
      if(length(id)>1) {id <- sample(id,1)} #resample the sample for any duplicated values
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
      synecop_knn<-syn_ecop_knn[seas,]
      syncop<-array(NA,c(length(seas),leads))
      id<-knn_lst[[i]]
      syncop<-synecop_knn[id,]
      syn_cop_knn[[i]]<-syncop
    }
    
    #initial set-up of 'app_mat' which is matrix of (t-3:t-1) previous month errors to append to current month
    #matrix to maintain VAR continuity
    app_mat<-syn_cop_knn[[mo_seq[1]]][1:3,]
    
    #proceed sequentially through 'yr_seq' index to maintain VAR continuity
    for(i in 1:length(yr_seq)){
      seas2<-which(ix2$mon==(mo_seq[i]-1))
      mat2<-syn_cop_knn[[mo_seq[i]]][which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
      cmn<-cmean[e,seas2,]
      cmean_mat<-cmn[which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
      coeff<-var_coefs[e,mo_seq[i],,]
      seas3<-which(ix2$mon==(mo_seq[i]-1) & ix2$year==yr_seq[i])
      ob<-cmean_mat
      #create empty residual matrix with 3 extra spaces for
      syn_resid_mat<-matrix(0,ncol=leads,nrow=(dim(mat2)[1]+ar))
      syn_resid_mat[1:3,]<-app_mat #add previous months errors to empty matrix
      #generate new errors using fitted VAR coefficients, sampled a_t, and heteroscedastic scaling
      for(j in 1:leads){
        for(k in (ar+1):(dim(mat2)[1]+ar)){
          #VAR linear model error estimation
          syn_resid_mat[k,j]<-t(matrix(c(syn_resid_mat[(k-1),],syn_resid_mat[(k-2),],syn_resid_mat[(k-3),]))) %*% matrix(coeff[j,]) + 
            #plus heteroscedastic scaled random error term
            sigma_t(gl_par_arr[e,mo_seq[i],j,1],gl_par_arr[e,mo_seq[i],j,2],ob[(k-ar),j])*mat2[(k-ar),j]
          #ensure residuals produce >= 0 output
          if(syn_resid_mat[k,j]>ob[(k-ar),j]) syn_resid_mat[k,j]<-ob[(k-ar),j]
        }
        syn_hefs_resid[m,e,seas3,j]<-syn_resid_mat[(ar+1):k,j]
        #synthetic forecasts are observation (in this case conditional mean) - errors
        syn_hefs_flow[m,e,seas3,j]<-ob[,j] - syn_resid_mat[(ar+1):k,j]
      }
      app_mat<-syn_resid_mat[(k-2):k,]
    }
  }
  print(paste(m,Sys.time())) #keep track of progress, verbose
  saveRDS(syn_hefs_flow, 'out/syn_hefs_flow_cm.rds')
  saveRDS(syn_hefs_resid, 'out/syn_hefs_resid_cm.rds')
}

#remove variables and clean environment
rm(list=ls());gc()

###################################END######################################