#Script to fit conditional mean, VAR, and GL model to forecast residuals and create
#Schaake shuffled empirical copula to sample from

#setwd('h:/firo_lamc/ensemble-will/')
library(fGarch)
library(BigVAR)

ens_num <- 68 #no of ensemble members
leads <- 14 #total number of lead times
ar<-3 #lags for VAR model
loc<-"ADOC1"

#1) Read in observed and simulated inflows
inf<-read.csv('data/ADOC_inflow.csv')
colnames(inf) <- c("GMT", "ADOC", "OCWD")
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),2] # col 2 is obsv ADOC flow

#set negative fnf observations to zero
obs[obs<0]<-0

#define matrices of observations that match forecast output for residual calculations
#e.g. repeat observations across 14 forecast leads, so this is a (n x 14) matrix
obs_mat<-cbind(matrix(rep(obs,leads),ncol=leads))

#date indices, fit is calculated starting 15 days after beginning of data since full forecast data for all leads required
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)
ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day')
hefs_idx<-which(ix3=='1985-10-15'):which(ix3=='2010-09-30')

#conditional mean function
cond_mean<-function(Q,f){
  cmean <- mean(f) + (cov(Q,f)/(sd(Q))^2) * (Q - mean(Q))
  return(cmean)
}

#forecast matrices
hefs_mat<-readRDS(paste0('data/', loc, '_hefs_ens_forc.rds'))

#scale by 1000 to convert from kcfs to cfs
hefs_mat<-hefs_mat[,hefs_idx,]*1000

# Drop members that create extreme outliers
dropMembers = c(14,60)
keepMembers = which(!(1:ens_num %in% dropMembers))
hefs_mat = hefs_mat[keepMembers,,]
ens_num = length(keepMembers)

#for loop to do across all `ens_num` ensembles, takes a while to fit all samples so 
#might want to comment this out to just play around with a single ensemble
print(Sys.time()) #start time 

#define matrices for raw residuals, conditional mean estimate, decorrelated residuals, and a_t
rresids<-array(NA,c(ens_num,length(ix),leads))
cmean<-array(NA,c(ens_num,length(ix),leads))
uc_resid<-array(NA,c(ens_num,length(ix),leads))
ats<-array(NA,c(ens_num,length(ix),leads))

gl_par_arr<-array(NA,c(ens_num,12,leads,4))
var_coefs<-array(NA,c(ens_num,12,leads*3,(leads*ar)))

for(e in 1:ens_num){

  #1a. Raw Resid Matrix
  #calculate raw residuals based on conditional mean (cmean) estimation (by month)
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    sim_inf<-hefs_mat[e,seas,]
    obs_inf<-obs_mat[seas,]
    cmean_inf<-array(NA,dim(obs_inf))
      for(j in 1:leads){
        cmean_inf[,j]<-cond_mean(obs_inf[,j],sim_inf[,j])
        cmean_inf[cmean_inf<0]<-0
      }
    rresids[e,seas,]<-cmean_inf - sim_inf
    cmean[e,seas,]<-cmean_inf
  }

  #2) BigVar model (LASSO penalized VAR model)
  #VAR model to create uncorrelated residuals

  #Calculate uncorrelated matrices
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    rresid_mat<-rresids[e,seas,]
    mc = list(intercept=F, MN=F) # model controls variable
    m1 = constructModel(rresid_mat, p = ar, struct = "Basic", gran = c(25, 10),IC = F,
                      verbose = F, VARX = list(), separate_lambdas = F, model.controls = mc)
    m1_res = cv.BigVAR(m1)
  
    var_coefs[e,i,,]<-m1_res@betaPred[,2:length(m1_res@betaPred[1,])]
    resids<-m1_res@fitted #calculate decorr residuals from 
  
    mat<-matrix(0,dim(rresid_mat)[1],dim(rresid_mat)[2])
    mat[(ar+1):dim(mat)[1],]<-resids #resids set to 0 for first columns = max lag
  
    uc_mat<-rresid_mat - mat
  
    uc_resid[e,seas,]<-uc_mat
  }

  #3) GL SGED Model
  #monthly model fits for each lead time
  source('GL_maineqs.R')
  st = c(0.5,0.5,0,1)
  lb = c(0.01,0,-0.99,0.1)
  ub = c(10,10,1,10)

  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    uc_sim<-uc_resid[e,seas,]
    cmean_mat<-cmean[e,seas,]

    at_arr<-array(NA,c(dim(uc_sim)[1],leads))
    for(j in 1:leads){

      gl_mle<-optim(par=st,fn=GL_fun_noscale_var,inflow=cmean_mat[,j],et=uc_sim[,j],
                  method = 'L-BFGS-B',lower = lb,upper = ub,
                  control = list(maxit=100000))
    
      at<-a_t_var(uc_sim[,j],sigma_t(gl_mle$par[1],gl_mle$par[2],cmean_mat[,j])) #cmean scaling
      gl_par_arr[e,i,j,]<-c(gl_mle$par[1],gl_mle$par[2],gl_mle$par[3],gl_mle$par[4])
      at_arr[,j]<-at
    }
  ats[e,seas,]<-at_arr
  }
  print(paste(e,Sys.time()))
} 


#save matrices in R data structure format
saveRDS(rresids,paste0('fit/', loc, '_rresids_cm.rds'))
saveRDS(cmean,paste0('fit/', loc, '_cmean.rds'))
saveRDS(ats,paste0('fit/', loc, '_ats_cm.rds'))
saveRDS(uc_resid,paste0('fit/', loc, '_uc_resid_cm.rds'))

saveRDS(gl_par_arr,paste0('fit/', loc, '_gl_par_arr_cm.rds'))
saveRDS(var_coefs,paste0('fit/', loc, '_var_coefs_cm.rds'))

#remove variables and clean environment
#rm(list=ls());gc()

############################################END#########################################################