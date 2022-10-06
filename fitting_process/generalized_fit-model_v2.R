#Script to fit conditional mean, VAR, and GL model to forecast residuals and create
#Schaake shuffled empirical copula to sample from

library(fGarch)
library(BigVAR)

#ens_num <- 68 #no of ensemble members
#leads <- 14 #total number of lead times
#ar<-3 #lags for VAR model

#date indices, fit is calculated starting 15 days after beginning of data since full forecast data for all leads required
ix<-seq(fit_start_date, fit_end_date,'day')
ix2<-as.POSIXlt(ix)

#required functions for normalization
lin_mod<-function(intcpt,slope,x){
  y = intcpt + slope * x
  return(y)
}

lin_fit<-function(pars,x,y){
  pred = pars[1] + pars[2] * x
  err = sum((y - pred)^2)
  return(err)
}


#for loop to do across all 61 ensembles, takes a while to fit all samples so 
#might want to comment this out to just play around with a single ensemble
print(Sys.time()) #start time 

#define matrices for normalized residuals, conditional mean estimate, decorrelated residuals, and a_t
nresids<-array(NA,c(ens_num,length(ix),leads))
ucresids<-array(NA,c(ens_num,length(ix),leads))
ats<-array(NA,c(ens_num,length(ix),leads))

loess_mat<-readRDS(paste0('fit/', location_name, '_loess_mat_v2.rds'))
rresids<-readRDS(paste0('fit/', location_name, '_rresids_v2.rds'))

param_lst<-vector('list',ens_num)

lbs<-c(0,0) #lower bounds for intcpt, slope
ubs<-c(1000,10) #upper bounds for intcpt, slope
sts<-c(0,1) #starting parameters for intcpt, slope

for(e in 1:ens_num){
  param_lst[[e]]<-vector('list',4)
  sd_arr<-array(NA,c(12,leads))
  norm_fit<-vector('list',12)
  norm_fit_sub<-vector('list',leads)
  
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    obmat<-loess_mat[seas,]
    rresid_mat<-rresids[e,seas,]
    rresid_res<-abs(rresid_mat)
    norm_fit[[i]]<-norm_fit_sub
    for(k in 1:leads){
      ob<-obmat[,k]
      res<-rresid_res[,k]
      res_sort<-sort(res)
      lbs[1]<-mean(res_sort[1:round(0.1*length(seas))])
      ubs[1]<-mean(res)
      sd_arr[i,k]<-mean(res)
      
      lfit<-optim(par=sts,lin_fit,x=ob,y=res,
                  method = 'L-BFGS-B',lower = lbs,upper = ubs,
                  control = list(fnscale=1,maxit=100000))
      
      norm_fit[[i]][[k]]<-lfit$par
    }
  }
  param_lst[[e]][[1]]<-norm_fit
  param_lst[[e]][[2]]<-sd_arr

  #2) Normalize

  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    obmat<-loess_mat[seas,]
    rresid_mat<-rresids[e,seas,]
    nresid_mat<-array(NA,dim(rresid_mat))
    for(k in 1:leads){
      ob<-obmat[,k]
      res<-rresid_mat[,k]
      norm_vec<-lin_mod(norm_fit[[i]][[k]][1],norm_fit[[i]][[k]][2],ob)
      norm_vec[norm_vec<=0]<-sd_arr[i,k]
      nresid_mat[,k]<-res / norm_vec
    }
    nresids[e,seas,]<-nresid_mat
  }

  var_coefs <- array(NA,c(12,leads,(leads*ar+1)))

  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    rresid_mat<-nresids[e,seas,]
    m1 = constructModel(rresid_mat, p = ar, struct = "BGR", gran = c(50, 10),IC = F,
                      verbose = T, VARX = list(),separate_lambdas = F,model.controls=list(intercept = F, MN=F))
  
    m1_res = cv.BigVAR(m1)
  
    var_coefs[i,,]<-m1_res@betaPred
    resids<-m1_res@fitted #calculate decorr residuals from 
  
    mat<-matrix(0,dim(rresid_mat)[1],dim(rresid_mat)[2])
    mat[(ar+1):dim(mat)[1],]<-resids #resids set to 0 for first columns = max lag
  
    uc_mat<-rresid_mat - mat
  
    ucresids[e,seas,]<-uc_mat
  }
  
  param_lst[[e]][[3]]<-var_coefs
  print(paste('varfit',e,'complete',Sys.time()))
  
  
  gl_par_arr<-array(NA,c(12,leads,4))
  at_lst<-vector('list',12)
  
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    uc_sim<-ucresids[e,seas,]
    at_arr<-uc_sim
    for(j in 1:leads){
      gl_mle<-sgedFit(uc_sim[,j])
      gl_par_arr[i,j,]<-gl_mle$par
    }
    ats[e,seas,]<-at_arr
  }
  param_lst[[e]][[4]]<-gl_par_arr
  
  print(paste('ensemble fit',e,'complete',Sys.time()))
} 


#save matrices in R data structure format
saveRDS(rresids,paste0('fit/', location_name, '_rresids_v2.rds'))
saveRDS(nresids,paste0('fit/', location_name, '_nresids_v2.rds'))
saveRDS(ucresids,paste0('fit/', location_name, '_ucresids_v2.rds'))
saveRDS(ats,paste0('fit/', location_name, '_ats_v2.rds'))

saveRDS(param_lst,paste0('fit/', location_name, '_param_lst_v2.rds'))

#remove variables and clean environment
#rm(list=ls());gc()

############################################END#########################################################
