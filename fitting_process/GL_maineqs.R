#Main equations for Generalized Likelihood Function

source('GL_subeqs.r') #bring in subequations in terms of parameters below

#1) GL function test that takes uncorrelated residuals(et) from VAR model
#pars in order: 
#par[1]: sigma_0, intercept for linearly scaled SD (heteroscedastic)
#par[2]: sigma_1, linear coefficient for scale SD
#par[3]: beta, kurtosis parameter (-1,1) for normalized SEP
#par[4]: xi, skewness parameter (0.1,10) for normalized SEP 
#for this version of function, but retained for more general version
GL_fun_noscale_var<-function(pars,inflow,et){
  n<-length(inflow)
  sig_xi<-sigma_xi(M1(pars[3]),M2,pars[4]) #Eq A8 from subequations
  om_b<-omega_beta(pars[3]) #Eq A2 from subequations
  sig_t<-sigma_t(pars[1],pars[2],inflow) #Eq 5 from subequations
  cb<-c_beta(pars[3]) #Eq A3 from subequations
  
  #Eq 6, derived residuals a_xt as a function of observed residuals modified by subequations
  a_xt<-a_xi_t(pars[4],mu_xi(M1(pars[3]),pars[4]),sig_xi,a_t_var(et,sig_t))
  
  #Eq 8, Generalized Log Likelihood Function
  gl_ll<-n*log((2*sig_xi*om_b)/(pars[4] + pars[4]^-1)) - sum(log(sig_t)) - cb * sum(abs(a_xt)^(2/(1+pars[3])))
  
  if (gl_ll == Inf|gl_ll == -Inf) ll<-0 #optimizer can't handle Inf or -Inf values
  else ll<-gl_ll
  return(-ll) #mult by -1 to enable maximization
}

#2) define SEP density as function of calculated xi and beta and random variable axt (a_xi_t)
SEP_dens<-function(xi,beta,axt){
  SEP_dens<-(2*sigma_xi(M1(beta),M2,xi))/(xi + xi^(-1))*omega_beta(beta)*
    exp(-c_beta(beta)*abs(a_xi_t(xi,mu_xi(M1(beta),xi),sigma_xi(M1(beta),M2,xi),axt))^(2/(1+beta)))
}

