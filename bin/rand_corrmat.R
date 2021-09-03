library(RcppArmadillo)
library(Rcpp)
sourceCpp("/sfs/qumulo/qhome/ycp6wm/Thresholding_Sims/bin/rand_corrmat.cpp")
rand_corr_wrapper <- function(d, alpha_func = NULL){
  if(is.null(alpha_func)){
  alpha <- c((d:2)/2,0)
  }else{
    alpha = alpha_func(d)
  }
  rands <- diag(d)
  
  for(i in 1:d-1){
    rands[i, (i+1):d] = (unlist(Map(rbeta,1, alpha[(-(d-i+1)):-d],alpha[(-(d-i+1)):-d]))-.5)*2
  }
  
  temp = rand_corr(d, rands)
  
  return(temp)
  
  
}