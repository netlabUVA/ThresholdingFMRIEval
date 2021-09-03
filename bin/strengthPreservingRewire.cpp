#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
sp_mat strength_preserve_rewire(sp_mat data, mat out_in, double mean_k, double mean_s, mat out_in_deg){
  
  sp_mat::iterator it = data.begin();
  sp_mat::iterator it_end = data.end();
  
  for(; it!=it_end; ++it){
    // Rcout << it.row() <<" " <<it.col() << std::endl;
    //Rcout << (mean_k/mean_s) << " " <<out_in(it.row(),0) << " " <<out_in(it.col(), 1) << " " << out_in_deg(it.col(),0) <<  " " <<out_in_deg(it.row(),1) << std::endl;
    
    *it = (mean_k/mean_s)*out_in(it.row(),0)*out_in(it.col(), 1)/(out_in_deg(it.col(),0)*out_in_deg(it.row(),1));
  }
  
  return data;
  
  
}