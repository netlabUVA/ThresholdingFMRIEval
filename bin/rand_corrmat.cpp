// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; using namespace arma;


//Implementation of the random correlation matrix algorithm by Joe 2006 "Generating random correlation matrices based on partial correlations"
// rands, to ensure that the space of random correlation matrices is uniform, takes a matrix of random betas calculated in a specific fashion.
// [[Rcpp::export]]
arma::mat rand_corr(int d, arma::mat rands) {
  
  
  arma::mat rand_corr = arma::eye(d, d);
  rand_corr.diag(1) = rands.diag(1);
  rand_corr.diag(-1) = rands.diag(1);
  for(int k = 2; k < (d); k++){
    for(int j = 0; j < (d-k); j++){
      // Rcout << j+1 << " " << j+k+1 <<std::endl;
      arma::rowvec r_1 = rand_corr(j, span(j+1, j+k-1));
      //  Rcout << r_1 << std::endl;
      arma::rowvec r_3 = rand_corr(j+k, span(j+1, j+k-1));
      //  Rcout << r_3 << std::endl;
      arma::mat R_2 = rand_corr.submat(j+1, j+1, j+k-1, j+k-1);
      //  Rcout << R_2 << std::endl;
      double rand = rands(j,j+k);
      
      arma::mat D_2jk = (1 - r_1*arma::inv(R_2)*arma::trans(r_1))*(1 - r_3*arma::inv(R_2)*arma::trans(r_3));
      rand_corr(j,j+k) = arma::as_scalar(r_1*arma::inv(R_2)*arma::trans(r_3) + rand*D_2jk);
      rand_corr(j+k,j) =rand_corr(j,j+k);
    }
  }
  return(rand_corr);
}