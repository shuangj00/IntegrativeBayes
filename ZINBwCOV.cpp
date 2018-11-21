#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double norm_rs(double a, double b);
double half_norm_rs(double a, double b);
double unif_rs(double a, double b);
double exp_rs(double a, double b);
double rnorm_trunc(double mu, double sigma, double lower, double upper);

// [[Rcpp::export]]
// main function
List zinb_w_cov(NumericMatrix Y_mat, NumericMatrix X_mat,
             NumericVector z_vec, 
             NumericVector s_vec,
             double mu0_mean,
             double tau_mukj = 1,
             int S = 20000, 
             double burn_rate = 0.5,
             double tau_mu0 = 1, 
             double tau_phi = 1, 
             double tau_beta = 1,
             double a_omega = 0.2, double b_omega = 1.8, // uniform prior 
             double a_pi = 1, double b_pi = 0.1, // control the sparsity of R matrix
             double a_p = 0.4, double b_p = 1.6, // control the sparsity of Delta matrix
             double a_mu = 2, double b_mu = 15, // spike and slab for mu_kj
             double a_phi = 10, double b_phi = 1, // dispersion parameter
             double a_beta = 2, double b_beta = 15, // spike and slab for beta_rj
             double phi_low = 1, 
             double beta_lim = 5
){
  // Sample statistics
  int n = Y_mat.nrow() ;
  int p = Y_mat.ncol() ;
  int R = X_mat.ncol() ;
  int K = 2 ;
  int count = 0;
  double pi = 3.14159265;
  
  // Initialize index vectors for variable selection
  // (1) Feature index vector:
  NumericVector J_vec(p), J_filter(p);
  for(int j =0; j < p; j++){
    J_vec[j] = j + 1; // feature selection index
    J_filter[j] = 1;
  }
  
  // In each iteration, we choose 5 of the fetures /5% of the covariates to perform between-model update
  int num_feature = ceil(0.05 * p);
  
  // Starting point of burn-in process
  int S_burn = ceil(S * burn_rate); 
  
  ////////////// pre-processing: remove the features with too many 0's in a group//////////////
  for(int j = 0; j < p; j++){
    // for each feature, we exam if it only has 0 or 1 none-zero elements:
    for(int k =1; k <= K; k++){
      int gp_size = 0;
      int zero_count = 0;
      for(int i = 0; i < n; i++){
        if( z_vec[i] == k){
          gp_size = gp_size + 1;
          if(Y_mat(i,j) == 0 ){zero_count = zero_count + 1; }
        }
      }
      if( zero_count > gp_size - 2 ){
        // we have 0 or only 1 none-zero Y_ij for group k
        J_filter[j] = 0;
        Rprintf("taxon %i need to be removed due to at most 2 observations for samples in group %i \n",j+1, k);
      }
    }
  }
  // maximum number of features could be included in the model after filtering:
  int p_filter = sum(J_filter); 
  
  ////////////////Initialize the matrices to store the MCMC samples////////////////////
  // store the results for phi_j:
  NumericVector phi_store(p);
  NumericVector phi_tmp(p);
  int accpt_phi = 0;
  
  // store the results for mu_0j
  NumericVector mu0_store(p);
  NumericMatrix mu0_mat(S, p);
  NumericVector mu0_tmp(p); 
  for(int j = 0; j < p; j++){
    if(J_filter[ j ] == 1){
      mu0_store[ j ] = mu0_tmp[ j ] = rnorm_trunc(mu0_mean, tau_mu0, 0, 40);
    }
  }
  mu0_mat(0, _) = mu0_store;
  int accpt_mu0 = 0;
  
  // store the results for gamma and mu_kj 
  NumericVector gamma_PPI(p), gamma_sum(S);
  NumericVector gamma_tmp(p);// could output the last one
  int row_muk = ceil(2 * S * burn_rate * 0.8 * p_filter);
  int col_muk = (K - 1)  + 3;
  int index_muk = 0;
  NumericMatrix mu_store(row_muk, col_muk); 
  NumericMatrix mu_mat_tmp(K, p);
  double omega = a_omega / (a_omega + b_omega);
  
  for(int j = 0; j < p; j++){
    if(J_filter[j] != 0)
    {
      gamma_tmp[j] = rbinom(1, 1, omega)[0];
      phi_tmp[j] = 10;
      //phi_tmp[j] = rgamma(1, a_phi, 1.0 / b_phi)[0];
    }
  }
  gamma_sum[0] = sum(gamma_tmp);
  for(int j = 0; j < p; j++){
    if(gamma_tmp[j] == 1){
      for(int k = 1; k < K; k++){ // k starts from 1 since we set all mu_0j = 0 as the reference group
        mu_mat_tmp(k, j) = rnorm(1, 0, tau_mukj)[0];
      }
    }
  }
  int accpt_mu_kj = 0, prop_mu_kj = 0;
  int accpt_gamma = 0, prop_gamma = 0;
  
  // store the results for R_ij:
  NumericMatrix R_sum(n, p), R_tmp(n, p), R_PPI(n, p);
  int rsize = n*p;
  for (int i = 0; i < rsize; i++) {
    if(Y_mat[i]!=0)
    { 
      R_sum[i] = R_tmp[i] = 0;
    }else{
      R_sum[i] = R_tmp[i] = 1;  
    }
  }
  
  // Clinical covariate index vector for variable selection:
  NumericVector R_vec(R);
  for(int r =0; r < R; r++){
    R_vec[r] = r + 1;
  }
  // update 5 of covariates in each between model step:
  int num_cov = 5;
  
  // store the results for Beta_ij:
  NumericMatrix beta_sum(R, p), delta_sum(R, p);
  NumericMatrix beta_tmp(R, p), delta_tmp(R, p);
  // the number of rows!!!
  int row_beta = 2 * S * burn_rate * p_filter * 0.8 * R;
  int col_beta = 4;
  int index_beta = 0;
  NumericMatrix beta_store(row_beta, col_beta), delta_store(R, p);
  double p_delta = a_p / (a_p + b_p); // set as the mean of the corresponding gamma distribution
  
  int msize = R * p;
  for(int i = 0; i < msize; i++){
    delta_tmp[i] = rbinom(1, 1, p_delta)[0];
    if(delta_tmp[i] != 0){
      beta_tmp[i] = rnorm(1, 0, tau_beta)[0];
    }
  }
  int accpt_beta = 0, prop_beta = 0;
  int accpt_delta = 0, prop_delta = 0; 
  
  ////////////////////////////////big loop (we start from s = 1) /////////////////////////////////////
  for (int s = 1; s < S; s++) {
    
    // Monitor the process
    if(s*100/S == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    // update phi_j vector:
    for(int j = 0; j < p; j++){
      // notice that we only need to estimate phi_j for J_filter[j] == 1
      if(J_filter[j] == 1) {
        double phi_old = phi_tmp[j];
        //propose a new phi_j*
        // double phi_star = rnorm_trunc(phi_old, tau_phi, phi_low, 1000);
        int count_2 = 0;
        double phi_star = 0;
        do {
          phi_star = rgamma(1, phi_old*phi_old/tau_phi, tau_phi/phi_old)(0);
          count_2++;
        } while (phi_star < phi_low && count_2 < 1000);
        if (count_2 == 1000)
        {
          phi_star = phi_low;
        }
        double log_mh_prior = (a_phi-1) * (log(phi_star) - log(phi_old)) + b_phi*(phi_old - phi_star);
        if(gamma_tmp[j] == 1){
          // case 1: discriminatory features:
          double log_mh = 0;
          for(int k = 1; k <= K; k++){
            for(int i = 0; i < n; i++){
              if(z_vec[i] == k && R_tmp(i, j)==0){
                double si = s_vec[i];
                double alpha_ij = 0;
                for(int r = 0; r < R; r++) {
                  alpha_ij = alpha_ij + X_mat(i, r) * beta_tmp(r, j);
                }
                alpha_ij = exp( mu0_tmp[j] + mu_mat_tmp(k-1, j) + alpha_ij );
                log_mh = log_mh + lgamma(Y_mat(i,j) + phi_star) + lgamma(phi_old) - lgamma(Y_mat(i,j) + phi_old) - lgamma(phi_star) + 
                  phi_star * (log(phi_star) - log(si * alpha_ij + phi_star)) - 
                  phi_old * (log(phi_old) - log(si * alpha_ij + phi_old) ) + 
                  Y_mat(i,j) * (log(si * alpha_ij + phi_old) - log(si * alpha_ij + phi_star));
              }
              
            }
          }
          log_mh = log_mh + log_mh_prior;
          if(log_mh > log(double(rand()%10001)/10000)){
            phi_tmp[j] = phi_star;
            accpt_phi = accpt_phi + 1;
            if(s >= S_burn) {
              phi_store[j] = phi_star + phi_store[j];
            }
          }else{
            if(s >= S_burn) {
              phi_store[j] = phi_old + phi_store[j];
            }
          }
          
        }else{
          // case 2: non-discriminatory features:
          double log_mh = 0;
          for(int i = 0; i < n; i++){
            if(R_tmp(i,j) == 0){
              double si = s_vec[i];
              double xb = 0;
              for(int r = 0; r < R; r++) {
                xb = xb + X_mat(i, r) * beta_tmp(r, j);
              }
              double alpha_ij = exp(mu0_tmp[j] + xb);
              log_mh = log_mh + lgamma(Y_mat(i,j) + phi_star) + lgamma(phi_old) - lgamma(Y_mat(i, j) + phi_old) - lgamma(phi_star) + 
                phi_star * (log(phi_star) - log(si * alpha_ij + phi_star)) - 
                phi_old * (log(phi_old) - log(si * alpha_ij + phi_old) ) + 
                Y_mat(i,j) * (log(si * alpha_ij + phi_old) - log(si * alpha_ij + phi_star));
            }
          }
          log_mh = log_mh + log_mh_prior;
          if(log_mh > log(double(rand()%10001)/10000)){
            phi_tmp[j] = phi_star;
            accpt_phi = accpt_phi + 1;
            if(s >= S_burn) {
              phi_store[j] = phi_star + phi_store[j];
            }
          }else{
            if(s >= S_burn) {
              phi_store[j] = phi_old + phi_store[j];
            }
          }
        }
      }
    }
    
    
    // Update mu_0j vector: 
    for(int j = 0; j < p; j++) {
      // notice that we only need to estimate mu_0j for J_filter[j] == 1
      if(J_filter[j] == 1) {
        double phi_j = phi_tmp[j];
        double mu0_old = mu0_tmp[j];
        // propose a new mu_0j* from N(mu_0_old,1) for each j:
        double mu0_star = rnorm_trunc(mu0_old, tau_mu0, 0, 40);
        // log(prior ratio):
        double log_mu0_prior = (mu0_old * mu0_old - mu0_star * mu0_star)/ 20;
         // (a_0 + 0.5) * (log(2 * b_0 + mu0_old * mu0_old) - log(2 * b_0 + mu0_star * mu0_star));
        // log(hasting ratio):
        double log_mh = 0;
        if(gamma_tmp[j] == 1){
          // discriminatory feature
          for(int k = 1; k <= K; k++) {
            double mu_kj = mu_mat_tmp(k - 1, j);
            for(int i = 0; i < n; i++) {
              if(z_vec[i] == k && R_tmp(i,j) == 0){
                double si = s_vec[i], y_ij = Y_mat(i,j);
                double xb = 0;
                for(int r = 0; r < R; r++) {
                  xb = xb + X_mat(i, r) * beta_tmp(r, j);
                }
                double exp_star = exp(mu0_star + mu_kj + xb); // mu0*
                double exp_old  = exp(mu0_old  + mu_kj + xb); // mu0_old
                log_mh = log_mh + phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + 
                  y_ij * (log(si * exp_old + phi_j) + log(exp_star) - 
                  log(s_vec[i] * exp_star + phi_j) - log(exp_old));     
              }
            }
          }
          // hasting ratio:
          log_mh = log_mh + log_mu0_prior;
          
          if(log_mh > log(double(rand()%10001)/10000)){
            mu0_tmp[j] = mu0_star;
            accpt_mu0 = accpt_mu0 + 1;
            mu0_mat(s, j) = mu0_star;
            if(s >= S_burn) {
              mu0_store[j] = mu0_star + mu0_store[j];
            }
          }else{
            mu0_mat(s, j) = mu0_mat(s - 1, j);
            if(s >= S_burn) {
              mu0_store[j] = mu0_old + mu0_store[j];
            }
          }
        }else{
          // non-discriminatory feature
          for(int i=0 ; i < n; i++){
            if(R_tmp(i,j) == 0){
              double si = s_vec[i], y_ij = Y_mat(i, j);
              double xb = 0;
              for(int r = 0; r < R; r++){
                xb = xb + X_mat(i, r) * beta_tmp(r, j);
              }
              double exp_star = exp(mu0_star + xb); // mu0*
              double exp_old  = exp(mu0_old + xb); // mu0_old
              log_mh = log_mh + phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + 
                y_ij * (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));       
            }
          }
          // hasting ratio:
          log_mh = log_mh + log_mu0_prior;
          if(log_mh > log(double(rand()%10001)/10000)){
            mu0_tmp[j] = mu0_star;
            accpt_mu0 = accpt_mu0 + 1;
            mu0_mat(s, j) = mu0_star;
            if(s >= S_burn) {
              mu0_store[j] = mu0_star + mu0_store[j];
            }
          }else{
            mu0_mat(s, j) = mu0_mat(s - 1, j);
            if(s >= S_burn) {
              mu0_store[j] = mu0_old + mu0_store[j];
            }
          }
        }
      }
      
    } // end of updating mu_0j 
    
    // Joint update mu_kj(k = 2,..., K) and gamma_j: 
    // (1) between model selection:
    for(int ss = 0; ss < num_feature; ss++){
      prop_mu_kj = prop_mu_kj + 1;
      prop_gamma = prop_gamma + 1;
      // get the current number of features included
      int num_ft_tmp = sum(gamma_tmp);
      ///// Variable selection /////
      // propose a feature to change 1 <-> 0 //:
      int j_cand = rand()%p;
      while(J_filter[j_cand] != 1) {
        j_cand = rand()%p;
      }
      // update current gamma vector:
      gamma_tmp[j_cand] = 1 - gamma_tmp[j_cand];
      // add - or - delete:
      if(gamma_tmp[j_cand] == 1) {
        // add:
        // propose mu_kj*
        NumericVector mu_star_vec = rnorm(K, 0, tau_mukj);
        // set 0 for the reference group:
        mu_star_vec[0] = 0; 
        
        // hasting ratio for between-model updating:
        double log_mh = 0, log_prior = 0, log_trans = 0;
        for(int k = 2; k <= K; k++){
          double mu_kj_star = mu_star_vec[k-1];
          // prior
          log_prior = log_prior - (a_mu + 0.5) * log(mu_kj_star * mu_kj_star / 2 + b_mu);
          
          // trans:
          log_trans = log_trans - R::dnorm(mu_kj_star, 0, tau_mukj, TRUE);
        }
        for(int k = 1; k <= K; k++){
          double mu_kj_star = mu_star_vec[k-1];
          // likelihood
          for(int i = 0; i < n; i++){
            if( z_vec[i] == k && R_tmp(i, j_cand) == 0){
              double mu_0j = mu0_tmp[j_cand], si = s_vec[i], phi_j = phi_tmp[j_cand], y_ij = Y_mat(i,j_cand);
              // calculate hasting ratio
              double xb = 0;
              for(int r = 0; r < R; r++) {
                xb = xb + X_mat(i, r) * beta_tmp(r, j_cand);
              }
              // calculate hasting ratio
              double exp_star = exp(mu_0j + mu_kj_star + xb); // mu_kj*
              double exp_old = exp(mu_0j + 0 + xb); // mu_kj_old
              log_mh = log_mh + phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + 
                y_ij * (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));  
            }
          }
        }
        // finish the calculate for log(prior)
        log_prior = log_prior + (K - 1) * (a_mu * log(b_mu) + lgamma(a_mu + 0.5) - 0.5 * log(2 * pi) - lgamma( a_mu ) ) + 
          lgamma(a_omega + num_ft_tmp + 1) + lgamma(p - num_ft_tmp - 1 + b_omega) - lgamma(a_omega + num_ft_tmp ) - lgamma(p - num_ft_tmp + b_omega);
        
        // get the hasting ratio:
        log_mh = log_mh + log_prior + log_trans;
        if(log_mh > log(double(rand()%10001)/10000))
        {
          accpt_mu_kj = accpt_mu_kj + 1;
          accpt_gamma = accpt_gamma + 1;
          
          // update mu_kj's 0 -> mu_star:
          for(int k = 1; k < K; k++){
            mu_mat_tmp(k, j_cand) = mu_star_vec[k];
          }
        }
        else{ gamma_tmp[j_cand] = 0; }
      } else {
        // delete
        double log_mh = 0, log_prior = 0, log_trans = 0;
        for(int k = 2; k <= K; k++){
          // prior
          double mu_kj_old = mu_mat_tmp(k - 1, j_cand);
          log_prior = log_prior + (a_mu + 0.5) * log(mu_kj_old * mu_kj_old / 2 + b_mu);
          
          // trans:
          log_trans = log_trans + R::dnorm(mu_kj_old, 0, tau_mukj, TRUE);
        }
        
        for(int k = 1; k <= K; k++){
          double mu_kj_old = mu_mat_tmp(k - 1, j_cand);
          
          // likelihood:
          for(int i = 0; i < n; i++){
            if( z_vec[i] == k && R_tmp(i, j_cand)==0){
              double mu_0j = mu0_tmp[j_cand], si = s_vec[i], phi_j = phi_tmp[j_cand], y_ij = Y_mat(i, j_cand);
              double xb = 0;
              for(int r = 0; r < R; r++) {
                xb = xb + X_mat(i, r) * beta_tmp(r, j_cand);
              }
              double exp_star = exp(mu_0j + 0 + xb); // mu_kj*
              double exp_old = exp(mu_0j + mu_kj_old + xb); // mu_kj_old
              log_mh = log_mh + phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + 
                y_ij * (log(si * exp_old + phi_j) +log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));  
            }
          }
        }
        
        log_prior = log_prior + (K - 1) * ( - a_mu * log(b_mu) - lgamma(a_mu + 0.5) + 0.5 * log(2 * pi) + lgamma(a_mu)) + 
          lgamma(a_omega + num_ft_tmp - 1) + lgamma(p - num_ft_tmp + 1 + b_omega) - lgamma(a_omega + num_ft_tmp) - lgamma(p - num_ft_tmp + b_omega);
        
        // get the hasting ratio:
        log_mh = log_mh + log_prior + log_trans;
        if(log_mh > log(double(rand()%10001)/10000))
        {
          accpt_mu_kj = accpt_mu_kj + 1;
          accpt_gamma = accpt_gamma + 1;
          
          // update mu_kj's mu_old -> 0:
          for(int k = 1; k < K; k++){
            mu_mat_tmp(k, j_cand) = 0;
          }
        } else {gamma_tmp[j_cand] = 1;}
      }
    }
    
    // add the current gamma vector to gamma_PPI after burn-in
    if( s >= S_burn ){
      gamma_PPI = gamma_PPI + gamma_tmp;
      for(int j = 0; j < p; j++){
        for(int k = 1; k < K; k++){
          if(mu_mat_tmp(k, j) > 0.3)
          {
            mu_store(index_muk, _) = NumericVector::create(s, j + 1, k + 1, mu_mat_tmp(k, j));
            index_muk = index_muk + 1;
          }
        }
      }
    }
    
    // (2) with-in model updating: mu_kj (k = 2,..., K) for all gamma_j != 0
    for(int j = 0; j < p; j++){
      if(gamma_tmp[j] == 1){
        for(int k = 2; k <= K; k++){
          prop_mu_kj = prop_mu_kj + 1;
          double mu_old = mu_mat_tmp(k - 1, j);
          double mu_prop = rnorm(1, mu_old, (tau_mukj / 2) )[0];  // random walk proposal
          
          // calculate hasting ratio:
          double lg_mh = 0;
          double lg_prior = (a_mu + 0.5) * (log(2 * b_mu + mu_old * mu_old ) - 
                             log(2 * b_mu + mu_prop * mu_prop));
          
          for(int i = 0; i < n; i++){
            if(z_vec[i] == k && R_tmp(i,j) == 0){
              double mu_0j = mu0_tmp[j], si = s_vec[i], phi_j = phi_tmp[j], y_ij = Y_mat(i,j);
              double xb = 0 ;
              for(int r = 0; r < R; r++) {
                xb = xb + X_mat(i, r) * beta_tmp(r, j);
              }
              double exp_star = exp(mu_0j + mu_prop + xb); // mu_kj*
              double exp_old = exp(mu_0j + mu_old + xb); // mu_kj_old
              lg_mh = lg_mh + phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + 
                y_ij * (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));  
            }
            lg_mh = lg_mh + lg_prior;
          }
          if(lg_mh > log(runif(1,0,1)[0]))
          {
            accpt_mu_kj = accpt_mu_kj + 1;
            mu_mat_tmp(k - 1, j) = mu_prop;
          }
        }
      }
    }
    // }
    if( s >= S_burn ){
      for(int j = 0; j < p; j++){
        for(int k = 1; k < K; k++){
          if(mu_mat_tmp(k, j) != 0)
          {
            mu_store(index_muk, _) = NumericVector::create(s, j + 1, k + 1, mu_mat_tmp(k, j));
            index_muk = index_muk + 1;
          }
        }
      }
    }
    gamma_sum[s] = sum(gamma_tmp);
    
    // Update R_ij matrix:
    for(int i = 0; i < n; i++){
      double si = s_vec[i];
      for(int j = 0; j < p; j++){
        // only need to update r_ij conditional on y_ij = 0
        if(Y_mat(i,j) == 0){
          double mu_0j = mu0_tmp[j], phi_j = phi_tmp[j];
          if(gamma_tmp[j] == 1){
            // discriminatory features
            for(int k = 1; k <= K; k++){
              if(z_vec[i]==k){
                double mu_kj = mu_mat_tmp(k-1, j);
                double xb = 0;
                for(int r = 0; r < R; r++) {
                  xb = xb + X_mat(i, r) * beta_tmp(r, j);
                }
                double alpha_ij = exp(mu_0j + mu_kj + xb);
                
                NumericVector lg_p_vec(2);
                // unormalized probability of r_ij = 1:
                lg_p_vec[0] = lgamma(a_pi + 1) + lgamma(b_pi) - lgamma(a_pi + 1 + b_pi);
                // unormalized probability of r_ij = 0:
                lg_p_vec[1] = phi_j * (log(phi_j) - log(si * alpha_ij + phi_j)) +
                  lgamma(a_pi) + lgamma(b_pi + 1) - lgamma(a_pi + 1 + b_pi);
                // normalize the probability:
                double max_tmp = max(lg_p_vec), sum_tmp = 0;
                for(int m = 0; m < 2; m++){
                  lg_p_vec[m] = lg_p_vec[m] - max_tmp;
                  lg_p_vec[m] = exp(lg_p_vec[m]);
                  sum_tmp = sum_tmp + lg_p_vec[m];
                }
                for(int m = 0; m < 2; m++){
                  lg_p_vec[m] = lg_p_vec[m] / sum_tmp;
                }
                int r_ij = rbinom(1, 1, lg_p_vec[0])[0];
                R_sum(i, j)  =  R_sum(i, j)+ r_ij;     
                R_tmp(i, j)  =  r_ij;    
                if( s >= S_burn ){
                  R_PPI(i, j) = R_PPI(i, j) + r_ij;
                }
              }
            }
          }else{
            // non-discriminatory features
            double xb = 0;
            for(int r = 0; r < R; r++) {
              xb = xb + X_mat(i, r) * beta_tmp(r, j);
            }
            double alpha_ij = exp(mu_0j + xb);
            
            NumericVector lg_p_vec(2);
            // unormalized probability of r_ij = 1:
            lg_p_vec[0] = lgamma(a_pi + 1) + lgamma(b_pi) - lgamma(a_pi + 1 + b_pi);
            // unormalized probability of r_ij = 0:
            lg_p_vec[1] = phi_j * (log(phi_j) - log(si * alpha_ij + phi_j)) +
              lgamma(a_pi) + lgamma(b_pi + 1) - lgamma(a_pi + 1 + b_pi);
            // normalize the probability:
            double max_tmp = max(lg_p_vec), sum_tmp = 0;
            for(int m = 0; m < 2; m++){
              lg_p_vec[m] = lg_p_vec[m] - max_tmp;
              lg_p_vec[m] = exp(lg_p_vec[m]);
              sum_tmp = sum_tmp + lg_p_vec[m];
            }
            for(int m = 0; m < 2; m++){
              lg_p_vec[m] = lg_p_vec[m] / sum_tmp;
            }
            int r_ij = rbinom(1,1,lg_p_vec[0])[0];
            R_sum(i, j)  =  R_sum(i, j) + r_ij;     
            R_tmp(i, j)  =  r_ij;
            if( s >= S_burn ){
              R_PPI(i, j) = R_PPI(i, j) + r_ij;
            }
          }
        }
      }
    }// end of updating R_ij matrix
    
    // Update Beta / Delta matrix (VS):
    for(int j = 0; j < p; j++ ){
      if(J_filter[j] != 0){
        IntegerVector delta_current(R);
        for(int r = 0; r < R; r++){
          delta_current[r] = delta_tmp(r, j);
        }
        // (1) between model update:
        for(int ss = 0; ss < num_cov; ss++){
          int num_cov_tmp = sum(delta_current);
          int r_cand = rand()%R;
          delta_current[r_cand] = 1 - delta_current[r_cand];
          if(delta_current[r_cand] == 1){
            // add:
            // propose beta_rj*
            double beta_star = rnorm_trunc(0, tau_beta, - beta_lim, beta_lim);
            //rnorm(1, 0, tau_beta)[0]; 
            
            // Joint update:
            double log_mh = 0;
            double phi_j = phi_tmp[j], mu_0j = mu0_tmp[j];
            double log_prior = lgamma(a_p + num_cov_tmp + 1) + lgamma(b_p + R - num_cov_tmp - 1) - lgamma(a_p + num_cov_tmp) - 
              lgamma(b_p + R - num_cov_tmp) + a_beta * log(b_beta) + lgamma(a_beta + 0.5) - 0.5 * log(2 * pi) - lgamma(a_beta) -
              (a_beta + 0.5) * log(beta_star * beta_star / 2 + b_beta);
            double log_trans = - R::dnorm(beta_star, 0, tau_beta, TRUE);
            
            // two cases(gamma_j = 0 or 1):
            if(gamma_tmp[j] == 1){
              // discriminatory feature
              for(int k = 0; k < K; k++){
                double mu_kj = mu_mat_tmp(k, j);
                for(int i = 0; i < n; i++){
                  if(z_vec[i] == k+1 && R_tmp(i,j) == 0){
                    double si = s_vec[i], y_ij = Y_mat(i, j), xb_old = 0;
                    for(int d = 0; d < R; d++){
                      xb_old = xb_old + X_mat(i, d) * beta_tmp(d, j);
                    }
                    double xb_star = xb_old + X_mat(i, r_cand) * beta_star;
                    double exp_old = exp(mu_0j + mu_kj + xb_old), exp_star = exp(mu_0j + mu_kj + xb_star);
                    log_mh = log_mh +  phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + y_ij * 
                      (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));
                  }
                }
              }
              // hasting ratio:
              log_mh = log_mh + log_prior + log_trans;
              
              if(log_mh > log(double(rand()%10001)/10000))
              {
                accpt_delta = accpt_delta + 1;
                accpt_beta = accpt_beta + 1;
                beta_tmp(r_cand, j) = beta_star;
              } else {
                delta_current[r_cand] = 0;}
            } else {
              // non discriminatory feature
              for(int i = 0; i < n; i++){
                if(R_tmp(i,j) == 0){
                  double si = s_vec[i], y_ij = Y_mat(i, j), xb_old = 0;
                  for(int d = 0; d < R; d++){
                    xb_old = xb_old + X_mat(i, d) * beta_tmp(d, j);
                  }
                  double xb_star = xb_old + X_mat(i, r_cand) * beta_star;
                  double exp_old = exp(mu_0j + xb_old), exp_star = exp(mu_0j + xb_star);
                  log_mh = log_mh +  phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + y_ij * 
                    (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));
                }
              }
              // hasting ratio:
              log_mh = log_mh + log_prior + log_trans;
              if(log_mh > log(double(rand()%10001)/10000))
              {
                accpt_delta = accpt_delta + 1;
                accpt_beta = accpt_beta + 1;
                beta_tmp(r_cand, j) = beta_star;
              } else {
                delta_current[r_cand] = 0;}
            }
          } else {
            // delete:
            prop_delta = prop_delta + 1;
            prop_beta = prop_beta + 1;
            // beta_rj* = 0
            double beta_old = beta_tmp(r_cand, j);
            // Joint update:
            double log_mh = 0;
            double phi_j = phi_tmp[j], mu_0j = mu0_tmp[j];
            
            double log_prior = - a_beta * log(b_beta) - lgamma(a_beta + 0.5) + 0.5 * log(2 * pi) + lgamma(a_beta) +
              (a_beta + 0.5) * log(beta_old * beta_old / 2 + b_beta) + lgamma(a_p + num_cov_tmp - 1) +
              lgamma(b_p + R - num_cov_tmp + 1) - lgamma(a_p + num_cov_tmp) - lgamma(b_p + R - num_cov_tmp);
            double log_trans = R::dnorm(beta_old, 0, tau_beta, TRUE);
            
            // two cases(gamma_j = 0 or 1):
            if(gamma_tmp[j] == 1){
              // discriminatory feature
              for(int k = 0; k < K; k++){
                double mu_kj = mu_mat_tmp(k, j);
                for(int i = 0; i < n; i++){
                  if(z_vec[i] == k+1 && R_tmp(i,j) == 0){
                    double si = s_vec[i], y_ij = Y_mat(i, j), xb_old = 0;
                    for(int d = 0; d < R; d++){
                      xb_old = xb_old + X_mat(i, d) * beta_tmp(d, j);
                    }
                    double xb_star = xb_old - X_mat(i, r_cand) * beta_old;
                    double exp_old = exp(mu_0j + mu_kj + xb_old), exp_star = exp(mu_0j + mu_kj + xb_star);
                    log_mh = log_mh +  phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + y_ij * 
                      (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));
                  }
                }
              }
              // hasting ratio:
              log_mh = log_mh + log_prior + log_trans;
              
              if(log_mh > log(double(rand()%10001)/10000))
              {
                accpt_delta = accpt_delta + 1;
                accpt_beta = accpt_beta + 1;
                beta_tmp(r_cand, j) = 0;
              } else { delta_current[r_cand] = 1;}
            } else {
              // non discriminatory feature
              for(int i = 0; i < n; i++){
                if(R_tmp(i,j) == 0){
                  double si = s_vec[i], y_ij = Y_mat(i, j), xb_old = 0;
                  for(int d = 0; d < R; d++){
                    xb_old = xb_old + X_mat(i, d) * beta_tmp(d, j);
                  }
                  double xb_star = xb_old - X_mat(i, r_cand) * beta_old;
                  double exp_old = exp(mu_0j + xb_old), exp_star = exp(mu_0j + xb_star);
                  log_mh = log_mh +  phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + y_ij * 
                    (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));
                }
              }
              // hasting ratio:
              log_mh = log_mh + log_prior + log_trans;
              if(log_mh > log(double(rand()%10001)/10000))
              {
                accpt_delta = accpt_delta + 1;
                accpt_beta = accpt_beta + 1;
                beta_tmp(r_cand, j) = 0;
              } else { delta_current[r_cand] = 1; }
            }
          }
        }
        // Finish the between-model step
        delta_tmp(_, j) = delta_current;
        // store the results after burn-in:
        if(s >= S_burn) {
          for(int r = 0; r < R; r++){
            delta_sum(r, j) = delta_sum(r, j) + delta_tmp(r, j);
            beta_sum(r, j) = beta_tmp(r, j) + beta_sum(r, j); 
            delta_store(r,j) = delta_store(r,j) + delta_tmp(r,j);
            if(delta_tmp(r, j) != 0){
              beta_store(index_beta, _) = NumericVector::create(s, j + 1, r + 1, beta_tmp(r, j));
              index_beta = index_beta + 1;
            }
          } 
        }
        
        // (2) With-in model update: 
        for(int r = 0; r < R; r++){
          if(delta_tmp(r, j) == 1){
            prop_beta = prop_beta + 1;
            // only update for significant clinical features
            double beta_old = beta_tmp(r, j);
            double beta_new = rnorm_trunc(beta_old, tau_beta * 0.5, - beta_lim, beta_lim);
            //rnorm(1, beta_old, (tau_beta / 1.0))[0]; // propose beta* 
            double lg_mh = 0;
            double lg_prior = (a_beta + 0.5) * (log(beta_old * beta_old + 2 * b_beta) - log(beta_new * beta_new + 2 * b_beta));      
            for(int i = 0; i < n; i++){
              double si = s_vec[i], phi_j = phi_tmp[j], mu_0j = mu0_tmp[j];
              double y_ij = Y_mat(i, j), xb_old = 0;
              if(gamma_tmp[j] == 1 && R_tmp(i,j) == 0){
                // discriminantory features
                for(int k = 0; k < K; k++){
                  if(z_vec[i] == k + 1){
                    double mu_kj = mu_mat_tmp(k, j);
                    for(int d = 0; d < R; d++){
                      xb_old = xb_old + X_mat(i, d) * beta_tmp(d, j);
                    }
                    double xb_star = xb_old - X_mat(i, r) * beta_tmp(r, j) + X_mat(i, r) * beta_new;
                    double exp_star = exp(mu_0j + mu_kj + xb_star), exp_old = exp(mu_0j + mu_kj + xb_old);
                    lg_mh = lg_mh + phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + y_ij * 
                      (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));
                  }
                }
              } else if(gamma_tmp[j] == 0 && R_tmp(i,j) == 0){
                // nondiscriminantory features
                for(int d = 0; d < R; d++){
                  xb_old = xb_old + X_mat(i, d) * beta_tmp(d, j);
                }
                double xb_star = xb_old - X_mat(i, r) * beta_tmp(r, j) + X_mat(i, r) * beta_new;
                double exp_star = exp(mu_0j + xb_star), exp_old = exp(mu_0j + xb_old);
                lg_mh = lg_mh + phi_j * (log(si * exp_old + phi_j) - log(si * exp_star + phi_j)) + y_ij * 
                  (log(si * exp_old + phi_j) + log(exp_star) - log(si * exp_star + phi_j) - log(exp_old));
              }
            }
            lg_mh = lg_mh + lg_prior;
            if(lg_mh > log(double(rand()%10001)/10000)){
              beta_tmp(r, j) = beta_new;
              accpt_beta = accpt_beta + 1;
            }
          }
        }
        
        // store the results after burn-in:
        if(s >= S_burn) {
          for(int r = 0; r < R; r++){
            beta_sum(r, j) = beta_tmp(r, j) + beta_sum(r, j); 
            if(beta_tmp(r, j) != 0){
              beta_store(index_beta, _) = NumericVector::create(s, j + 1, r + 1, beta_tmp(r, j));
              index_beta = index_beta + 1;
            }
          } 
        }
      }
    }
    
    
    
  } // end of big loop
  List result;
  // MCMC matrices:
  result["mukj store"] = mu_store; //result["mukj last"] = mu_mat_tmp;
  result["mu0 est"] = mu0_store / (S * (1.0 - burn_rate)); 
  result["mu0 store"] = mu0_mat;
  result["gamma sum store"] = gamma_sum; 
  result["gamma PPI"] = gamma_PPI / (S * (1.0 - burn_rate));
  //result["gamma last"] = gamma_tmp;
  result["phi est"] = phi_store / (S * (1.0 - burn_rate));
  // result["R sum"] = R_sum; 
  result["R PPI"] = R_PPI / (S * (1.0 - burn_rate));
  // result["R last"] = R_tmp;
  result["beta store"] = beta_store; 
  result["beta est"] = beta_sum / (2 * S * (1.0 - burn_rate)); 
  // result["Beta last"] = beta_tmp;
  // result["delta last"] = delta_tmp;
  result["delta PPI"] = delta_store /(S * (1.0 - burn_rate)); 
  
  // all the acceptance rate:
  // result["gamma accept rate"] = accpt_gamma / (double)prop_gamma;
  // result["mukj accept rate "] = accpt_mu_kj / (double)prop_mu_kj;
  // result["mu0 accept rate "] = accpt_mu0 / (double)(S * p_filter);
  // result["phi accept rate "] = accpt_phi / (double)(S * p_filter);
  // result["beta accept rate "] = accpt_beta / (double)prop_beta;
  // result["delta accept rate "] = accpt_delta / (double)prop_delta;
  return result;
}




// [[Rcpp::export]]
double rnorm_trunc(double mu, double sigma, double lower, double upper)
{
  int change;
  double a, b;
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double z, tmp, lograt;
  
  change = 0;
  a = (lower - mu)/sigma;
  b = (upper - mu)/sigma;
  
  // First scenario
  if( (a == R_NegInf)||(b == R_PosInf))
  {
    if(a == R_NegInf)
    {
      change = 1;
      a = -b;
      b = R_PosInf;
    }
    // The two possibilities for this scenario
    if(a <= 0.45) z = norm_rs(a, b);
    else z = exp_rs(a, b);
    if(change) z = -z;
  }
  
  // Second scenario
  else if((a*b) <= 0.0)
  {
    // The two possibilities for this scenario
    if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
    {
      z = norm_rs(a, b);
    }
    else z = unif_rs(a,b);
  }
  
  // Third scenario
  else
  {
    if(b < 0)
    {
      tmp = b; b = -a; a = -tmp; change = 1;
    }
    
    lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
    if(lograt <= logt2)
    {
      z = unif_rs(a,b);
    }
    else if((lograt > logt1)&&(a < t3))
    {
      z = half_norm_rs(a,b);
    }
    else
    {
      z = exp_rs(a,b);
    }
    if(change)
    {
      z = -z;
    }
  }
  double output;
  output = sigma*z + mu;
  return (output);
}

// [[Rcpp::export]]
double exp_rs(double a, double b)
{
  double  z, u, rate;
  rate = 1/a;
  
  // Generate a proposal on (0, b-a)
  z = R::rexp(rate);
  while(z > (b-a))
  {
    z = R::rexp(rate);
  }
  u = R::runif(0.0, 1.0);
  
  while( log(u) > (-0.5*z*z))
  {
    z = R::rexp(rate);
    while(z > (b-a))
    {
      z = R::rexp(rate);
    }
    u = R::runif(0.0,1.0);
  }
  return(z+a);
}

// [[Rcpp::export]]
double unif_rs(double a, double b)
{
  double xstar, logphixstar, x, logu;
  
  // Find the argmax (b is always >= 0)
  // This works because we want to sample from N(0,1)
  if(a <= 0.0)
  {
    xstar = 0.0;
  }
  else
  {
    xstar = a;
  }
  logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);
  
  x = R::runif(a, b);
  logu = log(R::runif(0.0, 1.0));
  while(logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
  {
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
  }
  return x;
}

// [[Rcpp::export]]
double half_norm_rs(double a, double b)
{
  double x;
  x = fabs(norm_rand());
  while((x<a)||(x>b))
  {
    x = fabs(norm_rand());
  }
  return x;
}

// [[Rcpp::export]]
double norm_rs(double a, double b)
{
  double x;
  x = Rf_rnorm(0.0, 1.0);
  while((x < a)||(x > b))
  {
    x = norm_rand();
  }
  return x;
}

