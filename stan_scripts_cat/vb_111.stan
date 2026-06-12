  /////////////////////////////////////////////////////////////////////////////
  //
  //  Von Bertalanffy Growth Model -                                        //
  //  Random site effects (MVN noncentered parametization)                  //
  //  Fixed Hydroperiod effects (Effect parameterization)                   //
  //
  ////////////////////////////////////////////////////////////////////////////

  // Global model:
  // all site effects (global)
  // all hydroperiod effects (global)

  // By William Annis - 12/15/2025

data {
  
  //////////////////
  /// Input data  //////////////////////////////////////////////////////////////
  /////////////////
  
  int<lower=1>                      N;  // number of fish
  int<lower=1>                N_SITES;  // number of sites
  int<lower=1>                 N_HYDR;  // number of hydroperiods
  row_vector[N]                LENGTH;  // fish lengths
  row_vector[N]                   AGE;  // fish ages
  array[N] int<lower = 0>          ID;  // site/year id
  array[N_SITES] int<lower = 0>  HYDR;  // hydroperiod

}
parameters {
  
  ///////////////////////
  /// Model parameters  ////////////////////////////////////////////////////////
  //////////////////////
  
  /// Hydroperiod effects ///
  sum_to_zero_vector[N_HYDR]  beta_Linf;  // Linf hydroperiod effect
  sum_to_zero_vector[N_HYDR]     beta_K;  // K hydroperiod effect
  sum_to_zero_vector[N_HYDR]    beta_t0;  // t0 hydroperiod effect
  
  /// Hyperparameters - means ///
  
  // Linf and K are log-transfomred to improve convergence 
  // and ensure postivevalues
  real                     mu_log_Linf;  // log-transformed Linf hyperparameter
  real                        mu_log_K;  // log-transformed K hyperparameter
  real                           mu_t0;  // t0 hyperparameter

  /// Hypermaters - Noncentered parametization error terms  ///
  cholesky_factor_corr[3]       L_omega;  // Cholesky transformed correlation matrix
  matrix[3,N_SITES]               alpha;  // site-error
  vector<lower=0>[3]                tau;  // error scaling term
  
  /// Fish-level error terms ///
  real<lower=0>            sigma_length;  // length error term

}	

transformed parameters {
  
  /////////////////////////////////////////////////////////////////
  /// Site-level MVN random effects - noncentered parametrization //////////////
  ////////////////////////////////////////////////////////////////
  
  row_vector[N_SITES]             Linf;  // site-level Linf's
  row_vector[N_SITES]                K;  // site-level K's
  row_vector[N_SITES]               t0;  // site-level t0's
  matrix[3,N_SITES]           beta_site;  // site-level error

  /// Correlated site-level error ///
  beta_site = diag_pre_multiply(tau,L_omega)*alpha;
  
  /// Site-specific parameters /// 
  
  // Hyperparameter plus site-level error and hydroperiod effect
  // Implies multi_normal(mu+beta_hydro, Sigma)
  Linf = exp(mu_log_Linf+beta_site[1]+to_row_vector(beta_Linf[HYDR]));
  K = exp(mu_log_K+beta_site[2]+to_row_vector(beta_K[HYDR]));
  t0 = mu_t0+beta_site[3]+to_row_vector(beta_t0[HYDR]);
}


model {
  
  ///////////////////
  /// Model priors /////////////////////////////////////////////////////////////
  //////////////////
  
  /// Hydroperiod effect priors ///
  beta_Linf ~ normal(0,10);
  beta_K  ~ normal(0,10);
  beta_t0 ~ normal(0,10);
  
  /// Hyperpriors - means  ///
  mu_log_Linf ~ normal(0,10);
  mu_log_K ~ normal(0,10);
  mu_t0 ~ normal(0,10);

  /// Hyperpriors - error  ///
  L_omega ~ lkj_corr_cholesky(2);  // cholesky correaltion
  to_vector(alpha) ~ normal(0,1);  // site error
  tau ~ student_t(3, 0, 2.5); // error scaling term
  
  /// Error prior ///
  sigma_length ~ student_t(3, 0, 2.5);
  
  ///////////////////////
  /// Model Likelihood  ////////////////////////////////////////////////////////
  ///////////////////////
  
  /// VBLF growth equation ///
  row_vector[N] length_hat;  // Vector containing predicted lengths based on model 
  length_hat = Linf[ID] .* (1 - exp(-K[ID] .*(AGE - t0[ID])));
  LENGTH ~ normal(length_hat,sigma_length);  // Likelihood
}	

generated quantities{
  
  ///////////////////////////
  /// Generated qunatities  ////////////////////////////////////////////////////
  ///////////////////////////
    
  row_vector[N_HYDR] hydr_Linf;  // hydroperiod specific Linf
  row_vector[N_HYDR]    hydr_K;  // hydroperiod specific K
  row_vector[N_HYDR]   hydr_t0;  // hydroperiod specific t0
  real                 mu_Linf;  // transformed Linf hyperparameter
  real                    mu_K;  // transformed K hyperparameter
  corr_matrix[3]       cor_mat;  // correlation matrix
  row_vector[N]        log_lik;  // log-likelihood vector - needed for LOO and WAIC
  
  /// Hydroperiod-specific hyperparameter means ///
  hydr_Linf = exp(mu_log_Linf+to_row_vector(beta_Linf));
  hydr_K = exp(mu_log_K+to_row_vector(beta_K));
  hydr_t0 = mu_t0+to_row_vector(beta_t0);
  
  /// Transform hyperparameter means out of log scale  ///
  mu_Linf = exp(mu_log_Linf);
  mu_K = exp(mu_log_K);
  
  /// retrieve correlation matrix from cholesky matrix  ///
  cor_mat = multiply_lower_tri_self_transpose(L_omega);
  
  /// Log pointwise predictive density  ///
  for(i in 1:N){
    log_lik[i] = normal_lpdf(LENGTH[i]|Linf[ID[i]] .* (1 - exp(-K[ID[i]] .*(AGE[i] - t0[ID[i]]))),sigma_length);
  }
}

