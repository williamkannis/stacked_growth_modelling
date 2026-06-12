  /////////////////////////////////////////////////////////////////////////////
  //
  //  Logistic Growth Model -                                               //
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
  sum_to_zero_vector[N_HYDR] beta_gninf;  // gninf hydroperiod effect
  sum_to_zero_vector[N_HYDR]    beta_ti;  // ti hydroperiod effect
  
  /// Hyperparameters - means ///
  
  // Linf and gninf are log-transfomred to improve convergence 
  // and ensure postivevalues
  real                     mu_log_Linf;  // log-transformed Linf hyperparameter
  real                    mu_log_gninf;  // log-transformed gninf hyperparameter
  real                       mu_log_ti;  // log-transformed ti hyperparameter

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
  row_vector[N_SITES]            gninf;  // site-level gninf's
  row_vector[N_SITES]               ti;  // site-level ti's
  matrix[3,N_SITES]           beta_site;  // site-level error

  /// Correlated site-level error ///
  beta_site = diag_pre_multiply(tau,L_omega)*alpha;
  
  /// Site-specific parameters /// 
  
  // Hyperparameter plus site-level error and hydroperiod effect
  // Implies multi_normal(mu+beta_hydro, Sigma)
  Linf = exp(mu_log_Linf+beta_site[1]+to_row_vector(beta_Linf[HYDR]));
  gninf = exp(mu_log_gninf+beta_site[2]+to_row_vector(beta_gninf[HYDR]));
  ti = exp(mu_log_ti+beta_site[3]+to_row_vector(beta_ti[HYDR]));
}


model {
  
  ///////////////////
  /// Model priors /////////////////////////////////////////////////////////////
  //////////////////
  
  /// Hydroperiod effect priors ///
  beta_Linf ~ normal(0,10);
  beta_gninf  ~ normal(0,10);
  beta_ti ~ normal(0,10);
  
  /// Hyperpriors - means  ///
  mu_log_Linf ~ normal(0,10);
  mu_log_gninf ~ normal(0,10);
  mu_log_ti ~ normal(0,10);

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
  length_hat = Linf[ID]./(1 + exp(-gninf[ID] .* (AGE - ti[ID])));
  LENGTH ~ normal(length_hat,sigma_length);  // Likelihood
}	

generated quantities{
  
  ///////////////////////////
  /// Generated qunatities  ////////////////////////////////////////////////////
  ///////////////////////////
    
  row_vector[N_HYDR]  hydr_Linf;  // hydroperiod specific Linf
  row_vector[N_HYDR] hydr_gninf;  // hydroperiod specific gninf
  row_vector[N_HYDR]    hydr_ti;  // hydroperiod specific ti
  real                  mu_Linf;  // transformed Linf hyperparameter
  real                 mu_gninf;  // transformed gninf hyperparameter
  real                    mu_ti;  // transformed ti hyperparameter
  corr_matrix[3]        cor_mat;  // correlation matrix
  row_vector[N]         log_lik;  // log-likelihood vector - needed for LOO and WAIC
  
  /// Hydroperiod-specific hyperparameter means ///
  hydr_Linf = exp(mu_log_Linf+to_row_vector(beta_Linf));
  hydr_gninf = exp(mu_log_gninf+to_row_vector(beta_gninf));
  hydr_ti = exp(mu_log_ti+to_row_vector(beta_ti));
  
  /// Transform hyperparameter means out of log scale  ///
  mu_Linf = exp(mu_log_Linf);
  mu_gninf = exp(mu_log_gninf);
  mu_ti = exp(mu_log_ti);
  
  /// retrieve correlation matrix from cholesky matrix  ///
  cor_mat = multiply_lower_tri_self_transpose(L_omega);
  
  /// Log pointwise predictive density  ///
  for(i in 1:N){
    log_lik[i] = normal_lpdf(LENGTH[i]|Linf[ID[i]]./(1 + exp(-gninf[ID[i]] .* (AGE[i] - ti[ID[i]]))),sigma_length);
  }
}

