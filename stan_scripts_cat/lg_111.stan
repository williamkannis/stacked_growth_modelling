  //----------------------------------------------------------------------------
  //
  //  Logistic Growth Model -                                       
  //  Random site effects (MVN noncentered parametization)                  
  //  Fixed categorical effects  (Effect parameterization)                   
  //
  //----------------------------------------------------------------------------

  // AUTHOR: William K. Annis
  // CREATED 12/15/2025

data {
  
  //----------------------------------------------------------------------------
  //  Input data
  //----------------------------------------------------------------------------
  
  int<lower=1>                      N;  // number of fish
  int<lower=1>                N_SITES;  // number of sites
  int<lower=1>                  N_CAT;  // number of groups
  row_vector[N]                LENGTH;  // fish lengths
  row_vector[N]                   AGE;  // fish ages
  int<lower=0>                     NU;  // degrees of freedom for students t errors. If zero, normal errors are estimated
  array[N] int<lower = 0>          ID;  // site/year id
  array[N_SITES] int<lower = 0>   CAT;  // groups

}
parameters {
  
  //----------------------------------------------------------------------------
  //  Model parameters
  //----------------------------------------------------------------------------
  
  // Hydroperiod effects
  sum_to_zero_vector[N_CAT]  beta_Linf;  // Linf group effect
  sum_to_zero_vector[N_CAT]    beta_g3;  // g3 group effect
  sum_to_zero_vector[N_CAT]    beta_ti;  // ti group effect
  
  // Hyperparameters - means
  // Linf and g3 are log-transfomred to improve convergence 
  // and ensure postivevalues
  real                     mu_log_Linf;  // log-transformed Linf hyperparameter
  real                       mu_log_g3;  // log-transformed g3 hyperparameter
  real                       mu_log_ti;  // log-transformed ti hyperparameter

  // Hyperparameters - Noncentered parametization error terms
  cholesky_factor_corr[3]       L_omega;  // Cholesky transformed correlation matrix
  matrix[3,N_SITES]               alpha;  // site-error
  vector<lower=0>[3]                tau;  // error scaling term
  
  // Fish-level error terms
  real<lower=0>            sigma_length;  // length error term

}	

transformed parameters {
  
  //----------------------------------------------------------------------------
  //  Site-level MVN random effects - noncentered parametrization
  //----------------------------------------------------------------------------
  
  row_vector[N_SITES]         site_Linf;  // site-level Linf's
  row_vector[N_SITES]           site_g3;  // site-level g3's
  row_vector[N_SITES]           site_ti;  // site-level ti's
  matrix[3,N_SITES]           beta_site;  // site-level error

  // Correlated site-level error
  beta_site = diag_pre_multiply(tau,L_omega)*alpha;
  
  // Hyperparameter plus site-level error and group effect
  // Implies multi_normal(mu+beta_cato, Sigma)
  site_Linf = exp(mu_log_Linf+beta_site[1]+to_row_vector(beta_Linf[CAT]));
  site_g3 = exp(mu_log_g3+beta_site[2]+to_row_vector(beta_g3[CAT]));
  site_ti = exp(mu_log_ti+beta_site[3]+to_row_vector(beta_ti[CAT]))-10;
}


model {
  
  //----------------------------------------------------------------------------
  //  Model priors
  //----------------------------------------------------------------------------
  
  // Hydroperiod effect priors
  beta_Linf ~ normal(0,10);
  beta_g3  ~ normal(0,10);
  beta_ti ~ normal(0,10);
  
  // Hyperpriors - means
  mu_log_Linf ~ normal(0,10);
  mu_log_g3 ~ normal(0,10);
  mu_log_ti ~ normal(0,10);

  // Hyperpriors - error
  L_omega ~ lkj_corr_cholesky(2);  // cholesky correaltion
  to_vector(alpha) ~ normal(0,1);  // site error
  tau ~ student_t(3, 0, 2.5); // error scaling term
  
  // Error prior
  sigma_length ~ student_t(3, 0, 2.5);
  
  //----------------------------------------------------------------------------
  //  Model Likelihood
  //----------------------------------------------------------------------------
  
  // growth equation
  row_vector[N] length_hat;  // Vector containing predicted lengths based on model 
  length_hat = site_Linf[ID]./(1 + exp(-site_g3[ID] .* (AGE - site_ti[ID])));
  
  // Likelihood
  if(NU == 0) LENGTH ~ normal(length_hat,sigma_length);  
  if(NU > 0) LENGTH ~ student_t(NU,length_hat,sigma_length); 
}	

generated quantities{
  
  //----------------------------------------------------------------------------
  //  Generated qunatities
  //----------------------------------------------------------------------------
    
  row_vector[N_CAT]    cat_Linf;  // group specific Linf
  row_vector[N_CAT]      cat_g3;  // group specific g3
  row_vector[N_CAT]      cat_ti;  // group specific ti
  real                  mu_Linf;  // transformed Linf hyperparameter
  real                    mu_g3;  // transformed g3 hyperparameter
  real                    mu_ti;  // transformed ti hyperparameter
  corr_matrix[3]        cor_mat;  // correlation matrix
  row_vector[N]         log_lik;  // log-likelihood vector - needed for LOO and WAIC
  
  // Hydroperiod-specific hyperparameter means
  cat_Linf = exp(mu_log_Linf+to_row_vector(beta_Linf));
  cat_g3 = exp(mu_log_g3+to_row_vector(beta_g3));
  cat_ti = exp(mu_log_ti+to_row_vector(beta_ti))-10;
  
  // Transform hyperparameter means out of log scale
  mu_Linf = exp(mu_log_Linf);
  mu_g3 = exp(mu_log_g3);
  mu_ti = exp(mu_log_ti)-10;
  
  // retrieve correlation matrix from cholesky matrix
  cor_mat = multiply_lower_tri_self_transpose(L_omega);
  
  // Log pointwise predictive density
  for(i in 1:N){
    if(NU==0) log_lik[i] = normal_lpdf(LENGTH[i]|site_Linf[ID[i]]./(1 + exp(-site_g3[ID[i]] .* (AGE[i] - site_ti[ID[i]]))),sigma_length);
    if(NU>0) log_lik[i] = student_t_lpdf(LENGTH[i]|NU,site_Linf[ID[i]]./(1 + exp(-site_g3[ID[i]] .* (AGE[i] - site_ti[ID[i]]))),sigma_length);
  }
}

