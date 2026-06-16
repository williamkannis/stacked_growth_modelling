  //----------------------------------------------------------------------------
  //
  //  Von Bertalanffy Growth Model -                                       
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
  
  int<lower=1>                       N;  // number of fish
  int<lower=1>                 N_SITES;  // number of sites
  int<lower=1>                   N_CAT;  // number of groups
  row_vector[N]                 LENGTH;  // fish lengths
  row_vector[N]                    AGE;  // fish ages
  int<lower=0>                      NU;  // degrees of freedom for students t errors. If zero, normal errors are estimated
  array[N] int<lower = 0>           ID;  // site/year id
  array[N_SITES] int<lower = 0>    CAT;  // groupings

}
parameters {
  
  //----------------------------------------------------------------------------
  //  Model parameters
  //----------------------------------------------------------------------------
  
  // Group effects
  sum_to_zero_vector[N_CAT]  beta_Linf;  // Linf group effect
  sum_to_zero_vector[N_CAT]    beta_g1;  // g1 group effect
  sum_to_zero_vector[N_CAT]    beta_t0;  // t0 group effect
  
  // Hyperparameters - means
  // Linf and g1 are log-transfomred to improve convergence 
  // and ensure postivevalues
  real                     mu_log_Linf;  // log-transformed Linf hyperparameter
  real                       mu_log_g1;  // log-transformed g1 hyperparameter
  real                           mu_t0;  // t0 hyperparameter

  // Hyperparamaters - Noncentered parametization error terms
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
  row_vector[N_SITES]           site_g1;  // site-level g1's
  row_vector[N_SITES]           site_t0;  // site-level t0's
  matrix[3,N_SITES]           beta_site;  // site-level error

  // Correlated site-level error
  beta_site = diag_pre_multiply(tau,L_omega)*alpha;
  
  // Hyperparameter plus site-level error and grouping effect
  // Implies multi_normal(mu+beta_hydro, Sigma)
  site_Linf = exp(mu_log_Linf+beta_site[1]+to_row_vector(beta_Linf[CAT]));
  site_g1 = exp(mu_log_g1+beta_site[2]+to_row_vector(beta_g1[CAT]));
  site_t0 = mu_t0+beta_site[3]+to_row_vector(beta_t0[CAT]);
}


model {
  
  //----------------------------------------------------------------------------
  //  Model priors
  //----------------------------------------------------------------------------
  
  // Grouping effect priors
  beta_Linf ~ normal(0,10);
  beta_g1  ~ normal(0,10);
  beta_t0 ~ normal(0,10);
  
  // Hyperpriors - means
  mu_log_Linf ~ normal(0,10);
  mu_log_g1 ~ normal(0,10);
  mu_t0 ~ normal(0,10);

  // Hyperpriors - error
  L_omega ~ lkj_corr_cholesky(2);  // cholesky correaltion
  to_vector(alpha) ~ normal(0,1);  // site error
  tau ~ student_t(3, 0, 2.5); // error scaling term
  
  // Error prior
  sigma_length ~ student_t(3, 0, 2.5);
  
  //----------------------------------------------------------------------------
  //  Model Likelihood
  //----------------------------------------------------------------------------
  
  // Growth equation
  row_vector[N] length_hat;  // Vector containing predicted lengths based on model 
  length_hat = site_Linf[ID] .* (1 - exp(-site_g1[ID] .*(AGE - site_t0[ID])));
  
  // Likelihood
  if(NU == 0) LENGTH ~ normal(length_hat,sigma_length);  // normal distribution
  if(NU > 0) LENGTH ~ student_t(NU,length_hat,sigma_length);  // student's T
  
}	

generated quantities{
  
  //----------------------------------------------------------------------------
  //  Generated qunatities
  //----------------------------------------------------------------------------
    
  row_vector[N_CAT]     cat_Linf;  // group specific Linf
  row_vector[N_CAT]       cat_g1;  // group specific g1
  row_vector[N_CAT]       cat_t0;  // group specific t0
  real                   mu_Linf;  // transformed Linf hyperparameter
  real                     mu_g1;  // transformed g1 hyperparameter
  corr_matrix[3]         cor_mat;  // correlation matrix
  row_vector[N]          log_lik;  // log-likelihood vector - needed for LOO and WAIC
  
  // Group-specific hyperparameter means
  cat_Linf = exp(mu_log_Linf+to_row_vector(beta_Linf));
  cat_g1 = exp(mu_log_g1+to_row_vector(beta_g1));
  cat_t0 = mu_t0+to_row_vector(beta_t0);
  
  // Transform hyperparameter means out of log scale
  mu_Linf = exp(mu_log_Linf);
  mu_g1 = exp(mu_log_g1);
  
  // retrieve correlation matrix from cholesky matrix
  cor_mat = multiply_lower_tri_self_transpose(L_omega);
  
  // Log pointwise predictive density
  for(i in 1:N){
    if(NU==0) log_lik[i] = normal_lpdf(LENGTH[i]|site_Linf[ID[i]] .* (1 - exp(-site_g1[ID[i]] .*(AGE[i] - site_t0[ID[i]]))),sigma_length);
    if(NU>0) log_lik[i] = student_t_lpdf(LENGTH[i]|NU,site_Linf[ID[i]] .* (1 - exp(-site_g1[ID[i]] .*(AGE[i] - site_t0[ID[i]]))),sigma_length);
  }
}

