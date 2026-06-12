  //----------------------------------------------------------------------------
  //
  //  Von Bertalanffy Growth Model -                                        
  //  Random site effects (MVN noncentered parametization)                  
  //  Fixed linear effects                                                  
  //
  //----------------------------------------------------------------------------

  // AUTHOR: William K. Annis
  // CREATED 4/9/2025


data {
  
  //----------------------------------------------------------------------------
  //  Input data
  //----------------------------------------------------------------------------
  
  int<lower=1>                      N;  // number of fish
  int<lower=1>                N_SITES;  // number of sites
  int<lower=1>                      K;  // number of growth predictors
  int<lower=1>                 N_PRED;  // number of input values for predicted inst. growth
  real<lower=1>              LENGTH_M;  // average length of fish, used for growth estimation
  vector[N]                    LENGTH;  // fish lengths
  vector[N]                       AGE;  // fish ages
  int<lower=0>                     NU;  // degrees of freedom for students t errors. If zero, normal errors are estimated
  array[N] int<lower = 0>          ID;  // site/year id
  matrix[N_SITES,K]                 X;  // predictor data
  matrix[N_PRED,K]             PRED_X;  // prediction input data

}

parameters {
  
  //----------------------------------------------------------------------------
  //  Model parameters
  //----------------------------------------------------------------------------
  
  // Second level effects
  vector[K]  beta_Linf;  // Linf effects
  vector[K]    beta_g1;  // g1 effects
  vector[K]    beta_t0;  // t0 effects
  
  // Hyperparameters - means
  // Linf and g1 are log-transfomred to improve convergence 
  // and ensure postivevalues
  real                     mu_log_Linf;  // log-transformed Linf hyperparameter
  real                       mu_log_g1;  // log-transformed g1 hyperparameter
  real                           mu_t0;  // t0 hyperparameter

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
  
  vector[N_SITES]             Linf;  // site-level Linf's
  vector[N_SITES]               g1;  // site-level g1's
  vector[N_SITES]               t0;  // site-level t0's
  matrix[3,N_SITES]      beta_site;  // site-level error

  // Correlated site-level error
  beta_site = diag_pre_multiply(tau,L_omega)*alpha;
  
  // Hyperparameter plus site-level error and hydroperiod effect
  // Implies multi_normal(mu+X*beta, Sigma)
  Linf = exp(X * beta_Linf + mu_log_Linf+to_vector(beta_site[1,]));
  g1 = exp(X * beta_g1 + mu_log_g1+to_vector(beta_site[2,]));
  t0 = X * beta_t0 + mu_t0+to_vector(beta_site[3,]);
}


model {

  //----------------------------------------------------------------------------
  //  Model priors
  //----------------------------------------------------------------------------
  
  // Hydroperiod effect priors
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
  
  // growth equation
  vector[N] length_hat;  // Vector containing predicted lengths based on model 
  length_hat = Linf[ID] .* (1 - exp(-g1[ID] .*(AGE - t0[ID])));
  
  // Likelihood
  if(NU == 0) LENGTH ~ normal(length_hat,sigma_length);  // normal distribution
  if(NU > 0) LENGTH ~ student_t(NU,length_hat,sigma_length);  // student's T
  
}	

generated quantities{
  
  //----------------------------------------------------------------------------
  //  Generated qunatities
  //----------------------------------------------------------------------------
    
  real                  mu_Linf;  // transformed Linf hyperparameter
  real                    mu_g1;  // transformed g1 hyperparameter
  corr_matrix[3]        cor_mat;  // correlation matrix
  matrix[N_PRED,K]    pred_Linf;  // predicted Linf based on range of predictor values
  matrix[N_PRED,K]      pred_g1;  // predicted g1 based on range of predictor values
  matrix[N_PRED,K]      pred_t0;  // predicted t0 based on range of predictor values
  matrix[N_PRED,K]      pred_ig;  // predicted inst. growth based on range of predictor values
  vector[N]             log_lik;  // log-likelihood vector - needed for LOO and WAIC

  // Transform hyperparameter means out of log scale
  mu_Linf = exp(mu_log_Linf);
  mu_g1 = exp(mu_log_g1);

  // retrieve correlation matrix from cholesky matrix
  cor_mat = multiply_lower_tri_self_transpose(L_omega);
  
  // Log pointwise predictive density
  for(i in 1:N){
    if(NU==0) log_lik[i] = normal_lpdf(LENGTH[i]|Linf[ID[i]] .* (1 - exp(-g1[ID[i]] .*(AGE[i] - t0[ID[i]]))),sigma_length);
    if(NU>0) log_lik[i] = student_t_lpdf(LENGTH[i]|NU,Linf[ID[i]] .* (1 - exp(-g1[ID[i]] .*(AGE[i] - t0[ID[i]]))),sigma_length);
  }
  
  
  //----------------------------------------------------------------------------
  //  Predictions across range of predictor values
  //----------------------------------------------------------------------------

  // predicted growth paramters
  for (i in 1:K) {
    pred_Linf[,i] = exp(mu_log_Linf + beta_Linf[i] .* PRED_X[,i]);
    pred_g1[,i] = exp(mu_log_g1 + beta_g1[i] .* PRED_X[,i]);
    pred_t0[,i] = mu_t0 + beta_t0[i] .* PRED_X[,i];
  }

  // Instantenous growth equation
  for (i in 1:K) {
    pred_ig[,i] = pred_g1[,i] .* (pred_Linf[,i]-LENGTH_M);
  }
 }
 
 
