// bernoulli_logistic transformed data function
data {
  
  int<lower=1> N;               // rows of data
  
  int<lower=0> d_t[N];          // Number of mosquitoes dying IRS HUTS
  int<lower=0> n_t[N];          // Total number of mosquitoes in IRS huts
  int<lower=0> f_t[N];          // Number of mosquitoes feeding in HUTS
  
  int<lower=0> n_det[N];      // Number deterred by spray
  int<lower=0> n_c[N];        // Denomimator for deterrence
  vector<lower=0>[N] x;       // predictor

  int<lower=1> N_IRS;           // IRS treatments
  int<lower=1, upper=N_IRS> IRS[N];
  
  }

parameters {

  real alpha1_tilde[N_IRS];
  real mu_alpha1;
  real sigma_alpha1;

  real alpha2_tilde[N_IRS];
  real mu_alpha2;
  real sigma_alpha2;

  real beta1_tilde[N_IRS];
  real mu_beta1;  // 
  real sigma_beta1;//

  real beta2_tilde[N_IRS];
  real mu_beta2;  // 
  real sigma_beta2;//

  real theta2_tilde[N_IRS];
  real mu_theta2;  // 
  real sigma_theta2;//

  real theta1_tilde[N_IRS];
  real mu_theta1;  // 
  real sigma_theta1;//


}

transformed parameters {
    real alpha1[N_IRS];
    real alpha2[N_IRS];

    real beta1[N_IRS];
    real beta2[N_IRS];

    real theta1[N_IRS];
    real theta2[N_IRS];

  for (v in 1:N_IRS) {
    alpha1[v] = mu_alpha1 + sigma_alpha1 * alpha1_tilde[v];
    alpha2[v] = mu_alpha2 + sigma_alpha2 * alpha2_tilde[v];

    beta1[v] = mu_beta1 + sigma_beta1 * beta1_tilde[v];
    beta2[v] = mu_beta2 + sigma_beta2 * beta2_tilde[v];

    theta1[v] = mu_theta1 + sigma_theta1 * theta1_tilde[v];
    theta2[v] = mu_theta2 + sigma_theta2 * theta2_tilde[v];

  }
}

model {
  real y_tilde[N];
  real k_tilde[N];
  real det_tilde[N];
  
  for(n in 1:N)
    y_tilde[n] = alpha2[IRS[n]] * x[n] + alpha1[IRS[n]];
    
  alpha1_tilde ~ normal(0, 1);
  mu_alpha1 ~ normal(0, 10);  
  sigma_alpha1 ~ normal(0, 5);


  alpha2_tilde ~ normal(0, 1);
  mu_alpha2 ~ normal(0, 10);  
  sigma_alpha2 ~ normal(0, 5);

  d_t ~ binomial_logit(n_t, y_tilde);


  for(n in 1:N)
    k_tilde[n] = beta2[IRS[n]] * x[n] + beta1[IRS[n]];

  beta2_tilde ~ normal(0, 1);
  mu_beta2 ~ normal(0, 10);  
  sigma_beta2 ~ normal(0, 5);

  beta1_tilde ~ normal(0, 1);
  mu_beta1 ~ normal(0, 10);  
  sigma_beta1 ~ normal(0, 5);

  f_t ~ binomial_logit(n_t, k_tilde);


 for(n in 1:N)
    det_tilde[n] = theta2[IRS[n]] * x[n] + theta1[IRS[n]];

  theta2_tilde ~ normal(0, 1);
  mu_theta2 ~ normal(0, 10);  
  sigma_theta2 ~ normal(0, 5);

  theta1_tilde ~ normal(0, 1);
  mu_theta1 ~ normal(0, 10);  
  sigma_theta1 ~ normal(0, 5);

  n_det ~ binomial_logit(n_c, det_tilde);

  }

generated quantities {
  real y_ppc[N_IRS, 100];
  real k_ppc[N_IRS, 100];
  real det_ppc[N_IRS, 100];

  for(v in 1:N_IRS) {
    for (t in 1:100) {
      y_ppc[v, t] = binomial_rng(100, inv_logit(alpha2[v] * t + alpha1[v])) / 100.0;
      k_ppc[v, t] = binomial_rng(100, inv_logit(beta2[v] * t + beta1[v])) / 100.0;
      det_ppc[v, t] = binomial_rng(100, inv_logit(theta2[v] * t + theta1[v])) / 100.0;
    }
  }
}
