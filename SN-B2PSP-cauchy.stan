data {
  
  int<lower=1> n1;  // # of sampled units
  int<lower=1> n2;  // # of non-sample units

  int<lower=1> p;
  int<lower=1> K;  // # of truncated polynomial bases
  
  vector[n1] y; // survey variable of units in the sample

  row_vector[p] X[n1];
  row_vector[K] Z[n1];  // truncated polynomials
  
  row_vector[p] predX[n2];
  row_vector[K] predZ[n2];  // truncated polynomials
}

parameters {
  
  vector[p] beta;
  vector[K] b;
  real alpha;
  vector[p] lambda;
  vector[K] nu;
  real<lower=0> sigma[n1];
  real<lower=0> sigmab;
  real<lower=0> sigmanu;
  
}

transformed parameters {
  
  real xi[n1];
  real<lower=0> omega[n1];
  real SPL2[n1];
  
  for (i in 1:n1) {
    xi[i] = X[i] * beta + Z[i] * b;
    omega[i] = sqrt( (pow(alpha, 2) + 1) ) * sigma[i];
    SPL2[i] = X[i] * lambda + Z[i] * nu;
  }
  
}

model {
  
  for (i in 1:n1) {
    y[i] ~ skew_normal(xi[i], omega[i], alpha);
    sigma[i] ~ lognormal(SPL2[i], 0.1);
  }
  
  for (l in 1:p) {
    beta[l] ~ normal(0, 1e3);
    lambda[l] ~ normal(0, 1e3);
  }  
  
  for (l in 1:K) {
    b[l] ~ normal(0, sigmab); 
    nu[l] ~ normal(0, sigmanu);
  }
  
  sigmab ~ cauchy(0, 1);
  sigmanu ~ cauchy(0, 1);
  
  alpha ~ normal(0, 10) T[0, ];
}

generated quantities {

  // real yhat[n2];

  vector[n2] predxi;
  vector[n2] predSPL2;
  real<lower=0> predomega[n2];

  for (i in 1:n2) {
       
    predxi[i] = predX[i] * beta + predZ[i] * b;
    predSPL2[i] = predX[i] * lambda + predZ[i] * nu;
    predomega[i] = sqrt( (pow(alpha, 2) + 1) ) * exp(predSPL2[i]);    
    // yhat[i] = skew_normal_rng(predxi[i], predomega[i], alpha);
    
  }
}