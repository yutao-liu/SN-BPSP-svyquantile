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
 
  real<lower=0> Pi[n1];
  real<lower=0> predPi[n2];
  
}

parameters {
  
  vector[p] beta;
  vector[K] b;
  real alpha;
  real<lower=0> sigma;
  real gamma;
  real<lower=0> sigmab;
  
}

transformed parameters {
  
  real xi[n1];
  real<lower=0> omega[n1];
  
  for (i in 1:n1) {
    xi[i] = X[i] * beta + Z[i] * b;
    omega[i] = sqrt( (pow(alpha, 2) + 1) ) * pow(Pi[i], gamma) * sigma;
  }
}

model {
  
  for (i in 1:n1) {
    y[i] ~ skew_normal(xi[i], omega[i], alpha);
  }
  
  for (l in 1:p) {
    beta[l] ~ normal(0, 1e3);
  }  
  
  for (l in 1:K) {
    b[l] ~ normal(0, sigmab); 
  }
  
  sigmab ~ cauchy(0, 1);
  sigma ~ cauchy(0, 1);
  
  alpha ~ normal(0, 10) T[0, ];
}

generated quantities{
  
  real yhat[n2];

  vector[n2] predxi;
  real<lower=0> predomega[n2];
  
  for (i in 1:n2) {
    
    predxi[i] = predX[i] * beta + predZ[i] * b;
    predomega[i] = sqrt( (pow(alpha, 2) + 1) ) * pow(predPi[i], gamma) * sigma ;
    // yhat[i] = skew_normal_rng(predxi[i], predomega[i], alpha);
    
  }
}