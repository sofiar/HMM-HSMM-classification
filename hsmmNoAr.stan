// HSMM with no autorregresive structure

data {
  int<lower=1> K;  // num categories
  int<lower=0> Nrep;   // num of ts
  int<lower=0> N[Nrep];   // num of NAS
  int<lower=0> NN;   // num of total NAS 
  int<lower=0> T[Nrep];  // num instances
  int<lower=0> TT;  // num total instances
  int<lower=1,upper=K> z[TT]; // behaviours
  int<lower=0> u[NN]; // sojuorn times
  int<lower=1,upper=K> Sest[NN]; // seq behaviours
  
  matrix[TT,3] y;
  matrix[K-1,K] counts;
  
  int inxIn[Nrep];
  int inxEn[Nrep];
  int nasIn[Nrep];
  //int nasEn[Nrep];
  
}

parameters {
  //sojourn times parameters
  
  real <lower=0> phi[K];
  real <lower=0> mu[K];
  
  //autorregresive models.
  real alphas[3,K];
  
  real<lower=0> sigmas[3,K];
  // K x K-1 tpm
  simplex[K-1] theta[K]; 
  
}

model {
  
  real lp;
  real lp_p1;
  real lp_px;
  real lp_py;
  real lp_pz;
  
  //////// PRIORS ////////
    
  //VedBA
  alphas[1] ~ normal(0, 3);  
  sigmas[1]~ normal(0,3);
  
  //Mean Pitch angle
  alphas[2,1] ~ normal(20, 5);  
  alphas[2,2] ~ normal(20, 5);  
  alphas[2,3] ~ normal(20, 5);  
  alphas[2,4] ~ normal(-20, 5);  
  sigmas[2]~ normal(0,10);

  //Sd Pitch angle
  alphas[3] ~ normal(0, 10);  
  sigmas[3]~ normal(0,3);
  
  
  //sojoun times
  mu[1]~normal(5, 5);
  mu[2]~normal(5, 5);
  mu[3]~normal(80, 30);
  mu[4]~normal(25, 20);
  
  
  phi[1]~normal(0, 5);
  phi[2]~normal(0, 5);
  phi[3]~normal(0, 5);
  phi[4]~normal(0, 5);
  
  
  //Transition probability matrix
  for (k in 1:K)
  {theta[k] ~ dirichlet(counts[,k]);}
  
  //Observation process
  for (l in 1:Nrep)
  {
    //First obs
    lp = log(0.25);
    
    for (t in (inxIn[l]):(inxEn[l])) { // looping over all observations pdfs
      
      //Acc x 
      lp_px = normal_lpdf(y[t,1] |alphas[1,z[t]], sigmas[1,z[t]]);
      //Acc y 
      lp_py =  normal_lpdf(y[t,2] |alphas[2,z[t]], sigmas[2,z[t]]); 
      //Acc z 
      lp_pz =  normal_lpdf(y[t,3] |alphas[3,z[t]], sigmas[3,z[t]]); 
      
      
      lp = lp+lp_px+lp_py+lp_pz;
      
    }
    
    //State Duration
    
    
    for (m in 1:N[l])
    {
      lp_p1 = neg_binomial_2_lpmf(u[nasIn[l]+m]|mu[Sest[nasIn[l]+m]],phi[Sest[nasIn[l]+m]]);
      lp = lp+lp_p1;  
    }
    target += lp;
  }
  
  
}

