//quiero modificarlo para tener diferentes longitudes de series temporales
// Binomial Negativa

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
  real betas1[3,K];
  real<lower=0> sigmas[3,K];
  // K x K-1 tpm
  simplex[K-1] theta[K]; 
  
  }

model {
  //vector[K] log_theta_tr[K];
  //real lpp;
  real lp;
  real lp_p1;
  real lp_px;
  real lp_py;
  real lp_pz;
  
  //////// PRIORS ////////
    
  //autorregresive parameters
 
  //VedBA
  alphas[1] ~ normal(0, 3);  
  
  //Pitch angle
  alphas[2,1] ~ normal(20, 5);  
  alphas[2,2] ~ normal(20, 5);  
  alphas[2,3] ~ normal(20, 5);  
  alphas[2,4] ~ normal(-20, 5);  
  // alphas[2,5] ~ normal(20, 5);  
  // alphas[2,6] ~ normal(-20, 5);  
  //alphas[2,7] ~ normal(-20, 5);  
  
  //
  alphas[3] ~ normal(0, 10);  
  
  for (i in 1:3)
  {
  betas1[i] ~ normal(0,5);
  }
  sigmas[1]~ normal(0,3);
  sigmas[2]~ normal(0,10);
  sigmas[3]~ normal(0,3);
  
  
  // for(i in 1:K)
  // {
  //   for(j in 1:3)
  //     sigmas[j,i] ~ normal(0,2);
  // }
  // 
  //sojoun times
  mu[1]~normal(5, 5);
  mu[2]~normal(5, 5);
  mu[3]~normal(80, 30);
  mu[4]~normal(25, 20);
  // mu[5]~normal(7, 5);
  // mu[6]~normal(10, 5);
  //mu[7]~normal(10, 5);
  
  phi[1]~normal(0, 5);
  phi[2]~normal(0, 5);
  phi[3]~normal(0, 5);
  phi[4]~normal(0, 5);
  // phi[5]~normal(0,1);
  // phi[6]~normal(0, 1);
 // phi[7]~normal(0, 1);

  //Transition probability matrix
    for (k in 1:K)
    {theta[k] ~ dirichlet(counts[,k]);}

  //Observation process
  for (l in 1:Nrep)
  {
    //First obs
    lp = log(0.25) +normal_lpdf(y[inxIn[l],1] | alphas[1,z[inxIn[l]]], sigmas[1,z[inxIn[l]]])+normal_lpdf(y[inxIn[l],2] | alphas[2,z[inxIn[l]]], sigmas[2,z[inxIn[l]]])+normal_lpdf(y[inxIn[l],3] | alphas[3,z[inxIn[l]]], sigmas[3,z[inxIn[l]]]);
    
    for (t in (inxIn[l]+1):(inxEn[l])) { // looping over all observations pdfs
      
      //Acc x --> AR(1)
      lp_px = normal_lpdf(y[t,1] |alphas[1,z[t]]+betas1[1,z[t]]*y[t-1,1], sigmas[1,z[t]]);
      //Acc y --> AR(1)
      lp_py =  normal_lpdf(y[t,2] |alphas[2,z[t]]+betas1[2,z[t]]*y[(t-1),2], sigmas[2,z[t]]); 
      //Acc z --> AR(1)
      lp_pz =  normal_lpdf(y[t,3] |alphas[3,z[t]]+betas1[3,z[t]]*y[(t-1),3], sigmas[3,z[t]]); 
      
      
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

