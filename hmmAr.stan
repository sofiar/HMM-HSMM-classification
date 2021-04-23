// HSMM with no autorregresive structure

data {
  int<lower=1> K;  // num categories
  int<lower=0> Nrep;   // num of ts
  int<lower=0> T[Nrep];  // num instances
  int<lower=0> TT;  // num total instances
  int<lower=1,upper=K> z[TT]; // behaviours
  
  matrix[TT,3] y;
  matrix[K,K] counts;
  
  int inxIn[Nrep];
  int inxEn[Nrep];

  
}

parameters {
  //sojourn times parameters
  
  
  //autorregresive models.
  real alphas[3,K];
  real betas1[3,K];
  
  
  real<lower=0> sigmas[3,K];
  // K x K-1 tpm
  simplex[K] theta[K]; 
  
}

model {
  
  real lp;
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
  
  //betas
  for (i in 1:3)
  {
  betas1[i] ~ normal(0,2);
  }

  

  //Transition probability matrix
  for (k in 1:K)
  {theta[k] ~ dirichlet(counts[,k]);}
  
  //Observation process
  for (l in 1:Nrep)
  {

 //First obs
    lp = log(0.25) +normal_lpdf(y[inxIn[l],1] | alphas[1,z[inxIn[l]]], sigmas[1,z[inxIn[l]]])+normal_lpdf(y[inxIn[l],2] | alphas[2,z[inxIn[l]]], sigmas[2,z[inxIn[l]]])+normal_lpdf(y[inxIn[l],3] | alphas[3,z[inxIn[l]]], sigmas[3,z[inxIn[l]]]);

    for (t in (inxIn[l]+1):(inxEn[l])) { // looping over all observations pdfs
      
      //Acc x 
      lp_px = normal_lpdf(y[t,1] |alphas[1,z[t]]+betas1[1,z[t]]*y[t-1,1], sigmas[1,z[t]]);
      //Acc y 
      lp_py =  normal_lpdf(y[t,2] |alphas[2,z[t]]+betas1[2,z[t]]*y[t-1,2], sigmas[2,z[t]]); 
      //Acc z 
      lp_pz =  normal_lpdf(y[t,3] |alphas[3,z[t]]+betas1[3,z[t]]*y[t-1,3], sigmas[3,z[t]]); 
      
      
      lp = lp+lp_px+lp_py+lp_pz;
      
    }
    

    
    target += lp;
  }
  
  
}

