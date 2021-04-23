FBAcpp <- function(pi, pdf, d, p){
  ################################################################
  ## input parameters:
  ##  
  ## pi: initial probs p(C_1=j) (J)
  ## pdf: probabilities of the obs of each state  (JxT)
  ## d: probabilities of sojourn times (JxM)
  ## p: tansition probability matrix (JxJ)
  ################################################################
  
  # source cpp functions 
  Rcpp::sourceCpp('Forwrd.cpp')
  Rcpp::sourceCpp('Backwrd.cpp')  
  
  # init params, vectors and matrix
  J = length(pi)
  tau = dim(pdf)[2]
  M = dim(d)[2] 

  ## Store D
  D=matrix(NA,ncol=tau,nrow=J)
  for (j in 1:J) {
    for (u in 1:(M)) {
      D[j,u] = sum(d[j,(u):(M)])
      #D[j,u] = 1-sum(d[j,(1):(u-1)])
      }
    u = M+1 
    while( u <= tau)
    {
      D[j,u] = 0
      u=u+1
    }
  }  
  
  N = rep(0, times = tau)
  Norm = matrix(NA, J , tau)
  StateIn = numeric(J*tau)    
  FF=numeric(J*tau)
  PDF=as.vector(t(pdf))
  Dd=as.vector(t(d))
  DD=as.vector(t(D))
  

  ### Recursion Forward
  
  a=Forwrd(FF,N,Norm,StateIn,J,tau,M,delta,PDF,Dd,DD,p)
  Ford=a$F
  dim(Ford)= c(tau,J)
  F=t(Ford)
  ST=a$StateIn
  dim(ST)= c(tau, J)
  StateIn=t(ST)
  N=a$N
  
  
  ### Recursion Backward
  LL = rep(NA, J*tau)
  GG = rep(NA, J*tau)
  LL1 = rep(NA, J*tau)
  FF=a$F
  NN=a$N
  SStatein=a$StateIn
  
  b=Backwrd(GG, LL, LL1, FF, NN, SStatein,J,tau, M, PDF, 
            Dd, DD,p)
  
  G=b$G
  L=b$L
  dim(G)= c(tau,J) # ?Es asi
  dim(L)= c(tau,J) # ?Es asi
  G=t(G)
  L=t(L)
  
  
  out=list('Forwrd'=F,'Backwrd'=G,'Gamma'=L)
  return(out)
}