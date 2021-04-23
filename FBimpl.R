##### Inicializo las variables

pi=delta
pdf=FB.pdf
dim(pdf) <- c(tau,J)
pdf=t(pdf)

d=(FB.d)
dim(d) <- c(M,J)
d=t(d)
p= (tpm)


F = matrix(NA, J , tau)
L = matrix(NA, J , tau)
G = matrix(0, J , tau)
L1 = matrix(NA, J , tau)
N = rep(0, times = tau)
Norm = matrix(NA, J , tau)
StateIn = matrix(NA, J , tau)
## Store D
D=matrix(NA,ncol=M,nrow=J)
for (j in 1:J) {
  for (u in 1:(M)) {
  D[j,u] = sum(d[j,(u):(M)])
  }
  u = M+1 
  while( u <= tau)
  {
    D[j,u] = 0
    u=u+1
  }
  }



### Iteration


## Recursion Forward
for (t in 1:tau)
{
  N[t] = 0
  for (j  in 1:J)
  {
    if (t == 1)
    {
      Norm[j,1] = pdf[j,1] * pi[j]
    }
    else
    {
      Norm[j,t] = pdf[j,t] * (StateIn[j,t] - F[j,t - 1] + Norm[j,t - 1])
    }
    N[t] = N[t]+ Norm[j,t]
  }
  if (N[t] <= 0)
  {
   print('N<0!')
  }
  
  for (j in 1:J)
  {
    Norm[j,t] =Norm[j,t]/N[t];
  }
  
  
  for (j in 1:J)
  {
    F[j,t] = 0
    Observ = 1
    
    if (t < tau )
    {
      for (u in 1:min(t, M)) 
      {							
        Observ =Observ * pdf[j,(t - u+1)] / N[t - u+1]
        
        if (u < t) 
        {
          F[j,t] =F[j,t] + Observ * d[j,u] * StateIn[j,t - u+1]
        }
        else
        {
          F[j,t] =F[j,t] + Observ * d[j,t] * pi[j]
        }
      }
    }
    
    else
    {
      for (u in 0:min(tau-1, M)) 
      {
        Observ =Observ * pdf[j,tau - u] / N[tau - u]
        if (u < tau-1) 
        {
          F[j,tau] =F[j,tau]+ Observ * D[j,(u+1)] * StateIn[j,tau - u]
        }
        # else
        # {
        #   F[j,tau] =F[j,tau] + Observ * D[j,tau] * pi[j]
        # }
      }
    }
    
    if (t < tau) 
    {
      for (j in 1:J) 
      {
        StateIn[j,t + 1] = 0;
        for (i in 1:J)
        {					
          StateIn[j,t + 1] =StateIn[j,t + 1]+ p[i,j] * F[i,t];
        }
      }
    }
  }				
  }


### Recursion Backward
for (j in 1:J)
{
  L[j,tau]=F[j,tau]
}


for (t in (tau-1):1) 
{
  H = array(NA,dim=c(J,tau,min(tau - t, M)))
  for (j in 1:J)
  {
    G[j,t + 1] = 0
    Observ = 1
    
    for (u in 1:min(tau - t, M)) 
    {
      Observ = Observ * pdf[j, t + u] / N[t + u]
      
       if (u < tau - t)
       {
         H[j,t + 1,u] = L1[j,t + u] * Observ * d[j, u] / F[j,t + u]
       }
      else
      {
        H[j,t + 1,u] = Observ * D[j,tau - t]
      }
    G[j,t + 1] = G[j,t + 1]+H[j,t + 1,u]
    }
  }

  for (j in 1:J) 
  {
    L1[j,t] = 0;
    for (k in 1:J) 
    {
      L1[j,t] = L1[j,t]+G[k,t + 1] * p[j,k]
    }
      L1[j,t] = L1[j,t] *F[j,t]
      L[j,t] = L1[j,t] + L[j,t + 1] - G[j,t + 1] * StateIn[j,t + 1]
  }  
  
}

### Prediccion 

pred=apply(L, 2, which.max)

mean(sim$path==pred)

