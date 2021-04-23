#include <Rcpp.h>
  using namespace Rcpp;

// [[Rcpp::export]]
DataFrame Backwrd(NumericVector GG, NumericVector LL, NumericVector LL1,
  NumericVector FF, NumericVector NN, NumericVector SStatein,
                 int JJ, int TAU, int MM,NumericVector PDF, 
                NumericVector Dd,NumericVector DDD,NumericMatrix TPM){
  NumericVector G(GG);
  NumericVector L (LL);
  NumericVector L1 (LL1); 
  NumericVector F(FF);
  NumericVector N(NN);
  NumericVector StateIn(SStatein);
  NumericVector pdf(PDF);
  NumericVector d(Dd);
  NumericVector D(DDD);
  NumericMatrix tpm(TPM);
  
  int J(JJ);
  int tau(TAU);
  int M(MM);
  double H;
  
  for (int j = 0; j <= J - 1; j++)
  {
    L[j*tau+tau - 1] = F[j*tau+tau - 1];
  }
  
  for (int t = tau - 2; t >=0; t--) 
  {
    for (int j = 0; j <= J - 1; j++)
    {
      G[j*tau+t + 1] = 0;
      
      double Observ = 1;
      for (int u = 1; u <=std::min(tau - 1 - t, M); u++) 
      {
        Observ *= pdf[j*tau+t + u]/N[t + u];
        //Rcout << Observ<< "\n";
        if (u < tau - 1 - t)
        {
          H = L1[j*tau+t + u] * Observ * d[j*M+u-1] / F[j*tau+t + u];
        //Rcout << d[j*M+u]<< "\n";
        }
        else
        {
          H = Observ * D[j*tau+tau -2- t];
        
        }
        
        G[j*tau+t + 1] += H;
        
      }
    }
    
    for (int j = 0; j <= J - 1; j++) 
    {
      L1[j*tau+t] = 0;
      for (int k = 0; k <= J - 1; k++) 
        L1[j*tau+t] += G[k*tau+t + 1] * tpm(j,k);
      L1[j*tau+t] *= F[j*tau+t];
      L[j*tau+t] = L1[j*tau+t] + L[j*tau+t + 1] - G[j*tau+t + 1] * StateIn[j*tau+t + 1];
    }
  }
  
  return DataFrame::create(Named("G")= G,Named("L")=L,Named("L1")=L1);
  
}