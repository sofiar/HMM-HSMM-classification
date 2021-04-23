#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
DataFrame Forwrd(NumericVector FF, NumericVector NN, NumericMatrix NNorm, NumericVector SStatein,
                 int JJ, int TAU, int MM, NumericVector delta,
                 NumericVector PDF, NumericVector Dd,NumericVector DDD,
                 NumericMatrix TPM){
  NumericVector F(FF);
  NumericVector N(NN);
  NumericMatrix Norm(NNorm);
  NumericVector StateIn(SStatein);
  NumericVector pi(delta);
  NumericVector pdf(PDF);
  NumericVector d(Dd);
  NumericVector D(DDD);
  NumericMatrix tpm(TPM);
  int J(JJ);
  int tau(TAU);
  int M(MM);
  
  
  
  for(int t=0; t<=tau-1; t++){
    //Rcout << t << "\n";
    N[t]=0;
    
     for (int j=0; j<J; j++)
     {
    if (t < 1)
    {
    Norm(j,0) = pdf[j*tau] * pi[j];
    }
    else
    {
    Norm(j,t) = pdf[j*tau+t] * (StateIn[tau*j+t] - F[j*tau+t-1] + Norm(j,t-1));
    }
    N[t] += Norm(j,t);
    }

   for (int j = 0; j < J; j++)
  {
    Norm(j,t) /= N[t];
  }
  
  for (int j = 0; j <J; j++)
  {
  F[j*tau+t] = 0;
  double Observ = 1;

    if (t < tau - 1)
    {
      for (int u = 1; u <= std::min(t + 1, M); u++)
      {
        Observ *= pdf[j*tau+t - u + 1] / N[t - u + 1];

        if (u < t + 1)
        {
          F[j*tau+t] += Observ * d[j*M+u-1] * StateIn[j*tau+t - u + 1];
         // Rcout << d[j*M+u]<< "\n";
          
          
        }
        else
        {
          F[j*tau+t] += Observ * d[j*M+t] * pi[j];
          //Rcout << d[j*M+t]<< "\n";
        }
      }
    }
   else
   {
    for (int u = 1; u <= std::min(tau, M); u++)
    {
    Observ *= pdf[j*tau+tau - u] / N[tau - u];
      if (u < tau)
      {
        F[j*tau+tau - 1] += Observ * D[j*tau+u-1] * StateIn[j*tau+tau - u];
      }
      else
      {
      F[j*tau+tau - 1] += Observ * D[j*tau+tau-1] * pi[j];
      }
    }
  }
  }
  if (t < tau - 1)
  {
    for (int j = 0; j < J; j++)
    {
      StateIn[j*tau+t + 1] = 0;
      //Rcout << j*tau+t + 1 << "\n";  
      for (int i = 0; i <J ; i++)
      {
        StateIn[j*tau+t + 1] += tpm(i,j)* F[i*tau+t];
      }

    }
  }

  }
  return DataFrame::create(Named("F")= F,Named("N")=N,Named("StateIn")=StateIn);
  
  }
  