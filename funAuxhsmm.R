########### Funciones auxiliares library(hsmm) ##########


get.pdf <- function(inputData, od, J, M, param){
  tau <- get.tau(inputData)
  pdf <- rep(0, times = J * tau)
  for (j in 1:J){
    # od = "Bernoulli"
    if (od == "bern"){
      for (t in 1:tau){
        if (inputData[t] == 1){
          pdf[t + (j - 1) * tau] <- param$b[j]
        }
        else {
          pdf[t + (j - 1) * tau] <- 1 - param$b[j]
        }
      }
    }
    # od = "Gaussian"
    if (od == "norm"){
      pdf[(1 + (j - 1) * tau):(tau + (j - 1) * tau)] <- dnorm(inputData[1:tau], mean = param$mean[j], sd = sqrt(param$var[j]))
    }
    # od = "Poisson"
    if (od == "pois"){
      pdf[(1 + (j - 1) * tau):(tau + (j - 1) * tau)] <- dpois(inputData[1:tau], lambda = param$lambda[j])
    }
    # od = "Student.t"
    if (od == "t"){
      pdf[(1 + (j - 1) * tau):(tau + (j - 1) * tau)] <- dtmod(inputData[1:tau], mu = param$mean[j], sigma = sqrt(param$var[j]), nu = param$df[j])
    }
    # od = "multivar.Gaussian"
    if (od == "mvnorm"){
      pdf[(1 + (j - 1) * tau):(tau + (j - 1) * tau)] <- dmvnorm(aperm(inputData[,1:tau]), mean = param$mean[,j], sigma = param$sigma[,,j])
    }
  }
  
  lower_bound <- 1e-300
  pdf[pdf < lower_bound] <- lower_bound
  return(pdf)
}


get.d <- function(rd, J, M, param){
  d <- c()
  for (j in 1:J){
    # rd = "non.parametric"
    if (rd == "nonp")
      d <- c(d, param$np[,j])
    # rd = "geometric"
    if (rd == "geom")
      d <- c(d, dgeom(c(0:(M - 1)), param$p[j]))
    # rd = "negative.binomial"
    if (rd == "nbinom")
      d <- c(d, dnbinom(c(0:(M - 1)), size = param$r[j], prob = param$pi[j]))
    # rd = "logarithmic.geometric"
    #    if (rd == "logarithmic.geometric")
    #      d <- c(d, dloggeom(c(0:(M - 1)), param$p.loggeom[j], param$theta[j]))
    # rd = "Poisson"
    if (rd == "pois")
      d <- c(d, dpois(c(0:(M - 1)), param$lambda[j]))
  }
  return(d)
}


get.tau <- function(inputData){
  if (length(dim(inputData)) != 2) {
    tau        <- as.integer(length(inputData))
  }
  if (length(dim(inputData)) == 2) {
    tau        <- as.integer(dim(inputData)[2])
  }
  return(tau)
}

FB <- function(censoring, tau, J, M, FB.d, FB.p.tpm, FB.pi.ini, FB.pdf, 
               F, L, G, L1, N, Norm, eta, xi, error){
  .C("FB", censoring, tau, J, M, FB.d, FB.p.tpm, FB.pi.ini, FB.pdf, 
     F, L, G, L1, N, Norm, eta, xi, error, PACKAGE="hsmm");
}

