## Viterbi, from Leos Barajas 2017

HMM.viterbi <- function(x, m, gamma, allprobs, delta=NULL, ...)
{
  if(is.null(delta)) delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  n <- length(x)
  xi <- matrix(0, n, m)
  foo <- delta*allprobs[1,]
  xi[1,] <- foo/sum(foo)
  for(i in 2:n)
  {
    foo <- apply(xi[i-1,]*gamma, 2, max)*allprobs[i,]
    xi[i,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for(i in (n-1):1)
    iv[i] <- which.max(gamma[, iv[i+1]]*xi[i,])
  return(iv)
}

# Local decoding: from Zucchini 2016

### A.1.7 Computing log(forward probabilities)

norm.HMM.lforward<-function(x,mod,allprobs)
{
  n <- length(x)
  lalpha <- matrix(NA,mod$m,n)
  foo <- mod$delta*allprobs[1,]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[,1] <- lscale+log(foo)
  for (i in 2:n)
  {
    foo          <- foo%*%mod$tpm*allprobs[i,]
    sumfoo       <- sum(foo)
    lscale       <- lscale+log(sumfoo)
    foo          <- foo/sumfoo
    lalpha[,i]   <- log(foo)+lscale
  }
  return(lalpha)
}



### A.1.8 Computing log(backward probabilities)

norm.HMM.lbackward<-function(x,mod,allprobs)
{
  n          <- length(x)
  m          <- mod$m
  lbeta      <- matrix(NA,m,n)
  lbeta[,n]  <- rep(0,m)
  foo        <- rep(1/m,m)
  lscale     <- log(m)
  for (i in (n-1):1)
  {
    foo        <- mod$tpm%*%(allprobs[i+1,]*foo)
    lbeta[,i]  <- log(foo)+lscale
    sumfoo     <- sum(foo)
    foo        <- foo/sumfoo
    lscale     <- lscale+log(sumfoo)
  }
  return(lbeta)
}


norm.HMM.state_probs <- function (x , mod ,allprobs)
{
  n <- length(x)
  la <- norm.HMM.lforward(x,mod,allprobs)
  lb <- norm.HMM.lbackward(x,mod,allprobs)
  c <- max(la[,n])
  llk <- c + log(sum(exp(la[,n]-c)))
  stateprobs <- matrix (NA,ncol=n,nrow = mod$m)
  for (i in 1:n) stateprobs [,i] <-exp(la[,i]+lb[,i]-llk)
  return (stateprobs)
}