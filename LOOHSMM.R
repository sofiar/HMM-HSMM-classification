######################################################################
################# LOO Prediction: HSMM Ar(1) ######################### 
######################################################################

library(rarhsmm)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(glmnet)
library(hsmm)
library(rstan)
library(parallel)

# 1.Load data  and tidy
n=10 # size window 
detach("package:signal", unload=TRUE)
source('tidyData.R')
source('UnifyandDelete.R') 

source('funAux.R')
source('FBcpp.R')
source('ViterbiImpl.R')


# 2. length of the series
Nrep=length(filesC)
nstates=4
Lars=numeric(Nrep)
nsamples=100

for (i in 1:Nrep)
{
  Lars[i]=length(filesC[[i]]$Acc_x)
}
rm(i)

# Loop LOO
Error01V=matrix(NA,nrow = Nrep,ncol=nsamples) 
CE=matrix(NA,nrow = Nrep,ncol=nsamples) 
Error01FB=matrix(NA,nrow = Nrep,ncol=nsamples)

PredsV=list()
PredsFB=list()

for (k in 19:Nrep)
{
  who=seq(1,Nrep)[-k]
  lars=Lars[-k]
  PredsFB[[k]]=matrix(NA,ncol=nsamples,nrow=Lars[k])
  PredsV[[k]]=matrix(NA,ncol=nsamples,nrow=Lars[k])
  
  ## count state changes
  counts=matrix(1,ncol=nstates,nrow=nstates)
  for (i in who)
  {
    seq.est=rle(filesC[[i]]$States)$values
    for (j in 1:length(seq.est))
    {
      counts[seq.est[j],seq.est[j+1]]=counts[seq.est[j],seq.est[j+1]]+1
    }  
  }  
  counts2=matrix(NA,ncol=nstates-1,nrow=nstates)
  for (i in 1:nstates)
  {counts2[i,]=counts[i,-i]}
  
  # input for stan
  
  data=list()
  data$K=nstates
  data$Nrep=Nrep-1
  data$T=lars
  data$TT=sum(lars)
  N=numeric(Nrep-1) #number of NAS of each ts
  z=c()
  u=c()
  y=matrix(NA,ncol=3,nrow = sum(lars))
  indx.inicial=c(1,cumsum(lars)+1)[1:(Nrep-1)]
  indx.final=cumsum(lars)
  S.est=c()
  ii=1
  for (i in who)
  {
    
  # States sequences
    z=c(z,filesC[[i]]$States)
    S.est=c(S.est,rle(filesC[[i]]$States)$values)
    
    # Sojourn times
    u=c(u,rle(filesC[[i]]$States)$lengths)
    N[ii]=length(rle(filesC[[i]]$States)$lengths)
    
    # Observations
    y[indx.inicial[ii]:indx.final[ii],]=cbind(log(filesC[[i]]$VeDBA),
                                            filesC[[i]]$Pitch.angle,
                                            log(filesC[[i]]$SD.Pitch))
 
    ii=ii+1
    }
  
  
  data$N=N
  data$NN=length(u)
  data$z=z
  data$u=u
  data$y=y
  data$inxIn=indx.inicial
  data$inxEn=indx.final
  nas.inicio=c(1,cumsum(N)+1)[1:(Nrep-1)]
  data$nasIn=nas.inicio-1
  data$Sest=S.est
  data$counts=t(counts2)
  
  # Fit stan
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  fit <- stan(file = 'fitRdataBN.stan',data = data,chain=2,cores=4)
  outs <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 
  
  ########################################################
  ###################### Prediction ######################
  ########################################################
  
  samp=sample(seq(1,length(outs$mu[,1])),nsamples)

  for (j in 1:nsamples)
  {
    # sample from the posterior
    mu.fit=outs$mu[samp[j],]
    phi.fit=outs$phi[samp[j],]
    alphas.fit=outs$alphas[samp[j],,]
    betas1.fit=outs$betas1[samp[j],,]
    sd.fit=outs$sigma[samp[j],,]
    tpm.fit=outs$theta[samp[j],,]
    tpm.change=matrix(NA,nstates,nstates)
    for (m in 1:nstates)
    {
      tpm.change[m,m]=0
      tpm.change[m,-m]=tpm.fit[m,]
    }
    rm(m)
  
    # probabilities 
    
    lar=Lars[k]
      
    y=as.matrix(cbind(log(filesC[[k]]$VeDBA),
                        filesC[[k]]$Pitch.angle,
                        log(filesC[[k]]$SD.Pitch)))
    
   M=min(1000,lar)
    
    ## Sojourn time (BN)
      FB.d=matrix(NA,M,nstates)
      for (m in 1:nstates)
      {
        FB.d[,m]= c(dnbinom(c(0:(M - 1)), mu=mu.fit[m],size=phi.fit[m]))
      }
      FB.d=t(FB.d)
      
      for (m in 1:nstates)
      {
        FB.d[m,M]=1-sum(FB.d[m,-M])
      }
      
      ## PDF 
      
      FB.pdf=matrix(NA,nstates,lar)
      for (m in 1:nstates)
      {
        FB.pdf[m,1] = dmvnorm(y[1,], mean = alphas.fit[,m], sigma = diag((sd.fit[,m])^2))
        for (s in (2):lar)
        {
          sd = diag((sd.fit[,m])^2)
          mu=alphas.fit[,m]+betas1.fit[,m]*y[s-1,]
          FB.pdf[m,s] <- dmvnorm(y[s,], mean = mu, sigma =sd)  
        }
      }
      
      
      lower_bound <- 1e-300
      FB.pdf[FB.pdf < lower_bound] <- lower_bound
      delta=rep(1,nstates)/nstates
      
      # Forwar backward    
      FB.Result=FBAcpp(pi=delta,pdf=FB.pdf,d=FB.d,p=tpm.change)
      pred=apply(FB.Result$Gamma,2,which.max)
      Error01FB[k,j]=mean(pred==filesC[[k]]$States)
      quienes=which(pred==filesC[[k]]$States)
      CE[k,j]=cross.entropy(FB.Result$Gamma,filesC[[k]]$States,quienes)
      PredsFB[[k]][,j]=pred
      
      # Viterbi
      resul=ViterbiImpl(delta, FB.pdf,FB.d, tpm.change)
      Error01V[k,j]=mean(filesC[[k]]$States==resul)
      PredsV[[k]][,j]=resul
    }
  save.image("hsmmAr.RData")
  print(k)
 
   }

