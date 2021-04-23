###########################################################################
############################### Loo HMM no Ar #############################
###########################################################################

library(rstan)
library(tidyverse)
library(mvtnorm)
source('predFun.R')
source('funAux.R')
    
#Load data and tidy
load("RealData.RData")   

#length of the series
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


for (k in 1:Nrep)
{
  who=seq(1,Nrep)[-k]
  lars=Lars[-k]
  PredsFB[[k]]=matrix(NA,ncol=nsamples,nrow=Lars[k])
  PredsV[[k]]=matrix(NA,ncol=nsamples,nrow=Lars[k])
  
  ## count state changes
  counts=matrix(1,ncol=4,nrow=4)
  
  for (i in who)
  {
    for (j in 1:(Lars[i]-1))
    {
      counts[filesC[[i]]$States[j],filesC[[i]]$States[j+1]]=counts[filesC[[i]]$States[j],filesC[[i]]$States[j+1]]+1
     }  
  }
  
  
  # stan input
  
  data=list()
  data$K=4
  data$Nrep=Nrep-1
  data$T=lars
  data$TT=sum(lars)
  z=c()
  y=matrix(NA,ncol=3,nrow = sum(lars))
  indx.inicial=c(1,cumsum(lars)+1)[1:(Nrep-1)]
  indx.final=cumsum(lars)
 
  ii=1
  for (i in who)
  {
    
    # States sequences
    z=c(z,filesC[[i]]$States)
    
    # Observations
    y[indx.inicial[ii]:indx.final[ii],]=cbind(log(filesC[[i]]$VeDBA),
                                              filesC[[i]]$Pitch.angle,
                                             log(filesC[[i]]$SD.Pitch))
   ii=ii+1
  }
  
  
  data$z=z
  data$y=y
  data$inxIn=indx.inicial
  data$inxEn=indx.final
  data$counts=t(counts)
  
  # fit with stan
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  # fit
  fit <- stan(file = 'hmmNoAr.stan',data = data,chain=2,cores=4)
  outs <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 
  
  ########################################################
  ###################### Prediction ######################
  ########################################################
  
  samp=sample(seq(1,length(outs$alphas[,1,1])),nsamples)
  x=cbind(log(filesC[[k]]$VeDBA),
          filesC[[k]]$Pitch.angle,
          log(filesC[[k]]$SD.Pitch))
  
  for (j in 1:nsamples)
  {
    # sample from the posterior
    
    mod=list()
    mod$m=4
    mod$mus=outs$alphas[samp[j],,]
    mod$sigmas=outs$sigmas[samp[j],,]
    mod$delta=rep(0.25,4)
    mod$tpm=outs$theta[samp[j],,]
    
    
    allprobs=matrix(NA,ncol=4,nrow=dim(x)[1])
    for (i in 1:(dim(x)[1]))
    {
      for (jj in 1:4)
        allprobs[i,jj]= dmvnorm(x[i,],mean=mod$mus[,jj],sigma=diag((mod$sigmas[,jj])^2))
    }
    
    
    
    # Forward backward    
    Pstates=norm.HMM.state_probs(x , mod ,allprobs)
    pred=apply(Pstates,2,which.max)
    Error01V[k,j]=mean(filesC[[k]]$States==pred)
    quienes=which(pred==filesC[[k]]$States)
    CE[k,j]=cross.entropy(Pstates,filesC[[k]]$States,quienes)
    PredsFB[[k]][,j]=pred
    
    # Viterbi
    predV=HMM.viterbi(x, mod$m, gamma=mod$tpm, allprobs, delta=mod$delta)
    Error01FB[k,j]=mean(filesC[[k]]$States==predV)
    PredsV[[k]][,j]=predV
    
    
  
  }
  
  
  print(k)
  save.image("hmmNoAr.RData")
}