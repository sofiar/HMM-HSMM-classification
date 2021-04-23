#### Predict Models: Sin LOO

### HSMM

nsamps=30

AccuracyFBHSMM=numeric(nsamps)
AccuracyVHSMM=numeric(nsamps)
CrossEHSMM=numeric(nsamps)


samp=sample(seq(1,length(outhsmm$mea[,1])),nsamps)

dd=TotalData %>% filter(track==ii) 
y=dd$Obs
lar=length(y)

M=min(1000,lar)



for (jj in 1:nsamps)
{
  mu.fit=outhsmm$mu[samp[jj],]
  phi.fit=outhsmm$phi[samp[jj],]
  mea.fit=outhsmm$mea[samp[jj],]
  sd.fit=outhsmm$sigma[samp[jj],]

    ## Sojourn time (BN)
    FB.d=matrix(NA,M,2)
    for (m in 1:2)
    {
      FB.d[,m]= c(dnbinom(c(0:(M - 1)), mu=mu.fit[m],size=phi.fit[m]))
    }
    FB.d=t(FB.d)
    rm(m)
    
    for (m in 1:2)
    {
      FB.d[m,M]=1-sum(FB.d[m,-M])
    }
    rm(m)
    
    ## PDF 
    
  FB.pdf=matrix(NA,2,lar)
  for (m in 1:2)
  {
  for (s in 1:lar)
      {
        sd = sd.fit[m]
        muu = mea.fit[m]
        FB.pdf[m,s] <- dnorm(y[s], mean = muu, sd =sd)  
      }
    }
    rm(m)
    
    lower_bound <- 1e-300
    FB.pdf[FB.pdf < lower_bound] <- lower_bound
    delta=rep(1,2)/2
    tpm=matrix(c(0,1,1,0),ncol=2)
    
    # Forwar backward    
    
    FB.Result=FBAcpp(pi=delta,pdf=FB.pdf,d=FB.d,p=tpm)
    pred=apply(FB.Result$Gamma,2,which.max)
    AccuracyFBHSMM[jj]=mean(pred==dd$States)
    quienes=which(pred==dd$States)
    CrossEHSMM[jj]=cross.entropy(FB.Result$Gamma,dd$States,quienes)
    
    # Viterbi   
    resul=ViterbiImpl(delta, FB.pdf,FB.d, tpm)
    AccuracyVHSMM[jj]=mean(dd$States==resul)

}
  



##### Prediccion HMM

source('predFun.R')

AccuracyFBHMM=numeric(nsamps)
AccuracyVHMM=numeric(nsamps)
CrossEHMM=numeric(nsamps)

dd=TotalData %>% filter(track==ii) 
x=dd$Obs
lar=length(x)

for (jj in 1:nsamps)
{
  mod=list()
  mod$m=2
  mod$mus=outhmm$mea[samp[jj],]
  mod$sigmas=outhmm$sigmas[samp[jj],]
  mod$delta=rep(0.5,2)
  mod$tpm=outhmm$theta[samp[jj],,]
  
  ## allprobs
  allprobs=matrix(NA,ncol=2,nrow=lar)
    for (i in 1:(lar))
    {
      for (j in 1:2)
        allprobs[i,j]= dnorm(x[i],mean=mod$mus[j],sd=mod$sigmas[j])
    }
    
    # FB
    Pstates=norm.HMM.state_probs(x , mod ,allprobs)
    pred=apply(Pstates,2,which.max)
    AccuracyFBHMM[jj]=mean(dd$States==pred)
    quienes=which(pred==dd$States)
    CrossEHMM[jj]=cross.entropy(Pstates,dd$States,quienes)
    
    # Viterbi
    predV=HMM.viterbi(x, mod$m, gamma=mod$tpm, allprobs, delta=mod$delta)
    AccuracyVHMM[jj]=mean(dd$States==predV)
   
  }


#### Create Data Frame 

index=rep(c(rep('FBAcc',nsamps),rep('VAcc',nsamps),rep('CE',nsamps)),2)
Model=c(rep('HSMM',nsamps*3),rep('HMM',nsamps*3))
values=c(as.vector(AccuracyFBHSMM),as.vector(AccuracyVHSMM),as.vector(CrossEHSMM),
         as.vector(AccuracyFBHMM),as.vector(AccuracyVHMM),as.vector(CrossEHMM))
#track=rep(rep(seq(1:Nrep),nsamps),6)
track=rep(ii,length(index))
sampl=rep(seq(1:nsamps),6)
