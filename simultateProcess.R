 
library(tidyverse)

### 1. Set parameters
lar=5000

# delta
delta=c(0.5,0.5)
Nrep=10

TotalObs=c()
TotalStates=c()
TotalTracks=c()

TotalStimes=c()
TotalStimesB=c()

### 2. Simulations
for (ll in 1:Nrep)
{
# Sojourn times
stimesBN=c()
statesBN=c()
sec.est=c()
q=sample(x=c(1,2),size=1)
ini=1
l=0
while (l<lar)
{
  samp=0
  samp= rnbinom(1,mu=mu[q],size=k[q])+1 
  statesBN[ini:(ini+samp)]=q
  sec.est=c(sec.est,q)
  ini= ini+samp
  
  if(q==1)
  {nq=2;}
  else  
  {nq=1;}
  q=nq
  l=sum(stimesBN)
  
}

TotalStimes=c(TotalStimes,stimesBN)
TotalStimesB=c(TotalStimesB,statesBN)

# Observaciones 
  
X=numeric(length(statesBN))


for (i in 1:length(X))
{
  X[i]=rnorm(1,mean=mean[statesBN[i]],sd=sigma[statesBN[i]])
}

TotalObs=c(TotalObs,X)
TotalStates=c(TotalStates,statesBN)
TotalTracks=c(TotalTracks,rep(ll,length(statesBN)))


}

TotalData=data.frame(States=TotalStates,Obs=TotalObs,track=TotalTracks)
TotalData$StatesF=as.factor(TotalData$States)


