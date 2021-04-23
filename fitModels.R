# Ajuste con stan 

# Filtramos datos para el LOO
Filt.Data=TotalData %>% filter(track!=ii)
qq=(1:Nrep)[-ii]
NrepLoo=Nrep-1
# 1. HSMM

# Preparo entrada para stan 
data=list()
data$Nrep=NrepLoo
lars=Filt.Data %>% group_by(track) %>% summarise(lar=length(Obs))
n.files=lars$lar
data$T=n.files
data$TT=sum(lars$lar)

N=numeric(NrepLoo) #number of NAS of each ts
z=c()
u=c()
nas=c()
y=c()
S.est=c()

indx.inicial=c(1,cumsum(n.files)+1)[1:NrepLoo]
indx.final=cumsum(n.files)


for (i in 1:NrepLoo)
{
  dd=Filt.Data %>% filter(track==qq[i])
  
  # States sequences
  z=c(z,dd$States)
  S.est=c(S.est,rle(dd$States)$values)
  
  # Sojourn times
  u=c(u,rle(dd$States)$lengths)
  N[i]=length(rle(dd$States)$lengths)
  
  # Observations
  y[indx.inicial[i]:indx.final[i]]= dd$Obs
}


data$N=N
data$NN=length(u)
data$z=z
data$u=u-1
data$y=y
data$inxIn=indx.inicial
data$inxEn=indx.final
nas.inicio=c(1,cumsum(N)+1)[1:NrepLoo]
data$nasIn=nas.inicio-1
data$Sest=S.est

# Ajuste con rstan
#stanc("hsmm.stan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit <- stan(file = 'hsmm.stan',data = data,chain=2,cores=3)
outhsmm <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 

#### Borramos cosas
rm(data,S.est,u,dd,fit)


### Chequeo el ajuste

# Obs parameters

# par(mfrow=c(1,2))
# for (j in 1:2)
# {
#   hist(outhsmm$mea[,j],main='Mean Obs')
#   abline(v=mean[j],col='red')
# }

# par(mfrow=c(1,2))
# for (j in 1:2)
# {
#   hist(outhsmm$sigmas[,j],main='Sigma Obs')
#   abline(v=sigma[j],col='red')
# }

# Sojourn time parameters

# par(mfrow=c(1,2))
# for (j in 1:2)
# {
#   hist(outhsmm$mu[,j],main='Mu')
#   abline(v=mu[j],col='red')
# }


#  par(mfrow=c(1,2))
#  for (j in 1:2)
# {
#   hist(outhsmm$phi[,j],main='k')
#   abline(v=k[j],col='red')
# }

###

# hist(rnbinom(1000,size=outhsmm$phi[100,1],mu=outhsmm$mu[100,1])+1)
# plot(dnbinom(seq(0:100),size=outhsmm$phi[100,1],mu=outhsmm$mu[100,1]))

# fit HMM

## contamos los cambios entre estados
counts=matrix(0,ncol=2,nrow=2)
for (i in 1:NrepLoo)
{
  dd=Filt.Data %>% filter(track==qq[i])
  long=length(dd$Obs)
  for (j in 1:(long-1))
  {
    counts[dd$States[j],dd$States[j+1]]=counts[dd$States[j],dd$States[j+1]]+1
  }  
}  

data=list()
data$Nrep=NrepLoo
data$T=n.files
data$TT=sum(lars$lar)
data$z=z
data$y=y
data$counts=counts
data$inxIn=indx.inicial
data$inxEn=indx.final


# Ajuste con rstan
#stanc("hmm.stan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit <- stan(file ='hmm.stan',data = data,chain=2,cores=3)
outhmm <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 

#### Borramos cosas
rm(data,fit,counts,z)
### Chequeo el ajuste

# Obs parameters
# 
# par(mfrow=c(1,2))
# for (j in 1:2)
# {
#   hist(outhmm$mea[,j],main='Mean Obs')
#   abline(v=mean[j],col='red')
# }
# 
# par(mfrow=c(1,2))
# for (j in 1:2)
# {
#   hist(outhmm$sigmas[,j],main='Sigma Obs')
#   abline(v=sigma[j],col='red')
# }


## tmp
# para el comportamiento que tiene distribucion geometrica

#par(mfrow=c(1,2))
# hist(outhmm$theta[,2,2])
# abline(v=1-ps[2],col='red')
#   
# hist(outhmm$theta[,2,1])
# abline(v=ps[2],col='red')

# para el comportamiento que tiene distribucion BN
# hist(outhmm$theta[,1,2])
# abline(v=1-ps[1],col='red')

# hist(outhmm$theta[,1,1])
# abline(v=ps[1],col='red')


# Diferencias entre los ajustes

# par(mfrow=c(1,2))
# hist(rgeom(1000,prob=outhmm$theta[100,2,1])+1)
# hist(rgeom(1000,prob=ps[2])+1)

# hist(rnbinom(1000,mu=40,size=45)+1)
# hist(rgeom(1000,prob=outhmm$theta[100,1,2])+1)





