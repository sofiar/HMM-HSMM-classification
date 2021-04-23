### Miro las formas de las distribuciones a simular
library(rstan)
source('FBcpp.R')
source('ViterbiImpl.R')
source('funAux.R')

# Caso mu1 = mu2
# Tiempos de espera
# hist(rnbinom(1000,size=50,mu=40)+1,col='red',freq=F) # Binomial negativo
# hist(rnbinom(1000,size=5,mu=40)+1, col=rgb(0, 1, 0, 0.5) ,add=T,freq=F)# Geometrico

# observations 
mean=c(0.3,0) #c(1,0) c(3,0)
sigma=c(1,1) 

# Sojourn Times
ks=list(c(80,100),c(30,50),c(10,30),c(1,10))#,c(80,85),c(30,35),c(10,15),c(1,5))
 # Mu=list(c(80,100),c(70,90),c(75,95), #1
 #         c(90,95),c(80,85),c(70,75), #2
 #         c(30,50),c(20,40),c(35,55), #3
 #         c(30,35),c(20,25),c(35,40), #4
 #         c(3,23),c(5,25),c(7,27), #5
 #         c(3,8),c(5,10),c(4,9)) #6
indxMu=rep(seq(1,6),each=3)

Mu=list()
Mu[[1]]=list(c(88.5,91.5),c(82.5,97.5),c(75,105))
Mu[[2]]=list(c(38.5,41.5),c(32.5,47.5),c(25,55))
Mu[[3]]=list(c(18.5,21.5),c(12.5,27.5),c(5,35))
mmean=c(90,40,20)
mdif=c(3,15,30)

ks=list(c(80,100),c(30,50),c(1,30),c(1,10))


Index=c()
Models=c()
Values=c()
Track=c()
Sampl=c()
KKs=c()
MuMean=c()
MuDif=c()
#MUs=c()


for (jnd in 1:length(Mu))
{
  for (knd in 1:length(Mu[[jnd]]))
  {
    for (ind in 3:length(ks))
    {
      k=ks[[ind]]
      mu=Mu[[jnd]][[knd]]
      
      # LOO
      
      source('simultateProcess.R')
      for (ii in 1:Nrep)
      {
        source('fitModels.R')  
        source('predictModels.R')
        Index=c(Index,index)
        Models=c(Models,Model)
        Values=c(Values,values)
        Track=c(Track,track)
        Sampl=c(Sampl,sampl)
        KKs=c(KKs,rep(ind,length(values)))
        MuMean=c(MuMean,rep(mmean[jnd],length(values)))
        MuDif=c(MuDif,rep(mdif[knd],length(values)))
       }
      save.image("E1.RData") 
}
   
    all.results=data.frame(Index,Models,Values,Track,Sampl,MuDif,MuMean,KKs)
    write.csv(all.results,paste('E1_jnd',jnd,'knd',knd,'.csv',sep=''))
    
   
    rm(Index,Models,Values,Track,Sampl,KKs,MuMean,MuDif,TotalData,Filt.Data)
      
    Index=c()
    Models=c()
    Values=c()
    Track=c()
    Sampl=c()
    KKs=c()
    MuMean=c()
    MuDif=c()
}
    
  
}





