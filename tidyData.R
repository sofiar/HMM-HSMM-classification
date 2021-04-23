###################################################################
################# Limpieza y promedio de los datos ################
###################################################################

library(dplyr)
library(lubridate)
library(ggplot2)
library(grid)
library(rstan)
library(mvtnorm)
library(rlist)
library(forecast)
library('TTR')

orden=c('Acc_x','Acc_y','Acc_z')

# 1.Load data
setwd("clasificados_sep")
temp = list.files(pattern="*.csv")
sep.files = lapply(temp, read.csv)

setwd("clasificados_dic")
temp = list.files(pattern="*.csv")
dic.files = lapply(temp, read.csv)

# combine all data
data=c(dic.files,sep.files)
ltot=length(data)

posibles.B=c('Walk', 'Fast.Walk','Resting','Eating','Search','Vigilance',
             'Vigilance.down')



filesC=list()
#n=10 # window size 

for (i in 1:ltot)
{ # Just the behavoirs of interest
  qn.si=which(data[[i]]$Behaviours  %in% posibles.B)  
  filesC[[i]]=data[[i]][qn.si,]
  

  filesC[[i]]$Acc_x=ma(filesC[[i]]$Acc_x,order=n)
  filesC[[i]]$Acc_y=ma(filesC[[i]]$Acc_y,order=n)
  filesC[[i]]$Acc_z=ma(filesC[[i]]$Acc_z,order=n)
  filesC[[i]]$VeDBA=ma(filesC[[i]]$VeDBA,order=n)
  filesC[[i]]$Pitch.angle=ma(filesC[[i]]$Pitch.angle,order=n)
  
  filesC[[i]]=filesC[[i]][!is.na(filesC[[i]]$Acc_x),]
 
  
  ####################################
  ########## Mean by second ##########
  ####################################
  
  filesC[[i]]=filesC[[i]] %>% group_by(DateTime) %>% summarise(Acc_x=mean(Acc_x),Acc_y=mean(Acc_y),
                                                               Acc_z=mean(Acc_z),Behaviours=unique(Behaviours),
                                                               SD.VeDBA=sd(VeDBA),SD.Pitch=sd(Pitch.angle),
                                                               VeDBA=mean(VeDBA),Pitch.angle=mean(Pitch.angle))
  #Behaviour to number
  beha.num=rep(NA,length(filesC[[i]]$Behaviours))
  beha.num[which(filesC[[i]]$Behaviours=='Walk')]=1
  beha.num[which(filesC[[i]]$Behaviours=='Fast.Walk')]=2
  beha.num[which(filesC[[i]]$Behaviours=='Resting')]=3
  beha.num[which(filesC[[i]]$Behaviours=='Eating')]=4
  beha.num[which(filesC[[i]]$Behaviours=='Search')]=5
  beha.num[which(filesC[[i]]$Behaviours=='Vigilance')]=6
  beha.num[which(filesC[[i]]$Behaviours=='Vigilance.down')]=7
  
  filesC[[i]]=filesC[[i]]%>% mutate(States=beha.num)
}






