
### unify Eating and Search
for (l in 1:ltot)
{
  filesC[[l]] = filesC[[l]] %>%
    mutate(Behaviours = recode(Behaviours, 
                               Eating = "Eating.Search",
                             Search = "Eating.Search"))
}

### unify Vigilances and resting as inactive
for (l in 1:ltot)
{
  filesC[[l]] = filesC[[l]] %>%
    mutate(Behaviours = recode(Behaviours, 
                               Vigilance = "Inactive",
                               Vigilance.down = "Inactive",
                               Resting='Inactive'))
}


# States
for (i in 1:ltot)
{
beha.num=rep(NA,length(filesC[[i]]$Behaviours))
beha.num[which(filesC[[i]]$Behaviours=='Walk')]=1
beha.num[which(filesC[[i]]$Behaviours=='Fast.Walk')]=2
beha.num[which(filesC[[i]]$Behaviours=='Inactive')]=3
beha.num[which(filesC[[i]]$Behaviours=='Eating.Search')]=4

filesC[[i]]=filesC[[i]]%>% mutate(States=beha.num)
}


#### Eliminate Resting when it las less than 1 minute



# for (i in 1: ltot)
# {
# changes=rle(filesC[[i]]$States)$values
# lengths=rle(filesC[[i]]$States)$length
# nuevodato=rep(rle(filesC[[i]]$States)$length,rle(filesC[[i]]$States)$length)
# filesC[[i]]=filesC[[i]] %>% mutate(stimes=nuevodato)
# 
# filesC[[i]]=filesC[[i]] %>% filter(!(Behaviours=='Resting' & stimes<300))
# 
# }

# saco la muestra de todo resting

#filesC[[40]]=NULL
