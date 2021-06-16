library("vegan")

Group_A<- read.csv("Group_ANOSIM_A.csv")

distance.bray<-vegdist(abund_table,method = 'bray')

anosim.result<-anosim(distance.bray,Group_A$Lake,permutations = 999)

summary(anosim.result)