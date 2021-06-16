####b-diversity####
library("vegan")
matrix <- read.csv("a.csv",header = T, row.names = 1)
braycurtis <- vegdist(matrix, "bray")
hist(braycurtis)
dis <- as.matrix(braycurtis)
write.csv(dis, "beta_vOTUs_gene_profile.reads_number3group.csv")

#continue with Kruskal and Wilcoxon test