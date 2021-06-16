comm<-read.table("comm.txt",head=T,row.names=1)
soil<-read.table("soil.txt",head=T,row.names=1)
geo<-read.table("geo.txt",head=T,row.names=1)

# Install and load (if you haven't already) the following packages

#install.packages("vegan")
library(vegan)

# Now, we will generate a distance matrix using "comm". We will do that in order to see how communities are dissimilar between each other in relation to their species. In our example, we will use "Bray-Curtis" distance.

commdist<-vegdist(comm, method="bray")

# Also, we will generate a distance matrix using "soil". We will do that in order to see how the sites are dissimilar between each other in relation to their soil variables (for example, C and N). In our example, we will use "Bray-Curtis" distance.

soildist<-vegdist(soil, method="bray")

# Then, we will generate a distance matrix using "geo". We will do that in order to see how communities are far from each other in space. In our example, we will use "Euclidean" distance.

geodist<-vegdist(geo, method="euclidean")

### Now, we are going to use a partial mantel test to test if communities have more different species compositions as the soil composition in the sites get increasingly different. We will do that while removing any possible spatial autocorrelation. So, the third distance matrix in the analysis will be the spatial distance between sites.

mantel.partial(commdist, soildist, geodist, method="pearson", permutations=999)
  