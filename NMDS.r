#####Modified: http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html

abund_table<-read.csv("SPE_A.csv",row.names=1,check.names=FALSE)
#Transpose the data to have sample names on rows
meta_table<-read.csv("Group_ANOSIM_A.csv",row.names=1,check.names=FALSE)
#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
meta_table<-meta_table[rownames(abund_table),]

#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))

#Load vegan library
library(vegan)
#Get MDS stats
sol<-metaMDS(abund_table,distance = "bray", k = 2, trymax = 50)

#Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],lakes=as.factor(grouping_info[,1]),GrowingSeason=as.factor(grouping_info[,2])
                ,filter=as.factor(grouping_info[,2]))

#Get spread of points based on countries###CHANGE X HERE!!!grouping_info[,x]
plot.new()
ord<-ordiellipse(sol, as.factor(grouping_info[,1]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()


#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Generate ellipse points
df_ell <- data.frame()
for(g in levels(NMDS$lakes)){
  if(g!="" && (g %in% names(ord))){
    
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$lakes==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                  ,lakes=g))
  }
}


#Generate mean values from NMDS plot grouped on Countries
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$lakes),mean)

#Now do the actual plotting
library(ggplot2)

shape_values<-seq(1:5)

p<-ggplot(data=NMDS,aes(x,y,colour=lakes))
p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)
p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
p<-p+geom_point(aes(shape=GrowingSeason))+scale_shape_manual(values=shape_values)+theme_bw() 
pdf("NMDS_lake_A_Aug.pdf")
print(p)
dev.off()