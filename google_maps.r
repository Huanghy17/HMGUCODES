#maps: all the maps
library(ggmap)
library(ggplot2)
LakeLoc <- geocode ("LakeName")#demo: MOALoc <- geocode ('Moaralmsee')
rownames(LakeLoc)<-c("LakeName")#demo: rownames(MOA)<-c("MOA")
LakesLocation <- rbind(LakeLoc1,LakeLoc2)#demo: LakesLocation <- rbind(MoALoc,GIGLoc)
Austria <- get_map(location = "Austria", zoom=7)
Austria <- ggmap(Sustria)
Austria
Geo <- Austria + geom_point(data=LakesLocation,aes(x=lon,y=lat),size=1,alpha=0.6,color="red",shape=19)
GEO
AusD <- get_map(location=c(lon=13.8,lat=47.5),zoom=10)
AusD <- ggmap(AusD)
AusD
GEOD <- AusD + geom_point(data=GIGLoc, aes(x=lon,y=lat),size=5,alpha=0.6,color="red",shape=19)+geom_text(data=GIGLoc,size=3,aes(x=lon,y=lat,label=rownames(GIGLoc)))
GEOD <- GEOD + geom_point(data=MOALoc, aes(x=lon,y=lat),size=5,alpha=0.6,color="red",shape=19)+geom_text(data=MOALoc,size=3,aes(x=lon,y=lat,label=rownames(MOALoc)))
GEOD <- GEOD + geom_point(data=OLALoc, aes(x=lon,y=lat),size=5,alpha=0.6,color="red",shape=19)+geom_text(data=OLALoc,size=3,aes(x=lon,y=lat,label=rownames(OLALoc)))
GEOD <- GEOD + geom_point(data=TWALoc, aes(x=lon,y=lat),size=5,alpha=0.6,color="red",shape=19)+geom_text(data=TWALoc,size=3,aes(x=lon,y=lat,label=rownames(TWALoc)))
GEOD <- GEOD + geom_point(data=WIRLoc, aes(x=lon,y=lat),size=5,alpha=0.6,color="red",shape=19)+geom_text(data=WIRLoc,size=3,aes(x=lon,y=lat,label=rownames(WIRLoc)))
Lake <- get_map(location = "LakeName",zoom=14)
Lake<-ggmap(LakeName)
Lake<-Lake+geom_text(x=13.938,y=47.765,label="Lake",size=4,color="darkblue")+geom_text(x=13.9388,y=47.763,label="(zoom=14)",size=3,color="black",fontface="italic")
Lake
#MOA <- get_map(location = "Moaralmsee",zoom=14)
#MOA<-ggmap(MOA)
#MOA<-MOA+geom_text(x=13.938,y=47.765,label="MOA",size=4,color="darkblue")+geom_text(x=13.9388,y=47.763,label="(zoom=14)",size=3,color="black",fontface="italic")
#MOA

#Rarefaction curves
library(vegan)
otu.raw <- read.delim("otu_table.txt")
otu2.raw<-as.matrix(otu.raw[,1:ncol(otu.raw)-1])
mapping <- read.delim("map_file.txt")
order<-cbind(colnames(otu2.raw), name="otus")
sort_mapping<-merge(order,mapping,by.x="V1", by.y="SampleID", sort=FALSE)
color <- as.character(sort_mapping$Site)
color[color=="GIG"] <- "red"
color[color=="MOA"] <- "blue"
color[color=="OLA"] <- "gray"
color[color=="TWA"] <- "orange"
color[color=="WIR"] <- "darkgreen"
pdf("Rarefaction curves.pdf", width=8, height=6, pointsize=12)
rarecurve(t(otu2.raw), col=color, step=100, label = FALSE, xlim=c(0,20000), xlab = "Sequence", ylab = "OTUs")
legend("bottomright", legend = c("GIG", "MOA", "OLA","TWA","WIR"), col=c("red", "blue", "gray","orange","darkgreen"), pch=15)
dev.off()