library(vegan)
library(OTUtable)
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

all_table<-read.csv("All_abundant_data.csv",row.names=1,check.names=FALSE)

table_Jul<-all_table[,1:5]

PGs_Jul_M <- lm(table_Jul$`1/kWTACiS` ~ table_Jul$metabolism, data = table_Jul)
summary(PGs_Jul_M)
AIC(PGs_Jul_M)
a <- summary(PGs_Jul_M)
sink("Chao1_Bshannon.txt")
a
sink()

p <- qplot(table_Jul$`1/kWTACiS`, table_Jul$metabolism, data=table_Jul)+
  geom_smooth(method = "lm")+
  labs(title="Metabolism r2=0.803 p<0.001*** AIC=-72.277")+
  theme(plot.title=element_text(size=14,hjust=0.5))+
  xlab("1/KWAS")+
  ylab("Abundance")

pdf("linear_Regression.pdf")
print(p)
dev.off()
