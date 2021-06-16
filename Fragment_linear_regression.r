####Draw fragment linear regression####
library(ggplot2)
library(vegan)
library(OTUtable)
library(plyr)

#change continus values in cyl to characters!IMPORTANT!
#Otherwise an error would be occured with "Error: A continuous variable can not be mapped to shape"

Alpha_ENV_C_Jul <- read.csv("Alpha_ENV_C_Jul.csv",header = T, row.names = 1)

Chao1 <- lm(Chao1 ~ ICD, data = Alpha_ENV_C_Jul)
summary(Chao1)
AIC(Chao1)
sink("C_Jul_Chao1.txt")
Chao1
AIC(Chao1)
sink()
p_Chao1_Whole <- qplot(Chao1, ICD, data=Alpha_ENV_C_Jul)+
  geom_smooth(method = "lm")+
  labs(title="ICD R2=0.17 p=0.036* AIC=255.21")+
  theme(plot.title=element_text(size=14,hjust=0.5))+
  xlab("Chao1")+
  ylab("ICD")
pdf("Alpha_ENV_C_Jul_Chao1_W.pdf")
print(p_Chao1_Whole)
dev.off()