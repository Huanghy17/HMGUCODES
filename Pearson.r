#### Pearson correlation ENV_alpha####

library("Hmisc")

Cor_ENV_ENV <- read.csv("a.csv",header = T, row.names = 1)

Cor_Pearson_ENV_ENV <- rcorr(as.matrix(Cor_ENV_ENV))

Cor_Pearson_ENV_ENV$r
Cor_Pearson_ENV_ENV$P
symnum(Cor_Pearson_ENV_ENV$P, cutpoints = c(0,  .001,.01,.05, .1, 1),
       symbols = c("***","**","*","."," "))

write.csv(Cor_Pearson_ENV_ENV$r,"Pearson_ENV_alpha_Species_alps_unique.csv")
write.csv(Cor_Pearson_ENV_ENV$P,"Pearson_p_ENV_alpha_Species_alps_unique.csv")
write.csv(symnum(Cor_Pearson_ENV_ENV$P, 
                 cutpoints = c(0,  .001,.01,.05, .1, 1), 
                 symbols = c("***","**","*","."," ")),
          "Pearson_p_ENV_alpha_l_Species_alps_unique.csv" )

####Draw figure
library(corrplot)

# Insignificant correlation are crossed, p<0.05
pdf("Pearson_ENV_alpha_Species_alps_unique.pdf")
p <- corrplot(Cor_Pearson_ENV_ENV$r, type="upper",
              p.mat = Cor_Pearson_ENV_ENV$P, sig.level = 0.05, insig = "blank")
dev.off()