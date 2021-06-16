####ANOVA####
a <- aov(Chao1 ~ Group, data = DATA)

##@@@@Turkey Test####
#if aov P<=0.05, continue the following test to find out the SD.
TukeyHSD(a, conf.level = 0.95)
