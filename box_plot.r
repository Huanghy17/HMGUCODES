Data <- read.csv("input.csv",header = T)

a <- ggplot(Data, aes(x=Group, y=Chao1, color=Group)) + 
  geom_boxplot()+
  geom_jitter(shape=16,size=3,  position=position_jitter(0.2))
  
pdf("box_plot.pdf", paper="special",height=10, width=15)
print(a)
dev.off()