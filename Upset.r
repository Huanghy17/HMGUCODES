#source("http://bioconductor.org/biocLite.R")
#BiocManager::install(c("UpSetR"))

library(UpSetR)

DATA = read.csv("a.csv",header=TRUE,row.names=1,check.names = FALSE)
knitr::kable(head(DATA[,1:7]))#预览矩阵前7行
p <- upset(DATA, mb.ratio = c(0.3, 0.7), order.by = "freq", 
      nsets = 20, number.angles = 0, point.size = 5, line.size = 1, mainbar.y.label = "Frequency",
      sets.x.label = "Frequency", text.scale = c(2, 2, 2, 2, 2, 2))
p

pdf("Upset.pdf", paper="special",height=10, width=15)
print(p)
dev.off()