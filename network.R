####This script is for co-occurrence network analysis
####YJ
####Modified from @Feng Ju @cite Ju F, Xia Y, Guo F, Wang ZP, Zhang T. 2014. 
###@Taxonomic relatedness shapes bacterial assembly in activated sludge of globally distributed wastewater treatment plants.
###@Environmental Microbiology. 16(8):2421-2432



####Co-occurrence Network
#setwd("E:\\Alpine lakes_algae\\STD\\network")

####Step 1. generate the function "co_occurrence_network"
#library(vegan)
#library(igraph)
#library(Hmisc)

co_occurrence_network<-function(matrix,cor.cutoff,p.cutoff){
  matrix1<-matrix
  matrix1[matrix1>0]<-1
  #correlation analysis based on spearman's co-efficient
  matrix.dist<-rcorr(t(matrix),type="spearman")
  ###matrix.dist<-rcorr(t(matrix),type="pearson")
  matrix.cor<-matrix.dist$r
  matrix.cor.p<-matrix.dist$P
  
  #Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
  matrix.cor.p <- p.adjust(matrix.cor.p, method="BH")
  
  #1.Consider positive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
  matrix.cor1<-matrix.cor
  matrix.cor1.p<-matrix.cor.p
  matrix.cor1[which(matrix.cor1 <= cor.cutoff)]=0
  matrix.cor1[which(matrix.cor1.p>p.cutoff)]=0
  # delete those rows and columns with sum = 0
  matrix.cor1<-matrix.cor1[which(rowSums(matrix.cor1)!=1),]
  matrix.cor1<-matrix.cor1[,which(colSums(matrix.cor1)!=0)]
  
  #2.Consider netagive cooccurence at given coefficient (-cor.cutoff) and p-value cutoffs
  ###matrix.cor2<-matrix.cor
  ###matrix.cor2.p<-matrix.cor.p
  ###matrix.cor2[which(matrix.cor2 > (-cor.cutoff))]=0
  ###matrix.cor2[which(matrix.cor2.p>p.cutoff)]=0
  # delete those rows and columns with sum = 0
  ###matrix.cor2<-matrix.cor2[which(rowSums(matrix.cor2)!=0),]
  ###matrix.cor2<-matrix.cor2[,which(colSums(matrix.cor2)!=0)]
  
  #3.Consider both positive and netagive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
  matrix.cor3<-matrix.cor
  matrix.cor3.p<-matrix.cor.p
  matrix.cor3[which(matrix.cor3>=(-cor.cutoff) & matrix.cor3 <= cor.cutoff)]=0
  matrix.cor3[which(matrix.cor3.p>p.cutoff)]=0
  
  # delete those rows and columns with sum = 0
  matrix.cor3<-matrix.cor3[which(rowSums(matrix.cor3)!=1),]
  matrix.cor3<-matrix.cor3[,which(colSums(matrix.cor3)!=0)]
  
  # generate graph using igraph
  g1<-graph.adjacency(matrix.cor1,weight=T,mode="undirected")
  g1<-simplify(g1)
  V(g1)$label <- V(g1)$name
  V(g1)$degree <- degree(g1)
  
  ###g2<-graph.adjacency(matrix.cor2,weight=T,mode="undirected")
  ###g2<-simplify(g2)
  ###V(g2)$label <- V(g2)$name
  ###V(g2)$degree <- degree(g2)
  
  g3<-graph.adjacency(matrix.cor3,weight=T,mode="undirected")
  g3<-simplify(g3)
  V(g3)$label <- V(g3)$name
  V(g3)$degree <- degree(g3)
  
  # append the output into results
  result<-list()
  result$matrix.cor<-matrix.cor
  result$matrix.cor.p<-matrix.cor.p
  
  result$matrix.cor1<-matrix.cor1
  result$graph1<-g1
  
  ###result$matrix.cor2<-matrix.cor2
  ###result$graph2<-g2
  
  result$matrix.cor3<-matrix.cor3
  result$graph3<-g3
  return(result)
}


####Step 2. generate network####
################## OTU filtering, network generation, topological analysis and export OTU table ###############################
library(igraph)
library(Hmisc)

Abu <- read.csv("Aug_C_2010_S.csv",header = T, row.names = 1)
Abu <- t(Abu)#如果不是列为sample，行为otu，需要此步转换
Abu<-as.matrix(Abu)

###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present)
table<-Abu
table[table>0]<-1
table.generalist<-Abu[which(rowSums(table)>=1),]
Abu<-table.generalist

###2. Creating gml files of network (to be visulized in Gephi or Cytoscape)
occor<-co_occurrence_network(Abu,0.6,0.05)  ## cutoffs for correlation coefficient and P-value

write.graph(occor$graph1,'Pos0.6-Aug_C_2010_S_a.gml',format='gml')    #network file for positive association
#write.graph(pattern$graph2,'Neg0.6-NW.gml',format='gml')   #network file for negative association (if any)
write.graph(occor$graph3,'PosNeg0.6-Aug_C_2010_S_a.gml',format='gml') #network file for all association

occor.r = occor$matrix.cor # ??????????????????R???
occor.p = occor$matrix.cor.p # ??????????????????p???

igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph
# NOTE:????????????weighted=NULL,????????????????????????????????????????????????????????????????????????,???????????????????????????????????????
# ?????????????????????????????????
# occor.r[occor.r!=0] = 1
# igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)

# ????????????????????????,????????????????????????
# remove isolated nodes,??????????????????otu??????????????????otu ?????????,????????????????????????
bad.vs = V(igraph)[degree(igraph) == 0]
igraph = delete.vertices(igraph, bad.vs)
igraph

# ???igraph weight???????????????igraph.weight
igraph.weight = E(igraph)$weight

# ???????????????igraph???weight??????,?????????????????????layout??????????????????
E(igraph)$weight = NA

# ????????????
# ?????????????????????,?????????????????????????????????????????????,?????????????????????????????????

#================================================

####Step 3. parameters####
###3. Calculating network topological properties
g<-occor$graph1   ###positive network
#g<-occor$graph1   ###negative network

c <- cluster_walktrap(g)
# Global toplogical features
modularity(c)
md <- modularity(g, membership(c), weights = NULL)
cc <- transitivity(g, vids = NULL,
                   weights = NULL)
spl <- average.path.length(g, directed=FALSE, unconnected=TRUE)
gd  <- graph.density(g, loops=FALSE)
nd  <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NA)

node.degree <- degree(g, v = V(g), mode="all")
ad  <- mean(node.degree)

e <- ecount(g)
v <- vcount(g)
global.topology <- data.frame(md, e,v,cc,spl,gd,nd,ad)
write.csv(global.topology, file="Pos0.6-NW-global.topology_Aug_C_2010_S_a.csv")

# Node toplogical features
betweenness.centrality <- betweenness(g, v=V(g), 
                                      directed = FALSE, weights = NA,
                                      nobigint = TRUE, normalized = FALSE)
closeness.centrality <- closeness(g, vids = V(g),
                                  weights = NA, normalized = FALSE)
node.transitivity <- transitivity(g, type = c("local"), vids = NULL,
                                  weights = NA)

node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
write.csv(node.topology, file="Pos0.6-NW-node.topology_Aug_C_2010_S_a.csv")

# Ploting node degreee distribution in a log-log plot
degree.df <- data.frame(table(degree=factor(node.degree, levels=seq_len(max(node.degree)))))
degree.df$degree <- as.numeric(as.character(degree.df$degree))

#4. Creating an abundance table for OTUs present in the positive and negative network
my.list1 <- row.names(occor$matrix.cor1)
###my.list2 <- row.names(pattern$matrix.cor2)

logical1 <- row.names(Abu)  %in% my.list1
###logical2 <- row.names(Abu)  %in% my.list2

tab.subset1 <- subset(Abu,logical1)
###tab.subset2 <- subset(Abu,logical2)

write.csv(tab.subset1,'Pos0.6-NW_Aug_C_2010_S_a.csv',row.names = F)
###write.table(tab.subset2,'Neg0.6-NW.txt',sep="\t")

####作图----------------####


# 读取otu-sample矩阵，行为sample，列为otu
otu = read.csv("Aug_C_2010_S.csv",header = T, row.names = 1)
#otu <- t(otu)#若不是，则转换

# 计算OTU间两两相关系数矩阵
# 数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue, 但p值需要借助其他函数矫正
occor = corr.test(otu,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occor.r = occor$r # 取相关性矩阵R值
occor.p = occor$p # 取相关性矩阵p值

# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 

# 构建igraph对象
igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph
# NOTE:可以设置weighted=NULL,但是此时要注意此函数只能识别相互作用矩阵内正整数，所以应用前请确保矩阵正确。
# 可以按下面命令转换数据
# occor.r[occor.r!=0] = 1
# igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)

# 是否去掉孤立顶点，根据自己实验而定
# remove isolated nodes，即去掉和所有otu均无相关性的otu 可省略，前期矩阵已处理过
bad.vs = V(igraph)[degree(igraph) == 0]
igraph = delete.vertices(igraph, bad.vs)
igraph

# 将igraph weight属性赋值到igraph.weight
igraph.weight = E(igraph)$weight

# 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
E(igraph)$weight = NA

# 简单出图
# 设定随机种子数，后续出图都从同一随机种子数出发，保证前后出图形状相对应
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# 如果构建网络时，weighted=NULL,此步骤不能统计
sum(igraph.weight>0)# number of postive correlation
sum(igraph.weight<0)# number of negative correlation

# set edge color，postive correlation 设定为red, negative correlation设定为blue
E.color = igraph.weight
E.color = ifelse(E.color>0, "red",ifelse(E.color<0, "blue","grey"))
E(igraph)$color = as.character(E.color)

# 改变edge颜色后出图
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

# 可以设定edge的宽 度set edge width，例如将相关系数与edge width关联
#E(igraph)$width = abs(igraph.weight)*4

# 改变edge宽度后出图
#set.seed(123)
#plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

# 模块性 modularity
fc = cluster_fast_greedy(igraph,weights =NULL)# cluster_walktrap cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity = modularity(igraph,membership(fc))
# 按照模块为节点配色
comps = membership(fc)
colbar = topo.colors(max(comps))
V(igraph)$color = colbar[comps] 

set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
#================================================
pdf("Aug_C_2010_S_a.pdf")
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
dev.off()

pdf("Aug_C_2010_Sname_a.pdf")
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
dev.off()

