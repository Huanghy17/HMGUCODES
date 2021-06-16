#Multiple groups
kruskal.test(dis~group, data = dat)

#if the overall difference is significant, use Wilcoxon and two-sided test for two groups test
wilcox.test(dis_env1, dis_env2, alternative = 'two.sided')
wilcox.test(dis_env1, dis_env3, alternative = 'two.sided')
