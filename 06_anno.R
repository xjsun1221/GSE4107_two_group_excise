rm(list = ls())  
load(file = 'step4output.Rdata')
#富集分析考验网速，因此给大家保存了Rdata
#上课运行示例数据无需修改，在做自己的数据时请注意把本行之后的load()去掉
library(clusterProfiler)
library(dplyr)
library(ggplot2)
source("kegg_plot_function.R")
#source表示运行整个kegg_plot_function.R脚本，里面是一个function
#以up_kegg和down_kegg为输入数据做图

#1.GO database analysis ----

#(1)输入数据
gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
#(2)GO分析，分三部分
#以下步骤耗时很长，实际运行时注意把if后面的括号里F改成T
if(F){
  #细胞组分
  ego_CC <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  #生物过程
  ego_BP <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  #分子功能：
  ego_MF <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  save(ego_CC,ego_BP,ego_MF,file = "ego_GSE42872.Rdata")
}
load(file = "ego_GSE42872.Rdata")
#(3)可视化
#条带图
barplot(ego_CC,showCategory=20)
#气泡图
dotplot(ego_CC)
#下面的图需要映射颜色，设置和示例数据一样的geneList
geneList = deg$logFC
names(geneList)=deg$ENTREZID
geneList = sort(geneList,decreasing = T)
#(3)展示top5通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(ego_CC, categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(ego_CC, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
#Enrichment Map
emapplot(ego_CC)
#(4)展示通路关系
goplot(ego_CC)
#(5)Heatmap-like functional classification
heatplot(ego_CC,foldChange = geneList)
#太多基因就会糊。可通过调整比例或者减少基因来控制。
pdf("heatplot.pdf",width = 14,height = 5)
heatplot(ego_CC,foldChange = geneList)
dev.off()

## 2.KEGG pathway analysis----
#上调、下调、差异、所有基因
#（1）输入数据
gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
#（2）对上调/下调/所有差异基因进行富集分析
#注意这里又有个F
if(F){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.9)
  save(kk.diff,kk.down,kk.up,file = "GSE42872kegg.Rdata")
}
load("GSE42872kegg.Rdata")
#(3)从富集结果中提取出结果数据框
kegg_diff_dt <- kk.diff@result

#(4)按照pvalue筛选通路
#在enrichkegg时没有设置pvaluecutoff，在此处筛选
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)
#(5)可视化
g_kegg <- kegg_plot(up_kegg,down_kegg)

#g_kegg +scale_y_continuous(labels = c(20,15,10,5,0,5))
ggsave(g_kegg,filename = 'kegg_up_down.png')

#gsea作kegg富集分析，可选----
#(1)查看示例数据
data(geneList, package="DOSE")
#(2)将我们的数据转换成示例数据的格式
geneList=deg$logFC
names(geneList)=deg$ENTREZID
geneList=sort(geneList,decreasing = T)
#(3)富集分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
#(4)可视化
kegg_plot(up_kegg,down_kegg)
ggsave('kegg_up_down_gsea.png')

