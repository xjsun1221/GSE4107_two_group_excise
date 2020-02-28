rm(list = ls())  
load(file = "step2output.Rdata")
#输入数据：exp和group_list

#Principal Component Analysis
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials

exp[1:4,1:4]
{
  dat=as.data.frame(t(exp))
  library(FactoMineR)
  library(factoextra) 
  dat.pca <- PCA(dat, graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  ggsave('all_samples_PCA.png')
}


#热图 
cg=names(tail(sort(apply(exp,1,sd)),1000))
n=exp[cg,]

#绘制热图
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
library(pheatmap)
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row")

dev.off()
