#此脚本无需修改，认清输入数据即可
rm(list = ls()) 
load(file = "step4output.Rdata")
#1.火山图,输入数据是deg
{
  library(dplyr)
  dat <- mutate(deg,v=-log10(P.Value))
  head(dat)
  library(ggpubr)
  ggscatter(dat, x = "logFC", y = "v", 
            color = "change",size = 0.5,
            label = "symbol", repel = T,
            #label.select = dat$symbol[1:30] ,
            label.select = c('CD36','DUSP6'), #挑选一些基因在图中显示出来
            palette = c("#00AFBB", "#999999", "#FC4E07"),
            ylab = "-log10p.value")
  ggsave("volcano.png")
}
dev.off()
#palette里面是RGB颜色编号，想换其他颜色可以尝试https://www.114la.com/other/rgb.htm

###### 取差异基因做热图
#输入数据是exp表达矩阵的子集
{
  load(file = 'step2output.Rdata')
  x=deg$logFC 
  names(x)=deg$probe_id 
  #上调下调基因各100个（可以自定义）
  # cg=c(names(head(sort(x),100)),names(tail(sort(x),100)))
  cg = deg$probe_id[deg$change!="stable"]
  n=exp[cg,]
}
#作热图
{
  library(pheatmap)
  annotation_col=data.frame(group=group_list)
  rownames(annotation_col)=colnames(n) 
  #保存
  pdf(file = "heatmap.pdf")
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           scale = "row",
           #cluster_cols = F, 
           annotation_col=annotation_col) 
  dev.off()
}
