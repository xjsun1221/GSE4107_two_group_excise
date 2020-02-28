rm(list = ls()) 
load(file = "step2output.Rdata")
#差异分析，用limma包来做
#需要表达矩阵和group_list，不需要改
{
  library(limma)
  design=model.matrix(~group_list)
  fit=lmFit(exp,design)
  fit=eBayes(fit)
  deg=topTable(fit,coef=2,number = Inf)
}

#为deg数据框添加几列
{
  #1.加probe_id列，把行名变成一列
  library(dplyr)
  deg <- mutate(deg,probe_id=rownames(deg))
  head(deg)
  #2.加symbol列，火山图要用
  deg <- inner_join(deg,ids,by="probe_id")
  head(deg)
  #按照symbol列去重复
  deg <- deg[!duplicated(deg$symbol),]
}

#3.加change列,下面两行是阈值，可修改
logFC_t=mean(deg$logFC)+2*sd(deg$logFC)
logP_t = 0.01

test1 = deg$P.Value < logP_t
test2 = deg$logFC < -logFC_t
test3 = deg$logFC > logFC_t

{
  change = ifelse(test1 & test2 , 
               "down" , 
               ifelse(test1 & test3 ,
                      "up", 
                      "stable"))
  deg <- mutate(deg,change)
  #4.加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  s2e <- bitr(deg$symbol, fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)#人类
  #其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
  deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
}

save(group_list,deg,file = "step4output.Rdata")
