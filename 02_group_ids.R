#group_list和芯片注释，每次都需要改
rm(list = ls())  
load(file = "step1output.Rdata")
#------------------

#第三类，ifelse
library(stringr)
group_list=ifelse(str_detect(pd$title,"healthy"),"normal","tumor")
#设置参考水平，对照在前，处理在后
group_list = factor(group_list,
                    levels = c("normal","tumor"))

#-----------------

#芯片注释，查找芯片平台对应的包,到此脚本中替换
if(F){
  #http://www.bio-info-trainee.com/1399.html
  #hgu133plus2
  if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
  library(hgu133plus2.db)
  ls("package:hgu133plus2.db")
  ids <- toTable(hgu133plus2SYMBOL)
  head(ids)
}else if(T){
  if(!file.exists(paste0(gpl,".soft"))) getGEO(gpl,destdir = ".")
  ids = data.table::fread(paste0(gpl,".soft"),header = T,skip = "ID",data.table = F)
  ids = ids[,c("ID","Gene Symbol")]
  colnames(ids) = c("probe_id","symbol")
  table(ids$symbol != "" & (!str_detect(ids$symbol," /// ")))
  ids = ids[ids$symbol != "" & (!str_detect(ids$symbol," /// ")),]
}else if(F){
  ids = idmap(gpl,type = "bioc")
}

save(exp,group_list,ids,file = "step2output.Rdata")


