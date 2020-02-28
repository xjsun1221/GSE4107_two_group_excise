#制作string的输入数据
load("step4output.Rdata")
gene_up= deg[deg$change == 'up','symbol'] 
gene_down=deg[deg$change == 'down','symbol'] 
write.table(gene_up,
            file="upgene.txt",
            row.names = F,
            col.names = F,
            quote = F)
write.table(gene_down,
            file="downgene.txt",
            row.names = F,
            col.names = F,
            quote = F)
write.table(deg$symbol[1:200],
            file="diffgene.txt",
            row.names = F,
            col.names = F,
            quote = F)
#--------------

#准备cytoscape的输入文件
tsv = read.table("string_interactionsdiff.tsv",comment.char = "!",header = T)
tsv2 = tsv[,c(1,2,ncol(tsv))]
head(tsv2)
write.table(tsv2,
            file = "cyto.txt",
            sep = "\t",
            quote = F,
            row.names = F)

p = deg[deg$change != "stable",c("symbol","logFC","P.Value")]
head(p)
write.table(p,
            file = "deg.txt",
            sep = "\t",
            quote = F,
            row.names = F)
