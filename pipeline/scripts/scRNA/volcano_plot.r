library(ggrepel)
DEG=read.delim("/public/scRNA/works/tangxuan/HT2019-6928-human/houxu20210725/2.sub_Fibroblasts_supply/Diffexp/clusters_5-vs-6-all_diffexp_genes_anno.xls",header = T,sep='\t')
rownames(DEG)=DEG[,1]
label_gene = read.delim("/public/scRNA/works/tangxuan/HT2019-6928-human/houxu20210725/2.sub_Fibroblasts_supply/Diffexp/top20_clusters_5-vs-6_genes.xls",header = T,sep='\t')
DEG$label <-''
label= unique(as.vector(label_gene$gene))
label = c(label,c("HLA-C","HLA-B","HLA-A","LGALS9","B2M"))
DEG[label,]$label = as.character(label) 

DEG$pvalue[which(DEG$pvalue < 5E-308 )] = 5E-300   ##极值处理

g = ggplot(data=DEG,aes(x=DEG$log2FoldChange, y=-log10(DEG$pvalue), color = DEG$log2FoldChange )) +
  geom_point( size=1.5) + guides(fill=guide_legend(title=expression(paste(log[2], " fold change")))) + scale_color_gradientn(colours=c("blue","yellow","red"), name = expression(paste(log[2], " fold change"))) +
  theme_set(theme_set(theme_bw(base_size=15)))+
  xlab(expression(paste(log[2], " fold change"))) + ylab(expression(paste("-", log[10], "pvalue"))) +
  geom_text_repel(aes(label = label), size = 4,vjust=-0.5, alpha=1,color="black")

ggsave("Volcano2.pdf",width=8)
ggsave("Volcano2.png",width=8)

g = ggplot(data=DEG,aes(x=DEG$log2FoldChange, y=-log10(DEG$pvalue), color = DEG$log2FoldChange )) +
  geom_point( size=1.5) + guides(fill=guide_legend(title=expression(paste(log[2], " fold change")))) + scale_color_gradientn(colours=c("blue","yellow","red"), name = expression(paste(log[2], " fold change"))) +
  theme_set(theme_set(theme_bw(base_size=15)))+
  xlab(expression(paste(log[2], " fold change"))) + ylab(expression(paste("-", log[10], "pvalue"))) 
  #+geom_text_repel(aes(label = label), size = 4,vjust=-0.5, alpha=1,color="black")

ggsave("Volcano3.pdf",width=8)
ggsave("Volcano3.png",width=8)