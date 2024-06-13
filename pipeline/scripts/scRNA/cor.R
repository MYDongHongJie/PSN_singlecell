# module load OESingleCell/2.0.0
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("magrittr"))
library(ggpubr)


setwd("/public/scRNA/works/liuhongyan/Project/scRNA/SA/DR2021/HT2021-15440-m/Further_analysis_2022-02-22/1.Transcriptome/")
##### Mas-AP_vs_WT-AP report #####
trome_kegg=read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/DR2021/HT2021-15440-m/Further_analysis_2022-02-22/file/trome_enrichment-kegg-Mas_AP-vs-WT_AP-Total.xls",row.names=1)

aa <- read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/DR2021/HT2021-15440-m/Further_analysis_2022-02-14/HT2021-15440_Report/8.enrichment/KEGG_enrichment/group_Mas_AP-vs-WT_AP/enrichment-kegg-group_Mas_AP-vs-WT_AP-Total.xls")

### 注意：以上内容，单细胞的数据所用的pathway，符合listhits>2，P<0.05h即可

aa <- aa %>% dplyr::filter(pval < 0.05 & ListHits > 2 )

write.table(aa,"report_8_31_sig_enrichment-kegg-group_Mas_AP-vs-WT_AP-Total.xls",sep="\t",quote=F,row.names=F)

inter_term <- intersect(aa$term, trome_kegg$Term)

aa_in_data <- aa[which(aa$term %in% inter_term),] %>% dplyr::rename(Term = term) %>% dplyr::select(Term,Enrichment_score)

trome_kegg_in_data <- trome_kegg[which(trome_kegg$Term %in% inter_term),] %>% dplyr::select(Term,Enrichment_score)


merge = dplyr::left_join( trome_kegg_in_data, aa_in_data,by="Term")
colnames(merge)=c("Term","bulk","scRNA")
# merge_qc = merge[which(merge[,2] >= 1), ]  


# x <- log10(merge_qc$bulk + 1)
# y <- log10(merge_qc$scRNA +1 )     # 取平均
x_lab <- "Bulk RNA-seq"
y_lab = "scRNA-seq"
suppressPackageStartupMessages(library(ggrepel)) 
d = merge
p <- ggplot(d,aes(bulk,scRNA))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),  axis.text.x = element_text(color="black") , axis.text.y = element_text(color="black") )+  xlab(x_lab)+ylab(y_lab) + labs(title = "Enrichment_score") + theme(plot.title = element_text(hjust = 0.5))
    #geom_label_repel(data=d,aes(d$bulk,d$scRNA),label = Term),size = 2.5, vjust=-0.5, force = 1,show.legend=F)
    

p+geom_point(cex=2,color="#7fc97f")+ geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)+stat_cor(data=d, method = "pearson") +
geom_text_repel(label = d$Term,size = 2.5, vjust=0.5, force = 1,show.legend=F)
# geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)  # 添加趋势线，se=FALSE不显示置信区间
ggsave("report_8_31_Mas-AP_vs_WT-AP_scatter.pdf")
#ggsave("report_8_31_Mas-AP_vs_WT-AP_scatter.png",dpi=1000)
system(paste0("convert ","-verbose -density 1500 -trim report_8_31_Mas-AP_vs_WT-AP_scatter.pdf -quality 100  -flatten ", "report_8_31_Mas-AP_vs_WT-AP_scatter.png")) 

write.table(d,"report_8_31_Mas-AP_vs_WT-AP_data.xls",sep="\t",quote=F,row.names=F)


#### WT-AP+91_vs_WT-AP
trome_kegg=read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/DR2021/HT2021-15440-m/Further_analysis_2022-02-22/file/trome_enrichment-kegg-WT_AP_91-vs-WT_AP-Total.xls",row.names=1)

aa <- read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/DR2021/HT2021-15440-m/Further_analysis_2022-02-14/HT2021-15440_Report/8.enrichment/KEGG_enrichment/group_WT_AP_91-vs-WT_AP/enrichment-kegg-group_WT_AP_91-vs-WT_AP-Total.xls")


aa <- aa %>% dplyr::filter(pval < 0.05 & ListHits > 2 )

write.table(aa,"report_8_31_sig_enrichment-kegg-group_WT_AP_91-vs-WT_AP-Total.xls",sep="\t",quote=F,row.names=F)

inter_term <- intersect(aa$term, trome_kegg$Term)

aa_in_data <- aa[which(aa$term %in% inter_term),] %>% dplyr::rename(Term = term) %>% dplyr::select(Term,Enrichment_score)

trome_kegg_in_data <- trome_kegg[which(trome_kegg$Term %in% inter_term),] %>% dplyr::select(Term,Enrichment_score)


merge = dplyr::left_join( trome_kegg_in_data, aa_in_data,by="Term")
colnames(merge)=c("Term","bulk","scRNA")
# merge_qc = merge[which(merge[,2] >= 1), ]  


# x <- log10(merge_qc$bulk + 1)
# y <- log10(merge_qc$scRNA +1 )     # 取平均
x_lab <- "Bulk RNA-seq"
y_lab = "scRNA-seq"
suppressPackageStartupMessages(library(ggrepel)) 
d = merge
p <- ggplot(d,aes(bulk,scRNA))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),  axis.text.x = element_text(color="black") , axis.text.y = element_text(color="black") )+  xlab(x_lab)+ylab(y_lab) + labs(title = "Enrichment_score") + theme(plot.title = element_text(hjust = 0.5))
    #geom_label_repel(data=d,aes(d$bulk,d$scRNA),label = Term),size = 2.5, vjust=-0.5, force = 1,show.legend=F)
    

p+geom_point(cex=2,color="#7fc97f")+ geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)+stat_cor(data=d, method = "pearson") +
geom_text_repel(label = d$Term,size = 3, vjust=0.5, force = 1,show.legend=F)
# geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)  # 添加趋势线，se=FALSE不显示置信区间
#ggsave("report_8_31_WT_AP_91-vs-WT_AP_scatter.png",dpi=1000)
ggsave("report_8_31_WT_AP_91-vs-WT_AP_scatter.pdf")
system(paste0("convert ","-verbose -density 1500 -trim report_8_31_WT_AP_91-vs-WT_AP_scatter.pdf -quality 100  -flatten ", "report_8_31_WT_AP_91-vs-WT_AP_scatter.png")) 

write.table(d,"report_8_31_WT_AP_91-vs-WT_AP_data.xls",sep="\t",quote=F,row.names=F)


#### 9-13

trome_kegg=read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/DR2021/HT2021-15440-m/Further_analysis_2022-02-22/file/trome_enrichment-kegg-Mas_AP-vs-WT_AP-Total.xls",row.names=1)

aa <- read.delim("/public/scRNA/works/zhangbing/houxuHT2021-15440_陈帅/HT2021-15440_陈帅老师__20210913/subresult_Endothelial_cells/5.enrichment/KEGG_enrichment/group_Mas_AP-vs-WT_AP/enrichment-kegg-group_Mas_AP-vs-WT_AP-Total.xls")

### 注意：以上内容，单细胞的数据所用的pathway，符合listhits>2，P<0.05h即可

aa <- aa %>% dplyr::filter(pval < 0.05 & ListHits > 2 )

write.table(aa,"SA_9_13_sig_enrichment-kegg-group_Mas_AP-vs-WT_AP-Total.xls",sep="\t",quote=F,row.names=F)

inter_term <- intersect(aa$term, trome_kegg$Term)

aa_in_data <- aa[which(aa$term %in% inter_term),] %>% dplyr::rename(Term = term) %>% dplyr::select(Term,Enrichment_score)

trome_kegg_in_data <- trome_kegg[which(trome_kegg$Term %in% inter_term),] %>% dplyr::select(Term,Enrichment_score)


merge = dplyr::left_join( trome_kegg_in_data, aa_in_data,by="Term")
colnames(merge)=c("Term","bulk","scRNA")
# merge_qc = merge[which(merge[,2] >= 1), ]  


# x <- log10(merge_qc$bulk + 1)
# y <- log10(merge_qc$scRNA +1 )     # 取平均
x_lab <- "Bulk RNA-seq"
y_lab = "scRNA-seq"
suppressPackageStartupMessages(library(ggrepel)) 
d = merge
p <- ggplot(d,aes(bulk,scRNA))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),  axis.text.x = element_text(color="black") , axis.text.y = element_text(color="black") )+  xlab(x_lab)+ylab(y_lab) + labs(title = "Enrichment_score") + theme(plot.title = element_text(hjust = 0.5))
    #geom_label_repel(data=d,aes(d$bulk,d$scRNA),label = Term),size = 2.5, vjust=-0.5, force = 1,show.legend=F)
    

p+geom_point(cex=2,color="#7fc97f")+ geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)+stat_cor(data=d, method = "pearson") +
geom_text_repel(label = d$Term,size = 2.5, vjust=0.5, force = 1,show.legend=F)
# geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)  # 添加趋势线，se=FALSE不显示置信区间
## ggsave("SA_9_13_Mas-AP_vs_WT-AP_scatter.png",dpi=1000)
ggsave("SA_9_13_Mas-AP_vs_WT-AP_scatter.pdf")
system(paste0("convert ","-verbose -density 1500 -trim SA_9_13_Mas-AP_vs_WT-AP_scatter.pdf -quality 100  -flatten ", "SA_9_13_Mas-AP_vs_WT-AP_scatter.png")) 

write.table(d,"SA_9_13_Mas-AP_vs_WT-AP_data.xls",sep="\t",quote=F,row.names=F)




#### WT-AP+91_vs_WT-AP
trome_kegg=read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/DR2021/HT2021-15440-m/Further_analysis_2022-02-22/file/trome_enrichment-kegg-WT_AP_91-vs-WT_AP-Total.xls",row.names=1)

aa <- read.delim("/public/scRNA/works/zhangbing/houxuHT2021-15440_陈帅/HT2021-15440_陈帅老师__20210913/subresult_Endothelial_cells/5.enrichment/KEGG_enrichment/group_WT_AP_91-vs-WT_AP/enrichment-kegg-group_WT_AP_91-vs-WT_AP-Total.xls")


aa <- aa %>% dplyr::filter(pval < 0.05 & ListHits > 2 )

write.table(aa,"SA_9_13_sig_enrichment-kegg-group_WT_AP_91-vs-WT_AP-Total.xls",sep="\t",quote=F,row.names=F)

inter_term <- intersect(aa$term, trome_kegg$Term)

aa_in_data <- aa[which(aa$term %in% inter_term),] %>% dplyr::rename(Term = term) %>% dplyr::select(Term,Enrichment_score)

trome_kegg_in_data <- trome_kegg[which(trome_kegg$Term %in% inter_term),] %>% dplyr::select(Term,Enrichment_score)


merge = dplyr::left_join( trome_kegg_in_data, aa_in_data,by="Term")
colnames(merge)=c("Term","bulk","scRNA")
# merge_qc = merge[which(merge[,2] >= 1), ]  


# x <- log10(merge_qc$bulk + 1)
# y <- log10(merge_qc$scRNA +1 )     # 取平均
x_lab <- "Bulk RNA-seq"
y_lab = "scRNA-seq"
suppressPackageStartupMessages(library(ggrepel)) 
d = merge
p <- ggplot(d,aes(bulk,scRNA))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),  axis.text.x = element_text(color="black") , axis.text.y = element_text(color="black") )+  xlab(x_lab)+ylab(y_lab) + labs(title = "Enrichment_score") + theme(plot.title = element_text(hjust = 0.5))
    #geom_label_repel(data=d,aes(d$bulk,d$scRNA),label = Term),size = 2.5, vjust=-0.5, force = 1,show.legend=F)
    

p+geom_point(cex=2,color="#7fc97f")+ geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)+stat_cor(data=d, method = "pearson") +
geom_text_repel(label = d$Term,size = 3, vjust=0.5, force = 1,show.legend=F)
# geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)  # 添加趋势线，se=FALSE不显示置信区间
#ggsave("SA_9_13_WT_AP_91-vs-WT_AP_scatter.png",dpi=1000)
ggsave("SA_9_13_WT_AP_91-vs-WT_AP_scatter.pdf")
system(paste0("convert ","-verbose -density 1500 -trim SA_9_13_WT_AP_91-vs-WT_AP_scatter.pdf -quality 100  -flatten ", "SA_9_13_WT_AP_91-vs-WT_AP_scatter.png")) 

write.table(d,"SA_9_13_WT_AP_91-vs-WT_AP_data.xls",sep="\t",quote=F,row.names=F)

