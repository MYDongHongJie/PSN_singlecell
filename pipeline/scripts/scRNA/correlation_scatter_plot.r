library(ggpubr)
bulk=read.delim("fpkm_symbol.xls",row.names=1)

fls2_CK_sc = as.data.frame(rowMeans(as.matrix(fls2_CK_ob$RNA@data)))  %>% rownames_to_column(var="gene")
fls1_flg22_bulk = as.data.frame(rowMeans(bulk[,c("fls1_flg22_1","fls1_flg22_2")])) %>% rownames_to_column(var="gene")
merge = left_join(fls1_flg22_bulk, fls2_CK_sc, by="gene")
colnames(merge)=c("gene","bulk","scRNA")
merge_qc = merge[which(merge[,2] >= 1), ]      # 过滤bulk表达小于1的基因

x <- log10(merge_qc$bulk + 1)
y <- log10(merge_qc$scRNA +1 )     # 取平均
x_lab <- "Bulk RNA-seq"
y_lab = "scRNA-seq"

d = merge_qc
p <- ggplot(d,aes(x,y))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),  axis.text.x = element_text(color="black") , axis.text.y = element_text(color="black") )+  xlab(x_lab)+ylab(y_lab)

p+geom_point(cex=1,color="grey")+ geom_smooth(se=FALSE,color = "black",size=0.5)+stat_cor(data=d, method = "pearson")
# geom_smooth(method = "lm",se=FALSE,color = "black",size=0.5)  # 添加趋势线，se=FALSE不显示置信区间
ggsave("fls2_CK_scRNA_fls1_flg22_bulk_cor.png")
ggsave("fls2_CK_scRNA_fls1_flg22_bulk_cor.pdf")
