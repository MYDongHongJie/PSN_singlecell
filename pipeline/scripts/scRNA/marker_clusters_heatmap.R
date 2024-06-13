suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("reticulate"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("OESingleCell"))

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}


opt=list()
seurat_ob = readRDSMC( "/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2018H2930S_HT2020-13912-h/HT2020-13912-human/houxu20201210_cluster5/singlecell_object.clustering_resolution0.4.rds" )

genelist="gene100.xls"
genelist <- read.table(genelist,header=T,sep=",")
regulonAUC_mat = seurat_ob[["RNA"]]@data
###R  grep  用法 去除 extended 
regulonAUC_mat_out = regulonAUC_mat[as.vector(genelist[,1]),]
#=================================================================================
# Step2. Inference of co-expression modules
#=================================================================================

cellInfo = seurat_ob@meta.data
col_anno = as.data.frame(seurat_ob@meta.data) %>% rownames_to_column(var="barcodes")
opt$groupby="clusters"


if ( dim(cellInfo)[1] > 60000 ){
    col_anno = col_anno[,c("barcodes",opt$groupby)] %>% group_by(get(opt$groupby)) %>% sample_frac(60000/length(col_anno$barcodes))
}else{
    col_anno = col_anno[,c("barcodes",opt$groupby)]
}


col_anno = col_anno %>% arrange(get(opt$groupby)) %>% column_to_rownames(var="barcodes")
regulonAUC_plotdata = regulonAUC_mat_out[,rownames(col_anno)]
bks <- unique(c(seq(-2.5,0, length=100),  seq(0,2.5, length=100)))
color_map = CustomCol2(1:length(unique(col_anno[[opt$groupby]])))

if (is.factor(col_anno[[opt$groupby]])) {
    names <- levels(col_anno[[opt$groupby]])
} else {
    names <- unique(col_anno[[opt$groupby]])
}
names(color_map) <- names

color_use = list()
for (i in names(color_map)){
    color_use[[opt$groupby]][[i]]=color_map[[i]]
}

library(dplyr)
# # 对前 7 列求和排序，以及对后 7 列求和排序
# aa <- as.data.frame(t(test_data)) %>% tibble::rownames_to_column("barcodes") %>% mutate(sum_before = rowSums(.[2:8]))  %>% mutate(sum_after = rowSums(.[9:15])) %>% mutate(mean_before = rowMeans(.[2:8]))  %>% mutate(mean_after = rowMeans(.[9:15])) %>% tibble::column_to_rownames("barcodes")
# aa <- as.matrix(aa)
# bb <- aa[order(-aa[,15],-aa[,16],aa[,17],aa[,18]),]
# bb<- as.data.frame(t(bb[,1:(ncol(bb)-4)]))
aa <- as.data.frame(t( regulonAUC_plotdata)) %>% tibble::rownames_to_column("barcodes") %>% mutate(diff = rowMeans(.[2:(1+nrow(regulonAUC_plotdata)/2)])-rowMeans(.[(2+nrow(regulonAUC_plotdata)/2):(1+nrow(regulonAUC_plotdata))]))   %>% tibble::column_to_rownames("barcodes")
aa <- as.matrix(aa)
## 降序
bb <- aa[order(-aa[,nrow(regulonAUC_plotdata)+1]),]

### 根据差值分成3个类群 ### 
library("ggplot")
data1 <- as.data.frame(bb)
diff <- data1$diff
data <- as.data.frame(diff)
data$type="diff"
p<-ggplot(data=data, aes(x=type,y=diff))+geom_boxplot(aes(fill="diff")) + xlab("") + ylab("diff")
ggsave("c5_diff.png", plot = p, dpi=1000,width = 7, height = 7)
p<-ggplot(data=data, aes(x=type,y=diff))+geom_point(aes(fill="diff")) + xlab("") + ylab("diff")
ggsave("c5_diff_p.png", plot = p, dpi=1000,width = 7, height = 7)

p<-ggplot(data=data, aes(x=type,y=diff))+geom_violin(aes(fill="diff")) + xlab("") + ylab("diff")
ggsave("c5_diff_v.png", plot = p, dpi=1000,width = 7, height = 7)

left <- subset(data1,diff> 0.2 )
right <- subset(data1,diff < (-0.2) )
union =union(rownames(left),rownames(right))
middle_c <- setdiff(rownames(data1),union)
middle <- data1[middle_c,]


bb<- as.data.frame(t(bb[,1:(ncol(bb)-1)]))

output_dir="markers100"
dir.create(output_dir)
p1 <- pheatmap::pheatmap( bb,
                    scale = "row",
                    cluster_cols=F,
                    cluster_rows=F,
                    show_colnames= F,
                    color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(200),
                    annotation_col = col_anno,
                    annotation_colors = color_use,
                    treeheight_col=10, 
                    border_color=NA,breaks=bks,fontsize_row=5)

pdf(file.path(output_dir, "marker_gene_heatmap.pdf"))
p1
dev.off()
png(file.path(output_dir, "marker_gene_heatmap.png"),res=300,height=2400,width=2400)
p1
dev.off()

# left<- as.data.frame(t(left[,1:(ncol(left)-1)]))
# dir.create(output_dir)
# p1 <- pheatmap::pheatmap(left,
#                     scale = "row",
#                     cluster_cols=F,
#                     cluster_rows=F,
#                     show_colnames= F,
#                     color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(200),
#                     annotation_col = col_anno,
#                     annotation_colors = color_use,
#                     treeheight_col=10, 
#                     border_color=NA,breaks=bks,fontsize_row=5)

# pdf(file.path(output_dir, "marker_gene_heatmap_left.pdf"))
# p1
# dev.off()
# png(file.path(output_dir, "marker_gene_heatmap_left.png"),res=300,height=2400,width=2400)
# p1
# dev.off()

# right<- as.data.frame(t(right[,1:(ncol(right)-1)]))
# p1 <- pheatmap::pheatmap(right,
#                     scale = "row",
#                     cluster_cols=F,
#                     cluster_rows=F,
#                     show_colnames= F,
#                     color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(200),
#                     annotation_col = col_anno,
#                     annotation_colors = color_use,
#                     treeheight_col=10, 
#                     border_color=NA,breaks=bks,fontsize_row=5)

# pdf(file.path(output_dir, "marker_gene_heatmap_right.pdf"))
# p1
# dev.off()
# png(file.path(output_dir, "marker_gene_heatmap_right.png"),res=300,height=2400,width=2400)
# p1
# dev.off()

# middle<- as.data.frame(t(middle[,1:(ncol(middle)-1)]))
# p1 <- pheatmap::pheatmap(middle,
#                     scale = "row",
#                     cluster_cols=F,
#                     cluster_rows=F,
#                     show_colnames= F,
#                     color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(200),
#                     annotation_col = col_anno,
#                     annotation_colors = color_use,
#                     treeheight_col=10, 
#                     border_color=NA,breaks=bks,fontsize_row=5)

# pdf(file.path(output_dir, "marker_gene_heatmap_middle.pdf"))
# p1
# dev.off()
# png(file.path(output_dir, "marker_gene_heatmap_middle.png"),res=300,height=2400,width=2400)
# p1
# dev.off()

p1 <- pheatmap::pheatmap( bb,
                    scale = "row",
                    cluster_cols=F,
                    cluster_rows=F,
                    show_colnames= F,
                    color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(200),
                    annotation_col = col_anno,
                    annotation_colors = color_use,
                    treeheight_col=10, 
                    gaps_col=c(nrow(left),(nrow(left)+nrow(middle))),
                    border_color=NA,breaks=bks,fontsize_row=5)

pdf(file.path(output_dir, "marker_gene_heatmap_col.pdf"))
p1
dev.off()
png(file.path(output_dir, "marker_gene_heatmap_col.png"),res=300,height=2400,width=2400)
p1
dev.off()

### 分组进行差异分析 left，right，
left_ob <- SubsetData(seurat_ob,cells=rownames(left))

right_ob <- SubsetData(seurat_ob,cells=rownames(right))

seurat_ob@meta.data$new_clusters="middle"
seurat_ob@meta.data[which(seurat_ob@meta.data$orig.ident %in% left_ob@meta.data$orig.ident),"new_clusters"] = "left"
seurat_ob@meta.data[which(seurat_ob@meta.data$orig.ident %in% right_ob@meta.data$orig.ident),"new_clusters"] = "right"

saveRDSMC(seurat_ob,file.path(output_dir, "new_clusters_seurat.rds"))