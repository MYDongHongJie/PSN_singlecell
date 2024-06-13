### hclust 计算距离是R自带的，ape包用来作图
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("dplyr"))

rm(list=ls())
### 调色板
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

seurat_ob <- readRDSMC("/public/scRNA/works/liuhongyan/Project/scRNA/done_report/HT2020-16492-chicken/Further_analysis_2021-11-01/new_celltype/new_celltype/new_celltype2021-11-28.rds")


subset_ob <- subset(seurat_ob ,subset=clusters %in% names(table(seurat_ob@meta.data$clusters))[-17] )

seurat_ob <- subset_ob
markers2vis=rownames(seurat_ob)
# count = as.matrix(seurat_ob@data[markers2vis,])
count = as.matrix(GetAssayData(seurat_ob, slot = "data")[markers2vis,])
meta.data = seurat_ob@meta.data
opt=list()
opt$collapseby="clusters"
collapseby = opt$collapseby

meta.data$id = rownames(meta.data)
collapsed_count = vector()
if ( !collapseby %in% colnames(meta.data) ){
    stop("NO specified column found!")
}

collapsed_group = meta.data %>% dplyr::group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))
if (collapseby == "clusters")  collapsed_group$clusters = paste("cluster",collapsed_group$clusters,sep="_")

for ( cells in collapsed_group$cellid ){
    samplex = unlist(strsplit(cells, ",", perl =T))
    collapsed_count= cbind(collapsed_count,rowMeans( count[,samplex,drop=F] ))
}
collapsed_count = as.matrix( collapsed_count )
collapsed_group = as.data.frame(collapsed_group)
colnames(collapsed_count) = as.matrix(collapsed_group[,1])

if (!is.null(opt$order)){
    order = unlist(strsplit( opt$order,",",perl = T))
    if (collapseby == "clusters") order = paste("cluster",order,sep="_")
    collapsed_count=collapsed_count[,order]
}

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

output_dir="./"
data <- as.data.frame(collapsed_count)
data = tibble::rownames_to_column( data,"GeneID")
write.table(data, file.path(output_dir,"clusters_mean_data.xls"),quote = F, row.names = F, sep = "\t")




data <- as.data.frame(collapsed_count)
# data <- data[colnames(data)[-7]]

res = hclust(dist(t(data)),method="ward.D")   # x为各cluster data平均表达量

######画图1
pdf("p.pdf")
plot(res)
dev.off()






## c( "#7FC97F","#BEAED4","#7FC97F","#FDC086","#7FC97F","#F0027F","#386CB0","#7FC97F","#BEAED4")
######画图2
library(ape)
pdf("plot.pdf",height=9,,width=8)
plot(as.phylo(res), type= "fan",  tip.color=CustomCol2( 1:length(colnames(data)) ) ,cex=1) 
dev.off()

# png("plot.png",height=4,,width=3,units = "px",res=96)
# plot(as.phylo(res), type= "fan",  tip.color=CustomCol2( 1:length(colnames(data)) ) ,cex=1) 
# dev.off()
