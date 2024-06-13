# module purge 
# module load OESingleCell/2.0.0

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

seurat_ob = readRDS( "/public/scRNA/works/liuhongyan/Project/scRNA/DR2021/HT2021-10023-4-h-nu/result/Further_analysis_2021-12-15/1.new_celltype/new_celltype/new_celltype2021-12-15.rds" )

if ( seurat_ob@version < 3 ){
    seurat_ob = UpdateSeuratObject(seurat_ob)
}

ident2use="clusters"
cluster_list="1,2,9,11,13,14,7,3,6,17,4,5,10,16,18,15,19,20"
cluster_list = unlist(strsplit( cluster_list,",",perl = T))
    seurat_ob = SubsetData( seurat_ob, subset.name = ident2use, accept.value = cluster_list)

    seurat_ob@meta.data[[ident2use]]=factor(seurat_ob@meta.data[[ident2use]],levels = sort(unique(seurat_ob@meta.data[[ident2use]])))


DefaultAssay(seurat_ob) = "RNA"

genelist = read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/DR2021/HT2021-10023-4-h-nu/result/Further_analysis_2022-03-24/2.heatmap/genelist.xls", sep="\t", header = T)


markers2vis = CaseMatch(search = as.vector(as.vector(genelist[,1])),match = rownames(seurat_ob))
filtered_gene = genelist[! genelist[,1] %in% names(markers2vis ),1]
if(length(filtered_gene)!=0){
    filtered_gene = as.data.frame(filtered_gene)
    colnames(filtered_gene) = "Gene"
    write.table(filtered_gene,file.path(output_dir,"filtered_gene.xls"),quote = F,row.names=F)
    print("There are some mismatched gene symbol, Please check filtered_gene.xls for the genename.")
}


# count = as.matrix(seurat_ob@data[markers2vis,])
count = as.matrix(GetAssayData(seurat_ob, slot = "data")[markers2vis,])
meta.data = seurat_ob@meta.data
collapseby = "clusters"

meta.data$id = rownames(meta.data)
collapsed_count = vector()
if ( !collapseby %in% colnames(meta.data) ){
    stop("NO specified column found!")
}

collapsed_group = meta.data %>% group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))
# if (collapseby == "clusters")  collapsed_group$clusters = paste("cluster",collapsed_group$clusters,sep="_")



for ( cells in collapsed_group$cellid ){
    samplex = unlist(strsplit(cells, ",", perl =T))
    collapsed_count= cbind(collapsed_count,rowMeans( count[,samplex,drop=F] ))
}
collapsed_count = as.matrix( collapsed_count )
collapsed_group = as.data.frame(collapsed_group)
colnames(collapsed_count) = as.matrix(collapsed_group[,1])
# if (!is.null(opt$order)){

    order="1,2,9,11,13,14,7,3,6,17,4,5,10,16,18,15,19,20"
    order = unlist(strsplit( order,",",perl = T))
    #if (collapseby == "clusters") order = paste("cluster",order,sep="_")
    collapsed_count=collapsed_count[,order]
# }


output_dir="./"
data = tibble::rownames_to_column(as.data.frame(collapsed_count),"GeneID")
write.table(data, file.path(output_dir,"heatmap_count.xls"),quote = F, row.names = F, sep = "\t")
if (dim(collapsed_count)[2]>2) {
    data_scaled = tibble::rownames_to_column(as.data.frame(pheatmap:::scale_rows(collapsed_count)),"GeneID")
    write.table(data_scaled, file.path(output_dir,"heatmap_count_scaled.xls"),quote = F, row.names = F, sep = "\t")
}


ind <- apply(collapsed_count, 1, mean) > 0
collapsed_count_filter <- collapsed_count[ind, ]


CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

## 注释 col 和 clusters 匹配
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

# > color_use
# $clusters
#         1         2         3         4         5         6         7         8
# "#7fc97f" "#beaed4" "#fdc086" "#386cb0" "#f0027f" "#a34e3b" "#666666" "#1b9e77"


palette <- colorRampPalette(c("#406AA8", "white", "#D91216"))(n=299)
annotation_row=NA
annotation_col = read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/DR2021/HT2021-10023-4-h-nu/result/Further_analysis_2022-03-24/2.heatmap/col_anno.xls", row.names=1)


## 注释 col 和 clusters 匹配
color_map = CustomCol2(1:length(unique(annotation_col[,1])))
# if (is.factor(col_anno[,1])) {
#     names <- levels(col_anno[,1])
# } else {
    names <- unique(annotation_col[,1])
# }
names(color_map) <- names
name=colnames(annotation_col)
color_use = list()
for (i in names(color_map)){
    color_use[[name]][[i]]=color_map[[i]]
}


opt=list()
gaps_row=NA
opt$showname
## parsing show rowname
if (is.null(opt$showname))  showname=ifelse(dim(collapsed_count_filter)[1]>100,FALSE,TRUE) else showname = opt$showname
cellwidth=36
cellheight=ifelse(showname,12,8) # if showname =F, set the panel height to 576/72point = 8 inches
opt$rowcluster=F
opt$colcluster=F
p = pheatmap(log2(collapsed_count_filter+0.0001),
        color=palette,
        cex=1,
        border=F,
        angle_col=45,
        treeheight_row=36, treeheight_col=36,
        lwd=1.1,
        cellheight=cellheight, cellwidth=cellwidth,
        scale=ifelse(dim(collapsed_count_filter)[2]==2,"none","row"),
        show_rownames=showname,
        #gaps_row=gaps_row,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = color_use,
        # display_numbers = display_number,
        cluster_rows=opt$rowcluster ,cluster_cols=opt$colcluster)
ggsave(file.path(output_dir,"heatmap.pdf"), plot= p, 
    width=(36*dim(collapsed_count_filter)[2]+144)/72+4, 
    height = ifelse(showname,(12*dim(collapsed_count_filter)[1]+108)/72,(8*dim(collapsed_count_filter)[1]+108)/72))
ggsave(file.path(output_dir,"heatmap.png"), plot= p, dpi=1000,
    width=(36*dim(collapsed_count_filter)[2]+144)/72+4,
    height = ifelse(showname,(12*dim(collapsed_count_filter)[1]+108)/72,(8*dim(collapsed_count_filter)[1]+108)/72))

