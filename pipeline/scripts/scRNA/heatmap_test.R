suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages( library(ComplexHeatmap) )

rds="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2021-12636-mouse/Further_analysis_2022-01-11/2.VEC/subtype/AEC/singlecell_object.clustering_resolution0.2.rds"

setwd("/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2021-12636-mouse/Further_analysis_2022-02-07/five_vecsubtype/")
seurat_ob = readRDS(rds)

if ( seurat_ob@version < 3 ){
    seurat_ob = UpdateSeuratObject(seurat_ob)
}

opt=list()
opt$outdir="heatmap"
if ( is.null(opt$outdir) ){
    output_dir = getwd()
} else {
    if ( file.exists(opt$outdir) ){
        output_dir = opt$outdir
    } else {
        output_dir = opt$outdir
        dir.create(output_dir,recursive=T)
    }
}

opt$assay="integrated"
if ( is.null(opt$assay) ){
    assay = "RNA"
} else {
    assay = opt$assay
}

opt$topn=100
if ( is.null( opt$topn )){
    topn = 25
} else {
    topn = opt$topn
}


opt$gaps_row=NULL
if ( !is.null( opt$gaps_row )){
    gaps_row = as.numeric(unlist(strsplit( opt$gaps_row,",",perl = T)))
} else {
    gaps_row = opt$gaps_row
}



DefaultAssay(seurat_ob) = assay

opt$genelist="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2021-12636-mouse/Further_analysis_2022-01-11/2.VEC/Marker/all_markers_for_each_cluster_anno.xls"
genelist = read.delim(opt$genelist, sep="\t", header = T)
opt$topn=100
opt$collapseby="new_celltypes"
opt$topby="avg_logFC"

## p_val,
## ,desc(gene_diff

if (dim(genelist)[2]>1){
    if ( "cluster" %in% colnames(genelist) ) {
        markers2vis = genelist
        #markers2vis[["cluster"]] = factor(markers2vis[["cluster"]] , levels = sort(unique(seurat_ob@meta.data[[opt$collapseby]])))
        markers2vis[["cluster"]] = factor(markers2vis[["cluster"]] , levels = as.vector(sort(unique(seurat_ob@meta.data[[opt$collapseby]]))) )
        topn_markers  = markers2vis %>% group_by(cluster) %>%
            arrange(desc(avg_logFC) ) %>%
            top_n(opt$topn,.data[[opt$topby]]) %>% arrange(cluster )

        write.table( as.data.frame(topn_markers), file.path(output_dir,paste0("order_",basename(opt$genelist)) ),quote = F,row.names=F,sep="\t")
        ## mark <-markers2vis %>% group_by(cluster)  %>% arrange(desc(avg_logFC) )%>% top_n(2,.data[[opt$topby]]) %>% arrange(cluster )
        topn_markers = topn_markers %>% mutate(folder_suffix = paste0("cluster",cluster)) %>% select(cluster,gene,folder_suffix)
        markers2vis = as.vector(topn_markers$gene)
        
    }else{
        up = filter(genelist,FoldChange > 1) %>% arrange(desc(log2FoldChange ))  %>% top_n(topn,log2FoldChange ) %>% select(gene)
        down = filter(genelist,FoldChange < 1) %>% arrange(log2FoldChange )  %>% top_n(as.numeric(paste0("-",topn)),log2FoldChange ) %>% select(gene)
        genelist = rbind(up, down)
        genelist[,1]=factor(as.character(genelist[,1]))
        markers2vis = genelist[,1]
    }
}else{
    markers2vis = CaseMatch(search = as.vector(as.vector(genelist[,1])),match = rownames(seurat_ob))
    filtered_gene = genelist[! genelist[,1] %in% names(markers2vis ),1]
    if(length(filtered_gene)!=0){
        filtered_gene = as.data.frame(filtered_gene)
        colnames(filtered_gene) = "Gene"
        write.table(filtered_gene,file.path(output_dir,"filtered_gene.xls"),quote = F,row.names=F)
        print("There are some mismatched gene symbol, Please check filtered_gene.xls for the genename.")
    }
}
if (length(markers2vis) <= 1) stop("The number of matching gene is smaller than 2, please check the input gene list.")


# count = as.matrix(seurat_ob@data[markers2vis,])
count = as.matrix(GetAssayData(seurat_ob, slot = "data")[markers2vis,])
meta.data = seurat_ob@meta.data
collapseby = opt$collapseby

meta.data$id = rownames(meta.data)
collapsed_count = vector()
if ( !collapseby %in% colnames(meta.data) ){
    stop("NO specified column found!")
}

collapsed_group = meta.data %>% group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))
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

# data = tibble::rownames_to_column(as.data.frame(collapsed_count),"GeneID")
# write.table(data, file.path(output_dir,"heatmap_count.xls"),quote = F, row.names = F, sep = "\t")

## unique_gene <- unique(rownames(collapsed_count) ) ## unique 去重会改变原有的顺序

index=c(1,2,opt$topn*1+1,opt$topn*1+2,opt$topn*2+1,opt$topn*2+2,opt$topn*3+1,opt$topn*3+2,opt$topn*4+1,opt$topn*4+2)
anno_gene <- rownames(collapsed_count)[index]
## index_all <- which(rownames(collapsed_count) %in% as.vector(anno_gene) )

# index_all <- c()
# for (gene in unique(rownames(collapsed_count)) ){
#     s <- which(rownames(collapsed_count) %in% gene )
#     index_all <- c(index_all,s[1])
# }

# subset_index <- setdiff(index_all,index )
# if ( length(subset_index) ){
#     collapsed_count <- collapsed_count[-subset_index,]
# }else{
#     print("pass")
# }


data = tibble::rownames_to_column(as.data.frame(collapsed_count),"GeneID")
write.table(data, file.path(output_dir,"heatmap_data.xls"),quote = F, row.names = F, sep = "\t")


## c("#406AA8", "white", "#D91216")
palette <- colorRampPalette(c("#406AA8", "white", "#D91216"))(n=299)
# palette <- colorRampPalette(c("blue", "white", "red"))(n=299)
#col_fun = colorRamp2(c(-1.6,0.6,2.8), c("greenyellow","white", "red"))
#palette <-colorRampPalette(c("navy", "white", "firebrick3"))(299)
ind <- apply(collapsed_count, 1, mean) > 0
collapsed_count_filter <- collapsed_count[ind, ]


if (dim(collapsed_count)[2]>2) {
    df <- t(scale(t(collapsed_count_filter)))
}else{
    df <- collapsed_count_filter
}

df_scale = tibble::rownames_to_column(as.data.frame(df),"GeneID")
write.table(df_scale, file.path(output_dir,"heatmap_data_scale.xls"),quote = F, row.names = F, sep = "\t")

# mark <- rownames(df)[1:7]
# mark
#[1] "ERMP1" "SEC31A" "LOC101795171" "ATOX1" "CORO2B" "ZEB1"
#[7] "PGAM1"
#anno_mark()至少需要两个参数，其中at是原始数据矩阵的索引，标签是对应的文本。

# genelist <- c("ADAP2","ANKRD44","TBCD","ZEB1","DHX8","GPX7","ATOX1","PCF11")
## index <- which(rownames(df) %in% as.vector(anno_gene) )


lab = rowAnnotation(ano = anno_mark(at = index,
labels = anno_gene,
labels_gp = gpar(fontsize = 8)))

pdf(file.path(output_dir,"heatmap.pdf") )
Heatmap(df,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    right_annotation = lab,
    col = palette,
    name=" ")
dev.off()

## -verbose -density 500 -trim  topmarker_gene_heatmap.pdf  -quality 100  -flatten  topmarker_gene_heatmap.png
system(paste0("convert -verbose -density 500 -trim ",file.path(output_dir, "heatmap.pdf -quality 100  -flatten "),file.path(output_dir, "heatmap.png")))

### all gene mean data
count = as.matrix(GetAssayData(seurat_ob, slot = "data") )
meta.data = seurat_ob@meta.data
collapseby = opt$collapseby

meta.data$id = rownames(meta.data)
collapsed_count = vector()
if ( !collapseby %in% colnames(meta.data) ){
    stop("NO specified column found!")
}

collapsed_group = meta.data %>% group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))
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

data = tibble::rownames_to_column(as.data.frame(collapsed_count),"GeneID")
write.table(data, file.path(output_dir,"all_genes_celltype_meandata.xls"),quote = F, row.names = F, sep = "\t")

