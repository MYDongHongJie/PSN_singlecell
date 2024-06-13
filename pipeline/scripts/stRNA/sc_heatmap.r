suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

option_list = list(
    make_option(c("-i", "--input"), type="character",
          help="The input exprssion matrix in RDS format."),
    make_option(c("-g", "--genelist"),type="character",
           help="The genelist used to generate the heatmap with header."),
    make_option(c("-c", "--collapseby"),type="character",
           help="The variable level of heatmap. e.g. clusters or sampleid or celltype"),
    make_option(c("-r", "--rowcluster"),type="logical", default = "F",
           help="boolean values determining if rows should be clustered or hclust object,"),
    make_option(c("-l", "--colcluster"),type="logical", default = "F",
           help="boolean values determining if columns should be clustered or hclust object"),
    make_option(c("-n", "--showname"),type="logical", default = NULL,
           help="Whether to display the row name"),
    make_option(c("-w", "--gaps_row"),type = "character", default = NULL,
           help="Whether to display gaps by row"),
    make_option(c("-p", "--gaps_col"),type = "character", default = NULL,
           help="Whether to display gaps by col"),
    make_option( c("--outdir","-o"),type="character", default = "./",
           help="the output directory of Clustering results." ),
    make_option( c("--assay","-a"),type="character", default = "RNA",
          help="the assay to use in case of multimodal data." ),
    make_option( c("--topn"), type="integer", default = 25,
                 help = "[OPTIONAL] the number of top DEGs to visualizse."),
    make_option( c("--topby"), type = "character", default = "avg_logFC",
                 help="the column used to pick top n marker gene to visulize.The
                 option can be one of the column in the input marker genes table."),
    make_option( c("--order"), type = "character", default = NULL,
                help = "[OPTIONAL] specify the order of colnames to present from left to right.
                        OR can also be used to show subset of levels for factors specifyed by --collapedby"),
    make_option( c("--var2use", "-q" ), type = "character", default = NULL,
                help = "[OPTIONAL]The column name in cell metadata used as identity
                        of each cell combined with levels4var."),
    make_option( c("--levels4var", "-u" ), type = "character", default = NULL,
                help = "[OPTIONAL] subset of factor levels for the specified factor by --var2use."),
    make_option(c("-x", "--rowanno"), type="character", default = NULL,
                help="A table with row annotation."),
    make_option(c("-y", "--colanno"), type="character", default = NULL,
                help="A table with col annotation.")
    # make_option(c("-s", "--sign"), type="logical", default = "F",
                # help="Boolean values determine whether or not to add significance markers.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#=================================================================================
# parse the command line parameters
#=================================================================================
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

if ( is.null(opt$assay) ){
    assay = "RNA"
} else {
    assay = opt$assay
}

if ( is.null( opt$topn )){
    topn = 25
} else {
    topn = opt$topn
}

if ( !is.null( opt$gaps_row )){
    gaps_row = as.numeric(unlist(strsplit( opt$gaps_row,",",perl = T)))
} else {
    gaps_row = opt$gaps_row
}

if ( !is.null( opt$gaps_col )){
    gaps_col = as.numeric(unlist(strsplit( opt$gaps_col,",",perl = T)))
} else {
    gaps_col = opt$gaps_col
}

#####################################################################
#read in the single cell expression matrix and the genelist
#####################################################################
seurat_ob = readRDS(opt$input)

if ( seurat_ob@version < 3 ){
    seurat_ob = UpdateSeuratObject(seurat_ob)
}

#get the subset of cells used for visualization if necessay
if ( !is.null(opt$levels4var)){
    if ( is.null(opt$var2use ) ){
        print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
        ident2use = "clusters"
    }else{
        ident2use = opt$var2use
    }
    cluster_list = unlist(strsplit( opt$levels4var,",",perl = T))
    seurat_ob = SubsetData( seurat_ob, subset.name = ident2use, accept.value = cluster_list)
    seurat_ob@meta.data[[ident2use]]=factor(seurat_ob@meta.data[[ident2use]],levels = sort(unique(seurat_ob@meta.data[[ident2use]])))
}

DefaultAssay(seurat_ob) = assay

genelist = read.delim(opt$genelist, sep="\t", header = T) #读取基因并进行筛选
if (dim(genelist)[2]>1){
    if ( "cluster" %in% colnames(genelist) ) {
        markers2vis = genelist
        markers2vis[["cluster"]] = factor(markers2vis[["cluster"]] , levels = sort(unique(seurat_ob@meta.data[[opt$collapseby]])))
        topn_markers  = markers2vis %>% group_by(cluster) %>%
            arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
            top_n(opt$topn,.data[[opt$topby]]) %>% arrange(cluster)
        write.table( as.data.frame(topn_markers), file.path(output_dir,paste0("order_",basename(opt$genelist)) ),quote = F,row.names=F,sep="\t")
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
count = as.matrix(GetAssayData(seurat_ob, slot = "data")[markers2vis,]) #基因表达数据
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
write.table(data, file.path(output_dir,"heatmap_count.xls"),quote = F, row.names = F, sep = "\t")
if (dim(collapsed_count)[2]>2) {
    data_scaled = tibble::rownames_to_column(as.data.frame(pheatmap:::scale_rows(collapsed_count)),"GeneID")
    write.table(data_scaled, file.path(output_dir,"heatmap_count_scaled.xls"),quote = F, row.names = F, sep = "\t")
}
#=================================================================================
# data transforamtion & heatmap plotting 
#=================================================================================
ind <- apply(collapsed_count, 1, mean) > 0
collapsed_count_filter <- collapsed_count[ind, ]

palette <- colorRampPalette(c("#406AA8", "white", "#D91216"))(n=299)
# palette <-colorRampPalette(c("navy", "white", "firebrick3"))(299)

if ( is.null(opt$rowanno) ){
    annotation_row = NA
} else {
    annotation_row = read.delim(opt$rowanno, row.names=1)
}

if( is.null(opt$colanno) ){
    annotation_col = NA
}else{
    annotation_col = read.delim(opt$colanno, row.names=1)
}

# if ( opt$sign && opt$rowanno == "up_down" ){
  # pmt <- exps$pvalue  #提取p值
  # #判断显著性
  # if (!is.null(pmt)){
    # ssmt <- pmt< 0.01
    # pmt[ssmt] <-'**'
    # smt <- pmt >0.01& pmt <0.05
    # pmt[smt] <- '*'
    # pmt[!ssmt&!smt]<- ''
  # } else {
    # pmt <- F
  # }
  # #构建显著性注释矩阵
  # display_number <- matrix(0,nrow(collapsed_count_filter),ncol(collapsed_count_filter))
  # for(i in 1:dim(collapsed_count_filter)[1]){
    # display_number[i,] <- pmt[i]
  # }
  
# }else{
  # display_number <- opt$sign
# }

## parsing show rowname
if (is.null(opt$showname))  showname=ifelse(dim(collapsed_count_filter)[1]>100,FALSE,TRUE) else showname = opt$showname
cellwidth=36
cellheight=ifelse(showname,12,8) # if showname =F, set the panel height to 576/72point = 8 inches
p = pheatmap(log2(collapsed_count_filter+0.0001),
        #color=palette,
        cex=1,
        border="gray",
        angle_col=45,
        treeheight_row=36, treeheight_col=36,
        lwd=1.1,
        cellheight=cellheight, cellwidth=cellwidth,
        scale=ifelse(dim(collapsed_count_filter)[2]==2,"none","row"),
        show_rownames=showname,
        gaps_row=gaps_row,
        gaps_col=gaps_col,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        # display_numbers = display_number,
        cluster_rows=opt$rowcluster ,cluster_cols=opt$colcluster)
ggsave(file.path(output_dir,"heatmap.pdf"), plot= p, 
    width=(36*dim(collapsed_count_filter)[2]+144)/72+2, 
    height = ifelse(showname,(12*dim(collapsed_count_filter)[1]+108)/72,(8*dim(collapsed_count_filter)[1]+108)/72))
ggsave(file.path(output_dir,"heatmap.png"), plot= p, dpi=1000,
    width=(36*dim(collapsed_count_filter)[2]+144)/72+2,
    height = ifelse(showname,(12*dim(collapsed_count_filter)[1]+108)/72,(8*dim(collapsed_count_filter)[1]+108)/72))

