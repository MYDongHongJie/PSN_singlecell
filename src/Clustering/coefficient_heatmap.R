#!/usr/bin/env Rscript
# Created by: tangxuan
# Created on: 2020-06-11
# Updated on: -
# this script is used to export the gene expression in different group sepecified by the user.

rm(list=ls())

#=================================================================================
# customized function definition
#=================================================================================

suppressWarnings({
	suppressPackageStartupMessages( library("Seurat") )
	suppressPackageStartupMessages( library("optparse") )
	suppressPackageStartupMessages(library("dplyr"))
	suppressPackageStartupMessages(library("ggplot2"))
	suppressPackageStartupMessages(library("tibble"))
	suppressPackageStartupMessages(library("pheatmap"))
})

#=command line parameters setting=============================
option_list = list(
	make_option( c("--RDS", "-v"), type = "character", default = "TRUE",
		help = "the seurat object saved as R object in RDS format."),
	make_option( c("--output","-o"),type="character", default = "./",
		help="the output directory of results.", metavar="character"),
	make_option( c("--groupby", "-g"), type = "character", default = "clusters",
		help = "[OPTIONAL]The grouppinig variable in the metadata for
				separate the cells to calculate the correlation."),
	make_option( c("--sample_ratio"), type = "double", default = NULL,
		help = "[OPTIONAL]the ratio of random subsample for each group when calculating the correlation."),
	make_option( c("--slot"), type = "character", default = "data",
		help = "[OPTIONAL]the slot in the assay to use."),
	make_option( c("--assay", "-e"), type = "character", default = "RNA",
		help = "[OPTIONAL]the array result to use in case of multimodal analysis."),
	make_option( c("--var2use", "-q" ), type = "character", default = "clusters",
		help = "[OPTIONAL]The column name in cell metadata used as identity
				of each cell combined with levels4var."),
	make_option( c("--levels4var", "-u" ), type = "character", default = NULL,
		help = "[OPTIONAL] subset of factor levels for the specified factor by --var2use."),
    make_option( c("--regu_heatmap"), type= "logical",default = FALSE,
		help = "[OPTIONAL]visulize gene expression by regular heatmap."),
	make_option( c("--genelist"), type= "character",default = NULL,
		help = "[OPTIONAL]the gene list used to plot regular heatmap.")
	);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$RDS) ){
	stop("the seurat object is NOT AVAILABLE!")
}else{
	seurat_ob = readRDS(opt$RDS)
	if ( seurat_ob@version < 3){
		seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
	}

	# change the default assay for reduction if necessary
	if ( !is.null( opt$assay) ){
		DefaultAssay(seurat_ob) = opt$assay
	}else{
		DefaultAssay(seurat_ob) = "RNA"
	}
	metadata = seurat_ob@meta.data
	if ( is.null(metadata$clusters) ){
		seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
	}else{
		seurat_ob = SetIdent( seurat_ob, value = "clusters")
	}
	seurat_ob@meta.data$clusters = factor(seurat_ob@meta.data$clusters, levels = sort( unique(as.numeric( seurat_ob@meta.data$clusters))))
}
# output directory setting
if ( is.null(opt$output) ){
	print("NO output directory specified,the current directory will be used!")
	output_dir = getwd()
}else{
	if ( file.exists(opt$output)){
		output_dir = opt$output
	}else{
		output_dir = opt$output
		dir.create(output_dir)
	}
}

if ( is.null( opt$groupby ) ){
	print( "NO groupping variable AVAILABLE for cell groupping! The default cell clusters id will be used!")
	groupby = "clusters"
}else{
	groupby = opt$groupby
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
}


if ( !is.null(opt$sample_ratio) ){
	sampled_cellmeta = seurat_ob@meta.data %>% rownames_to_column() %>%
						group_by( .dots= groupby ) %>%
						sample_frac( size = opt$sample_ratio,replace = F) %>% column_to_rownames()
	seurat_ob = SubsetData(seurat_ob, cells = rownames(sampled_cellmeta))
		}

# export the gene expression matrix
groupby_data = vector()
for (i in names(table(seurat_ob[[groupby]]))  ) {
    sub_ob = SubsetData(seurat_ob, subset.name= groupby,accept.value=i)
    normalized_data = as.matrix(sub_ob[[opt$assay]]@data)
    meta.data = sub_ob@meta.data %>% tibble::rownames_to_column(var = "id")
    groupby_data = cbind(groupby_data,rowMeans(normalized_data))
}
colnames(groupby_data) = names(table(seurat_ob[[groupby]]))     
data = rownames_to_column(as.data.frame(groupby_data),var="GeneID")
write.table(data, file.path(output_dir,paste0("normalized_data_groupby_",groupby,".xls")),quote = F, row.names = F, sep = "\t")

#=================================================================================
# calculate the correlation and visulize by heatmap
#=================================================================================
colnames(groupby_data) = gsub('^',paste0(groupby,"_"),colnames(groupby_data))
matrix<-cor(groupby_data,method="pearson")
wid<-5+1.5*log2(length(colnames(data)))
hig<-5+1.5*log2(length(colnames(data)))
coefficient = pheatmap(matrix,
                       display_numbers = F,
                       border_color = "white",
                       scale = "none",
                       fontsize_number=(10.0+0.0001*log2(length(colnames(data)))),
                       number_format = "%.1f",
                       fontsize_row = (10.0+0.0001*log2(length(colnames(data)))),
                       fontsize_col = (10.0+0.0001*log2(length(colnames(data)))),
                       number_color="black",
                       angle_col=45)
ggsave(file.path(output_dir,"coefficient_heatmap.pdf"), plot = coefficient, height = hig, width = wid)
ggsave(file.path(output_dir,"coefficient_heatmap.png"), plot = coefficient, height = hig, width = wid, dpi = 1000)
matrix =rownames_to_column(as.data.frame(matrix),var="person")
write.table(matrix, file.path(output_dir,paste0("pearson_data_groupby_",groupby,".xls")),quote = F, row.names = F, sep = "\t")

#=================================================================================
# visulize gene expression by regular heatmap 
#=================================================================================
if (opt$regu_heatmap == T){
	if (!is.null(opt$genelist)){
		genelist = read.table(opt$genelist, header=T, sep = '\t')
	}else{
		stop("Please provide a genelist!")
	}
	
	filtered_data = groupby_data[as.vector(genelist[,1]),]
	ind <- apply(filtered_data, 1, mean) > 0
	heatmap_data <- filtered_data[ind, ]
	regu_heatmap = pheatmap(log2(heatmap_data+0.0001), cex=1, border=F,
							scale="row", fontsize_row=6, angle_col=45, treeheight_row=80, 
							treeheight_col=50, margin=c(18,8,8,5), show_rownames=T, col.names=T, cluster_rows=T, cluster_cols=T)
	ggsave(file.path(output_dir,"regular_heatmap.pdf"), plot = regu_heatmap, height = hig+2, width = wid)
	ggsave(file.path(output_dir,"regular_heatmap.png"), plot = regu_heatmap, height = hig+2, width = wid, dpi = 1000)
	gene_data = rownames_to_column(as.data.frame(filtered_data),var="GeneID")
	write.table(gene_data, file.path(output_dir,paste0("heatmap_data_groupby_",groupby,".xls")),quote = F, row.names = F, sep = "\t")
}
