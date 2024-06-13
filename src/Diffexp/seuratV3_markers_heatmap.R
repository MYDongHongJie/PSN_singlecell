#!/usr/bin/env Rscript
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages(library("dplyr"))

#=command line parameters setting=============================
option_list = list(
    make_option( c("--RDS", "-v"), type = "character", default = "TRUE",
            help = "the seurat object saved as R object in RDS format."),
    make_option( c("--markers","-l"), type ="character",
            help="the list file of marker genes to be visulized."),
    make_option( c("--output","-o"),type="character", default = "./",
            help="the output directory of results.", metavar="character"),
            help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user."),
    make_option( c("--topn", "-n"), type="integer", default = 1,
            help = "the number of top markers for each cluster to visualizse."),
    make_option( c("--topby", "-c"), type = "character", default = "avg_LogFC",
            help="the column used to pick top n marker gene to visulize.The
                 option can be one of the column in the input marker genes table."),
    make_option( c("--groupby", "-g"), type = "character", default = "clusters",
            help = "[OPTIONAL]The grouppinig variable in the metadata for
                            separate the cells to visulize marker genes."),
    make_option( c("--ident2use", "-q" ), type = "character", default = NULL,
            help = "[OPTIONAL]The column name in cell metadata used as identity
                            of each cell combined with which_cell."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
            help = "[OPTIONAL] subset of cluster ids used for subtyping.")
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
    seurat_ob = UpdateSeuratObject(seurat_ob)
}

#get the subset of cells used for visualization if necessay
if ( !is.null(opt$which_cells)){
    if ( is.null(opt$ident2use ) ){
        print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
        ident2use = "clusters"
    }else{
        ident2use = opt$ident2use
    }
    cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
    seurat_ob = SubsetData(seurat_ob,cells.use = FastWhichCells(seurat_ob,
    group.by = ident2use, subset.value = cluster_list ),subset.raw = T)
}

if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    root_dir = getwd()
}else{
    if ( file.exists(opt$output)){
        root_dir = opt$output
    }else{
        root_dir = opt$output
        dir.create(root_dir)
    }
}

if ( is.null( opt$groupby ) ){
    print( "NO groupping variable AVAILABLE for cell groupping! The default cell clusters id will be used!")
    groupby = "clusters"
}else{
    groupby = opt$groupby
}

if ( is.null( opt$markers ) & is.null(opt$extraGene)){
    stop("NO marker genes is AVAILABLE!")
}

topn_markers = data.frame()
if ( !is.null(opt$markers) ){
    markers2vis = read.table(  opt$markers, sep="\t", header = T)
    topn_markers  = markers2vis %>% group_by(cluster) %>%
        arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
        top_n(opt$topn,opt$topby) %>% arrange(cluster) %>% select(cluster,gene)
}
#=================================================================================
# visualize the markers in heatmap
#=================================================================================
seurat_ob = SetAllIdent( seurat_ob, id = groupby )
markers2vis4heatmap = unique(as.vector(topn_markers$gene))
seurat_ob = SetAllIdent( seurat_ob, id = groupby )
ggheat = DoHeatmap( object = SubsetData(object = seurat_ob,max.cells.per.ident = 100),
                    features=markers2vis4heatmap,size=4) +
                    theme(axis.text.y = element_text(size = 4))
ggsave(file.path(root_dir,paste0("top","marker_gene_heatmap.pdf", collapse = "_")))
ggsave(file.path(root_dir, paste0("top", "marker_gene_heatmap.png", collapse = "_")), dpi = 1000 ,limitsize = F)
