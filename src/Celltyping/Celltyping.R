#!/usr/bin/env Rscript
# This is the script used to identify the cell types
# in the sequenced data after data clustering and predicting gene activity with Cicero
# the package used here is developed by ...
# the cell typing method used here is based on the reference data, which
# compare all the cells in the sample with all the cell expression in
# the referenced expression matrix.
# The reference matrix usually come from the previous microarray results,
# RNA-seq and single cell
# RNA-seq data from GEO or any other place. The reference data must
# have pure cell types annotation for each sample which will be used
# to make a singler object. So this method used here is limited by
# the cell type annotation data in databases.
# to run this program, the following things are necessary:
# 1. the seurat object/filtered gene activity matrix
# 2. the reference expression matrix with sce format

.libPaths(c("/home/luyao/R/x86_64-conda_cos6-linux-gnu-library/3.6/"))

rm(list=ls())
suppressWarnings({
    suppressPackageStartupMessages(library(SingleR))
    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(Seurat))
    suppressPackageStartupMessages(library(Matrix))
    suppressPackageStartupMessages(library(scran))
    suppressPackageStartupMessages(library(scater))
    suppressPackageStartupMessages(library(pheatmap))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(RColorBrewer))
    suppressPackageStartupMessages(library(future))
    suppressPackageStartupMessages(library(OESingleCell))
})

####################################################################
#function definition
####################################################################
CustomCol =   function (n) {
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual", ]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_palette <- col_vector
    return(col_palette[n]) 
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

# characterizing the cell types based on the reference data
celltyping <- function(test.sce,ref.sce,nmakers,ncpu,types){
    # screen for overlapping genes
    common <- intersect(rownames(test.sce), rownames(ref.sce))
    test.sce <- test.sce[common,]
    ref.sce <- ref.sce[common,]

    print("[TOPN mode] Find markers")
    out <- pairwiseTTests(logcounts(ref.sce), ref.sce[[types]], direction="up")
    markers <- getTopMarkers(out$statistics, out$pairs, n=nmakers)
    trained <- trainSingleR(ref.sce, factor(ref.sce[[types]]),genes=markers)
    print("[TOPN mode] Identifying cell types")
    pred <- classifySingleR(test.sce, trained,BPPARAM = MulticoreParam(workers = ncpu))
    return(pred)
}

# plotScoreHeatmap
plotScoreHeatmap = function (results, cells.use = NULL, labels.use = NULL, clusters = NULL, 
    show.labels = FALSE, show.pruned = FALSE, max.labels = 40, 
    normalize = TRUE, cells.order = NULL, order.by.clusters = FALSE, 
    annotation_col = NULL, annotation_colors = "CustomCol2" ,show_colnames = FALSE, group_by="clusters",...) 
{
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    if (is.null(annotation_col)) {
        annotation_col <- data.frame(row.names = rownames(results))
    }
    if (show.pruned) {
        prune.calls <- results$pruned.labels
        if (!is.null(prune.calls)) {
            names(prune.calls) <- rownames(results)
            Pruned <- data.frame(Pruned = as.character(is.na(prune.calls)), 
                row.names = rownames(results))[rownames(annotation_col), 
                ]
            annotation_col <- cbind(Pruned, annotation_col)
        }
    }
    if (show.labels) {
        labels <- results$labels
        names(labels) <- rownames(results)
        annotation_col$Labels <- labels[rownames(annotation_col)]
    }
    if (!is.null(clusters)) {
        names(clusters) <- rownames(results)
        annotation_col$Clusters <- clusters[rownames(annotation_col)]
    }
    scores <- results$scores
    rownames(scores) <- rownames(results)
    if (!is.null(cells.use)) {
        scores <- scores[cells.use, , drop = FALSE]
        clusters <- clusters[cells.use]
        cells.order <- cells.order[cells.use]
    }
    if (!is.null(labels.use)) {
        scores <- scores[, labels.use, drop = FALSE]
    }
    cluster_cols <- FALSE
    if (order.by.clusters) {
        order <- order(clusters)
    } else if (!is.null(cells.order)) {
        order <- order(cells.order)
    } else {
        order <- seq_len(nrow(scores))
        cluster_cols <- TRUE
    }
    m <- rowMaxs(scale(t(scores)))
    to.keep <- head(order(m, decreasing = TRUE), max.labels)
    if (normalize) {
        mmax <- rowMaxs(scores)
        mmin <- rowMins(scores)
        scores <- (scores - mmin)/(mmax - mmin)
        scores <- scores^3
        breaks <- seq(0, 1, length.out = 101)
    } else {
        breaks <- seq(-1, 1, length.out = 101)
    }
    scores <- scores[, seq_len(ncol(scores)) %in% to.keep, drop = FALSE]
    scores <- t(scores)
    args <- list(mat = scores[, order, drop = FALSE], border_color = NA, 
        show_colnames = show_colnames, clustering_method = "ward.D2", 
        cluster_cols = cluster_cols, breaks = breaks,...)
    if (normalize) {
        # args$color <- (grDevices::colorRampPalette(c("blue","red")))(100)
        # args$color <- (grDevices::colorRampPalette(c("royalblue","yellow","red")))(100)
        # args$color  = wes_palette("Zissou1", 100, type = "continuous")
        args$legend <- TRUE
    }

    if (ncol(annotation_col) > 0) {
        args$annotation_col <- annotation_col
    }
    if (annotation_colors == "CustomCol2") {
        args$annotation_colors <-  list(Clusters=CustomCol2(as.numeric(levels(clusters))))
        names(args$annotation_colors$Clusters) <- as.numeric(levels(clusters))
    } else if (annotation_colors == "CustomCol") {
        args$annotation_colors <-  list(Clusters=CustomCol(as.numeric(levels(clusters))))
        names(args$annotation_colors$Clusters) <- as.numeric(levels(clusters))
    }  else args$annotation_colors <- SingleR:::.make_heatmap_annotation_colors(args, show.pruned)
    do.call(pheatmap::pheatmap, args)
}

####################################################################
#command line parameters setting
####################################################################
option_list = list(
  make_option( c( "--INEXPRESS","-i"), type = "character",
               help = "the input expression matrix in seurat format to be analyzed (rds file)." ),
  make_option( c( "--ref_dataset","-r"), type = "character",
               help = "the reference data set in SingleCellExperiment format (rds file)."),
  make_option( c( "--REFDATA"), type = "character",
               help = "[OPTIONAL] the file path to the non-built-in reference dataset in 1.Seurat format (.rds file) 2. SingleCellExperiment (.rds file) or
                       3.matrix (.txt file) Note: the gene name should be in the first column of hte matrix file
                       ; the number of column in the matrix file should equal to the number of nrow in the metadata file"),
  make_option( c( "--REFMETA"), type = "character",
               help = "[OPTIONAL] the file path to the non-built-in reference metadata .txt file (with barcode in the 1st, label in the 2nd column)."),
  make_option( c("--singler" ), type = "character",default = NULL,
               help = "[OPTIONAL] the previous celltyping results stored in singler object for reuse.[default:NULL]" ),
  make_option( c("--genes", "-g" ), type = "character",default = "de",
               help = "SingleR calculate default parameters for marker genes.[default:de]" ),
  make_option( c("--mode","-m" ), type = "character",default = "separate",
               help = "celltyping mode when mulitple references are provided. can be either combined or separate .[default:separate]" ),
  make_option( c("--nmakers", "-n" ), type = "integer",default = 25,
               help = "[OPTIONAL]the number of markers each cell. [default:25]" ),
  make_option( c("--assay"), type = "character", default = "RNA",
              help = "[OPTIONAL]the default assay to use  for the this run.[default:RNA]"),
  make_option( c("--CORES", "-j" ), type = "integer", default = 10,
               help = "[OPTIONAL]the core number used to run this script.[default:10]" ),
  make_option( c("--pointsize", "-v"), type = "double", default = NULL,
               help = "[OPTIONAL]the point size in the plot."),
  make_option( c("--color_scheme", "-c"), type = "character", default = "CustomCol2",
               help = "[OPTIONAL]the color scheme used in the plot. Could one of:
                       CustomCol (report v2), CustomCol2 (v3), [default:CustomCol2]"),
  make_option( c("--reduct"), type = "character", default = "tsne",
               help = "[OPTIONAL]the reduction coordination to visualize the  results.[default:tsne]"),
  make_option( c("--TOPN", "-N" ), type = "integer",default = 25,
               help = "the top number of cell typing results to display on heatmap.[default:25]" ),
  make_option( c("--SPECIES", "-s" ), type = "character",
               help = "the specie where the samples come from.Currently,
                        Human and Mouse are supported as default."),
  make_option( c("--LEVEL", "-l" ), type = "character",default = "main",
               help = "the cell typing level for this run.
                       1. For non-built-in dataset this is the column name of the celltype labels.
                       2. For built-in refdata there are two levels supported for built-in refdata: main and single.
                       'single' means the subtype level and 'main' means the main cell type level.
                       If 'single'is specified, please make sure that your reference expression
                       matrix has the subtype information in 'pdata' file. [default:main]"),
  make_option( c("--resolution"), type = "double",default = NULL,
                  help = "Clustering chart resolution.[default:NULL]"),
  make_option( c("--OUTDIR", "-o" ), type = "character",
               help = "the output directory of this program." ),
  make_option( c("--WHICH_GROUP", "-q" ), type = "character",
                help = "[Optional]select the groupping column in metadata used to do subsettig of cell if necessary."),
  make_option( c("--WHICH_CELLS", "-u" ), type = "character",
                help = "[Optional] the level id list in selected cell groupping used to celltyping.For some cell 
                       clusters with high hertergensity this may be useful for sub-celltyping
                       combined with the option -l/--LEVEL with single mode.
                       If not specified with cell clusters's ID, all cells will be used.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

####################################################################
# parse the command line parameters
#####################################################################
assay2use = opt$assay

if ( !is.null(opt$CORES) ){
    nCPUs = opt$CORES
} else {
    nCPUs = 10
}

if ( !is.null(opt$nmakers) ){
    nmakers = opt$nmakers
} else {
    nmakers = 25
}

if ( is.null( opt$reduct ) ){
    reduct.use = "tsne"
} else {
    reduct.use = opt$reduct
}

if ( is.null( opt$LEVEL ) ){
    LEVEL = "main"
} else {
    LEVEL = opt$LEVEL
}

if ( is.null(opt$color_scheme) ){
    colorscheme = "CustomCol2"
} else {
    colorscheme = opt$color_scheme
}

if ( is.null(opt$OUTDIR) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
} else {
    if ( file.exists(opt$OUTDIR) ){
        output_dir = opt$OUTDIR
    } else {
        output_dir = opt$OUTDIR
        dir.create(output_dir, recursive = T)
    }
}
output_dir = normalizePath(output_dir)

if ( is.null(opt$SPECIES) ){
    stop("NO specie specified!")
} else {
    species = opt$SPECIES
}

#####################################################################
# prepare the refrence datasets
#####################################################################
ref.sce = list()
if ( !is.null(opt$ref_dataset) ){
    builtindatasets = unlist(strsplit(opt$ref_dataset, ",", perl = T))
    for (refdata in builtindatasets) {
        refrds = paste0(refdata,".rds")
        ref.sce[[refdata]] = readRDSMC(file.path("/data/database/celltype_refdata/logNorm_rds", refrds), cores = availableCores() )
    }
}
if (!is.null(opt$REFDATA)) {
    if ( !is.null(opt$REFMETA) ){ ## 1.refdataset in [[matrix table + meta table]] format
        refdata = gsub("\\..[^\\.]*$","",make.names(basename(opt$REFDATA)))
        ref.metadata = read.table(opt$REFMETA,header = TRUE,row.names = 1,sep = "\t")
        ref.matrix = read.table(opt$REFDATA,header = TRUE,row.names = 1,sep = "\t")
        ref.sce[[refdata]] = SingleCellExperiment(assays = list(counts = ref.matrix), colData = ref.metadata)
        ref.sce[[refdata]] <- ref.sce[[refdata]][, colSums(counts(ref.sce[[refdata]])) > 0]
    } else {
        refdata = gsub("\\..[^\\.]*$","",make.names(basename(opt$REFDATA)))
        REF_OBJ = readRDSMC(opt$REFDATA, cores = availableCores() )
        # if the input seurat object version is less than 3, upgrade it to version 3
        if ( REF_OBJ@version < 3 ){
            REF_OBJ = UpdateSeuratObject(REF_OBJ) #make sure the seurat object match with the latest seurat package
        }
        if (class(REF_OBJ) == "Seurat") { ## 2.refdataset in seurat [[rds]] format
            ref.metadata = REF_OBJ@meta.data %>% select(!!LEVEL)
            ref.matrix = GetAssayData(REF_OBJ,assay = "RNA",slot = "counts")
            ref.sce[[refdata]] = SingleCellExperiment(assays = list(counts = ref.matrix), colData = ref.metadata)
        } else if (class(REF_OBJ) == "SingleCellExperiment") { ## 3.refdataset in sce [[rds]] format
            ref.sce[[refdata]] = REF_OBJ
        }
    }
    ref.sce[[refdata]] = logNormCounts(ref.sce[[refdata]])
    # saveRDS(ref.sce[[refdata]],paste0("sce_logNorm_",refdata,".rds")) ## if you want to save the reference for another time
}
if (length(ref.sce) == 0) stop("NO reference dataset specified!")

#####################################################################
#read in the single cell expression matrix in seruat format
#####################################################################
if ( !is.null(opt$INEXPRESS) ){
    # the input is a seurat object which may contain more than one sample
    seurat_ob = readRDSMC( opt$INEXPRESS, cores = availableCores()  )
    # if the input seurat object version is less than 3, upgrade it to version 3
    if ( seurat_ob@version < 3 ){
        seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
    }
    DefaultAssay(seurat_ob) = opt$assay
    # parse the resolution. Note: it only serve as part of the ouput file name.
    if ( is.null(opt$resolution) ){
        findcluster_record = Command(seurat_ob, command = "FindClusters")
        resolution = findcluster_record$resolution
        if ( is.null(resolution)) stop("please specify the resolution with --resolution") 
    } else {
        resolution = opt$resolution
    }
    if ( is.null(opt$pointsize) ){
        if (dim(seurat_ob)[2] < 500){
            pointsize = 1.5
        } else pointsize = 0.5
    } else {
        pointsize = opt$pointsize
    }
}

#####################################################################
# identifying cell types
#####################################################################
if ( is.null(opt$singler) ){
    ## subset. Curretnly not supporting subset singler object. 
    if ( !is.null(opt$WHICH_CELLS)){
        cluster_list = unlist(strsplit( opt$WHICH_CELLS,",",perl = T))
        seurat_ob = SubsetData(seurat_ob,  subset.name = opt$WHICH_GROUP, accept.value = cluster_list)
    }
    # seurat to sce
    test.m = GetAssayData(seurat_ob,assay = assay2use,slot = "counts")
    test.sce = SingleCellExperiment(assays = list(counts = test.m))
    test.sce = logNormCounts(test.sce)

    pred=list()
    if (length(ref.sce) > 1 & opt$mode == "combined" ) {
         print("using combined result from multiple datasets")
         for (refdata in names(ref.sce)){
             ref.sce[[refdata]][[LEVEL]] = paste(refdata,ref.sce[[refdata]][[LEVEL]],sep = "_")
         }
         refdata.combine = paste0( names(ref.sce),collapse="_")
         pred[[refdata.combine]] = SingleR(test.sce,ref.sce,labels = lapply(ref.sce, FUN = function(x) factor(x[[LEVEL]])) ,BPPARAM = MulticoreParam(workers = nCPUs))
         pred[[refdata.combine]]$scores = do.call(cbind,lapply(pred[[refdata.combine]]$orig.results,FUN = function(x) x$scores))
    } else {
        for (refdata in names(ref.sce)) {
            if (opt$genes == "de"){
                pred[[refdata]] = SingleR(test.sce,ref.sce[[refdata]],labels = factor(ref.sce[[refdata]][[LEVEL]]),BPPARAM = MulticoreParam(workers = nCPUs))
            } else {
                pred[[refdata]] = celltyping(test.sce,ref.sce[[refdata]],nmakers,nCPUs,LEVEL)
            }
        }
    }
    if (!is.null(opt$REFDATA)) {
        print("1.Saving singler object")
        saveRDSMC(pred,file.path(output_dir,"singler.rds"))
    }
} else pred = readRDSMC(opt$singler, cores = availableCores() )

#####################################################################
# add the celltyping result as the metadata & plotting
#####################################################################
for (refdata in names(pred)) {
    # add the celltyping result as the metadata
    raw.metaname<-paste0("raw_",LEVEL,"_celltype")
    seurat_ob <- AddMetaData(object = seurat_ob,
                             metadata = pred[[refdata]]$labels, col.name = raw.metaname) 
    seurat_ob@meta.data$clusters = factor(seurat_ob@meta.data$clusters,levels = sort(as.numeric(unique(seurat_ob@meta.data$clusters))))
    #count the cell number for each cell type
    celltyping_stat = seurat_ob@meta.data %>% select( clusters, !!raw.metaname)%>%
                           group_by(clusters)%>% dplyr::count(.data[[raw.metaname]])  %>% dplyr::rename(cell_num=n)
    write.table(celltyping_stat, quote = F,
                file.path(output_dir,paste0(species,"ref_",refdata,"_",LEVEL,"_celltyping_statistics.xls",collapse="")),
                sep = "\t",row.names = F)
    ## heatmap 
    print("2.Plot celltyping heatmap")
    if (length(Cells(seurat_ob)) <= 40000){
        ggheat = plotScoreHeatmap(pred[[refdata]], clusters = seurat_ob@meta.data$clusters, max.labels = opt$TOPN, 
                                  annotation_colors = colorscheme, order.by.clusters = T)
    } else {
        print("Screening of 50,000 cells in proportion...")
        meta.data = seurat_ob@meta.data
        meta.data$cells = Cells(seurat_ob)
        meta.data = meta.data %>% group_by(clusters) %>% sample_frac(40000/length(Cells(seurat_ob)))
        ggheat = plotScoreHeatmap(pred[[refdata]], cells.use = meta.data$cells, clusters = seurat_ob@meta.data$clusters, max.labels = opt$TOPN, 
                                  annotation_colors = colorscheme, order.by.clusters = T)
    }
    ggsave(file.path(output_dir,paste0(species,"ref_",refdata,"_",LEVEL,"_celltyping_heatmap.pdf",collapse=".")),plot = ggheat,
           width = 6.5+max(nchar(names(head(sort(table(as.factor(pred[[refdata]]$labels)),decreasing = T),opt$TOPN))))/15)
    ggsave(file.path(output_dir,paste0(species,"ref_",refdata,"_",LEVEL,"_celltyping_heatmap.png",collapse = ".")),plot = ggheat,
           width = 6.5+max(nchar(names(head(sort(table(as.factor(pred[[refdata]]$labels)),decreasing = T),opt$TOPN))))/15,bg="white")
    ## plot raw celltype
    print("3.Plot celltyping tsne")
    seurat_ob = SetIdent( seurat_ob, value = raw.metaname)
    nlevel = length(unique(seurat_ob@meta.data[,raw.metaname]))
    if (nlevel <= 30){
        if (colorscheme == "CustomCol2" ){
            ggtsne = DimPlot(object = seurat_ob, reduction = reduct.use,pt.size = pointsize ) +
              theme( plot.title = element_text(hjust = 0.5),legend.key.size = unit(0.9,"lines"),legend.text = element_text(size = 10)) +
              scale_colour_manual( values = CustomCol2(1:nlevel)) + 
              guides(colour = guide_legend(ncol = 1,override.aes = list(size=2)))
        } else if (colorscheme == "CustomCol"){
            ggtsne = DimPlot(object = seurat_ob, reduction = reduct.use,pt.size = pointsize ) +
              theme( plot.title = element_text(hjust = 0.5),legend.key.size = unit(0.9,"lines"),legend.text = element_text(size = 10)) +
              scale_colour_manual( values = CustomCol(1:nlevel)) +
              guides(colour = guide_legend(ncol = 1,override.aes = list(size=2)))
        }
        ggsave(file.path(output_dir,paste0(species,"ref_",refdata,"_",LEVEL,"_celltyping_on_",reduct.use,"_resolution",resolution,"_plot.pdf")),
               width = (6.5+ max(nchar(levels(factor(as.character(seurat_ob@meta.data[[raw.metaname]])))))/13))
        ggsave(file.path(output_dir,paste0(species,"ref_",refdata,"_",LEVEL,"_celltyping_on_",reduct.use,"_resolution",resolution,"_plot.png")),dpi=1000,
               width = (6.5+ max(nchar(levels(factor(as.character(seurat_ob@meta.data[[raw.metaname]])))))/13),bg="white")
    } else {
        suppressPackageStartupMessages(library(RColorBrewer))
        n = as.integer(nlevel/50)+1
        ggtsne = DimPlot(object = seurat_ob, reduction = reduct.use,pt.size = pointsize ) +  
          theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(0.65,"lines"),legend.text = element_text(size = 8)) +  
          scale_colour_manual(values = colorRampPalette(brewer.pal(8,"Set1"))(nlevel) )+  
          guides(colour = guide_legend(ncol = n,override.aes = list(size = 1.5)))
        ggsave(file.path(output_dir,paste0(species,"ref_",refdata,"_",LEVEL,"_celltyping_on_",reduct.use,"_resolution",resolution,"_plot.pdf",collapse="")),width = 7+2*n)
        ggsave(file.path(output_dir,paste0(species,"ref_",refdata,"_",LEVEL,"_celltyping_on_",reduct.use,"_resolution",resolution,"_plot.png")),dpi=1000,width = 7+2*n,bg="white")
    } 
    #####################################################################
    ## TOPN result
    #####################################################################
    seurat_ob = SetIdent( seurat_ob, value = "clusters")
    top_celltype = celltyping_stat %>% group_by(clusters) %>% top_n(1,cell_num)
    from.id = as.vector(top_celltype$clusters)
    to.id = as.vector(top_celltype[[raw.metaname]])
    seurat_ob = SetIdent(seurat_ob, value = plyr::mapvalues(x = Idents(seurat_ob),from = from.id,to = to.id))
    seurat_ob[["celltype"]] <- Idents(object = seurat_ob)
    seurat_ob[[paste0(refdata,".",LEVEL,".celltype")]] <- Idents(object = seurat_ob)
    if ('rawbc' %in% colnames(seurat_ob@meta.data)) {
        Barcode_content = 'rawbc'
    }else{
        Barcode_content = 'orig.ident'
    }   
    full_celltyping_df = seurat_ob@meta.data %>%
                         tibble::rownames_to_column(var = "barcode_inuse") %>%
                         dplyr::rename( "Barcode" = Barcode_content ) %>% select( Barcode, everything())
    write.table(full_celltyping_df, quote = F,
                file.path(output_dir,paste0(species,"ref_",refdata,"_resolution",resolution,
                          "_",LEVEL,"_celltyping_results.xls",collapse = "")), sep = "\t", row.names = F)
    simplified_celltyping_df = seurat_ob@meta.data %>%  
                               dplyr::rename( "Barcode" = Barcode_content ) %>%  
                               select( Barcode, sampleid, celltype, clusters,group)
    write.table(simplified_celltyping_df, quote = F,sep =",",row.names = F,
                file.path(output_dir,paste0(species,"ref_",refdata,"_resolution",resolution, 
                          "_",LEVEL,"_simplified_celltyping_results.csv",collapse = "")))
    
    print("4.plot celltyping tsne(top)")
    ## obtain the top cluster number for each celltype (for the purpose of coloring)
    top_celltype$clusters = as.integer(top_celltype$clusters)
    top_cluster_for_celltype = top_celltype %>% group_by(.data[[raw.metaname]]) %>%
                               top_n(-1,clusters) %>% select(-cell_num)
    top_cluster_for_celltype  = as.data.frame(top_cluster_for_celltype)
    rownames(top_cluster_for_celltype) = top_cluster_for_celltype[[raw.metaname]]
    color_level = as.numeric(top_cluster_for_celltype[levels(seurat_ob@meta.data$celltype),"clusters"])
    if (colorscheme == "CustomCol2" ){
        ggtsne = DimPlot(object = seurat_ob, reduction = reduct.use,pt.size = pointsize ) +
          theme( plot.title = element_text(hjust = 0.5)) + scale_colour_manual( values = CustomCol2(color_level))
    } else if (colorscheme == "CustomCol"){
        ggtsne = DimPlot(object = seurat_ob, reduction = reduct.use,pt.size = pointsize ) +
          theme( plot.title = element_text(hjust = 0.5)) + scale_colour_manual( values = CustomCol(color_level))
    }
    ggsave(file.path(output_dir,paste0(species,"ref_",refdata,"_top.",LEVEL,"_celltyping_on_",reduct.use,"_resolution",resolution,"_plot.pdf",collapse="")),plot = ggtsne,
           width = (6.5 + max(nchar(levels(factor(as.character(seurat_ob@meta.data$celltype)))))/10))
    ggsave(file.path(output_dir,paste0(species,"ref_",refdata,"_top.",LEVEL,"_celltyping_on_",reduct.use,"_resolution",resolution,"_plot.png")),dpi=1000, plot = ggtsne,
           width = (6.5 + max(nchar(levels(factor(as.character(seurat_ob@meta.data$celltype)))))/10),bg="white")
    seurat_ob = SetIdent(seurat_ob, value = "clusters")

}

#####################################################################
## save seurat object
#####################################################################
if ( !is.null(opt$WHICH_CELLS)){
  saveRDSMC(seurat_ob, file.path(output_dir, paste0("subseted_", opt$WHICH_GROUP, "_", paste0(cluster_list, collapse = "-"), "_seurat.rds")))
} else{
  saveRDSMC(seurat_ob, opt$INEXPRESS)
}

print("all done!")
