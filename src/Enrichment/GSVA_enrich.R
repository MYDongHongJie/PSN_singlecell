#!/usr/bin/env Rscript
#this is the script used to realize the Function Class type gene set enrichment
#analysis of results from the single cell RNA-seq sequencing expriements

#=================================================================================
#function definition
#=================================================================================

#count the quantification matrix of each gene for each sample/cell
#meta.data the annotation of each sample/cell
#collapseby the vector of header names of the meta.data data.frame
collapse_count <- function(count, meta.data, collapseby,use ){

    suppressPackageStartupMessages(library(dplyr))
    meta.data$id = rownames(meta.data)
    collapsed_count = vector()
    if ( !collapseby %in% colnames(meta.data) ){
        stop("NO specified column found!")
    }

    collapsed_group = meta.data %>% group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))

    for ( cells in collapsed_group$cellid ){
        samplex = unlist(strsplit(cells, ",", perl =T))
        if ( use == "sum" ){
            collapsed_count= cbind(collapsed_count,rowSums( count[,samplex] ))
        }else if ( use == "mean"){
            collapsed_count= cbind(collapsed_count,rowMeans( count[,samplex] ))
        }else{
            stop("Neither mean or sum is supplied.")
        }
    }
    collapsed_count = as.matrix( collapsed_count )
    collapsed_group = as.data.frame(collapsed_group)
    colnames(collapsed_count) = collapsed_group$clusters
    return(collapsed_count)
}


#=================================================================================
#main workflow
#=================================================================================
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(methods))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

#command line interface definition using docopt package
option_list = list(
make_option(c("--input", "-i"), type = "character",
            help = "the gene expression matrix. It may come in different
                    format specified by the informat parameter. The normoalized
                    data will be used if it comes to seurat object."),
make_option( c( "--informat", "-f"), type = "character",default = "seurat",
            help = "[Required]the format of the input file.Options can be:seurat,raw."),
make_option( c( "--metadata", "-m"), type = "character",
            help = "[Optional]the metadata of each sample/cell in the input file."),
make_option( c( "--method", "-M"), type = "character",default= "gsva",
            help = "Method to employ in the estimation of gene-set
                    enrichment scores per sample. By default this is set
                    to gsva (Hanzelmann et al, 2013) and other options
                    are ssgsea (Barbie et al, 2009), zscore (Lee et
                    al, 2008) or plage (Tomfohr et al, 2005)."),
make_option( c("--gmt", "-g" ), type = "character",
            help = "the gene sets in gmt format from the MsigDB/KEGG/GO database etc."),
# make_option( c("--group", "-l"), type = "character", default = "cluster",
#             help = "[OPTIONAL]the groupping field of cells used to collapse the single cell count matrix."),
make_option( c( "--splitby", "-n" ), type = "character",
        help = "[OPTIONAL]Parallelized processing by spliting the cell size using one groupping factors in the metadata
                        in case of too many cell for calculation"),
make_option( c( "--chunkby", "-c" ), type = "integer",
        help = "[OPTIONAL]Parallelized processing by even chunked cell numbers for improving performace
                        in case of too many cell for calculation"),
make_option( c("--OUTDIR", "-o" ), type = "character",
            help = "the output directory of the analysis results."),
make_option( c("--kcdf", "-k"), type = "character",default= "Gaussian",
            help = "Character string denoting the kernel to use during the
                    non-parametric estimation of the cumulative distribution
                    function of expression levels across samples when method='gsva'.
                    The option can be Guassian, Possion and none.
                    By default, kcdf='Gaussian' which is suitable when input
                    expression values are continuous, such as microarray
                    fluorescent units in logarithmic scale, RNA-seq log-
                    CPMs, log-RPKMs or log-TPMs. When input expression
                    values are integer counts, such as those derived from
                    RNA-seq experiments, then this argument should be set
                    to kcdf='Poisson'. This argument supersedes arguments
                    rnaseq and kernel, which are deprecated and will be
                    removed in the next release."),
make_option( c("--abs_rank", "-a"), type = "logical",default= F,
            help = "Flag used only when mx_diff=TRUE. When
                    abs_ranking=FALSE (default) a modified Kuiper
                    statistic is used to calculate enrichment scores,
                    taking the magnitude difference between the largest
                    positive and negative random walk deviations. When
                    abs.ranking=TRUE the original Kuiper statistic that
                    sums the largest positive and negative random walk
                    deviations, is used. In this latter case, gene sets
                    with genes enriched on either extreme (high or low)
                    will be regarded as'highly' activated."),
make_option( c("--min_sz", "-s"), type = "integer",  default= 2,
            help = "Minimum size of the resulting gene sets."),
make_option( c("--max_sz", "-S"), type = "integer", default= 99999,
            help = "Maximum size of the resulting gene sets."),
make_option( c("--parallel_sz", "-j"), type = "integer", default = 8,
            help = "Number of processors to use when doing the"),
make_option( c( "--mx_diff", "-x" ), type = "logical", default = T,
            help = "Offers two approaches to calculate the enrichment (ES)
                    from the KS random walk statistic.mx_diff=FALSE: ES is
                    calculated as the maximum distance of the random walk
                    from 0. mx_diff=TRUE (default): ES is calculated as the
                    magnitude difference between the largest positive and
                    negative random walk deviations."),
make_option( c("--tau", "-t" ), type = "double",
            help = "Exponent defining the weight of the tail in the random
                    walk performed by both the gsva (Hanzelmann et al.,
                    2013) and the ssgsea (Barbie et al., 2009) methods. By
                    default, this tau=1 when method='gsva' and tau=0.25
                    when method='ssgsea' just as specified by Barbie et
                    al. (2009) where this parameter is called alpha."),
make_option( c("--WHICH_GROUP", "-q" ), type = "character",
            help = "[Optional]select the groupping column in metadata used to do subsettig of cell if necessary."),                                                                                          
make_option( c("--WHICH_CELLS", "-u" ), type = "character",
            help = "[Optional] the level id list in selected cell groupping used to celltyping.For some cell 
                   clusters with high hertergensity this may be useful for sub-celltyping
                   combined with the option -l/--LEVEL with single mode.
                   If not specified with cell clusters's ID, all cells will be used."),
make_option( c("--downsample", "-d"),type = "character", default = "40000",
            help = "the downsample number of cells ")

);
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);
#=================================================================================
#check the command line parameters
#=================================================================================
if ( is.null(opts$OUTDIR) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opts$OUTDIR) ){
        output_dir = opts$OUTDIR
    }else{
        output_dir = opts$OUTDIR
        dir.create(output_dir, recursive = T)
    }
}
output_dir = normalizePath(output_dir )


if (! is.null(opts$metadata) ){
    assay_metadata = read.csv(opts$metadata,sep=",",header =T, row.names=1)
}

#parameters for gsva
abs.ranking=as.logical(toupper(opts$abs_rank))
if ( is.null(opts$min_sz) ){
    min.sz = 2
}else{
    min.sz=as.numeric(opts$min_sz)
}

if ( !is.null(opts$method )){
    method = opts$method
}else{
    print("NO enrichment algrithms is AVAILABLE!GSVA will be used as default!")
    method = "gsva"
}

if( is.null(opts$max_sz)) {
    max.sz = 100000
} else {
    max.sz = as.numeric(opts$max_sz)
}

if ( is.null(opts$parallel_sz) ){
    parallel.sz= 6
}else{
    parallel.sz=as.numeric(opts$parallel_sz)
}

if ( is.null(opts$mx_dff) ){
    if(method == "gsva") {
        mx.diff= T
    } else {
        mx.diff= F
    }
}else{
    mx.diff=as.logical(toupper(opts$mx_diff))
}

if( is.null(opts$tau) ) {
    if(method == "gsva") {
        tau= 1
    } else {
        tau = 0.25
    }
}else{
    tau = as.numeric(tau)
}

if ( is.null( opts$kcdf) ){
    kcdf = "Poisson"
}else{
    kcdf=opts$kcdf
}
if (is.null(opts$downsample)) {
  downsample <- 40000
} else {
  downsample <- opts$downsample
}
#=================================================================================
#main steps of this script
#=================================================================================
#get the expression matrix and metadata to construct the
#expressionset object
if ( opts$informat == "seurat" ){
    suppressPackageStartupMessages(library(Seurat))
    seurat_ob = readRDS(opts$input)
    if ( seurat_ob@version < 3 ){
        seurat_ob = UpdateSeuratObject(seurat_ob)
    }
    if ( is.null(seurat_ob@meta.data$clusters) ){
        seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
    }else{
        seurat_ob = SetIdent( seurat_ob, value = "clusters")
    }
    
    if ( !is.null(opts$WHICH_CELLS)){
        cluster_list = unlist(strsplit( opts$WHICH_CELLS,",",perl = T))
        seurat_ob = SubsetData(seurat_ob,  subset.name = opts$WHICH_GROUP, accept.value = cluster_list)
    }
    if (ncol(seurat_ob) > 40000) {
            ratio <- as.numeric(opts$downsample) / ncol(seurat_ob)
            metadata_temp <- as.data.frame(seurat_ob@meta.data)
            cells_sample <- c()
            for (i in unique(seurat_ob$clusters)) {
              cells_temp <- rownames(metadata_temp)[which(metadata_temp$clusters == i)]
              cells_temp_sample <- sample(cells_temp, ceiling(length(cells_temp) * ratio), replace = FALSE, prob = NULL)
              cells_sample <- append(cells_sample, cells_temp_sample)
            }
            seurat_ob <- subset(seurat_ob, cells = cells_sample)
            print(dim(seurat_ob))
    }

    cellnames = Cells(seurat_ob)
    if ( "assay_metadata" %in% ls()) {
        #the index order is the same as the row index of the assay metadata
        sampleidx =  gsub("(_|-)[ATGC]{16,}.*","",cellnames,perl=T)
        #integrate the additional metadata from the assay design
        additional_cell_meta = vector( )
        for ( colidx in colnames(assay_metadata) ){
            additional_cell_meta = cbind(additional_cell_meta, as.vector(assay_metadata[sampleidx, colidx]))
        }
        colnames(additional_cell_meta) = colnames(assay_metadata)
        rownames(additional_cell_meta) = cellnames
        additional_cell_meta = as.data.frame(additional_cell_meta)
        seurat_ob = AddMetaData( seurat_ob, additional_cell_meta)

        assay_metadata = seurat_ob@meta.data %>%
                             tibble::rownames_to_column(var = "cellbarcode") %>%
                             select( cellbarcode, everything() )
        write.table(assay_metadata,file = file.path(output_dir,"metadata4each_cell.xls"),quote =F, col.names = T,row.names = F,sep = ",")
        # seurat_ob@raw.data = seurat_ob@raw.data[,colnames(seurat_ob@data)]
    }
    #make the geneSet object from the command line supplied gmt file for customize use
    gene_sets = getGmt(opts$gmt)
    suppressPackageStartupMessages(library("doParallel"))
    suppressPackageStartupMessages(library("foreach"))

    if ( !is.null( opts$chunkby) ){
        cell_num = length(Cells(seurat_ob))
        if ( cell_num < opts$chunkby){
            chunkby = cell_num
        }else{
            chunkby = opts$chunkby
        }
        intervals = seq(1,cell_num, by = chunkby)
        registerDoParallel(cores=length(intervals))
        gsva_enrichment = foreach( x = intervals ) %dopar%  gsva(
                        as.matrix( GetAssayData(seurat_ob, slot = "counts")[,cellnames[seq(x, min(x+chunkby-1, cell_num ))]]),
                        gene_sets, method=method, kcdf=kcdf, min.sz=min.sz, max.sz=max.sz,
                        parallel.sz=parallel.sz, mx.diff=mx.diff )
        doParallel::stopImplicitCluster()
        gsva_scores = as.data.frame(gsva_enrichment[[1]],check.names = F)
        for ( indx in 2:length(gsva_enrichment) ){
            gsva_scores = transform(merge(gsva_scores,as.data.frame(gsva_enrichment[[indx]]),by=0,all=F),
                                            row.names=Row.names, Row.names=NULL)
        }
    }else if ( !is.null(opts$splitby) ){
        seurat_ob_list = SplitObject(seurat_ob, attribute.1 = opts$splitby, subset.raw =T )
        registerDoParallel(cores=length(seurat_ob_list))
        gsva_enrichment = foreach( seurat_subset = seurat_ob_list )%dopar% gsva(
                                    as.matrix(GetAssayData(seurat_ob, slot = "counts")),gene_sets, method=method,
                                    kcdf=kcdf, min.sz=min.sz, max.sz=max.sz, parallel.sz=parallel.sz,
                                    mx.diff=mx.diff )
        doParallel::stopImplicitCluster()
        gsva_scores = as.data.frame(gsva_enrichment[[1]],check.names = F)
        for ( indx in 2:length(gsva_enrichment) ){
            gsva_scores = transform(merge(gsva_scores,as.data.frame(gsva_enrichment[[indx]]),by=0,all=F), row.names=Row.names, Row.names=NULL)
        }
    }else{
        gsva_scores = gsva( as.matrix(as.matrix(GetAssayData(seurat_ob, slot = "counts"))),gene_sets, method=method, kcdf=kcdf,
                            min.sz=min.sz, max.sz=max.sz, parallel.sz=parallel.sz,
                            mx.diff=mx.diff )
    }
    gsva_scores = as.data.frame(gsva_scores)
    gsva_scores$geneset = rownames(gsva_scores)
    gsva_scores = gsva_scores %>% select(geneset,everything())
    write.table(gsva_scores,file = file.path(output_dir,"GSVA_enrichment_results.xls"), col.names = T,row.names = F,sep = "\t")
    # eset <- varFilter(eset,var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
    # gsva_scores = gsva( exprs(eset),gene_sets, method=method, kcdf=kcdf,
    # min.sz=min.sz, max.sz=max.sz, parallel.sz=parallel.sz, mx.diff=mx.diff )
    rm(seurat_ob)
}else if ( opts$informat == "raw" ){ #for cell/sample raw expression matrix
    if ( !"assay_metadata" %in% ls()) { stop("Please provide metadata file for raw matrix using -m")}
    exprs_matrix = read.table(opts$input, sep ="\t", header = T, row.names = 1)
    exprs_matrix = as.matrix(exprs_matrix)
    eset <- ExpressionSet(assayData=expers_matrix,phenoData = AnnotatedDataFrame(assay_metadata))
    if (ncol(eset) > 40000) {
        ratio <- as.numeric(opts$downsample) / ncol(eset)
        metadata_temp <- as.data.frame(pData(eset))
        cells_sample <- c()
        for (i in unique(metadata_temp$clusters)) {
              cells_temp <- rownames(metadata_temp)[which(metadata_temp$clusters == i)]
              cells_temp_sample <- sample(cells_temp, ceiling(length(cells_temp) * ratio), replace = FALSE, prob = NULL)
              cells_sample <- append(cells_sample, cells_temp_sample)
            }
    eset <- eset[, cells_sample]
    }
    #make the geneSet object from the command line supplied gmt file for customize use
    gene_sets = getGmt(opts$gmt)
    gsva_scores = gsva( exprs(eset),gene_sets, method=method, kcdf=kcdf,
                        min.sz=min.sz, max.sz=max.sz, parallel.sz=parallel.sz,
                        mx.diff=mx.diff )
    gsva_scores = as.data.frame(gsva_scores)
    gsva_scores$geneset = rownames(gsva_scores)
    gsva_scores = gsva_scores %>% select(geneset,everything())
    #write.table(gsva_scores,file = file.path(output_dir,"GSVA_enrichment_results.xls"), col.names = T,row.names = F,sep = "\t")
    fwrite(gsva_scores,file.path(output_dir,"GSVA_enrichment_results.xls"),col.names=T,row.names=F,sep="\t",quote=F)
}
# else if ( opts$informat == "gsva" ){
#     gsva_scores = read.table(opts$input, sep ="\t", header = T, row.names = 1,check.names = F)
# }else{
#     stop("NO format information of input is AVAILABLE!")
# }

#here the gsva score act as the gene expression matrix with pathway as row names.
#So we can apply the general differential expression to the gsva score matrix to
# find the singificant expressed pathway.
# if ( !is.null(opts$group) ){
#     collapse_list = unlist(strsplit(opts$group,",", perl=T))
#     for ( group2collapse in collapse_list ){
#         gsva_collapsed_scores = collapse_count(gsva_scores,assay_metadata,collapseby=group2collapse,use="mean")
#         colnames(gsva_collapsed_scores) = paste0(group2collapse,colnames(gsva_collapsed_scores))
#         pdf(file.path(output_dir,paste0("GSVA_enrichment_score_by_",group2collapse,"_heatmap.pdf")))
#         ComplexHeatmap::draw( ComplexHeatmap::Heatmap(rbind(head(gsva_collapsed_scores,20),tail(gsva_collapsed_scores,20)),
#         row_names_gp = grid::gpar(fontsize = 8), name = "GSVA\nscore"),
#         heatmap_legend_side = "right", annotation_legend_side = "bottom" )
#         dev.off()
#     }
# }

