#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_PYTHON = "/home/lipeng/miniconda3/envs/Scenic/bin/python3")
# construct the gene regulation network from the scRNA-seq
# expression matrix with TF annotations
# this work is based on the manuscript:
# Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering.
# Nature Methods. doi: 10.1038/nmeth.4463.

#' fast parallel block-based cor function for large gene expression matrix adapted from
#' https://rmazing.wordpress.com/2013/02/22/bigcor-large-correlation-matrices-in-r/
#'
#' @param x the gene expression matrix, where sample/cell is row and gene in on the column.
#' @param y NULL (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).
#' @param size chunck size for chunck calculation and parallization.
#' @param cores the number of threads used to run
#' @param fun the function used to call.
#' @param method the method for cor, options can be peason, spearman,kendrall. pearson as default.
#' @param verbose
#'
bigcor <- function(
  x,
  y = NULL,
  # fun = c("cor", "cov"),
  size = 2000,
  cores = 8,
  verbose = TRUE,
  ...
) {
  tictoc::tic()
  # if (fun == "cor") FUN <- cor else FUN <- cov
  # if (fun == "cor") STR <- "Correlation" else STR <- "Covariance"
  if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")

  NCOL <- ncol(x)
  if (!is.null(y)) YCOL <- NCOL(y)

  ## calculate remainder, largest 'size'-divisible integer and block size
  REST <- NCOL %% size
  LARGE <- NCOL - REST
  NBLOCKS <- NCOL %/% size

  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  if (is.null(y)) {
    resMAT <- ff::ff(vmode = "double", dim = c(NCOL, NCOL))
  }else{
    resMAT <- ff::ff(vmode = "double", dim = c(NCOL, YCOL))
  }

  ## split column numbers into 'nblocks' groups + remaining block
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0){
    GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  }
  SPLIT <- split(1:NCOL, GROUP)

  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))

  require(doMC)
  ncore = min(future::availableCores(), cores)
  doMC::registerDoMC(cores = ncore)
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  results <- foreach(i = 1:nrow(COMBS)) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    ## if y = NULL
    if (is.null(y)) {
      if (verbose) message("bigcor: ", sprintf("#%d:Block %s and Block %s (%s x %s) ... ",
                                               i, COMB[1], COMB[2], length(G1),  length(G2)))
      flush.console()
      RES<- do.call("cor", list(x = x[, G1], y = x[, G2], ... ))
      # RES <- FUN(x[, G1], x[, G2], ...)
      resMAT[G1, G2] <- RES
      resMAT[G2, G1] <- t(RES)
    } else {## if y = smaller matrix or vector
      if (verbose) message("bigcor: ", sprintf("#%d:Block %s and 'y' (%s x %s) ... ",
                                               i, COMB[1], length(G1),  YCOL))
      flush.console()
      RES<- do.call("cor", list(x = x[, G1], y = y, ... ))
      # RES <- FUN(x[, G1], y, ...)
      resMAT[G1, ] <- RES
    }
  }

  if ( is.null(y) ){
    resMAT <- resMAT[1:ncol(x),1:ncol(x)]
    colnames(resMAT) <- colnames(x)
    rownames(resMAT) <- colnames(x)
  }else{
    resMAT <- resMAT[1:ncol(x),1:ncol(y)]
    colnames(resMAT) <- colnames(x)
    rownames(resMAT) <- colnames(y)
  }
  tictoc::toc()
  return(resMAT)
}

write.gmt <- function(geneSet=kegg2symbol_list,gmt_file='kegg2symbol.gmt'){
   sink( gmt_file )
   for (i in 1:length(geneSet)){
     cat(names(geneSet)[i])
     #cat('\tNA\t')
     cat('\t')
     cat(paste(geneSet[[i]],collapse = '\t'))
     cat('\n')
   }
   sink()
}


FilterGenes <- function (object, min.value=1, min.cells = 0, filter.genes = NULL ) {
  genes.use <- rownames(object)

  if (min.cells > 0) {
    num.cells <- Matrix::rowSums( GetAssayData(object,slot="counts") > min.value)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
      object = subset( object, features = genes.use)
    # object = SetAssayData(object, new.data = GetAssayData(object, slot ="data")[genes.use,])
    # object@data <- object@data[genes.use, ] # Seurat V2.x
  }
  if (!is.null(filter.genes)) {
    filter.genes = CaseMatch(search = filter.genes, match = rownames(seurat_ob))
    genes.use <- setdiff(genes.use, filter.genes) #keep genes not in filter.genes
      object = subset( object, features = genes.use)
    # object = SetAssayData(object, new.data = GetAssayData(object, slot ="data")[genes.use,])
    # object[["RNA"]]@data = object[["RNA"]]@data[genes.use,]
    # object@data <- object@data[genes.use, ] #seurat V2.x
  }
  object <- LogSeuratCommand(object)
  return(object)
}

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("SCENIC"))
suppressPackageStartupMessages(library("reticulate"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("dplyr"))
#suppressPackageStartupMessages(library("OESingleCell"))


#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
        help = "The filtered exprssion matrix in prefered format."),
    make_option( c("--informat", "-f" ), type = "character", default = "seurat",
        help = "The indication of type of input expression matrix, the possible type can be:
                            seurat: the seurat object from the clustering results."),
    make_option( c("--cisTargetdb", "-d"), type = "character",
        help = "[REQUIRED]the glob of downloaded official motif annotation database for specified spieces.
                The files are usually suffixed with .feather, hg19-500bp-upstream-7species.mc9nr.feather as example."),
    make_option( c("--species", "-s"), type = "character",
        help = "[REQUIRED]the spieces abstraction, the current options can be 'hgnc' for human and 'mgi' for mouse."),
    make_option( c("--minCell4gene","-x" ),type="double", default = 0.01,metavar = "minimium proportion",
        help="the minimium cell number one gene detected.If the value is less than 1,
                it will be interpreted as a proportion."),
    make_option( c("--tfs"), type = "character",
        help = "the TF list of interested species in file"),
    make_option( c("--coexMethod"), type = "character",
        help = "the co-expression method for finding regulons. The current supported methods are
                w001,w005,top50,top50perTarget,top10perTarget,top5perTarget."),
    make_option( c("--ncores", "-j"), type = "integer", default = 10,
        help = "[OPTIONAL]the CPUs used to run this job, the more the better for project with more than 10k cells."),
    make_option(c("--downsample", "-e"),type = "character", default = "30000",
        help = "the downsample number of cells "),
    make_option( c("--hvg","-v"),type="logical", default = F,
        help="use the high variable features to do analysis", metavar="high variable features"),
    make_option( c("--extended"),type="logical", default = F,
        help="whether to use the extended regulons for calculation and visualization"),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory.", metavar="outputdir")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# setting the output directory
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir,recursive = T)
    }
}
output_dir = normalizePath(output_dir )

if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){# the input is a seurat object which may contain more than one sample
        seurat_ob = readRDS( opt$input )
        # in case of seurat pbject derived from different version, update it to be consensus with the
        # current version in this script
        if ( seurat_ob@version < 3){
            seurat_ob = UpdateSeuratObject(seurat_ob)
        }
    }
}

if ( opt$extended == "F" ){
   non_extended = FALSE
}else{
    non_extended = TRUE
}

if ( is.null( opt$downsample) ){
  downsample = 30000
}else{
  downsample = opt$downsample
}
if (ncol(seurat_ob) > 40000) {
   # library(sampling) need install,instead with manual function
    ratio <- as.numeric(downsample) / ncol(seurat_ob)
    metadata_temp <- as.data.frame(seurat_ob@meta.data)
    # strata(metadata_temp,stratanames="clusters",ratio,description=FALSE)
    cells_sample <- c()
    for (i in unique(seurat_ob$clusters)) {
        cells_temp <- rownames(metadata_temp)[which(metadata_temp$clusters == i)]
        cells_temp_sample <- sample(cells_temp, ceiling(length(cells_temp) * ratio), replace = FALSE, prob = NULL)
        cells_sample <- append(cells_sample, cells_temp_sample)
       }
    seurat_ob <- subset(seurat_ob, cells = cells_sample)
	saveRDS(seurat_ob,"sub_cells.rds")
   }
print(dim(seurat_ob))


if ( is.null( opt$species ) ){
    stop("NO motif annotation for your specified species.")
}else{
    # human = "hgnc", mouse = "mgi"
    org = tolower(opt$species)
}
scenicOptions <- initializeScenic(org=org, dbDir= opt$cisTargetdb, nCores=opt$ncores)
#scenicOptions@settings$seed <- 123

# do primary filtering of genes if it was not carried out before.
# filter genes by using the minimium cell number one gene is detected
# this step shoud be run after cell filtering
if ( opt$minCell4gene < 1 ){ #the parameter is a percentage
    minCell4gene = round(opt$minCell4gene * ncol(seurat_ob))
}else{ #the parameter is a integer
    minCell4gene = opt$minCell4gene
}

# determine the final expression matrix according to the usage of variable features
# filte genes
# seurat_ob = FilterGenes(seurat_ob, min.cells = minCell4gene, filter.genes = NULL )
if ( opt$hvg ){
    exprMat = GetAssayData(seurat_ob, slot = "counts")[VariableFeatures(seurat_ob),]
}else{
    exprMat = GetAssayData(seurat_ob, slot = "counts")
}

genesKept <- geneFiltering(as.matrix(exprMat),
                    scenicOptions=scenicOptions,
                    minCountsPerGene=1,
                    minSamples=minCell4gene)
exprMat_filtered <- exprMat[genesKept, ] # now the final matrix used for GRN inference

#=================================================================================
# Step1. Inference of co-expression modules
#=================================================================================
# load the list of known TFs for specified species from the command line file or annotated database
# R version:run Genies3 to do co-expression network analysis
# but for the high performance we change to python version here
# runGenie3(exprMat_filtered, scenicOptions)
# python version for this step using GRNBoost from arboreto
if ( !is.null(opt$tfs) ){
    arb.util = import("arboreto.utils")
    tf_names = arb.util$load_tf_names(opt$tfs)
}else{
    tf_names = getDbTfs( scenicOptions )
}
tf_names = CaseMatch(search=tf_names, match = rownames(seurat_ob)) # all the tf should be in the genes of matrix
arb.algo = import("arboreto.algo")
#adjacencies = arb.algo$grnboost2(as.data.frame(as.matrix(t(exprMat_filtered))), tf_names=tf_names, verbose=T, seed=123L)
adjacencies = arb.algo$grnboost2(as.data.frame(t(as.matrix(exprMat_filtered))), tf_names=tf_names, verbose=T, seed=123L)
colnames(adjacencies) = c( "TF", "Target", "weight" )
saveRDS( adjacencies, file = getIntName(scenicOptions, "genie3ll") )

# correlation analysis to distinguish positive from negative activity for each TF target
    # detect the block size according to the number of genes in count matrix
corrMat = bigcor(t(as.matrix(exprMat_filtered)),size = 2000, cores = opt$ncores, method = "spearman")
saveRDS(corrMat, file= getIntName(scenicOptions, "corrMat"))

### Build and score the GRN
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod = opt$coexMethod)

scenicOptions <- initializeScenic(org=org, dbDir= opt$cisTargetdb, nCores=1)
runSCENIC_3_scoreCells(scenicOptions,log2(as.matrix(exprMat_filtered)+1))

regulonAUC  = loadInt(scenicOptions, "aucell_regulonAUC")
if ( non_extended ){
    regulonAUC = regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
}
regulonAUC_mat = getAUC(regulonAUC)
rownames(regulonAUC_mat) = gsub("_", "-", rownames(regulonAUC_mat))
regulonAUC_mat_out = regulonAUC_mat[-grep(pattern="-extended",rownames(regulonAUC_mat)),]
write.table(as.data.frame(regulonAUC_mat_out) %>% tibble::rownames_to_column(var = "regulon"),
            file.path(output_dir,"regulon_activity.xls"),
            sep = "\t", col.names =T, row.names =F)

seurat_ob[["SCENIC"]] = CreateAssayObject(counts = regulonAUC_mat)
seurat_ob = ScaleData(seurat_ob, assay = "SCENIC")
seurat_ob@tools$RunAUCell = regulonAUC
# Tool(seurat_ob) = regulonAUC
# saveRDSMC(seurat_ob,"SCENIC_seurat.rds")

regulonTargetsInfo = loadInt(scenicOptions, "regulonTargetsInfo")
write.table(regulonTargetsInfo,
            file.path(output_dir, "0.1.TF_target_enrichment_annotation.xls"),
            sep = "\t", col.names =T, row.names =F, quote =F)

regulons <- loadInt(scenicOptions, "regulons")
sub_regulons = gsub(" .*","",rownames(regulonAUC_mat_out))
regulons = regulons[sub_regulons]
write.gmt(regulons, gmt_file = file.path(output_dir, "0.2.regulon_annotation.xls"))

if(!file.exists(file.path(output_dir, "Regulon调控子分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/Regulon调控子分析说明.docx",
  file.path(output_dir, "Regulon调控子分析说明.docx"))
}
if(!file.exists(file.path(output_dir, "0.3.MotifEnrichment_preview.html"))){
  file.copy("./output/Step2_MotifEnrichment_preview.html",
  file.path(output_dir, "0.3.MotifEnrichment_preview.html"))
}
