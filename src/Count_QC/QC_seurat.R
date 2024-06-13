#!/usr/bin/env Rscript
# This is the script used to carry out quality control of the gene expression matrix.
# Because of the variablity of single cell RNA-seq expression ,there are
# many confound factors to be denoised. The factors to be considered as follows:
#   1. filter out the cells with low quality, which is identified by setting a
#       threshold of the total number of genes or UMIs of the cell.
#   2. filter the genes with the threshold of counts across all the cells
#   3. filter the cells with high proportion of genes/UMIs mapped to the mitochondrion
#   4. filter the empty cells and multiplets which contains more than one cell in the droplet,
#      if the input matrix come from the cellranger count pipeline, the empty cell are filtered
#      as default since cellranger V3.0.
# All the QC process are carried out on the normalized matrix but raw counts. So you can
# rerun this script with different parametes for your purpose to find the best results.

#========import packages=====================================
rm(list=ls())
suppressWarnings({
    suppressPackageStartupMessages(library("Seurat"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("Matrix"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("future"))
    suppressPackageStartupMessages(library("future.apply"))
    suppressPackageStartupMessages(library("OESingleCell"))
    library(tictoc)
    library(DoubletFinder)
    library(Seurat)
})

RemoveDoublets <-function(
  object,
  identity,
  doublet.rate,
  pN=0.25,
  PCs=1:30,
  use.SCT=FALSE,
  num.cores=1,
  quietly=TRUE
){

  tic("get sweep parameters")
  # calculate parameters
  if (quietly==TRUE){
    invisible(capture.output(sweep.res.list <- paramSweep_v3(object, PCs = PCs, sct=use.SCT, num.cores=num.cores)))
    invisible(capture.output(sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)))
    ff <- tempfile()
    png(filename=ff)
    invisible(capture.output(bcmvn <- 
find.pK(sweep.stats)))
    dev.off()
    unlink(ff)
  }else{
    sweep.res.list <- paramSweep_v3(object, PCs = PCs, sct=use.SCT, num.cores=num.cores)
    sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
    ff <- tempfile()
    png(filename=ff)
    bcmvn <- 
find.pK(sweep.stats)
    dev.off()
    unlink(ff)
  }
  toc()

  # choose parameters
  maxBCmetric    <- max(bcmvn$BCmetric, na.rm = TRUE)
  pK <- as.numeric(as.character(bcmvn[bcmvn$BCmetric==maxBCmetric, ]$pK))

  # compute doublet scores
  tic("Removing doublets")
  annotations    <- object@meta.data$identity  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi       <- round(doublet.rate*length(colnames(x = object)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))
  seu.scored     <- doubletFinder_v3(object, PCs =PCs, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = use.SCT)
  toc()

  # pick out doublets
  cname <-colnames(seu.scored[[]])
  DF<-cname[grep('^DF',cname)]
  seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")

  # remove doublets
  seu.removed <- subset(seu.scored, subset = doublet != 1)
  return(list(removed=seu.removed, original=seu.scored))
}
rate=function(num){
	if (num >10000) rate= 0.1 else rate= num*0.0000077
	return(rate)
}

#=================================================================================
#function definition
#=================================================================================
# filter genes specified by user or according to the cell number where one gene can be detected
# min.value the minium
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

# filter the low quality cells with specified thresholds for different parameters
# keep the QC parameter in the seurat object for provinciable
CellRemover = function( object, parameters,lower_limit, high_limit, sdfold ){
    lower_threshold = vector()
    upper_threshold = vector()
    param_list = list()
    for ( paramx in parameters ){
        param_value = object@meta.data[,paramx]
        param_vector = param_value[param_value>0]
        mean4paramx = mean(log10(param_vector))
        sd4paramx = sd(log10(param_vector))
        upper_bound <- 10^(mean4paramx + sdfold*sd4paramx)
        lower_bound <- 10^(mean4paramx - sdfold*sd4paramx)
        if ( is.null( lower_limit[paramx]) | is.na( lower_limit[paramx]) |lower_limit[paramx] %in% c("NA","NULL") ){
            # lower_threshold = append( lower_threshold, lower_bound)
            lower_threshold = lower_bound
        }else{
            lower_threshold = as.numeric(lower_limit[paramx])
            # lower_threshold = append( lower_threshold, as.numeric(lower_limit[paramx]) )
        }
        if ( is.null( high_limit[paramx]) | is.na( high_limit[paramx]) | high_limit[paramx] %in% c("NA","NULL") ){
            # upper_threshold = append( upper_threshold, upper_bound)
            upper_threshold = upper_bound
        }else{
            # upper_threshold = append( upper_threshold, as.numeric(high_limit[paramx]) )
            upper_threshold = as.numeric(high_limit[paramx])
        }
        suppressWarnings({
            object = SubsetData(object, subset.name = paramx,
            low.threshold = lower_threshold, high.threshold = upper_threshold)
        })
        # param_list[[unique(seurat_ob@meta.data$sampleid)]][[paramx]] = c( lower = lower_threshold, upper = upper_threshold)
    }
    # object = SubsetData(object, subset.name = parameters,
    #           low.threshold = lower_threshold, high.threshold = upper_threshold)
    return( object )
}

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
                 help = "The input exprssion matrix in several possible format."),
    make_option( c("--informat", "-f" ), type = "character", default = "tenx",
                 help = "The indication of type of input expression matrix, the possible type can be:
                        seurat: the seurat object from the clustering results."),
    make_option( c("--minCell4gene","-x" ),type="double", default = 0.01,metavar = "minCell4gene",
        help="the minimium cell number one gene detected.If the value is less than 1,
             it is a proportion, otherwise a integer cell number.[default: 10]"),
    make_option( c("--output","-o"),type="character", default = "./Count_QC",
        help="the output directory of QC results.[default:./Count_QC]", metavar="outputdir"),
    make_option( c("--empty2remove", "-y" ), type = "logical", default = F, metavar = "logical",
                help = "[OPTIONAL]Wether to remove empty droplets for this project.[default:F]"),
    make_option( c("--doublets2remove", "-e" ), type = "logical", default = F, metavar = "logical",
                help = "[OPTIONAL] to remove doublets for this project.[default:F]"),
    make_option( c("--use_rlm","-u"), type = "logical", default = T, metavar = "logical",
                 help = "[OPTIONAL]using the robust linear model to remove the outliers of cells.[default:T]"),
    make_option( c("--removebatch","-b"), type = "character", default = NULL,
                 help = "use the integration method the remove the batch effect.[default:NULL]"),
    make_option( c("--filters","-c"), type = "character", default = "nFeature_RNA,nCount_RNA,percent.mito",
                help = "the filtering variables used to remove the outliers of cells.[default:nFeature_RNA,nCount_RNA,percent.mito];
                        The other options can be log10GenesPerUMI, percent.ribo etc."),
    make_option( c("--lower_threshold","-l"), type="character",
                help="the minimium value for the parameters specified by --filters.[default:200,1000]" ),
    make_option( c("--upper_threshold", "-L"), type="character",
                help="the maximium value list for the parameters specified by --filters.If not specified,
                     will be determined by the parameter fold2sd automatically." ),
    make_option( c("--fold2sd","-t" ),type="double", default = 2.0, metavar="double",
                help="the fold of standard variance for parameters specified by --filters determine the outlier threshold by mean +/- n*sd.[default:2.0]"),
    make_option( c("--normmeth", "-m" ), type="character", default = "LogNormalize", metavar="character",
                help="the normalization method to produce the data slot in the specified assay by --assay.
                       Options can be sctransform, LogNormalize, scran, CR, CLR etc. For feature-barcoding/cite-seq,
                       CLR is recommended. For scRNA-seq, sctransform is recommended. Notice that if sctransform used,
                       the default should be changed to SCT automatically in the later analysis."),
    make_option( c("--assay" ), type="character", default = "RNA",
                help="[OPTIONAL]the assay used to calulation in case of multimodal data."),
    make_option( c("--vars2regress", "-r" ), type = "character", default = "nCount_RNA,percent.mito",
                help = "unwanted source of variables to regress out from the single cell expression matrix.
                        The other option can be percent.ribo, cellcycle etc.[default:nCount_RNA,percent.mito]"),
    make_option( c("--genes2filter", "-z" ), type = "character", default = NULL,
                help = "[OPTIONAL]the specified gene list file with gene as header to remove from the gene-cell matrix."),
    make_option( c("--ncores", "-p"), type="integer", default = 10,
            help="the number of CPUs used to improve the performace."),
    make_option( c("--rmdoublets"), type="logical", default = TRUE,
            help="to Remove doublets useing DoubletFinder."),
    make_option( c("--mito", "-M" ), type = "character", default = NULL,
                  help = "mitochondrial gene list.") );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
# parse the command line parameters
#=================================================================================
# if ( !file.exists( opt$metadata) ){
#     assay_metadata = read.csv(opt$metadata,sep=",",header =T )
# }

if ( is.null( opt$fold2sd) ){
    print( "NO fold2sd AVAILABLE!The default will be used!")
    fold2sd = 2
}else{
    fold2sd = opt$fold2sd
}

if ( is.null( opt$filters) ){
    print("NO filtering parameters AVAILABLE! The default variables will be selected!")
    filter_params = c( "nFeature_RNA", "nCount_RNA", "percent.mito")
}else{
    filter_params = unlist(strsplit(opt$filters, ",", perl =T))
}

#determine the minimium number of cells one gene should be detected
if ( is.null(opt$minCell4gene) ){
    print( "NO minimium cell number threshold for a expressed gene PROVIDED,the default will be used!")
    min_cell4gene = 0.01
}else{
    min_cell4gene = opt$minCell4gene
}

if ( !is.null(opt$genes2filter) ){
    if ( file.exists( opt$genes2filter) ){
        genes2filter = read.csv(opt$genes2filter,sep=",",header =T )
    }
    genes2filter = as.vector(genes2filter$gene)
}else{
    genes2filter = opt$genes2filter
}

if ( is.null(opt$lower_threshold) ){
    print( "NO lower threshold for each variable AVAILABLE!The default will be used!")
    lower_threshold = c( nFeature_RNA = "NULL", nCount_RNA = "NULL",percent.mito = 0)
}else{
    lower_threshold = unlist(strsplit(opt$lower_threshold, ",", perl=T))
    if ( length(lower_threshold) != length(filter_params) ){
        stop("The lower threshold setting is not consistant with the parameters in --filters!")
    }
    names(lower_threshold) = filter_params
}

if ( is.null(opt$upper_threshold) ){
    print( "NO upper threshold for each variable is AVAILABLE!The default will be calculated according to their distribution!")
    upper_threshold = c( nFeature_RNA = "NULL", nCount_RNA = "NULL", percent.mito = 0.05)
}else{
    upper_threshold = unlist(strsplit(opt$upper_threshold, ",", perl=T))
    if ( length(upper_threshold) != length(filter_params) ){
        stop("The upper threshold setting is not consistant with the parameters in --filters!")
    }
    names(upper_threshold) = filter_params
}

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

if ( is.null(opt$vars2regress )){
    #the default regressed cofounders
    vars2regress = c("nCount_RNA", "percent.mito")
}else{
    vars2regress = unlist(strsplit(opt$vars2regress, ",", perl =T))
}

if ( is.null(opt$ncores) ) {
    nCores = 10
}else{
    nCores = opt$ncores
}

# #####################################################################
#====== Setup the Seurat Object (Required)=========================
# ####################################################################
options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
plan("multicore", workers = min(detectCores(), nCores)) # parallization using specified CPUs start from here

if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){# the input is a seurat object which may contain more than one sample
        seurat_ob = readRDSMC( opt$input , cores = availableCores())
        seurat_ob = SetIdent(seurat_ob, value = "sampleid")
        # in case of seurat pbject derived from lower versions, update it to be consensus with the
        # current version in this script
        if ( seurat_ob@version < 3){
            seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
        }

        # set the default assay to use
        DefaultAssay(object = seurat_ob) = opt$assay
    }
}

# The number of genes and UMIs are automatically calculated for every object by Seurat.
# For non-UMI data, nUMI represents the sum of the non-normalized values within a cell.
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using AddMetaData.
# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
# ,since this represents non-transformed and non-log-normalized counts.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
#calculate the mitochondrial gene derived transcript proportion and add to the metadata
raw.counts = GetAssayData(seurat_ob, slot = "counts")

if ( !opt$assay %in% c("ADT", "CRISPR") ){
    if ( "percent.mito" %in% vars2regress ){
        if ( !is.null(opt$mito) ){
            mito.genes <- read.delim(opt$mito)
            percent.mito <- Matrix::colSums(raw.counts[mito.genes[,1], ])/Matrix::colSums(raw.counts)
        }else{
            mito.genes <- grep(pattern = "^(MT|mt)(-|_)", x = rownames(seurat_ob), value = T,perl=T)
            percent.mito <- Matrix::colSums(raw.counts[mito.genes, ])/Matrix::colSums(raw.counts)
        }
        # percent.mito <- PercentageFeatureSet(seurat_ob, features = mito.genes )
        seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.mito, col.name = "percent.mito")
        vars2regress = append(vars2regress,"percent.mito")
    }

    #whether to regress out the ribosome genes
    if ( "percent.ribo" %in% vars2regress ){
        #Calculate percent ribosomal genes.
        ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]|^[mM]Rp[sl][[:digit:]]",
        rownames(seurat_ob), value = TRUE,ignore.case = T)
        percent.ribo <- Matrix::colSums(raw.counts[ribo.genes, ])/Matrix::colSums(raw.counts)
        # percent.ribo <- PercentageFeatureSet(seurat_ob, features = ribo.genes )
        seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.ribo, col.name = "percent.ribo")
        vars2regress = append(vars2regress,"percent.ribo")
    }

    if ( "top50" %in% filter_params ){
        top50_gene = apply( raw.counts, 2, function(x) sum(x[order(x, decreasing = T)][1:50])/sum(x))
        seurat_ob <- AddMetaData(object = seurat_ob, metadata = top50_gene, col.name = "top50")
    }

    if ( "UMI.per.gene" %in% filter_params ){
        tmp = seurat_ob@meta.data$nCount_RNA/seurat_ob@meta.data$nFeature_RNA
        names(tmp) = row.names(seurat_ob@meta.data)
        seurat_ob = AddMetaData(seurat_ob, metadata = tmp, col.name = "UMI.per.Gene")
    }

	if ( "log10GenesPerUMI" %in% filter_params ){
		seurat_ob@meta.data$log10GenesPerUMI <- log10(seurat_ob@meta.data$nFeature_RNA)/log10(seurat_ob@meta.data$nCount_RNA)  
	}

    #statistics of metadat before QC
    statistics_beforeQC = seurat_ob@meta.data %>%
        group_by(sampleid) %>%
        do(data.frame( Mean_nUMI_beforeQC = mean(.$nCount_RNA),
        Mean_nGene_beforeQC = mean(.$nFeature_RNA),
        Mean_mito.percent_beforeQC = mean(.$percent.mito),
        Total_cells_beforeQC = nrow(.) ))
    # rename(sampleid=orig.ident)

    #=================================================================================
    # remove the empty droplets using dropletuitls if necessary
    if ( opt$empty2remove == T ){
        seurat_ob = RemoveEmptyDrops(seurat_ob, emptyDropNIters = 1000)
    }

    if ( opt$doublets2remove == T ){
        seurat_ob = RemoveDoublets(seurat_ob)
    }

    #remove the outliers cells using the linear model
    if ( opt$use_rlm == TRUE ){
        UMIs_per_cell = Matrix::colSums(raw.counts)
        genes_per_cell = Matrix::colSums(raw.counts>0) # count gene only if it has non-zero reads mapped.
        df = data.frame(UMIs_per_cell=UMIs_per_cell, genes_per_cell=genes_per_cell)
        suppressPackageStartupMessages(library(MASS))
        df = df[order(df$UMIs_per_cell),] # order by UMIs_per_cell
        pdf(file.path(output_dir,"outliers.pdf"))
        plot(df, log='xy')
        m <- rlm(genes_per_cell~UMIs_per_cell,data=df) # robust linear model, not sens to outliers
        p.level = 1e-3
        # predict genes_per_cell based on observed UMIs_per_cell
        suppressWarnings(pb <- data.frame(predict(m, interval='prediction',
        level = 1-p.level, # define conf interval
        type="response")))
        polygon(c(df$UMIs_per_cell, rev(df$UMIs_per_cell)), c(pb$lwr, rev(pb$upr)), col=adjustcolor(2,alpha=0.1), border = NA)
        # identifier outliers as having observed genes_per_cell outside the prediction confidence interval
        outliers <- rownames(df)[df$genes_per_cell > pb$upr | df$genes_per_cell < pb$lwr]
        points(df[outliers,],col=2, bg="red",pch = 16,cex=0.6)
        title("Outlier of Cells")
        dev.off()
        png(file.path(output_dir,"outliers.png"))
        plot(df, log='xy')
        m <- rlm(genes_per_cell~UMIs_per_cell,data=df) # robust linear model, not sens to outliers
        p.level = 1e-3
        # predict genes_per_cell based on observed UMIs_per_cell
        suppressWarnings(pb <- data.frame(predict(m, interval='prediction',
        level = 1-p.level, # define conf interval
        type="response")))
        polygon(c(df$UMIs_per_cell, rev(df$UMIs_per_cell)), c(pb$lwr, rev(pb$upr)), col=adjustcolor(2,alpha=0.1), border = NA)
        # identifier outliers as having observed genes_per_cell outside the prediction confidence interval
        outliers <- rownames(df)[df$genes_per_cell > pb$upr | df$genes_per_cell < pb$lwr]
        points(df[outliers,],col=2, bg="red",pch = 16,cex=0.6)
        title("Outlier of Cells")
        dev.off()


        seurat_ob = subset(seurat_ob, cells = colnames(raw.counts)[!colnames(raw.counts) %in% outliers])
    }


    # We filter out cells that have unique gene counts over 2,500 or less than
    # 200. Note that low.thresholds and high.thresholds are used to define a
    # 'gate'. -Inf and Inf should be used if you don't want a lower or upper
    # threshold.
    # if each filter's lower or upper threshold is specified from the command line
    # the user defined value will be used. Otherwise,the default threshold will
    # be calculated using the distribution significant level specified automatically
    # All the thresholds are specific for each sample, so when the threshold and the
    # standard variance fold are both
    # available from the command line, we will compare the value with the one derived
    # from the distribution lower outliner
    if ( length(unique(as.vector(seurat_ob@meta.data$sampleid))) > 1 ){
        seurat_by_sample = future_lapply(SplitObject(seurat_ob, split.by = "sampleid" ),
        function(x)CellRemover(object = x, parameters = filter_params,lower_limit = lower_threshold,
        high_limit = upper_threshold, sdfold = fold2sd ))
        merged_seurat = seurat_by_sample[[1]]
        for( idx in 2:length(seurat_by_sample) ){
            merged_seurat = merge(x = merged_seurat, y = seurat_by_sample[[idx]],
            do.scale = F, do.center = F, do.normalize = F)
        }
        seurat_ob = merged_seurat
        # if ( dim( seurat_ob@data )[2]/dim( seurat_ob@raw.data)[2] < 0.7) {
        #     warning("The cells after filtering is less than 70% of the previous, Which means possible over-filtering.")
        # }
    }else{
        seurat_ob = CellRemover(object = seurat_ob,parameters = filter_params,lower_limit = lower_threshold,
        high_limit = upper_threshold, sdfold = fold2sd )
        # if ( dim( seurat_ob@data )[2]/dim( seurat_ob@raw.data)[2] < 0.7) {
        #     warning("The cells after filtering is less than 70% of the previous, Which means possible over-filtering.")
        # }
    }


    # filter genes by using the minimium cell number one gene is detected
    # this step shoud be run after cell filtering, because there may be
    # satisfied cell number for one gene before filtering
    # but fails after cell filtering
    if ( min_cell4gene < 1 ){ #the parameter is a percentage
        min.cell_N = round(min_cell4gene * ncol(seurat_ob))
    }else{ #the parameter is a integer
        min.cell_N = min_cell4gene
    }
    seurat_ob = FilterGenes(seurat_ob, min.cells = min.cell_N, filter.genes = genes2filter )

    # statistics of metadat after QC
    statistics_afterQC = seurat_ob@meta.data %>%
        group_by(sampleid) %>%
        do(data.frame(
        Mean_nUMI_afterQC = mean(.$nCount_RNA),
        Mean_nGene_afterQC = mean(.$nFeature_RNA),
        Mean_mito.percent_afterQC = mean(.$percent.mito),
        Total_cells_afterQC = nrow(.) ))
    #     #merge the statitics by columns
    #cell_statitics = right_join(statistics_beforeQC,statistics_afterQC,by="sampleid")
    #write.table(cell_statitics,file.path(output_dir,"cell_statitics_before_after_QConly.xls"),sep="\t",col.names=T,row.names=F)


    ########## DoubletFinder ##########
    if ( opt$rmdoublets == TRUE ){
        print("removing doublet using DoubletFinder")
        seurat_ob = SetIdent( seurat_ob, value = "sampleid")
        obj = SplitObject(seurat_ob,split.by = "sampleid")
        obj_rm=list()
        for( i in names(obj)){
            obj[[i]] <- NormalizeData(obj[[i]])
            obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
            obj[[i]] <- ScaleData(obj[[i]])
            obj[[i]] <- RunPCA(obj[[i]])
            obj[[i]] <- RunUMAP(obj[[i]], dims = 1:10)
            obj_rm[[i]] = RemoveDoublets(obj[[i]], doublet.rate=rate(dim(obj[[i]])[2]),  num.cores=4)
        }
        # metadata + doublet pANN DF.classifications
        removed= lapply(obj_rm,FUN=function(x) x = x$removed)
        if ( length(removed) > 1 ) {
            seurat_ob = merge(removed[[1]],do.call(c,removed[-1]))
        } else {
            seurat_ob = removed[[1]]
        }
        #saveRDS(seurat_ob,"doublet_seurat.rds")

        # statistics of metadat after QC
        statistics_afterDF = seurat_ob@meta.data %>%
            group_by(sampleid) %>%
            do(data.frame(
            Mean_nUMI_afterQC = mean(.$nCount_RNA),
            Mean_nGene_afterQC = mean(.$nFeature_RNA),
            Mean_mito.percent_afterQC = mean(.$percent.mito),
            Total_cells_afterQC = nrow(.) ))
        #     #merge the statitics by columns
        cell_statitics = right_join(statistics_beforeQC,statistics_afterDF,by="sampleid")
        write.table(cell_statitics,file.path(output_dir,"cell_statitics_before_after_QC.xls"),sep="\t",col.names=T,row.names=F)

    } else {
        print("regular QC only.")
        cell_statitics = right_join(statistics_beforeQC,statistics_afterQC,by="sampleid")
        write.table(cell_statitics,file.path(output_dir,"cell_statitics_before_after_QC.xls"),sep="\t",col.names=T,row.names=F)
    }

    #remove the batch effect using the integration method
    if ( !is.null(opt$removebatch) ){
        anchors = FindIntegrationAnchors( object.list = SplitObject(seurat_ob, split.by = opt$removebatch ), dims = 1:20)
        seurat_ob = IntegrateData(anchorset = anchors, dims = 1:20)
    }
}


if ( tolower(opt$normmeth ) == "sctransform" ){
    seurat_ob = SCTransform(seurat_ob, vars.to.regress = vars2regress,
                            verbose = FALSE,return.only.var.genes = FALSE)
}else{
    # normalize the UMI count data in assay using the
    # specified method by normalization.method parameter
    seurat_ob <- NormalizeData(object = seurat_ob,
                            normalization.method = opt$normmeth,scale.factor = 10000)
    # How to choose top variable features. Choose one of :
    # vst: First, fits a line to the relationship of log(variance)
    #       and log(mean) using local polynomial regression (loess).
    #       Then standardizes the feature values using the observed
    #       mean and expected variance (given by the fitted line).
    #       Feature variance is then calculated on the standardized
    #       values after clipping to a maximum (see clip.max parameter).
    # mean.var.plot: First, uses a function to calculate average expression
    #       (mean.function) and dispersion (dispersion.function) for
    #       each feature. Next, divides features into num.bin (deafult 20)
    #       bins based on their average expression, and calculates z-scores
    #       for dispersion within each bin. The purpose of this is to
    #       identify variable features while controlling for the strong
    #       relationship between variability and average expression.
    # dispersion: selects the genes with the highest dispersion values

    # Cell Cycle effect
    # use the CaseMatch implemented in seurat package which has the same function of following functions
    # s.genes = gene4data[unlist(lapply( cc.genes$s.genes,
    #                                     function(x)grep(paste0("^",x,"$"),gene4data, perl = T, ignore.case = T)))]
    # g2m.genes = gene4data[unlist(lapply(cc.genes$g2m.genes,
    #                                     function(x)grep(paste0("^",x,"$"),gene4data,gene4data, perl = T, ignore.case = T)))]
    if (seurat_ob@meta.data$specie == "GRCh38" || seurat_ob@meta.data$specie == "mm10" ){
        genes.inuse = rownames(GetAssayData(seurat_ob, slot="counts"))
        s.genes = CaseMatch(search = cc.genes$s.genes, match = genes.inuse)
        g2m.genes = CaseMatch(search = cc.genes$g2m.genes, match = genes.inuse)
        bin=48
        flag = NA
        while( (class(flag)=="try-error" || is.na(flag[1]) ) && bin > 1  ) {
            bin = round(bin/2)
            flag = try( CellCycleScoring(object = seurat_ob, s.features = s.genes, g2m.features = g2m.genes,set.ident = F,nbin=bin))
        } ; print(paste0("number of bin used: ",bin))
        seurat_ob <- flag
        if ( "cellcycle" %in% vars2regress ){
            CC.Difference <- seurat_ob@meta.data$S.Score - seurat_ob@meta.data$G2M.Score
            seurat_ob <- AddMetaData(object = seurat_ob, metadata = CC.Difference, col.name = "CC.Difference")
            # vars2regress = c(vars2regress,c("S.Score", "G2M.Score"))
            vars2regress = c(vars2regress, "CC.Difference")
            vars2regress = vars2regress[-which(vars2regress=="cellcycle")]
        }
    }

    if ( opt$assay != "RNA" ){
        seurat_ob <- ScaleData(object = seurat_ob, features = rownames(seurat_ob), verbose = T )
    }else{
        #remove unwanted source of variance
        vars2regress = unique(vars2regress)
        seurat_ob = FindVariableFeatures(object= seurat_ob, loess.span = 0.3,
                            clip.max = "auto", mean.function = "FastExpMean",
                            dispersion.function = "FastLogVMR", num.bin = 20,
                            nfeature = 4000, binning.method = "equal_width" )
        # regress out all the specified the varibales
        seurat_ob <- ScaleData(object = seurat_ob, features = rownames(seurat_ob),
                            vars.to.regress = vars2regress, verbose = T ) #takes some time
    }
}

saveRDSMC(seurat_ob, file.path(dirname(output_dir), "filtered_seurat.rds"), threads = nCores)
