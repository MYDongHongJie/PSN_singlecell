#!/usr/bin/env Rscript
#this script is used to find all the differential expressed gene for single cell UMI counts
#and pick the marker genes with specified restrictions

rm(list=ls())
suppressWarnings({
    suppressPackageStartupMessages(library("future"))
    suppressPackageStartupMessages(library("future.apply"))
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("Matrix") )
    suppressPackageStartupMessages( library("optparse") )
    suppressPackageStartupMessages( library("OESingleCell") )
    suppressPackageStartupMessages(library("dplyr"))
})

# Calculate ROC score for each gene in seperating clusters
# data DEG matrix from FindAllMarkers
addROC <- function( object, data, assay = "RNA" ){
    DefaultAssay(object) <- assay
    expr = GetAssayData( object, slot = "data")
    data$auc <- 0
    lapply( unique(data$cluster), function(i){
        row.idx <- which(data$cluster == i)
        for( j in row.idx ){
            this.gene <- data$gene[j]
            this.exp <- expr[this.gene,]
            this.class <- as.character(Idents(object) )
            this.class[this.class == i ] <- 1
            this.class[this.class != 1 ] <- 0
            data$auc[j] <- as.numeric( pROC::roc(this.class, this.exp)$auc )
        }
    })
    return(data)
}

#=command line parameters setting=============================
option_list = list(
    make_option( c("--RDS", "-v"), type = "character", default = "TRUE",
                 help = "the seurat object saved as R object in RDS format."),
    make_option( c("--min_pct1", "-t"), type = "double", default = 0.5,metavar = "Filter option",
                 help ="the minimium ratio of cells expressing one specific gene in a cluster."),
    make_option( c("--max_pct2", "-T"), type = "double", default = 0.5,metavar = "Filter option",
                 help ="the maximiium ratio of cells expressing one specific gene in all other clusters."),
    make_option( c("--pct_fold", "-c"), type = "integer", default = 2,metavar = "Filter option",
                 help ="the minimiu fold of pct1 for gene in a specific cluster against pct2 for all other cluster."),
    make_option( c("--topn_marker","-N"), type = "integer",default = 10,metavar = "Filter option",
                 help = "the maximium number of ranked marker genes on the top for each cluster "),
    make_option( c("--avg_logFC","-k"), type = "double", default = 1,metavar = "Filter option",
                 help = "The average logFC of the gene UMI count in its cluster against the all other clusters."),
    make_option( c("--pvalue","-p"), type = "double", default = 0.05,
                 help = "the P-value of the gene differential expression.",metavar = "Filter option"),
    make_option( c("--cluster_name", "-n"), type = "character",
                 help = "the name of groupping column from clustering result metadata used to find markers. The example
                 can be tsne.2D.res.1."),
    make_option( c("--FDR","-q"), type = "double", default = 0.05,
                 help = "the FDR of the gene differential expression.",metavar = "Filter option"),
    make_option( c("--output","-o"),type="character", default = "./",
                 help="the output directory of QC results.", metavar="character"),
    make_option( c("--ncores","-j"),type="integer", default = "8",
                 help="the number of CPUs used to parallize this job."),
    make_option( c("--strict","-s"), type = "logical", default = F,
                 help = "whether to use strict mode, which will use all the filtering options, to find the markers. Notice that,
                 this may result in no markers remain for some clusters."),
    make_option( c("--assay" ), type="character", default = "RNA",
                 help="[OPTIONAL]the assay used to calulation in case of multimodal data."),
    make_option( c("--DEGtest","-e"), type = "character", default = "bimod",
                 help = "the test methods used to find differential expressed genes.Options are:
                   'wilcox' : Wilcoxon rank sum test (default).
                   'presto' : performace improved Wilcoxon rank sum test for big data.
                   'roc' : Standard AUC classifier.
                   't' : Student\'s t-test.
                   'bimod' : Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013).
                   'tobit' : Tobit-test for differential gene expression (Trapnell et al., Nature Biotech, 2014).
                   'poisson' : Likelihood ratio test assuming an underlying poisson distribution. Use only for UMI-based datasets.
                   'negbinom' : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets.
                   'MAST' : GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015).
                   'DESeq2' : DE based on a model using the negative binomial distribution (Love et al, Genome Biology, 2014).
                   'Venice' : package from bioTuning with function VeniceAllMarkers imtating the FindAllMarkers in Seurat but much faster.")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$RDS) ){
  stop("the seurat object is NOT AVAILABLE!")
}else{
  seurat_ob = readRDSMC(opt$RDS)
  if ( seurat_ob@version < 3){
    seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
  }
  # change the default assay for reduction if necessary
  if ( !is.null( opt$assay) ){
      DefaultAssay(seurat_ob) = opt$assay
  }else{
      DefaultAssay(seurat_ob) = "RNA"
  }
}

if ( is.null(opt$min_pct1) ){
  pct1_cutoff = 0.5
}else{
  pct1_cutoff = opt$min_pct1
}

if ( is.null(opt$max_pct2) ){
  pct2_cutoff = 0.5
}else{
  pct2_cutoff = opt$max_pct2
}

if ( is.null(opt$pct_fold) ){
  pct_fold_cutoff = 2
}else{
  pct_fold_cutoff = opt$pct_fold
}
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output)){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir, recursive = T)
    }
}

if ( is.null(opt$DEGtest) ){
  print("NO differential expression test method supplied,using the bimod as default,")
  dftest = "bimod"
}else{
  dftest = opt$DEGtest
}

if ( is.null(opt$topn_marker) ){
  print("The top n ranked marker genes for each cluster will be preserved!")
  topn = opt$topn_marker
}else{
  topn = opt$topn_marker
}

#in default, the FindAllMarkers() function will use the default identity class
#In order to select the different clustering result to find markers, change this
#class value to one of the column of the stored clustering in seurat_object@meta.data
suppressWarnings(
    if ( !opt$cluster_name %in% colnames(seurat_ob@meta.data) ){
        seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
    }else{
        seurat_ob = SetIdent( seurat_ob, value = opt$cluster_name)
    }
)

# setting the cores for parallization
options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
plan("multicore", workers = min(availableCores(), opt$ncores)) # parallization using specified CPUs start from here

#find the differential expressed genes for each clusters
#FindAllMarkers() is primarily used to find markers, but here it was also used
#to find differentially expressed genes. The default logfc threshold was set to 0.25
#,here we set to 0 as no prefiltering and then manually filter the genes to find markers
# note that if test.use is "negbinom", "poisson", or "DESeq2", slot will be set to "counts" automatically
if ( tolower(dftest) == "presto" ){
    suppressPackageStartupMessages( library("presto") )
    global_DEGs = wilcoxauc(seurat_ob, group_by  = opt$cluster_name, assay = "data")
    data = GetAssayData(seurat_ob, slot = "counts")
    results = future_lapply(unique(Idents(seurat_ob)), function(idx){
        cells.1 = Cells(subset(seurat_ob, idents = idx))
        cells.2 = Cells(subset(seurat_ob, idents = idx, invert = T))
        deg4idx = global_DEGs %>% filter(group == idx)
        pct.1 = round( Matrix::rowSums(data[deg4idx$feature, cells.1] > 0)/length(cells.1), digits = 3)
        pct.2 = round( Matrix::rowSums(data[deg4idx$feature, cells.2] > 0)/length(cells.2), digits = 3)
        data.alpha <- cbind(pct.1, pct.2)
        colnames(x = data.alpha) <- c("pct.1", "pct.2")
        xe = cbind(deg4idx, data.alpha)
        alpha.max <- apply( data.alpha, MARGIN = 1, FUN = max)
        names(x = alpha.max) <- rownames(data.alpha)
        # alpha.diff <- alpha.max - apply(X = data.alpha, MARGIN = 1, FUN = min)
        features <- names( which(alpha.max > 0.25 ) )
        if (length(x = features) == 0) { stop("No features pass min.diff.pct threshold") }
        alpha = xe %>% filter( feature %in% features)
    })
    global_DEGs =  do.call(rbind, results)
    global_DEGs = global_DEGs %>% filter( logFC > 0 ) %>%
                    select( feature, group, logFC, pval,padj,pct.1, pct.2, auc) %>%
                    rename(gene = feature, cluster = group,avg_logFC = logFC, p_val_adj = padj, p_val = pval )
} else if ( tolower(dftest) == "venice" ){
    suppressPackageStartupMessages( library("Venice") )
    global_DEGs = VeniceAllMarkers(object = seurat_ob,only.pos = T, pvalue = +Inf, logfc.threshold = 0)
    data = GetAssayData(seurat_ob, slot = "counts")
    results = future_lapply(unique(Idents(seurat_ob)), function(idx){
        cells.1 = Cells(subset(seurat_ob, idents = idx))
        cells.2 = Cells(subset(seurat_ob, idents = idx, invert = T))
        deg4idx = global_DEGs %>% filter(cluster == idx)
        pct.1 = round( Matrix::rowSums(data[deg4idx$Gene.Name, cells.1] > 0)/length(cells.1), digits = 3)
        pct.2 = round( Matrix::rowSums(data[deg4idx$Gene.Name, cells.2] > 0)/length(cells.2), digits = 3)
        data.alpha <- cbind(pct.1, pct.2)
        colnames(x = data.alpha) <- c("pct.1", "pct.2")
        xe = cbind(deg4idx, data.alpha)
        alpha.max <- apply( data.alpha, MARGIN = 1, FUN = max)
        names(x = alpha.max) <- rownames(data.alpha)
        # alpha.diff <- alpha.max - apply(X = data.alpha, MARGIN = 1, FUN = min)
        features <- names( which(alpha.max > 0.25 ) )
        if (length(x = features) == 0) { stop("No features pass min.diff.pct threshold") }
        alpha = xe %>% filter( Gene.Name %in% features)
    })
    global_DEGs =  do.call(rbind, results)
    global_DEGs$Log10.adjusted.p.value = 10^global_DEGs$Log10.adjusted.p.value
    global_DEGs$Log10.p.value = 10^global_DEGs$Log10.p.value
    global_DEGs = global_DEGs %>% rename(gene = Gene.Name,
                        avg_logFC = Log2.fold.change,
                        p_val_adj = Log10.adjusted.p.value,
                        p_val = Log10.p.value )
}else{
    global_DEGs = FindAllMarkers(object = seurat_ob,only.pos = T,test.use = dftest, logfc.threshold = 0, min.pct = 0.25)
}

global_DEGs = global_DEGs %>% mutate( gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>% select( gene, everything())
write.table(global_DEGs,file = file.path(output_dir,"all_markers_for_each_cluster.xls"),
    col.names =T,row.names = F,sep = "\t",quote=F)
#find the significant differential expressed genes for each clusters against all other clusters
#feature plot of potential marker gene for each cluster
#for each gene to be a potential marker,in add to be a significant expressed gene, the gene should account for 
#large proportion in the interested cluster but as small as possiable in the other clusters, that means the pct.1 should
# be bigger than the pct.2.
if ( opt$strict == T ){
    if ( !is.null(opt$pvalue) ){
        topn_markers  = global_DEGs %>% group_by(cluster) %>% 
		      filter(avg_logFC>=opt$avg_logFC & p_val < opt$pvalue&pct.1 > pct1_cutoff & pct.2 < pct2_cutoff & gene_diff > pct_fold_cutoff)  %>%
		      arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
		      filter(gene_diff > pct_fold_cutoff)  %>% 
		      top_n(topn,gene_diff)
    }else{
        topn_markers  = global_DEGs %>% group_by(cluster) %>% 
		      filter(avg_logFC>=opt$avg_logFC &p_val_adj < opt$FDR & pct.1 > pct1_cutoff & pct.2 < pct2_cutoff & gene_diff > pct_fold_cutoff)  %>%
		      arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
		      filter(gene_diff > pct_fold_cutoff)  %>% 
		      top_n(topn,gene_diff)
    }
}else{
  topn_markers  = global_DEGs %>% group_by(cluster) %>% 
              #filter(p_val < opt$pvalue ) %>%
              arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
              # filter(gene_diff > pct_fold_cutoff)  %>% 
              top_n(topn,gene_diff)
}
write.table(topn_markers,file = file.path(output_dir,paste0("top", topn, "_markers_for_each_cluster.xls", collapse = "")),
            col.names =T,row.names = F,sep = "\t",quote=F)
