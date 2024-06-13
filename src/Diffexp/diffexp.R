#!/usr/bin/env Rscript
#this script is used to find all the differential expressed gene for single cell UMI counts
#and pick the differential expressed  genes with specified restrictions.
#The differential expression analysis is somewhat different from the
#the marker gene analysis in that the marker genes are only genes up-regulated.
rm(list=ls())
##########################################################################################
##########function definition ####################
##########################################################################################
##differential gene expression ananlysis using DESeq2#####
diff_gene_analysis <- function( 
  object,
  features,
  defunction,
  contrast,
  fdr,
  fc_threshold,
  pvalue_threshold,
  outputdir,
  ...
){
    #get the differential expressed genes and the intermediated result if necessary
    contrasts = unlist( strsplit(contrast,":",perl = T) )
    numerator_subset = unlist( strsplit(contrasts[2],",",perl = T) )
    denominator_subset =unlist( strsplit(contrasts[3],",",perl = T) )
    numerator = OldWhichCells( object,  subset.name= contrasts[1], accept.value = numerator_subset) #cells used as interested cells
    denominator = OldWhichCells( object,  subset.name= contrasts[1], accept.value = denominator_subset) #cells used as interested cells
    dftest4seurat = c( "wilcox", "bimod", "roc", "t", "tobit", "poisson", "negbinom", "MAST")
    zinbwave_based = c( "DESeq2", "zinbwave-edgeR" )
    if ( defunction %in% dftest4seurat ){
        if ( defunction == "MAST" ){
            suppressPackageStartupMessages( library("MAST") )
        }
        res = Seurat::FindMarkers(object, ident.1 = numerator_subset,
                          ident.2 = denominator_subset,
                          group.by = contrasts[1],
                          test.use = defunction,
                          only.pos = F )
        res = res %>% tibble::rownames_to_column(var = "gene") %>%
              dplyr::rename( pvalue = p_val, padj = p_val_adj)
    }else if ( defunction %in% zinbwave_based ){
        if ( defunction == "DESeq2"){
            suppressPackageStartupMessages( library("BiocParallel") )
            res = Seurat::FindMarkers(object, ident.1 = numerator_subset,
                                    ident.2 = denominator_subset,
                                    group.by = contrasts[1],
                                    test.use = defunction,
                                    only.pos = F,
                                    BPPARAM = MulticoreParam(detectCores()) )
            res = res %>% tibble::rownames_to_column(var = "gene") %>%
                  dplyr::rename( pvalue = p_val, padj = p_val_adj)
            # res = diffExp(
            #             object = object,
            #             numerator = numerator,
            #             denominator = denominator,
            #             caller = "DESeq2" )
            # res = as.data.frame(res)
        }else{
            suppressPackageStartupMessages( library("pointillism") )
            object = as.SingleCellExperiment( object )
            # object = as(object, "SingleCellExperiment")
            suppressPackageStartupMessages( library("edgeR") )
            res_edgeR = diffExp(
                              object = object,
                              numerator = numerator,
                              denominator = denominator,
                                  caller = "edgeR" )
            res <- topTags(res_edgeR, n=nrow(res_edgeR))
            res = as.data.frame(res)
            res = res %>% dplyr::rename( pvalue = PValue, padj = FDR )
        }
    }else if ( defunction == "scde" ){
        suppressPackageStartupMessages( library("scde") )
        groups = slot(object,name = "meta.data")[c(numerator,denominator),contrast[1]]
        res = RunSCDE( object = object,
                       groups = groups,
                        ident.1 = numerator,
                       ident.2 = denominator )
    }else if ( defunction == "venice"){
        suppressPackageStartupMessages( library("Venice") )
        res = VeniceFindMarkers(object, ident.1 = numerator_subset ,
                                ident.2 = denominator_subset ,
                                group.by = contrasts[1],
                                logfc.threshold = -Inf)
        res$Log10.adjusted.p.value = 10^res$Log10.adjusted.p.value
        res$Log10.p.value = 10^res$Log10.p.value
        res = res %>% rename(gene = Gene.Name,
                            avg_logFC = Log2.fold.change,
                            padj = Log10.adjusted.p.value,
                            pvalue = Log10.p.value )
    }

    if ( nrow(res) == 0 ){
        warning(paste0("NO Differential Genes Identified for ", contrast,sep =""))
        return(0)
    }
    numerator_means = Matrix::rowMeans(GetAssayData(object,slot = "data")[res$gene,numerator])
    denominator_means = Matrix::rowMeans(GetAssayData(object,slot = "data")[res$gene,denominator])
    # res$M <- log2( exp(numerator_means) / exp(denominator_means) )
    res$baseMean <- 1/2 * (log2(numerator_means) + log2(denominator_means))
    
    res = res %>% dplyr::select(gene,everything())  #put the Gene ID column to the first
    # the seurat in default use e the log function base, so we need to change it back
    # to 2-based log function
    res$FoldChange = exp(1)^res$avg_logFC #add the FoldChange column
    res = res %>% select(-avg_logFC)
    res$log2FoldChange = log2(res$FoldChange)
    base_prefix = paste0(contrasts[1],"_",contrasts[2],"-vs-",contrasts[3])
    write.table(res, file.path(outputdir,paste0(base_prefix,"-all_diffexp_genes",".xls")),
                sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="")
    #filter the genes by the threshold specified by the user
    if(! is.null(pvalue_threshold)){
        res_Significant = dplyr::filter(res, pvalue < pvalue_threshold , abs(log2FoldChange) > log2(fc_threshold))
    }
    if(! is.null(fdr)){
        res_Significant = dplyr::filter(res, padj < fdr , abs(log2FoldChange) > log2(fc_threshold))
    }
    if ( nrow(res_Significant) == 0 ){
        warning(paste0("NO Significant Differential Genes Identified for ", contrast,sep =""))
        return(0)
    }
    res_Significant[which(res_Significant$log2FoldChange  > 0), "up_down"] <- "Up"
    res_Significant[which(res_Significant$log2FoldChange  < 0), "up_down"] <- "Down"
    #write out the desired differential gene expression analysis results
    if(! is.null(pvalue_threshold)){
        write.table(res_Significant,
                    file.path(outputdir,
                    paste0(base_prefix,"-diff-","pval-",pvalue_threshold,"-FC-",fc_threshold,".xls")),
                    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="")
        stat<-matrix(c("Case","Control","Up_diff","Down_diff",
        paste("Total_diff(","pvalue<",pvalue_threshold,"&FoldChange>", fc_threshold,")",sep="")),
        ncol=5)
    }
    if(! is.null(fdr)){
        write.table(res_Significant,
                    file.path(outputdir,
                    paste0(base_prefix,"-diff-","padj-",fdr,"-FC-",fc_threshold,".xls")),
                    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="")
        stat<-matrix(c("Case","Control","Up_diff","Down_diff",
        paste("Total_diff(","padj<",fdr,"&FoldChange>",fc_threshold,")",sep="")),
                    ncol=5)
    }
    up_num<-length(which(res_Significant$up_down=="Up"))
    down_num<-length(which(res_Significant$up_down=="Down"))
    total<-sum(up_num,down_num)
    #need to change here case_name and control name is the level name,the
    #case and control includes all the sample id in each level
    stat1<-c(paste(contrasts[2],"(",contrasts[1],")",sep=""),paste(contrasts[3],"(",contrasts[1],")",sep=""),up_num,down_num,total)
    if(file.exists(file.path(outputdir,"diffexp_results_stat.xls"))){
        stat<-rbind(stat1)
    }else{
        stat<-rbind(stat,stat1)
    }
    write.table(stat,file.path(outputdir,"diffexp_results_stat.xls"),
                append=file.exists(file.path(outputdir,"diffexp_results_stat.xls")),
                quote=F,row.names=F,col.names=F,sep="\t",na="")
}

#=================================================================================
#mian workflow
#=================================================================================
suppressPackageStartupMessages(library("future"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages( library("Matrix") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages( library("dplyr") )

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character", 
                 help = "The input exprssion matrix in seurat object format."),
    make_option( c("--design", "-d"), type = "character", default = NULL,
                 help = "The group design for cell clusters or samples to make differential expression analysis."),
    make_option( c("--addition_metadata", "-M" ), type="character", default = NULL,
                help="[Optional]additional metadata for each sample which includes sample id and additional sample groupping info.",
                metavar="character"),
    make_option( c("--which_group", "-v" ), type = "character",
                help = "[Optional]select the groupping column in metadata used to do subsettig of cell if necessary."
                ),
    make_option( c("--WHICH_CELLS", "-u" ), type = "character",
                help = "[Optional]the level id list in selected cell groupping used to do differential expression analysis."
                ),
    make_option( c("-c","--contrast"),type = "character",default = NULL,
                 help = "[Optional]levels of a factor used to compare with for final differenetial results.
                        The format is Factor:interesting_levle:reference_level."),
    make_option( c("--FC","-k"), type = "double", default = 1,
                 help = "The average FC of the gene UMI count in its cluster against the all other clusters."),
    make_option( c("--pvalue","-p"), type = "double", default = 0.05,
                 help = "the P-value of the gene differential expression.",metavar = "P-value"),
    make_option( c("--fdr","-q"), type = "double", default = NULL,
                 help = "the FDR of the gene differential expression.",metavar = "FDR"),
    make_option( c("--ncores","-j"),type="integer", default = "4",
                 help="the number of CPUs used to parallize this job."),
    make_option( c("--assay"), type = "character", default = "RNA",
                help = "[OPTIONAL]the default assay to use  for the this run.[default:RNA]"),
    make_option( c("--output","-o"),type="character", default = "./",
                 help="the output directory of QC results.", metavar="character"),
    make_option( c("--DEGtest","-e"), type = "character", default = "wilcox",
                 help = "the test methods used to find differential expressed genes.Options are:wilcox,roc,t,bimod,tobit,poission
                        negbinom,MAST,DESeq2,zinbwave-DESeq2,zinbwave-edgeR."));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# "wilcox" : Wilcoxon rank sum test (default).
# "roc" : Standard AUC classifier.
# "t" : Student\'s t-test.
# "bimod" : Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013).
# "tobit" : Tobit-test for differential gene expression (Trapnell et al., Nature Biotech, 2014).
# "poisson" : Likelihood ratio test assuming an underlying poisson distribution. Use only for UMI-based datasets.
# "negbinom" : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets.
# "MAST" : GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015).
# "DESeq2" : DE based on a model using the negative binomial distribution (Love et al, Genome Biology, 2014).
# "zinbwave-DESeq2": Perform pairwise differential expression across groups of cells by fitting to a
#          DESeq2-based zero-inflated negative binomial (ZINB) model using the zinbwave package.
# "zinbwave-edgeR": Perform pairwise differential expression across groups of cells by fitting to a
#          edgeR-based zero-inflated negative binomial (ZINB) model using the zinbwave package. '

#=================================================================================
#parse the command line parameters
#=================================================================================
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
  print("NO differential expression test method supplied,using the DESeq2as default,")
  dftest = "zinbwave_DESeq2"
}else{
  dftest = opt$DEGtest
}

#parse the design
if ( is.null( opt$design ) ){
    print_help(opt_parser)
    stop("NO assay design is PROVIDED\n")
}else{
    design = opt$design
}

if ( is.null(opt$pvalue) && is.null(opt$fdr) ){
    print_help(opt_parser)
    stop("None of P-value or FDR is AVAILABLE! Filtering can be performed using any one of (-p), (-f) at one time.", call.=FALSE)
}else if ( !is.null(opt$fdr)){
    fdr = as.numeric(opt$fdr)
    pvalue = NULL
}else{
    pvalue = as.numeric(opt$pvalue)
    fdr = NULL
}

#read in the clustering results in seurat object format
if ( is.null(opt$input) ){
    stop("NO seurat object is FOUND!")
}else{
    seurat_ob = readRDS( opt$input )
    if ( seurat_ob@version < 3 ){
        seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
    }
    DefaultAssay(seurat_ob) = opt$assay
}

#check the metadata of the assay design,which describe the expriement groupping of each
#sample in the assay
if ( !is.null( opt$addition_metadata)  ){ #the additional metadata for each sample
    add_assay_metadata = read.table(opt$addition_metadata,sep=",",header =T )
    # cellnames = seurat_ob@cell.names
    cellnames = colnames(seurat_ob) #seurat v3 style
    sampleidx =  gsub("_.*","",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
    #integrate the additional metadata from the assay design
    additional_cell_meta = vector( )
    for ( colidx in colnames(add_assay_metadata) ){
        additional_cell_meta = cbind(additional_cell_meta, as.vector(add_assay_metadata[sampleidx, colidx]))
    }
    colnames(additional_cell_meta) = colnames(assay_metadata)
    rownames(additional_cell_meta) = cellnames
    additional_cell_meta = as.data.frame(additional_cell_meta)
    seurat_ob = AddMetaData( seurat_ob, additional_cell_meta)

    # assay_metadata = seurat_ob@meta.data
    assay_metadata = seurat_ob[[]] # seurat v3 style
}else{
    assay_metadata = seurat_ob[[]] # seurat v3 style
}

if ( is.null(assay_metadata$clusters) ){
    seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
}else{
    seurat_ob = SetIdent( seurat_ob, value = "clusters")
}

##subset the cell for differential expression analysis
if ( !is.null(opt$WHICH_CELLS)){
    cluster_list = unlist(strsplit( opt$WHICH_CELLS,",",perl = T))
    seurat_ob = SubsetData( seurat_ob, subset.name = opt$which_group, accept.value = cluster_list)
}


#=================================================================================
#parse the contrast for this assay
#=================================================================================
#if the constrast is not specified by the user explicitly from the command line,
#the final differential result will use last level of the last variable in the
#design formula over the first level of this variable. The levels for each
#factor is aphabetly ordered by default,which can be reordered by the user.
#if available,the contrast string from the user must obey the right format:
#the_interested_factor:the_interested_levels_in_this_factor:the_reference_level_in_this_factor
if ( is.null(opt$contrast ) ){ #no contrast is provided
    factors_indesign = strsplit(opt$design,"[~+ ]+",perl = T)
    last_factor_indesign = factors_indesign[length(factors_indesign)]
    if ( is.null(assay_metadata[,last_factor_indesign]) ){
        stop("The factor in design formula does not exist in assay metadata.")
    }
    variable_levels = levels(assay_metadata[,last_factor_indesign])
    contrast = paste( last_factor_indesign,variable_levels[length(variable_levels)],variable_levels[1],sep = ":" )
}else{
    contrast = opt$contrast
}
contrasts = unlist( strsplit(contrast,":",perl = T) )
all_levels = as.vector(unique(assay_metadata[,contrasts[1]]))

if ( contrasts[2] == "all" & contrasts[3] != "all" ){
    all_levels = all_levels[-which(all_levels==contrasts[3])] #delete the reference level
    all_comparisions = paste(contrasts[1],all_levels,contrasts[3],sep = ":")
}else if( contrasts[2] == "all" & contrasts[3] == "all" ){

    # combine_of2 = combn(all_levels,2) #random combination of two group
    # all_comparisions = c( paste(contrasts[1],combine_of2[2,],combine_of2[1,],sep = ":"),
    #                       paste(contrasts[1],combine_of2[1,],combine_of2[2,],sep = ":"))
    all_comparisions = lapply(all_levels,
                    function(x) paste(contrasts[1],x,paste0(all_levels[-which(all_levels==x)],collapse = ","),sep = ":"))
    all_comparisions = unlist(all_comparisions)
}else if ( contrasts[2] != "all" & contrasts[3] == "all" ){
    #delete the interested level in the  reference level
    ref_levels = paste0(all_levels[-which(all_levels==contrasts[2])],collapse = ",")
    all_comparisions = paste(contrasts[1],contrasts[2],ref_levels,sep = ":")
}else if (contrasts[2] == "each" & contrasts[3] != "each"){
    all_levels = all_levels[-which(all_levels==contrasts[3])]


}else{
    all_comparisions = contrast
}


# Differential expression analysis
#parse the contrast string to a list for later use
# suppressPackageStartupMessages(library("doParallel"))
# suppressPackageStartupMessages(library("foreach"))
#
# registerDoParallel(cores=length(all_comparisions))
# foreach( contrast = all_comparisions )%dopar% diff_gene_analysis( object = seurat_ob,defunction = dftest,fdr = fdr, fc_threshold = opt$FC,pvalue_threshold = pvalue,contrast = contrast, outputdir = output_dir)
# doParallel::stopImplicitCluster()
#Differential expression analysis
#parse the contrast string to a list for later use

# setting the cores for parallization
options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
plan("multicore", workers = min(availableCores(), opt$ncores)) # parallization using specified CPUs start from here
future_lapply( all_comparisions, function( contrast){
    diff_gene_analysis( object = seurat_ob,
            defunction = dftest,fdr = fdr,
            fc_threshold = opt$FC,
            pvalue_threshold = pvalue,
            contrast = contrast,
            outputdir = output_dir)
})

