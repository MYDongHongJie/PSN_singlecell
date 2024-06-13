#!/usr/bin/env Rscript
# Created by: mhan
# Created on: 18-8-29
# Updated on: 19-6-25
# this script is used to visualize the gene expression in different plot sepecified by the user.
# the supported visualization types can be as follows:
# tSNE plot: as featureplot on the command line
# violin plot: as marker_gene_violin_plot
# heatmap:
# featuregroup
#
rm(list=ls())
autotranslate= function (genes, targetGenes, possibleOrigins = NULL, possibleTargets = NULL, 
    returnAllPossible = FALSE, db = homologene::homologeneData) 
{
    pairwise = db$Taxonomy %>% unique %>% utils::combn(2) %>% 
        {
            cbind(., .[c(2, 1), ], rbind(db$Taxonomy %>% unique, 
                db$Taxonomy %>% unique))
        }
    if (!is.null(possibleOrigins)) {
        possibleOrigins[possibleOrigins == "human"] = 9606
        possibleOrigins[possibleOrigins == "mouse"] = 10090
        pairwise = pairwise[, pairwise[1, ] %in% possibleOrigins, 
            drop = FALSE]
    }    else {
        possibleOrigins = db$Taxonomy %>% unique
    }
    if (!is.null(possibleTargets)) {
        possibleTargets[possibleTargets == "human"] = 9606
        possibleTargets[possibleTargets == "mouse"] = 10090
        pairwise = pairwise[, pairwise[2, ] %in% possibleTargets, 
            drop = FALSE]
    }    else {
        possibleTargets = db$Taxonomy %>% unique
    }
    possibleOriginData = db %>% dplyr::filter(Taxonomy %in% possibleOrigins & 
        (Gene.Symbol %in% genes | Gene.ID %in% genes)) %>% dplyr::group_by(Taxonomy)
    possibleOriginCounts = possibleOriginData %>% dplyr::summarise(n = dplyr::n())
    possibleTargetData = db %>% dplyr::filter(Taxonomy %in% possibleTargets & 
        (Gene.Symbol %in% targetGenes | Gene.ID %in% targetGenes)) %>% 
        dplyr::group_by(Taxonomy)
    possibleTargetCounts = possibleTargetData %>% dplyr::summarise(n = dplyr::n())
    pairwise = pairwise[, pairwise[1, ] %in% possibleOriginCounts$Taxonomy, 
        drop = FALSE]
    pairwise = pairwise[, pairwise[2, ] %in% possibleTargetCounts$Taxonomy, 
        drop = FALSE]
   # genes=intersect(genes,possibleTargetData$Gene.Symbol)
    possibleTranslations <- pairwise %>% apply(2, function(taxes) {
        homologene(genes, inTax = taxes[1], outTax = taxes[2])
    }) %>% {
        .[purrr::map_int(., nrow) > 0]
    }
    translationCounts <- possibleTranslations %>% sapply(function(trans) {
        sum(c(trans[, 2], trans[, 4]) %in% targetGenes)
    })
    if (any( translationCounts > 0) ){
        if (!returnAllPossible) {
            possibleTranslations <- translationCounts %>% which.max %>% 
                {
                    possibleTranslations[[.]]
                }
            if (sum(translationCounts > 0) > 1) {
                bestMatch = translationCounts %>% which.max
                nextBest = max(translationCounts[-bestMatch])
                warning("There are other pairings, best of which has ", 
                    nextBest, " matching genes")
            }
        }    else {
            possibleTranslations = possibleTranslations[translationCounts != 
                0]
        }
    }
    return(possibleTranslations)
}

casematch = function (search, match=rownames(seurat_ob),possibleOrigins=c(9606,10090),possibleTargets)
{
    search.match <- sapply(X = search, FUN = function(s) {
        s=sub(" ","",s)
        res= grep(pattern = paste0("^", s, "$"), x = match,
            ignore.case = TRUE, perl = TRUE, value = TRUE)
        if (length(res) != 0) {
            return(res)
        } else {
            res= try(autotranslate(genes=s, targetGenes=match, possibleOrigins = possibleOrigins , returnAllPossible = FALSE, db = homologene::homologeneData)[,2],TRUE)
            if (class(res)!="try-error") {
                return(res)
            } else {
                return(character())
            }
        }
    })
    return(unlist(x = search.match))
}



#=================================================================================
# customized function definition
#=================================================================================

suppressWarnings({
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("optparse") )
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("gridExtra"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("tibble"))
    suppressPackageStartupMessages(library("RColorBrewer"))
    suppressPackageStartupMessages(library("OESingleCell"))
})

#=command line parameters setting=============================
option_list = list(
    make_option( c("--RDS", "-v"), type = "character", default = "TRUE",
                 help = "the seurat object saved as R object in RDS format."),
    make_option( c("--markers","-l"), type ="character",
                help="the list file of marker genes to be visulized."),
    make_option( c("--output","-o"),type="character", default = "./",
                help="the output directory of results.", metavar="character"),
    make_option( c("--topn", "-n"), type="integer", default = 25,
                 help = "the number of top markers for each cluster to visualizse."),
    make_option( c("--topby", "-c"), type = "character", default = "gene_diff",
                 help="the column used to pick top n marker gene to visulize.The
                 option can be one of the column in the input marker genes table."),
    make_option( c("--extraGene", "-x"), type = "character", default = NULL,
                 help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user."),
    make_option( c("--diffGene", "-d"), type = "character", default = NULL,
                 help = "[OPTIONAL]The diff gene list of interest to visualize specified by the user."),
    make_option( c("--groupby", "-g"), type = "character", default = NULL,
                help = "[OPTIONAL]The grouppinig variable in the metadata for
                        separate the cells to visulize marker genes."),
    make_option( c("--splitby", "-y"), type = "character", default = NULL,
                help = "[OPTIONAL]The variable in the metadata used to split the graph by the variable levels to
                        comparing the gene expression difference in different levels."),
    make_option( c("--pointsize", "-s"), type = "double", default = NULL,
        help = "[OPTIONAL]the point size in the plot."),
    make_option( c("--sample_ratio"), type = "double", default = 0.6,
        help = "[OPTIONAL]the ratio of random subsample for each group when drawing heatmap."),
    make_option( c("--reduct"), type = "character", default = "tsne",
        help = "[OPTIONAL]the previous calculated reduction result used in the featureplot,."),
    make_option( c("--alpha2use", "-a"), type = "double", default = 0,
        help = "[OPTIONAL]the opacity of the pooints on the violin plot."),
    make_option( c("--slot"), type = "character", default = "scale.data",
        help = "[OPTIONAL]the slot in the assay to use."),
    make_option( c("--assay", "-e"), type = "character", default = "RNA",
        help = "[OPTIONAL]the array result to use in case of multimodal analysis."),
    make_option( c("--vismethod","-m"), type= "character",default="vlnplot,featureplot",
                 help = "the visulization methods for the marker genes of each cell cluster.
                 he methods can be ridgeplot,vlnplot,flip_vlnplot,splitby_featureplot,dotplot,featureplot,heatmap,featurebygroup,diff_heatmap,geneset,corscatter,ggstatsplot"),
    make_option( c("--var2use", "-q" ), type = "character", default = NULL,
                help = "[OPTIONAL]The column name in cell metadata used as identity
                        of each cell combined with levels4var."),
    make_option(c("-w", "--order"),type = "character", default = NULL,
                help="Whether to display the Ordering of group"),
    make_option( c("--levels4var", "-u" ), type = "character", default = NULL,
                help = "[OPTIONAL] subset of factor levels for the specified factor by --var2use."),
    make_option( c("--colors" ), type = "character", default = "Spectral",
                help = "[OPTIONAL] DOTPLOT picture color selection"),
    make_option( c("--pvalue", "-p"),type= "character",default = NULL,
                help = "[OPTIONAL]use like GROUP1:GROUP2+GROUP2:GROUP3 to add pvalue on vlnplot or boxplot."),
    make_option( c("--scoredata"),type= "character",default = NULL,
                help = "[OPTIONAL] the file of geneset.addmodeulescore_plot.xls.")
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

    metadata = seurat_ob@meta.data
    if ( is.null(metadata$clusters) ){
        seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
    }else{
        # if it is the first time to run this script on a seurat object, the
        # clusters here is actually the sample index but cell cluster ids.
        # After running this script, the cluster id will be overwrite with
        # the actual cell cluster id.
        seurat_ob = SetIdent( seurat_ob, value = "clusters")
    }
	seurat_ob@meta.data$clusters = factor(seurat_ob@meta.data$clusters, levels = sort( unique(as.numeric(as.vector( seurat_ob@meta.data$clusters)))))
}

if ( is.null( opt$groupby ) ){
    print( "NO groupping variable AVAILABLE for cell groupping! The default cell clusters id will be used!")
    groupby = "clusters"
}else{
    groupby = opt$groupby
}

if ( !is.null( opt$order )){
    print("Useing customize factor from ordering")
    group_order = as.character(unlist(strsplit( opt$order,",",perl = T)))
    seurat_ob@meta.data[[groupby]] = factor(seurat_ob@meta.data[[groupby]],levels = group_order)
} else {
    print(paste0("Useing default factor for ",groupby ))
}

if ( !is.null( opt$splitby ) ){
    splitby = opt$splitby
    #facetbyx = unlist(strsplit(splitby,",",perl=T))
    if ( length(table(seurat_ob@meta.data[,splitby])) >2 ){
        print( "如果组数大于2 , 不建议使用splitby模式绘制小提琴图！")
    }
    if ( splitby == groupby ){
        stop( "The variable specified by --splitby conflicted with the --groupby parameter, NULL will be used! ")
        splitby = NULL
    }
}else{
    splitby = NULL
}

#get the subset of cells used for visualization if necessay
if ( !is.null(opt$levels4var)){
    if ( is.null(opt$var2use ) ){
        print("NO cell identity column name AVAILABLE! The groupby parameter will be used as default.")
        ident2use = groupby
    }else{
        ident2use = opt$var2use
    }
    cluster_list = unlist(strsplit( opt$levels4var,",",perl = T))
    seurat_ob = SubsetData( seurat_ob, subset.name = ident2use, accept.value = cluster_list)
    seurat_ob@meta.data[[ident2use]]=factor(seurat_ob@meta.data[[ident2use]],levels = sort(unique(seurat_ob@meta.data[[ident2use]])))
}

# output directory setting
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    root_dir = getwd()
}else{
    if ( file.exists(opt$output)){
        root_dir = opt$output
    }else{
        root_dir = opt$output
        dir.create(root_dir, recursive = T)
    }
}



if (is.null( opt$markers ) & is.null(opt$extraGene) & is.null(opt$diffGene) & is.null(opt$scoredata) ){
     stop("NO marker genes is AVAILABLE!")
}    

topn_markers = data.frame()
if ( !is.null(opt$markers) ){
    markers2vis = read.table(  opt$markers, sep="\t", header = T,quote="")
    markers2vis[["cluster"]] = factor(markers2vis[["cluster"]] , levels = sort(unique(seurat_ob@meta.data[[groupby]])))
    topn_markers  = markers2vis %>% group_by(cluster) %>%
        arrange(desc(gene_diff)) %>%
        top_n(opt$topn,.data[[opt$topby]]) %>% arrange(cluster) %>% mutate(folder_suffix = paste0("cluster",cluster)) %>% select(cluster,gene,folder_suffix)
}
if ( !is.null(opt$extraGene) ){
    if (sub(".*\\.","",opt$extraGene,perl=T)=="xlsx" ){
        library(readxl);extra_gene = read_xlsx(opt$extraGene)
    } else extra_gene = read.table(opt$extraGene, sep="\t", header = T)
    if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene"){colnames(extra_gene)[1]="extra"}
    formated_extra_gene = as.data.frame(tidyr::gather(extra_gene,key = "cluster",value = "GENE"))
    match = casematch(search = as.vector(formated_extra_gene$GENE),match = rownames(seurat_ob))
    filtered_gene = formated_extra_gene$GENE[!formated_extra_gene$GENE %in% names(match )& formated_extra_gene$GENE != ""]
    if(length(filtered_gene)!=0){
        filtered_gene = as.data.frame(filtered_gene)
        colnames(filtered_gene) = "Gene"
        write.table(filtered_gene,file.path(root_dir,"genes_not_matched.xls"),quote = F,row.names=F)
        print("There are some mismatched gene symbol, Please check genes_not_matched.xls for the genename.")}
    formated_extra_gene = formated_extra_gene %>% dplyr::filter(GENE %in% names(match)) %>% dplyr::mutate(gene = match,folder_suffix = cluster) %>% dplyr::select(cluster, gene,folder_suffix)
    topn_markers = rbind(topn_markers, formated_extra_gene)
}
if ( !is.null(opt$diffGene) ){
    diffGene_path = normalizePath(opt$diffGene)
    diffGene_name = gsub(".*/","",diffGene_path)
    if(grepl("-diff-",diffGene_name)){
        name = gsub("-diff-.*","",diffGene_name)
    }else if(grepl("-all_",diffGene_name)){
        name = gsub("-all_.*","",diffGene_name)
    }

    markers2vis = read.delim( opt$diffGene, sep="\t", header = T,quote="")
    markers2vis = subset(markers2vis,!grepl("^(mt-|Rps|Rpl|MT-|RPS|RPL)",markers2vis$gene))
    up = filter(markers2vis,FoldChange > 1) %>% arrange(desc(log2FoldChange ))  %>% top_n(opt$topn,log2FoldChange )
    down = filter(markers2vis,FoldChange < 1) %>% arrange(log2FoldChange )  %>% top_n(as.numeric(paste0("-",opt$topn)),log2FoldChange )
    topn_markers = rbind(up, down)
    write.table(topn_markers,file.path(root_dir,paste0("top", opt$topn, "_", name, "_genes.xls", collapse = "")),quote = F,row.names = FALSE,sep="\t")
}

if ( is.null(opt$vismethod) ){
    print("NO marker gene visulization method provided,the default method vlnplot and featureplot will be used!")
    vismethods = c("vlnplot","featureplot")
}else if( opt$vismethod == "all" ){
    vismethods = c("vlnplot","featureplot","flip_vlnplot","splitby_featureplot","ridgeplot","dotplot","heatmap", "violinEnsemble")
}else{
    vismethods = unlist(strsplit(opt$vismethod,","))
}

if ( !is.null(opt$colors) ){
    colors = opt$colors
}else{
    colors = c("Spectral")
}
if ( !is.null(opt$pvalue) ){
    pvalue = opt$pvalue
}else{
    pvalue = NULL
}

# 兼容老版本脚本未存储rawbc情况
if (!"rawbc" %in% colnames(seurat_ob@meta.data)) { seurat_ob[["rawbc"]] <- seurat_ob[["orig.ident"]] }

# Palette define
Palette = colorRampPalette(rev(brewer.pal(11, colors)))
sc <- scale_colour_gradientn(colours = Palette(100))

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

if ( is.null(opt$pointsize) ){
    if (dim(seurat_ob)[2] < 500){
        pointsize = 1.5
    } else pointsize = 1
} else {
    pointsize = opt$pointsize
}



#=================================================================================
# visualize the markers in different ways
#=================================================================================
for ( vismethod in vismethods ){
    if ( vismethod == "flip_vlnplot" ){
        .libPaths("/home/ziqingzhen/R/x86_64-conda_cos6-linux-gnu-library/3.6")
        suppressPackageStartupMessages(library("patchwork"))

        for ( clusterx in unique(topn_markers$folder_suffix) ){
            topn = topn_markers %>% filter( folder_suffix == clusterx) %>% select(cluster,gene,folder_suffix)

            path4vis = file.path(root_dir, paste0("markers_vis4",clusterx,collapse = ""))
            if ( file.exists( path4vis ) ){
            output_dir = path4vis
            }else{
            output_dir = path4vis
            dir.create(output_dir, recursive = TRUE)
            }

            modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
            p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab("") + ggtitle(feature) +
               coord_flip() +
               theme(legend.position = "none",
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
            #    axis.ticks.x = element_line(),
                axis.text.x = element_text(size = 8,angle = 0, hjust = 0.5),
               axis.title.x = element_text(size = 12, angle = 0 ),
               plot.margin = plot.margin )
            return(p)
            }

            StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
            plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[1]]  <- plot_list[[1]] + 
            xlab(groupby) + 
            theme(axis.text.y=element_text(size = 12), axis.ticks.y = element_line() )
            p <- patchwork::wrap_plots(plotlist = plot_list, ncol = length(plot_list) )
            return(p)
            }
            topn_markers2vis = as.vector(topn$gene)
            colors2use = CustomCol2(1:length(unique(seurat_ob@meta.data[,groupby])))
            flip <- StackedVlnPlot(seurat_ob, topn_markers2vis, pt.size=0, cols=colors2use, group.by = groupby)     
            pdf(file.path(output_dir,paste0("top", "marker_gene_flip_vlnplot_plot",".pdf", collapse = "_")),
                width = length(topn_markers2vis)*2 + max(nchar(as.character(seurat_ob@meta.data[,groupby])))/4 , height = 8)
            print(flip)
            dev.off()
            png(file.path(output_dir,paste0("top","marker_gene_flip_vlnplot_plot",".png",collapse = "_")),
                width = length(topn_markers2vis)*2+max(nchar(as.character(seurat_ob@meta.data[,groupby])))/4, height = 8 , res = 96, units = "in")
            print(flip)
            dev.off()
        }
    }

    if ( vismethod == "vlnplot" ){
        # Draws a violin plot of single cell data (gene expression, metrics, PC scores, etc.)

        for ( clusterx in unique(topn_markers$folder_suffix) ){
            topn = topn_markers %>% filter( folder_suffix == clusterx) %>% select(cluster,gene,folder_suffix)
            topn_markers2vis = as.vector(topn$gene)

            path4vis = file.path(root_dir, paste0("markers_vis4",clusterx,collapse = ""))
            if ( file.exists( path4vis ) ){
                output_dir = path4vis
            }else{
                output_dir = path4vis
                dir.create(output_dir, recursive = TRUE)
            }

            colors2use = CustomCol2(1:length(unique(seurat_ob@meta.data[,groupby])))
            gs = lapply(topn_markers2vis, function(x) VlnPlot(seurat_ob, features = x, cols = colors2use,x.lab.rot = T,
                                                            pt.size = pointsize, alpha = opt$alpha2use,
                                                            group.by = groupby, split.by = splitby)+
                                                      labs(title = "",y = x) + 
                                                      theme(legend.position = "none", 
                                                            # panel.spacing = unit(.05, "lines"),
                                                            axis.title.x = element_text(size = 0),
                                                            axis.title.y = element_text(size = 12),
                                                            axis.text.y=element_text(size = 8)))
            if(any(nchar(levels(seurat_ob@meta.data[[groupby]])) >20)){char_length_excess = TRUE}else{char_length_excess = FALSE}
	    if (char_length_excess == TRUE) {
            i=1
            while (i <=  length(gs)){
              gs_sub = gs[i:(i+9)]
              #print(paste0("gs_sub ",gs_sub," i = ", i))
              gs_sub = gs_sub[which(!sapply(gs_sub, is.null))]

              pdf(file.path(output_dir,paste0("top", "marker_gene_violin_plot",ifelse(i==1,"",i%/%10+1),".pdf", collapse = "_")),
                              width = 8, height = length(gs_sub)*2 + max(nchar(as.character(seurat_ob@meta.data[,groupby])))/2)
              grid.arrange(grobs = gs_sub, ncol=1)
              dev.off()
              png(file.path(output_dir,paste0("top","marker_gene_violin_plot",ifelse(i==1,"",i%/%10+1),".png",collapse = "_")),
                              width = 8, height = length(gs_sub)*2+max(nchar(as.character(seurat_ob@meta.data[,groupby])))/2 , res = 96, units = "in")
              grid.arrange(grobs = gs_sub, ncol=1)
              dev.off()
              i=i+10
            }}else{
	    i=1
            while (i <=  length(gs)){
              gs_sub = gs[i:(i+49)]
              gs_sub = gs_sub[which(!sapply(gs_sub, is.null))]

              pdf(file.path(output_dir,paste0("top", "marker_gene_violin_plot",ifelse(i==1,"",i%/%50+1),".pdf", collapse = "_")),
                              width = 8, height = length(gs_sub)*2 + max(nchar(as.character(seurat_ob@meta.data[,groupby])))/4)
              grid.arrange(grobs = gs_sub, ncol=1)
              dev.off()
              png(file.path(output_dir,paste0("top","marker_gene_violin_plot",ifelse(i==1,"",i%/%50+1),".png",collapse = "_")),
                              width = 8, height = length(gs_sub)*2+max(nchar(as.character(seurat_ob@meta.data[,groupby])))/4 , res = 96, units = "in")
              grid.arrange(grobs = gs_sub, ncol=1)
              dev.off()
              i=i+50}
            }
        }
    }
    if ( vismethod == "featureplot" ){
        for ( clusterx in unique(topn_markers$folder_suffix) ){

            if ( !opt$reduct %in% names(Key(seurat_ob)) ){
                stop( "NO specified reduction found in the object!")
            }else{
                reduct = opt$reduct
            }
            topn = topn_markers %>% filter( folder_suffix == clusterx) %>% select(cluster,gene,folder_suffix)
            topn_markers2vis = as.vector(topn$gene)

            path4vis = file.path(root_dir, paste0("markers_vis4",clusterx,collapse = ""))
            if ( file.exists( path4vis ) ){
                output_dir = path4vis
            }else{
                output_dir = path4vis
                dir.create(output_dir, recursive = TRUE)
            }
            suppressMessages({
                ggfeature = lapply(topn_markers2vis,
                                    function(x) FeaturePlot(seurat_ob,
                                            features = x,cols = c("grey","red"),
                                            reduction= reduct,ncol = 2,
                                            pt.size = pointsize, order= T ) +
                                            theme( plot.title = element_text(hjust = 0.5)) + sc)
            })
        #     nrow = ceiling(length(topn_markers2vis)/2)
        #     pdf(file.path(output_dir,paste0("top", "marker_gene_featureplot.pdf", collapse = "_")),
        #         width =  ifelse(length(topn_markers2vis)==1,yes = 6, no =12),
        #         height = nrow*6)
	    # grid.arrange(grobs = ggfeature, nrow = nrow, ncol= ifelse(length(topn_markers2vis)==1,yes = 1, no =2))
        #     dev.off()
        #     png(file.path(output_dir,paste0("top","marker_gene_featureplot.png", collapse = "_")),
        #         width = ifelse(length(topn_markers2vis)==1,yes = 6, no = 12),
        #         height = nrow*6, res = 96, units = "in")
        #     grid.arrange(grobs = ggfeature, nrow = nrow, ncol=ifelse(length(topn_markers2vis)==1,yes = 1, no =2))
        #     dev.off()
            i=1
            while (i <=  length(ggfeature)){
              ggfeature_sub = ggfeature[i:(i+49)]
              ggfeature_sub = ggfeature_sub[which(!sapply(ggfeature_sub, is.null))]
              nrow = ceiling(length(ggfeature_sub)/2)
              pdf(file.path(output_dir,paste0("top", "marker_gene_featureplot",ifelse(i==1,"",i%/%50+1),".pdf", collapse = "_")),
                width =  ifelse(length(ggfeature_sub)==1,yes = 6, no =12),
                height = nrow*6)
              grid.arrange(grobs = ggfeature_sub, ncol= ifelse(length(ggfeature_sub)==1,yes = 1, no =2))
              dev.off()
              png(file.path(output_dir,paste0("top","marker_gene_featureplot",ifelse(i==1,"",i%/%50+1),".png",collapse = "_")),
                width = ifelse(length(ggfeature_sub)==1,yes = 6, no = 12),
                height = nrow*6, res = 96, units = "in")
              grid.arrange(grobs = ggfeature_sub, ncol= ifelse(length(ggfeature_sub)==1,yes = 1, no =2))
              dev.off()
              i=i+50
            }
        }
    }
    if ( vismethod == "splitby_featureplot" ){
        for ( clusterx in unique(topn_markers$folder_suffix) ){

            if ( !opt$reduct %in% names(Key(seurat_ob)) ){
                stop( "NO specified reduction found in the object!")
            }else{
                reduct = opt$reduct
            }
            topn = topn_markers %>% filter( folder_suffix == clusterx) %>% select(cluster,gene,folder_suffix)
            topn_markers2vis = as.vector(topn$gene)

            path4vis = file.path(root_dir, paste0("markers_vis4",clusterx,collapse = ""))
            if ( file.exists( path4vis ) ){
                output_dir = path4vis
            }else{
                output_dir = path4vis
                dir.create(output_dir, recursive = TRUE)
            }
            split <- SplitObject(seurat_ob,split.by=splitby)
            if( length(levels(seurat_ob@meta.data[,splitby])) == length(unique(seurat_ob@meta.data[,splitby])) ) {
            order=levels(seurat_ob@meta.data[,splitby]) 
            }else{
            order=sort(unique(seurat_ob@meta.data[,splitby])) 
            } 
        for (i in topn_markers2vis){ 
            suppressMessages({
            gene_range <- range(FetchData(object = seurat_ob, vars = i))
            min.cutoff <- gene_range[1]
            max.cutoff <- gene_range[2]
            ggfeature = lapply( order,
                    function(x)p =  FeaturePlot(split[[ x ]],
                        features = i,
                        reduction= reduct,max.cutoff=max.cutoff,
                        ncol = 2, pt.size = pointsize, order= T  ) +  
                        theme( plot.title = element_text(hjust = 0.5) ) + 
                        labs(title=x ) )
            ggfeature <- lapply(ggfeature, function(x) x
                        + scale_colour_gradientn(colors= c("grey","red"),limits = c(min.cutoff,max.cutoff)
                                ) )
            plot <- do.call(ggpubr::ggarrange,
                c(ggfeature, list(ncol = 2,
                            nrow = ceiling(length(ggfeature) / 2),
                            common.legend = TRUE,
                            legend = "right",
                            align = "none")))
            g = plot+ labs(title=i)+
                theme(plot.title = element_text( hjust = 0.5,face ="bold",size=20),
				plot.margin = margin(t = 10,r = 10,b = 10,l = 10))
			}) 
            nrow=ceiling(length(ggfeature) / 2)       
            ggsave(file.path(output_dir,paste0(i,"_",reduct,ifelse(is.null(splitby),"",paste0("_splitby_",splitby)),".pdf",collapse = "")),plot=g,width = ifelse(length(ggfeature)==1,yes = 5, no = 10),height = nrow*5,limitsize = FALSE)
            ggsave(file.path(output_dir,paste0(i,"_",reduct,ifelse(is.null(splitby),"",paste0("_splitby_",splitby)),".png",collapse = "")),plot=g, width = ifelse(length(ggfeature)==1,yes = 5, no = 10),height = nrow*5,limitsize = FALSE ,dpi = 1000)
            }
        }
    }

    if ( vismethod == "ggstatsplot" ){
        seurat_ob = SetIdent( seurat_ob, value = groupby )
        for ( clusterx in unique(topn_markers$folder_suffix) ){
            topn = topn_markers %>% filter( folder_suffix == clusterx) %>% select(cluster,gene,folder_suffix)
            topn_markers2vis = as.vector(topn$gene)

            path4vis = file.path(root_dir, paste0("ggstatsplot_",clusterx,collapse = ""))
            if ( file.exists( path4vis ) ){
                output_dir = path4vis
            }else{
                output_dir = path4vis
                dir.create(output_dir, recursive = TRUE)
            }
        df1 = FetchData(seurat_ob, vars=topn_markers2vis)
        df1 = df1[order(rownames(df1)),]
        data = seurat_ob@meta.data %>%
                               dplyr::rename( "Barcode" = "rawbc")%>%
                               dplyr::select( Barcode, groupby)
        data = data[order(rownames(data)),]
        new_data = cbind(data,df1) %>% dplyr::arrange(!!sym(groupby))
        if (!is.null(new_data$df1)){
            colnames(new_data)[colnames(new_data) == "df1"] <- topn_markers2vis
        }
        write.table(new_data, quote = F,sep ="\t",row.names = F,
             file.path(output_dir,paste0("ggstatsplot_plot.xls",collapse = "")))
        #小提琴
        geneset_R = "/public/scRNA_works/pipeline/scRNA-seq_further_analysis/ggstatsplot.R"
		if ( is.null(opt$pvalue) ){
		system(glue::glue("module purge && source  /home/liuhongyan/miniconda3/bin/activate scvdj  &&
        Rscript {geneset_R}   --input {output_dir}/ggstatsplot_plot.xls -g {groupby} -m ggstatsplot  -o {output_dir}  && module purge && module load OESingleCell/2.0.0 ")) 
        }else{
        pvalue = opt$pvalue
        system(glue::glue("module purge && source  /home/liuhongyan/miniconda3/bin/activate scvdj  &&
        Rscript {geneset_R}   --input {output_dir}/ggstatsplot_plot.xls -g {groupby} -m ggstatsplot -o {output_dir}  -p {pvalue} && module purge && module load OESingleCell/2.0.0 ")) 
         }
      }
   }

    if ( vismethod == "geneset" ){
        path4vis = file.path(root_dir,"geneset_visualization")
        if ( file.exists( path4vis ) ){
            output_dir = path4vis
        }else{
            output_dir = path4vis
            dir.create(output_dir, recursive = T)
        }
        seurat_ob = SetIdent( seurat_ob, value = groupby )

        if (!is.null(opt$scoredata)){
            scoredata <- read.delim(opt$scoredata,sep='\t',row.name=1)
            if(all(seurat_ob@meta.data$rawbc %in% rownames(scoredata))) {       
                scoredata = scoredata[which(rownames(scoredata) %in% seurat_ob@meta.data$rawbc),-1]
                seurat_ob@meta.data = cbind(seurat_ob@meta.data,scoredata)
                print(head(seurat_ob@meta.data))
            }else{
                print("rds里的barcode和打分文件中的barcode不匹配")
            }
           
            metadata = seurat_ob@meta.data
            metadata= metadata[,c(groupby,colnames(scoredata))]
            data = seurat_ob@meta.data %>%
                               dplyr::rename( "Barcode" = "rawbc") %>%
                               dplyr::select( Barcode, sampleid, clusters,group,!!colnames(scoredata))
            write.table(data, quote = F,sep ="\t",row.names = F,
                 file.path(output_dir,paste0("geneset",".addmodeulescore.xls",collapse = "")))
            matrix = seurat_ob@meta.data %>%
                               dplyr::rename( "Barcode" = "rawbc") %>%
                               dplyr::select( Barcode,groupby,!!colnames(scoredata)) %>%
                               dplyr::arrange(!!sym(groupby))
            write.table(matrix, quote = F,sep ="\t",row.names = F,
                 file.path(output_dir,paste0("geneset",".addmodeulescore_plot.xls",collapse = "")))
 
            plot_term=colnames(scoredata)
        }else{
            topn_markers2vis=list()
            for ( clusterx in unique(topn_markers$folder_suffix) ){
                topn_markers2vis[[clusterx]] = subset(topn_markers,folder_suffix == clusterx)$gene
            }
            seurat_ob = AddModuleScore(seurat_ob,features=topn_markers2vis,name=names(topn_markers2vis),seed=1)
            colnames(seurat_ob@meta.data)[(dim(seurat_ob[[]])[2]-length(topn_markers2vis)+1):dim(seurat_ob[[]])[2]] = names(topn_markers2vis)
            metadata = seurat_ob@meta.data
            metadata= metadata[,c(groupby,names(topn_markers2vis))]
            data = seurat_ob@meta.data %>%
                               dplyr::rename( "Barcode" = "rawbc") %>%
                               dplyr::select( Barcode, sampleid, clusters,group,!!names(topn_markers2vis) )
            write.table(data, quote = F,sep ="\t",row.names = F,
                 file.path(output_dir,paste0("geneset",".addmodeulescore.xls",collapse = "")))
            matrix = seurat_ob@meta.data %>%
                               dplyr::rename( "Barcode" = "rawbc") %>%
                               dplyr::select( Barcode,groupby,!!names(topn_markers2vis) ) %>%
                               dplyr::arrange(!!sym(groupby))
            write.table(matrix, quote = F,sep ="\t",row.names = F,
                 file.path(output_dir,paste0("geneset",".addmodeulescore_plot.xls",collapse = "")))

            plot_term=names(topn_markers2vis)
        }
        #绘图
        geneset_R = "/public/scRNA_works/pipeline/scRNA-seq_further_analysis/ggstatsplot.R"
		if ( is.null(opt$pvalue) ){
		system(glue::glue("module purge && source  /home/liuhongyan/miniconda3/bin/activate scvdj  &&
        Rscript {geneset_R}   --input {output_dir}/geneset.addmodeulescore_plot.xls -g {groupby} -o {output_dir}  && module purge && module load OESingleCell/2.0.0 ")) 
        }else{
        pvalue = opt$pvalue
        system(glue::glue("module purge && source  /home/liuhongyan/miniconda3/bin/activate scvdj  &&
        Rscript {geneset_R}   --input {output_dir}/geneset.addmodeulescore_plot.xls -g {groupby} -o {output_dir}  -p {pvalue} && module purge && module load OESingleCell/2.0.0 ")) 
        }
        feat=list()
        for (geneset in plot_term){
          if ( !opt$reduct %in% names(Key(seurat_ob)) ){
              stop( "NO specified reduction found in the object!")
          }else{
              reduct = opt$reduct
          }
          feat[[geneset]]=FeaturePlot(seurat_ob,features = geneset,
                                      cols = c("grey","red"),split.by = NULL, reduction = reduct, ncol = 2, pt.size = 0.4, order = T) +
            theme( plot.title = element_text(hjust = 0.5)) + sc
          ggsave(file.path(output_dir,paste0(geneset,"_score_featureplot.pdf")),plot=feat[[geneset]])
          ggsave(file.path(output_dir,paste0(geneset,"_score_featureplot.png")),plot=feat[[geneset]])
        }
    }

if ( vismethod == "corscatter" ){
        path4vis = file.path(root_dir,"geneset_corscatter")
        if ( file.exists( path4vis ) ){
            output_dir = path4vis
        }else{
            output_dir = path4vis
            dir.create(output_dir, recursive = T)
        }
        if ( length(unique(topn_markers$folder_suffix))!=2 ){stop("NO geneset genes list is applicable, it must be two column")}
        topn_markers2vis=list()
        for ( clusterx in unique(topn_markers$folder_suffix) ){
            topn_markers2vis[[clusterx]] = subset(topn_markers,folder_suffix == clusterx)$gene
        }
        seurat_ob = AddModuleScore(seurat_ob,features=topn_markers2vis,name=names(topn_markers2vis),seed=1)
        ### rename 
        colnames(seurat_ob@meta.data)[(dim(seurat_ob[[]])[2]-length(topn_markers2vis)+1):dim(seurat_ob[[]])[2]] = names(topn_markers2vis)
        # data1 addmodule score
        data=seurat_ob@meta.data[,c(names(topn_markers2vis),groupby)]

        # data2 ave
        nlevel=length(table(seurat_ob[[groupby]]))
        groupby_data = vector()
        for (i in names(table(seurat_ob[[groupby]]))  ) {
            sub_ob = SubsetData(seurat_ob, subset.name= groupby,accept.value=i)
            normalized_data = as.matrix(sub_ob@meta.data[,names(topn_markers2vis)])
            meta.data = sub_ob@meta.data %>% tibble::rownames_to_column(var = "id")
            groupby_data = cbind(groupby_data,colMeans(normalized_data))
        }
        colnames(groupby_data) = names(table(seurat_ob[[groupby]]))
        groupby_data = t(groupby_data) %>% as.data.frame %>% tibble::rownames_to_column(groupby)
        p=ggplot(data) + 
        geom_point(alpha=0.6,size=pointsize,aes_string(names(topn_markers2vis)[1],names(topn_markers2vis)[2],col=groupby))+
        scale_color_manual(values=CustomCol2(1:nlevel))+
        geom_point(data=groupby_data,shape=17,size=5, aes_string(names(topn_markers2vis)[1],names(topn_markers2vis)[2],col=groupby))+
        theme_bw() + 
        labs(x=paste0(names(topn_markers2vis)[1],"_score"),y=paste0(names(topn_markers2vis)[2],"_score")) + 
        theme(strip.background = element_rect(fill = NA,color = NA), 
        plot.background= element_rect(fill = NA,color = NA),
        panel.background = element_rect(fill = NA,color = NA))+
        # theme(panel.grid = element_blank()) +
        theme( axis.text.x= element_text(
                size=9,angle=45,hjust=1,lineheight=1.2),
                plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))

        # p = p+ geom_hline(yintercept = 0, linetype = 5, color = c("#ADB6B6FF"), size = 0.1) + geom_vline(xintercept = 0, linetype = 5, color = c("#ADB6B6FF"), size = 0.1)
        ggsave(file.path(output_dir,paste0("geneset_corscatter.pdf")))
        ggsave(file.path(output_dir,paste0("geneset_corscatter.png")),dpi = 1000)
}

    if ( vismethod == "dotplot" ){
        for ( clusterx in unique(topn_markers$folder_suffix) ){
            topn = topn_markers %>% filter( folder_suffix == clusterx) %>% select(cluster,gene,folder_suffix)
            topn_markers2vis = unique(as.vector(topn$gene))
            path4vis = file.path(root_dir,paste0("markers_vis4",clusterx,collapse = ""))
            if ( file.exists( path4vis ) ){
                output_dir = path4vis
            }else{
                output_dir = path4vis
                dir.create(output_dir, recursive = T)
            }
            seurat_ob = SetIdent( seurat_ob, value = groupby )
            dot_plot = DotPlot(object = seurat_ob, features = topn_markers2vis ) + RotatedAxis() +sc
            dot_plot + guides(color = guide_colorbar(order = 1, title = "Average Expression"))
            ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.pdf", collapse = "_")),limitsize = FALSE,width=(0.3*length(topn_markers2vis)+2.5+max(nchar(names(table(seurat_ob@meta.data[,groupby]))))/10))
            ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.png", collapse = "_")),dpi = 300 ,limitsize = F ,width=(0.3*length(topn_markers2vis)+2.5+max(nchar(names(table(seurat_ob@meta.data[,groupby]))))/10))
        }
    }

    if ( vismethod == "ridgeplot"){
        for ( clusterx in unique(topn_markers$folder_suffix) ){
            topn = topn_markers %>% filter( folder_suffix == clusterx) %>% select(cluster,gene,folder_suffix)
            topn_markers2vis = as.vector(topn$gene)
            colors2use = CustomCol2(1:length(unique(seurat_ob@meta.data[,groupby])))

            path4vis = file.path(root_dir,paste0("markers_vis4",clusterx,collapse = ""))
            if ( file.exists( path4vis ) ){
                output_dir = path4vis
            }else{
                output_dir = path4vis
                dir.create(output_dir, recursive = T)
            }
            seurat_ob = SetIdent( seurat_ob, value = groupby )
            nrow = ceiling(length(topn_markers2vis)/2)
            RidgePlot( object = seurat_ob, features = topn_markers2vis,
                        cols = colors2use, ncol = 2)
            ggsave(file.path(output_dir,paste0( "top", "marker_gene_ridgeplot.pdf",collapse = "_" )),width =  ifelse(length(topn_markers2vis)==1,yes = 6, no =12),height = nrow*5)
            ggsave(file.path(output_dir,paste0( "top", "marker_gene_ridgeplot.png", collapse = "_")),dpi = 1000 ,limitsize = F,width =  ifelse(length(topn_markers2vis)==1,yes = 6, no =12),height = nrow*5)
        }
    }

    if ( vismethod == "violinensemble"){
        # colors2use = CustomCol2(1:length(unique(seurat_ob@meta.data[,groupby])))
        # ggensemble = ViolinEnsemble( object =  seurat_ob, features = topn_markers2vis, cols = colors2use,
                                        # group.by = groupby, show_point = T, pt_size = pointsize)
        # ggsave(file.path(output_dir,paste0( "all_top", "marker_gene_enhanced_violin_plot.pdf",collapse = "_" )))
        # ggsave(file.path(output_dir,paste0( "all_top", "marker_gene_enhanced_violin_plot.png", collapse = "_")),
                    # dpi = 1000 ,limitsize = F)
        for ( clusterx in unique(topn_markers$cluster) ){
            topn = topn_markers %>% filter( cluster == clusterx) %>% select(cluster,gene)
            topn_markers2vis = as.vector(topn$gene)
        
            path4vis = file.path(root_dir,paste0("markers_vis4cluster",clusterx,collapse = ""))
            if ( file.exists( path4vis ) ){
                output_dir = path4vis
            }else{
                output_dir = path4vis
                dir.create(output_dir, recursive = T)
            }
			colors2use = CustomCol2(1:length(unique(seurat_ob@meta.data[,groupby])))
            ggensemble = ViolinEnsemble( object =  seurat_ob, features = topn_markers2vis, cols = colors2use,
                                         show_point = F)
            ggsave(file.path(output_dir,paste0( "top", "marker_gene_enhanced_violin_plot.pdf",collapse = "_" )))
            ggsave(file.path(output_dir,paste0( "top", "marker_gene_enhanced_violin_plot.png", collapse = "_")),
                            dpi = 1000 ,limitsize = F)
        }
    }

    if ( vismethod == "heatmap" ){
        markers2vis4heatmap = unique(as.vector(topn_markers$gene))
        if ( is.null(opt$sample_ratio) ){
            subseted_seurat = seurat_ob
        }else{
            sampled_cellmeta = seurat_ob@meta.data %>% rownames_to_column() %>%
                                group_by( .dots= groupby ) %>%
                                sample_frac( size = opt$sample_ratio,replace = F) %>% column_to_rownames()
            subseted_seurat = SubsetData(seurat_ob, cells = rownames(sampled_cellmeta))
        }
        subseted_seurat@meta.data[,groupby] = as.factor(subseted_seurat@meta.data[,groupby])
        colors2use = CustomCol2(1:length(unique(subseted_seurat@meta.data[,groupby])))
        if (length(markers2vis4heatmap) > 135){
            sz = 4-log(length(markers2vis4heatmap)/100)
            heig = 5+log2(length(markers2vis4heatmap)/10)
            wid = 7.5
        }else if (length(markers2vis4heatmap) < 75){
            sz = 6-log2(length(markers2vis4heatmap)/80);heig = 7;wid = 7
        }else{
            sz = 4-log2(length(markers2vis4heatmap)/120);heig = 7;wid = 7
        }

        ggheat = DoHeatmap( object = subseted_seurat,
                            features = markers2vis4heatmap,
                            group.colors = colors2use,
                            group.by = groupby, group.bar = T, label = F) +
                            theme(axis.text.y = element_text(size = sz, face = "bold"))
                            # group.cex = 10, cex.row = 4,
                            # slim.col.label = T, group.label.rot = F)
        ggheat + guides(fill = guide_colorbar( title.position = "top", order = 1), color = guide_legend(order = 2, override.aes = list(alpha = 1)))
        ggsave(file.path(root_dir,paste0("top","marker_gene_heatmap.pdf", collapse = "_")),height = heig, width = wid)
        # ggsave(file.path(root_dir, paste0("top", "marker_gene_heatmap.png", collapse = "_")),height = heig, width = wid, dpi = 1000 ,limitsize = F)
        setwd(root_dir)
        print("Convert pdf to png...")
        system(" /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 800 -trim  topmarker_gene_heatmap.pdf  -quality 500  -flatten  topmarker_gene_heatmap.png ")
    }

    if ( vismethod == "diff_heatmap" ){
        markers2vis4heatmap = unique(as.vector(topn_markers$gene))
        if ( is.null(opt$sample_ratio) ){
            subseted_seurat = seurat_ob
        }else{
            sampled_cellmeta = seurat_ob@meta.data %>% rownames_to_column() %>%
                                group_by( .dots= groupby ) %>%
                                sample_frac( size = opt$sample_ratio,replace = F) %>% column_to_rownames()
            subseted_seurat = SubsetData(seurat_ob, cells = rownames(sampled_cellmeta))
        }
        subseted_seurat@meta.data[,groupby] = as.factor(subseted_seurat@meta.data[,groupby])
        colors2use = CustomCol2(1:length(unique(subseted_seurat@meta.data[,groupby])))
        ggheat = DoHeatmap( object = subseted_seurat,
                            features = markers2vis4heatmap,
                            group.colors = colors2use,
                            group.by = groupby, group.bar = T, label = F) +
                            theme(axis.text.y = element_text(size = 6-log2(length(markers2vis4heatmap)/80), face = "bold"))
                            # group.cex = 10, cex.row = 4,
                            # slim.col.label = T, group.label.rot = F)
        ggheat +  guides(fill = guide_colorbar( title.position = "top", order = 1), color = guide_legend(order = 2, override.aes = list(alpha = 1)))
        ggsave(file.path(root_dir, paste0("top", opt$topn, "_", name, "_heatmap.pdf", collapse = "")))
        # ggsave(file.path(root_dir, paste0("top", opt$topn, "_", name, "_heatmap.png", collapse = "")), dpi = 1000 ,limitsize = F)
        setwd(root_dir)
        print("Convert pdf to png...")
        system("for i in  `ls top*_heatmap.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")
    }


    if ( vismethod == "velocity" ){
        markers2vis4velo = unique(as.vector(topn_markers$gene))
        emb = Embeddings( seurat_ob, reduction = reduct)
        cell.dist <- as.dist(1-armaCor(t( emb )))
        emat = seurat_ob[["spliced"]]
        nmat = seurat_ob[["unspliced"]]
        colors2use = CustomCol2(1:length(unique(seurat_ob@meta.data[,groupby])))
        rvel.cd = Tool(object = seurat_ob, slot = "RunVelocity")
        gene.relative.velocity.estimates(emat,nmat,deltaT=1,
                    kCells = 25,kGenes=1,fit.quantile=0.2,
                    cell.emb=emb,cell.colors=colors2use,cell.dist=cell.dist,
                    show.gene=gene,old.fit=rvel.cd,do.par=T)
    }
}

