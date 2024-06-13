#!/usr/bin/env Rscript
# Objective : Processing and visualization of SCENIC results
# Created by: hanmin
# Created on: 2020/5/20
#
rm(list=ls())

BinaryCount <- function(
  object,
  method = c("kmeans", "aucell" ),
  auc = NULL,
  nCores = 10,
  do.tfidf = FALSE,
  ...
){
  if ( !is.null(auc) ){
    if ( class(auc) == "matrix" ){
      AUC_mat <- auc
    }else if ( class(auc) == "aucellResults"){
      regulonAUC <- auc
      AUC_mat <- AUCell::getAUC(regulonAUC)
    }
  }else{
    regulonAUC <- Tool(object, slot = "RunAUCell")
    if ( is.null(regulonAUC) ){
      stop("NO regulon AUC matrix supplied or found in the object!")
    }else{
      AUC_mat <- AUCell::getAUC(regulonAUC)
    }
  }

  binary_mat <- switch (tolower(method),
    "aucell" = {
      cells_AUCellThresholds <- AUCell::AUCell_exploreThresholds(regulonAUC,
                                                                 smallestPopPercent=0.25,
                                                                 assignCells=TRUE, plotHist=FALSE,
                                                                 verbose=FALSE, nCores=nCores)
      # Get cells assigned to each regulon
      cellsAssigned <- AUCell::getAssignments(cells_AUCellThresholds)
      # cellsAssigned   <- lapply(cells_AUCellThresholds, function(x) x$assignment)
      assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
      colnames(assignmentTable)[2] <- "geneSet"
      assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
      binary_mat <- tidyr::spread(as.data.frame.table(assignmentMat),key = "Var2","Freq" ) %>%
                      tibble::column_to_rownames(var = "Var1") %>% as.matrix()
    },
    "kmeans"= {
    # Iterate over each regulon in the AUC matrix
    AUC_df <- as.data.frame(AUC_mat) %>%
          tibble::rownames_to_column(var = "regulon") %>%
          tidyr::gather("cell", "auc", -regulon) %>%
          dplyr::filter( auc > 0) %>%
          dplyr::group_by( regulon ) %>%
          dplyr::mutate( cluster = as.factor(kmeans(auc, centers = 2)$cluster)) %>%
          dplyr::ungroup()
    kmeans_thresholds <- lapply(split(AUC_df,as.factor(AUC_df$regulon),drop = F), function(df){
      cluster1_max <- max(subset(df,cluster == 1)$auc)
      cluster2_max <- max(subset(df,cluster == 2)$auc)

      if(cluster1_max > cluster2_max){
        df <- df %>% mutate("cluster" = gsub(2,3,cluster)) %>%
          mutate("cluster" = gsub(1,2,cluster)) %>%
          mutate("cluster" = gsub(3,1,cluster))
      }

      df <- df %>% arrange(desc(auc))
      df_sub <- df %>% subset(cluster == 1)
      auc_thresholds <- df_sub[1,]$auc
    })
    binary_mat <- as.data.frame(AUC_mat) %>%
          tibble::rownames_to_column(var = "regulon") %>%
          tidyr::gather("cells", "auc", -regulon) %>%
          dplyr::group_by( regulon ) %>%
          mutate("values"= if_else(auc >= kmeans_thresholds[regulon],1,0)) %>%
          dplyr::ungroup() %>% select(-auc) %>%
          tidyr::spread("cells", "values") %>%
          tibble::column_to_rownames(var = "regulon") %>% as.matrix()
  })

  object <- subset(object, cells = colnames(binary_mat) )
  row.names(binary_mat) = gsub("_", "-",rownames(binary_mat) )
  object[["SCENIC"]] <- CreateAssayObject(data = binary_mat)
  if ( do.tfidf ){
    object[["SCENIC"]] <- SetAssayData(object[["SCENIC"]], slot = "data",
                                     new.data =as(TF.IDF(binary_mat), 'dgCMatrix') )
  }

  object = LogSeuratCommand(object)
  return( object )
}


##' Calculates Regulon specificity score (RSS) from binary regulon activity.
##' Iterate over all cell types and perform jensen shannon divergence test
##' using binary regulon activity and genotype
##'
##' @param object the Seurat object with binarized assay "SCENIC".
##' @param assay the SCENIC assay as default.
##' @param slot "data" as default.
##' @param metadata Dataframe containing metadata about cells. Has to create a column named cell_type that assigns groupings to cells.
##' Can be the meta.data slot from a Seurat object.
##' @param binary_regulons Data frame with binary regulons, where regulons are rows and columns are cells. Can be created from output of binarize_regulons().
##' @param group.by the groupping factor for RSS calculaion.
##' @import  philentropy JSD
##' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
##' @export
##' @examples
RunRSS <- function(
  object,
  assay = "SCENIC",
  slot = "data",
  binary_regulons = NULL,
  metadata = NULL,
  group.by = "cell_type"
){
  if ( !is.null( binary_regulons) ){
    regulons <- rownames(binary_regulons)
  }else{
    binary_regulons <- Seurat::GetAssayData(object, assay = assay, slot = slot)
    regulons <- rownames(binary_regulons)
  }

  if ( is.null(metadata) ){
    metadata <- object@meta.data
  }
  cell_types <- unique(metadata[,group.by])
  jsd_matrix_ct <- data.frame("regulon" = c(), "cell_type" = c(), "jsd" = c())

  cell_type_counter <- 0
  for(ct in unique(cell_types)) {
    cell_type_counter <- cell_type_counter + 1
    print(paste("Processing cell type:",cell_type_counter,ct,sep=" "))
    for(regulon_no in 1:length(regulons)) {
      regulon <- regulons[regulon_no]
      regulon_vec <- binary_regulons[regulon,]
      regulon_vec_sum <- sum(regulon_vec)
      ## Check that there are cells with binary activity > 0 for this regulon
      if(regulon_vec_sum > 0){
        #progress(regulon_no)
        regulon_norm <- regulon_vec/regulon_vec_sum
        genotype_vec <- metadata[colnames(binary_regulons),]
        genotype_vec <- genotype_vec %>%
          mutate("cell_class" = if_else(get(group.by) == ct,1,0))
        genotype_vec <- genotype_vec$cell_class
        genotype_norm <- genotype_vec/sum(genotype_vec)
        dist_df <- rbind(regulon_norm,genotype_norm)
        ## Calculate the Jensen-Shannon divergence
        jsd_divergence <- suppressMessages(philentropy::JSD(dist_df))
        ## Calculate Jensen-Shannon distance
        rss <- 1-sqrt(jsd_divergence)
        regulon_jsd <- data.frame("regulon" = regulon, "cell_type" = ct, "RSS" = rss[1])
        jsd_matrix_ct <- rbind(jsd_matrix_ct,regulon_jsd)
      }else if(regulon_vec_sum == 0){
        print(paste("Filtered out:",regulon,". No cells with binary activity > 0 identified. Please check your threshold for this regulon!",sep=""))
      }
    }
  }

  jsd_matrix_ct <- jsd_matrix_ct %>% dplyr::arrange(desc(RSS))
  jsd_matrix_ct[,group.by] <- jsd_matrix_ct$cell_type
  jsd_matrix_ct <- jsd_matrix_ct %>% dplyr::select(-cell_type)
  return(jsd_matrix_ct)
}

##' Calculates Regulon specificity score (RSS) from binary regulon activity.
##'
##' @param rrs_df Data frame containing RSS scores for all regulons over all cell types. Can be created with calculate_rrs.
##' @param cell_type Cell type for which to plot jsd ranking. Select "all" to plot a facet plot over all cell types.
##' @param ggrepel_force same as the force parameter for geom_text_repel.
##' @param ggrepel_point_padding same as the force parameter for geom_text_repel.
##' @param top_genes Number of top genes to label in plot using ggrepel.
##' @keywords SCENIC, regulons, RRS, cell type classification
##' @export
##' @examples
##'
## Plot JSD enrichment plot for specific cell type
RSSRanking <- function(
  rrs_df,
   group.by,
   ggrepel_force = 1,
   ggrepel_point_padding = 0.2,
   top_genes = 4,
   plot_extended = FALSE
){
  require(ggrepel)
  require(cowplot)

  if(plot_extended == TRUE){
    rrs_df <- rrs_df %>%
      subset(grepl("extended",regulon))
  }else if(plot_extended == FALSE){
    rrs_df <- rrs_df %>%
      subset(!grepl("extended",regulon))
  }

  rrs_df_sub <- rrs_df %>% dplyr::group_by(.dots = group.by) %>%
    mutate("rank" = order(order(RSS, decreasing = TRUE)))

  #jsd_matrix_sub$regulon <- factor(jsd_matrix_sub$regulon,levels = unique(jsd_matrix_sub$regulon))
  rrs_ranking_plot <- ggplot(rrs_df_sub,aes(rank,RSS,label = regulon)) +
    geom_point(color = "grey20",size = 2) +
    geom_point(data = subset(rrs_df_sub,rank < top_genes),
               color = "red",size = 2) +
    geom_text_repel(data = subset(rrs_df_sub,rank < top_genes),
                    force = ggrepel_force,point.padding = ggrepel_point_padding) +
    theme_bw() + theme(panel.grid =element_blank()) +
    labs(x = "Rank", y = "RSS", title = group.by) +
    facet_wrap(eval(expr(~!!ensym(group.by))), ncol = 2, scales = "free_y" )
  return(rrs_ranking_plot)
}


#================================================================================
# package loading
#================================================================================
suppressWarnings({
    #========import packages=====================================
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages( library("OESingleCell") )
    suppressPackageStartupMessages( library("tibble") )
})

#=command line parameters setting=============================
option_list = list(
make_option( c("--input", "-i"), type = "character",
    help = "the seurat object saved as R object in RDS format." ),
make_option( c("--auc", "-v"), type = "character",  default = NULL, 
   help = "The regulon activity scores results from the SCENIC output." ),
make_option( c("--aucformat", "-f" ), type = "character", default = "rds",
    help = "The indication of type for input AUC, the choices can be:
                            rds: the aucellResults object in RDS format,
                                 usually the 3.4_regulonAUC.Rds from SCENIC;
                            xsv: the delimited table from AUC."),
make_option( c("--topGenes", "-t"), type = "integer", default = NULL,
    help = "Number of top genes to label in plot of rrs ranking." ),
make_option( c("--groupby", "-c"), type = "character",
    help = "the groupping column of cells in the metadata of seurat object." ),
make_option( c("--binmethod", "-m"), type = "character", default = "aucell",
    help = "the binary methods used to binarize the regulon activity matrix element to 0/1.
            Options can be aucell and kmeans. Note that 'aucell' is only avaiable for aucellResults object" ),
make_option( c("--ncores", "-j" ), type="integer", default = 10,
    help="the number of CPUs used to improve the performace."),
make_option( c("--threshold", "-s" ), type="double",
    help="subset the regulon according to the threshold of RSS"),
make_option( c("--extended"),type="logical", default = FALSE,
    help="whether to use the extended regulons for calculation and visualization"),
make_option( c("--output", "-o"), type = "character", default = "./",
    help = "the output directory of QC results.", metavar = "outdir" )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if ( is.null(opt$topGenes) ){
    topGenes = 4
}else{
    topGenes = opt$topGenes + 1
}

if ( opt$extended ){
    non_extended = TRUE 
}else{
    non_extended = FALSE 
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

seurat_ob = readRDSMC( opt$input, cores = 6 )
if ( seurat_ob@version < 3 ){
    seurat_ob = UpdateSeuratObject(seurat_ob)
}

if ( !is.null(opt[["auc"]])){ # using the user supplied regulon activity score matrix
    if ( tolower(opt$aucformat) == "rds" ){
        auc = readRDS(opt$auc)
        auc = auc[,rownames(seurat_ob@meta.data)]
        auc = auc[which(rowSums(getAUC(auc))!=0),]
    }else if ( tolower(opt$aucformat) == "xsv" ){
        auc = vroom::vroom(opt$auc, col_names = T) %>% as.data.frame()
        rownames(auc) = auc[,1]
        auc = as.matrix(auc[,-1])
    }
}else { # use the AUC result integrated in the Seurat object by running RunAUCell
    auc = Tool(seurat_ob, slot = "RunAUCell")
    auc = auc[,rownames(seurat_ob@meta.data)]
    auc = auc[which(rowSums(getAUC(auc))!=0),]
}

seurat_ob = BinaryCount(seurat_ob, method = opt$binmethod, auc = auc, nCores = opt$ncores)
#if (!is.null(opt[["auc"]])) saveRDSMC(seurat_ob, opt$input)
rss_df = RunRSS(seurat_ob, group.by = opt$groupby)
rss_df_out = rss_df %>% subset(!grepl("extended",regulon))
write.table(rss_df_out, file.path(output_dir, "regulon_RSS_annotation.xls"),
            sep = "\t", col.names = T, row.names =F, quote =F)
gg_rss_rank = RSSRanking(rss_df, group.by = opt$groupby, top_genes = topGenes, plot_extended = non_extended )
ggsave(file.path(output_dir, "RSS_ranking_plot.pdf" ), plot =gg_rss_rank, width = 6,
        height = ceiling(length(unique(rss_df[,opt$groupby]))/2) * 8)

rss_df_wide <- rss_df_out %>% tidyr::spread_( opt$groupby,"RSS")
rownames(rss_df_wide) <- rss_df_wide$regulon
rss_df_wide <- rss_df_wide[,2:ncol(rss_df_wide)]
## Subset all regulons that don't have at least an RSS of 0.4 for one cell type
if ( is.null(opt$threshold) ){
    rss_threshold = 0
}else{
    rss_threshold = as.numeric(opt$threshold)
}

rss_df_wide_specific <- rss_df_wide[apply(rss_df_wide,MARGIN = 1 ,FUN = function(x) any(x > rss_threshold)),]
pdf(file.path(output_dir, "RSS_heatmap.pdf"),
    width = dim(rss_df_wide_specific)[2]*1.5+1,
    height = dim(rss_df_wide_specific)[1]*0.2)
pheatmap::pheatmap( rss_df_wide_specific,
                    color=colorRampPalette(c("blue","white","red"))(100),
                    angle_col = 45,treeheight_col=10, border_color=NA)

dev.off()
