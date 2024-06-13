#!/usr/bin/env Rscript
# Title     : CSI analysis of reuglons
# Objective : Processing and visualization of SCENIC results
# Created by: hanmin
# Created on: 2020/5/20

#' Creates a list of unique color values used for plotting
#'
#' @return A named vector of unique hexedecimal color values, either generated from a preselected
#'         vector of 20 unique colors, or from a sequence of colors in hsv colorspace.
#'
#' @param seurat.obj A singular preprocessed Seurat object
#' @param gradient Setting to TRUE will use a sequence of hsv colors instead of 20 unique colors,
#'                 useful for comparisons of more than 20 cell types.
#' @param value The Seurat metadata slot to generate colors for. Defaults to "celltype".
#'
#' @import SingleCellExperiment
#' @import Seurat
#'
#' @seealso \code{\link{as.SingleCellExperimentList}}
#' @seealso \code{\link{ExtractGenes}}
#' @seealso \code{\link{DecoderVariance}}
#' @seealso \code{\link{MeanDecoderVariance}}
#' @seealso \code{\link{GetCharMetadata}}
#'
#' @examples
#' DimPlot(object = seurat.obj,
#'         reduction = "tsne",
#'         cols = SelectColors(seurat.obj),
#'         group.by = "celltype",
#'         label = TRUE,
#'         repel = TRUE)
#'
#' @export

save_ggplots <- function( filename = NULL,
    plot = NULL,
    width = 6,
    height = 8,
    dpi = 300,
    to.pdf = TRUE,
    to.png = TRUE,
    ...
) {
    if ( to.png ) {
        ggplot2::ggsave(paste(filename, "png", sep = "."),
            plot,
            device = "png",
            width = width,
            height = height,
            units = "in",
            dpi = dpi,
            ...
        )
    }
    if ( to.pdf ) {
        ggplot2::ggsave(paste(filename, "pdf", sep = "."),
            plot,
            device = "pdf",
            width = width,
            height = height,
            units = "in",
            dpi = dpi,
            ...
        )
	}
}

SelectColors <- function(
    object = NULL,
    palette = "blindless",
    value = "celltype",
    n = NULL
){
    if ( !is.null(object) ){
        if ( class(object) == "data.frame" ){
            colid <- ifelse( is.null(value), colnames(object)[1], value )
            if (is.factor(object[[colid]])) {
                names= levels(object[[colid]])
            } else {
                names <- unique(object[[colid]])
            }
        }
        if ( is.factor(object) ){
            names <- levels(object)
        }
        if ( is.vector(object) ){
            names <- unique(object)
        }
        n = length(names)
    }else if ( !is.null(n) ) {
        names = NULL
    }

    colors2pick = switch(palette,
    ##ref: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
    ditto = c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
    "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
    "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
    "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
    "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"),
    CustomCol2 = c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17"),
    seurat = hcl( h = seq(15, 375, length = n+1), l = 65, c = 100),
    col50 =  c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
    "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
    "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
    "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
    "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
    "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
    "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
    "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
    "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
    "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c"),
    paired = brewer.pal(n = n, 'Paired'),
    colx22 = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
    '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
    '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
    '#000075', '#808080', '#4f34ff', '#f340F0'),
    jet = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
    "#FF7F00", "red", "#7F0000" ),
    tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
    "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
    "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
    "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
    tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
    "#CDCC5D", "#6DCCDA"),
    colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
    "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
    "#FFBC79", "#CFCFCF"),
    trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
    "#FFC156", "#69B764", "#F26C64", "#FFDD71",
    "#9FCD99"),
    purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
    "#5F5A41", "#B4B19B", "#995688", "#D898BA",
    "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
    bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
    "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
    "#BD0A36", "#F4737A"),
    greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
    "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
    "#39737C", "#86B4A9", "#82853B", "#CCC94D"),
    cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
    "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
    "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
    "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"),
    CustomCol = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
    blindless = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
    "#BF5B17", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4",
    "#B3DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3", "#BC80BD",
    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
    "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
    "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8",
    "#984EA3",  "#FFFF33", "#A65628", "#F781BF", "#999999","#FFED6F",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9","#666666"),
    )

    if ( is.null(n) ){
        colors_use <- colors2pick
    }else{
        colors_use <- colors2pick[1:n]
    }
    if ( !is.null(names) ){ names(colors_use) <- names }
    return(colors_use)
}


my_color <- function(n){
  my_palette=pal_d3("category20")(10)
  return(my_palette[n])
}

#' Calculates CSI values for all regulon pairs
#'
#' @param regulonAUC_mat The AUC  matrix for all regulons as calculated by SCENIC.
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @export
#' @examples
calculate_csi <- function(
  regulonAUC_mat = AUCell::getAUC(regulonAUC),
  calc_extended = FALSE,
  verbose = FALSE
){
  if (calc_extended == FALSE){
    regulonAUC_mat <- subset(regulonAUC_mat,!grepl("extended",rownames(regulonAUC_mat)))
  }

  pearson_cor <- cor(t(regulonAUC_mat))
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    tidyr::gather(regulon_2,pcc,-regulon_1) %>%
    dplyr::mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))

  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names)*length(regulon_names)
  csi_regulons <- data.frame(matrix(nrow=num_of_calculations,ncol = 3))
  colnames(csi_regulons) <- c("regulon_1", "regulon_2", "CSI")
  num_regulons <- length(regulon_names)

  f <- 0
  for(reg in regulon_names){
    for(reg2 in regulon_names){
      f <- f + 1
      # fraction_lower <- calc_csi(reg,reg2,pearson_cor)
      test_cor <- pearson_cor[reg,reg2]
      total_n <- ncol(pearson_cor)
      pearson_cor_sub <- subset(pearson_cor,rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)

      # sums <- apply(pearson_cor_sub,MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
      sums <- apply(pearson_cor_sub, 2,function(m) ifelse( length(m[m>test_cor]) == length(m), 0, length(m)) )
      fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)]) / total_n
      csi_regulons[f,] <- c(reg,reg2,fraction_lower)
    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}

#' Calculate CSI module activity over all cell types
#'
#' @param clusters_df -
#' @param regulonAUC_mat The AUC  matrix for all regulons as calculated by SCENIC.
#' @param metadata -
#' @keywords SCENIC, regulons, CSI activity
#' @export
#' @examples
calc_csi_module_activity <- function(
  clusters_df,
  regulonAUC_mat = AUCell::getAUC(regulonAUC),
  metadata
){
  cell_types<- unique(metadata$cell_type)
  regulons <- unique(clusters_df$regulon)

  regulonAUC_mat <- regulonAUC_mat[as.character(regulons),]

  csi_activity_matrix_list <- list()
  csi_cluster_activity <- data.frame("csi_module" = c(),
                                    "mean_activity" = c(),
                                    "cell_type" = c())

  cell_type_counter <- 0
  groupby_1 =c()   # select group with only one cell
  for(i in names(table(metadata$cell_type)) ){
	if( table(metadata$cell_type)[i]==1 ){
		groupby_1= i
	}
  }
  cell_types <- setdiff(cell_types, groupby_1)
  
  regulon_counter <-
    for(ct in cell_types) {
      cell_type_counter <- cell_type_counter + 1
      cell_type_aucs <- rowMeans(regulonAUC_mat[,rownames(subset(metadata,cell_type == ct))])
	  
      cell_type_aucs_df <- data.frame("regulon" = names(cell_type_aucs),
                                      "activtiy"= cell_type_aucs,
                                      "cell_type" = ct)
      csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
    }
  if(!is.null(groupby_1)){
    regulonAUC_mat_1 = regulonAUC_mat[,rownames(subset(metadata,cell_type == groupby_1))]
    csi_activity_matrix_list_1=data.frame("regulon" = names(regulonAUC_mat_1),
                                      "activtiy"= as.vector(regulonAUC_mat_1),
                                      "cell_type" = groupby_1)
    csi_activity_matrix_list[[groupby_1]] <- csi_activity_matrix_list_1
  }
  for(ct in names(csi_activity_matrix_list)){
    for(cluster in unique(clusters_df$csi_module)){
      csi_regulon <- subset(clusters_df,csi_module == cluster)
      csi_regulon_activtiy <- subset(csi_activity_matrix_list[[ct]],regulon %in% csi_regulon$regulon)
      csi_activtiy_mean <- mean(csi_regulon_activtiy$activtiy)
      this_cluster_ct_activity <- data.frame("csi_module" = cluster,
                                             "mean_activity" = csi_activtiy_mean,
                                             "cell_type" = ct)
      csi_cluster_activity <- rbind(csi_cluster_activity,this_cluster_ct_activity)
    }
  }

  csi_cluster_activity[is.na(csi_cluster_activity)] <- 0

  csi_cluster_activity_wide <- csi_cluster_activity %>%
    spread(cell_type,mean_activity)

  rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_module
  csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])

  return(csi_cluster_activity_wide)
}

#' Plots a heatmap for the connection specificty indices for all regulons.
#'
#' @param csi_df Data frame containing CSI values for all pairwise regulons.
#' @param nclust Number of clusters to divide the heatmap into
#' @param font_size_regulons Font size for regulon names.
#' @keywords SCENIC, regulons, CSI
#' @export
#' @examples
#' @import tidyr
#' @import pheatmap
#' @import viridis
plot_csi_modules <- function(
  csi_df,
  nclust = 10,
  row_anno = NULL,
  font_size_regulons = 6
){
  ## subset csi data frame based on threshold
  csi_test_mat <- csi_df %>% tidyr::spread(regulon_2,CSI)

  future_rownames <- csi_test_mat$regulon_1
  csi_test_mat <- as.matrix(csi_test_mat[,2:ncol(csi_test_mat)])
  rownames(csi_test_mat) <- future_rownames

  color_map = SelectColors(object=NULL, n = nclust, palette = "ditto")
  names(color_map) <- 1:nclust
  color_use = list()
  #for (i in names(color_map)){
  #    color_use[["csi_module"]][[i]]=color_map[[i]]
  #}
  color_use[[opt$groupby]] <- color_map  ##different with local script 
  
  pheatmap::pheatmap(csi_test_mat,
           show_colnames = TRUE,
           #border_color = NA,
           color = viridis::viridis(n = 10),
           fontsize_row = font_size_regulons,
           fontsize_col = font_size_regulons,
           angle_col = 90,
           cutree_cols = nclust,
           cutree_rows = nclust,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           annotation_row = row_anno,
           annotation_colors = color_use,
           treeheight_row = 20,
           treeheight_col = 20,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           width = 2000,
           height = 3200)
}



#================================================================================
# package loading
#================================================================================
suppressWarnings({
    suppressPackageStartupMessages(library("optparse"))
    #suppressPackageStartupMessages(library("OESingleCell"))
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("viridis"))
    suppressPackageStartupMessages(library("Seurat"))
})

#=command line parameters setting=============================
option_list = list(
make_option( c("--input", "-i"), type = "character",
    help = "the seurat object saved as R object in RDS format." ),
make_option( c("--auc", "-v"), type = "character", 
    help = "The regulon activity scores results from the SCENIC output,
		 usually the 3.4_regulonAUC.Rds from SCENIC" ),
make_option( c("--aucformat", "-f" ), type = "character", default = "rds",
    help = "The indication of type of input AUC, the possible type can be:
              rds: the aucellResults object in RDS format,
	      xsv: the delimited table from AUC.
	      [default \"%default\"]"),
make_option( c("--groupby", "-c"), type = "character",
    help = "the groupping column of cells in the metadata of seurat object." ),
make_option( c("--nclust", "-n" ), type="integer", default = 5,
    help="the number of csi modules to use [default \"%default\"]"),
make_option( c("--extended"),type="logical", default = FALSE,
    help="whether to use the extended regulons for calculation and visualization. [default \"%default\"]"),
make_option( c("--predicate" ), type = "character", default = NULL, 
             help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell."),
make_option( c("--output", "-o"), type = "character", default = "./",
    help = "the output directory of QC results.[default \"%default\"]", metavar = "outdir" )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if ( opt$extended){
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

seurat_ob = readRDS( opt$input)

## subset 
if (!is.null(opt$predicate)) {
    df <- seurat_ob@meta.data
    desired_cells <- subset(df, eval(parse(text = opt$predicate)))
    seurat_ob <- seurat_ob[, rownames(desired_cells)]
}

if ( seurat_ob@version < 3 ){
    seurat_ob = UpdateSeuratObject(seurat_ob)
}

if ( !is.null(opt[["auc"]])){ # using the user supplied regulon activity score matrix
    if ( tolower(opt$aucformat) == "rds" ){
        auc = readRDS(opt$auc)
        seurat_ob=seurat_ob[,colnames(auc)]
        auc = auc[,rownames(seurat_ob@meta.data)]
        auc = auc[which(rowSums(getAUC(auc))!=0),]
        auc = AUCell::getAUC(auc)
    }else if ( tolower(opt$aucformat) == "xsv" ){
        auc = vroom::vroom(opt$auc, col_names = T) %>% as.data.frame()
        rownames(auc) = auc[,1]
        auc = as.matrix(auc[,-1])
    }
}else {
    auc = Tool(seurat_ob, slot = "RunAUCell")
    auc = auc[,rownames(seurat_ob@meta.data)]
    auc = auc[which(rowSums(getAUC(auc))!=0),]
    auc = AUCell::getAUC(auc)
    #stop("NO regulon activity results supplied!")
}

regulons_csi <- calculate_csi(auc,calc_extended = non_extended)

csi_csi_wide <- regulons_csi %>% tidyr::spread(regulon_2,CSI)
future_rownames <- csi_csi_wide$regulon_1
csi_csi_wide <- as.matrix(csi_csi_wide[,2:ncol(csi_csi_wide)])
rownames(csi_csi_wide) <- future_rownames
regulons_hclust <- hclust(dist(csi_csi_wide,method = "euclidean"))

clusters <- cutree(regulons_hclust,k= opt$nclust)
clusters_df <- data.frame("regulon" = names(clusters),"csi_module" = clusters)
write.table(clusters_df, file.path(output_dir, "3.1.csi_module_annotation.xls"),
            sep = "\t", col.names = T, row.names =F, quote =F)

row_anno <- clusters_df
rownames(row_anno) <- NULL
row_anno$csi_module <- as.character(row_anno$csi_module)
plot = plot_csi_modules(regulons_csi,
                  nclust = opt$nclust,
                  row_anno = row_anno %>% column_to_rownames(var="regulon"),
                  font_size_regulons = 8)
save_ggplots(file.path(output_dir, paste0("3.2.regulons_csi_correlation_heatmap", collapse = "")),
    plot = plot,
    width = length(unique(regulons_csi[,1])) * 0.15+2,
    height = length(unique(regulons_csi[,1])) * 0.15,
    dpi = 1000,
    limitsize = F
)

cellinfo = seurat_ob@meta.data
cellinfo[,"cell_type"] = cellinfo[,opt$groupby]
csi_cluster_activity_wide <- calc_csi_module_activity(clusters_df,auc, cellinfo)

plot = pheatmap::pheatmap(csi_cluster_activity_wide,
                    show_colnames = TRUE,
                    color = viridis::viridis(n = 10),
                    cellwidth = 24,
                    cellheight = 24,
                    cluster_cols = TRUE,
                    cluster_rows = TRUE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean")
save_ggplots(file.path(output_dir, paste0("3.3.csi_module_activity_heatmap", collapse = "")),
    plot = plot,
    width = 8,
    height = 8,
    dpi = 1000,
    limitsize = F
)
if(!file.exists(file.path(output_dir, "Regulon调控子分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/Regulon调控子分析说明.docx",
  file.path(output_dir, "Regulon调控子分析说明.docx"))
}
