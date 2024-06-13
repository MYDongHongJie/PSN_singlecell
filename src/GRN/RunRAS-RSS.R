#!/usr/bin/env Rscript
# Objective : Processing and visualization of SCENIC results
# Created by: hanmin
# Created on: 2020/5/20
#
rm(list=ls())

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

## RAS
RemoveOutlier <- function(
  metric,
  nmads = 5,
  type = c("both", "lower", "higher"),
  log = FALSE,
  subset = NULL,
  batch = NULL,
  min_diff = NA
) {
    if (log) {
        metric <- log10(metric)
    }
    if (any(is.na(metric))) {
        warning("missing values ignored during outlier detection")
    }

    if (!is.null(batch)) {
        N <- length(metric)
        if (length(batch) != N) {
            stop("length of 'batch' must equal length of 'metric'")
        }

        # Coercing non-NULL subset into a logical vector.
        if (!is.null(subset)) {
            new.subset <- logical(N)
            names(new.subset) <- names(metric)
            new.subset[subset] <- TRUE
            subset <- new.subset
        }

        # Computing QC metrics for each batch.
        by.batch <- split(seq_len(N), batch)
        collected <- logical(N)
        all.threshold <- vector("list", length(by.batch))
        for (b in seq_along(by.batch)) {
            bdx <- by.batch[[b]]
            current <- Recall(metric[bdx], nmads = nmads, type = type, log = FALSE, subset = subset[bdx], batch = NULL, min_diff = min_diff)
            all.threshold[[b]] <- attr(x, "thresholds")
            collected[bdx] <- current
        }

        all.threshold <- do.call(cbind, all.threshold)
        colnames(all.threshold) <- names(by.batch)
        # return(.store_thresholds(collected, all.threshold, logged=log))
        if ( log ){ val <- 10^all.threshold }
        attr(collected, "thresholds") <- val
        return( collected )
    }
    # Computing median/MAD (possibly based on subset of the data).
    if (!is.null(subset)) {
        submetric <- metric[subset]
        if (length(submetric) == 0L) {
            warning("no observations remaining after subsetting")
        }
    } else {
        submetric <- metric
    }
    cur.med <- median(submetric, na.rm = TRUE)
    cur.mad <- mad(submetric, center = cur.med, na.rm = TRUE)

    diff.val <- max(min_diff, nmads * cur.mad, na.rm = TRUE)
    upper.limit <- cur.med + diff.val
    lower.limit <- cur.med - diff.val

    type <- match.arg(type)
    if (type == "lower") {
        upper.limit <- Inf
    } else if (type == "higher") {
        lower.limit <- -Inf
    }

    kx = metric < lower.limit | upper.limit < metric
    val = c(lower=lower.limit, higher=upper.limit)
    if ( log ){
      val <- 10^val
    }
    attr(kx, "thresholds") <- val
    return( kx )

    # .store_thresholds <- function(x, val, logged=FALSE) {
    #     if (logged) val <- 10^val
    #     attr(x, "thresholds") <- val
    #     x
    # }
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
## RSS
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
    #suppressPackageStartupMessages( library("OESingleCell") )
    suppressPackageStartupMessages( library("dplyr") )
    suppressPackageStartupMessages( library("SCENIC") )
    suppressPackageStartupMessages( library("ggplot2") )
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
make_option( c("--predicate" ), type = "character", default = NULL, 
             help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell."),
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

if ( length(levels(seurat_ob@meta.data[[opt$groupby]])) != length(unique(seurat_ob@meta.data[[opt$groupby]])) ) {
    seurat_ob@meta.data[[opt$groupby]] = factor(seurat_ob@meta.data[[opt$groupby]], levels= sort(unique(seurat_ob@meta.data[[opt$groupby]])) )
}

if ( !is.null(opt[["auc"]])){ # using the user supplied regulon activity score matrix
    if ( tolower(opt$aucformat) == "rds" ){
        regulonAUC = readRDS(opt$auc)
        seurat_ob = seurat_ob[,colnames(regulonAUC)]
        regulonAUC = regulonAUC[,rownames(seurat_ob@meta.data)]
        regulonAUC = regulonAUC[which(rowSums(AUCell::getAUC(regulonAUC))!=0),]
    }else if ( tolower(opt$aucformat) == "xsv" ){
        regulonAUC = vroom::vroom(opt$auc, col_names = T) %>% as.data.frame()
        rownames(regulonAUC) = regulonAUC[,1]
        regulonAUC = as.matrix(regulonAUC[,-1])
    }
}else { # use the AUC result integrated in the Seurat object by running RunAUCell
    regulonAUC = Tool(seurat_ob, slot = "RunAUCell")
    regulonAUC = regulonAUC[,rownames(seurat_ob@meta.data)]
    regulonAUC = regulonAUC[which(rowSums(AUCell::getAUC(regulonAUC))!=0),]
}


#####################################
# RAS
#####################################
# matrix prep
regulonAUC= regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC_mat = AUCell::getAUC(regulonAUC)
rownames(regulonAUC_mat) = gsub("_", "-", rownames(regulonAUC_mat))
regulonAUC_mat_out = regulonAUC_mat[grep(pattern="-extended",rownames(regulonAUC_mat),invert=T),]
## groupby 
cellInfo = seurat_ob@meta.data
col_anno = as.data.frame(seurat_ob@meta.data) %>% rownames_to_column(var="barcodes")
if ( dim(cellInfo)[1] > 65536 ){
    col_anno = col_anno[,c("barcodes",opt$groupby)] %>% group_by(get(opt$groupby)) %>% sample_frac(65536/length(col_anno$barcodes))
}else{
    col_anno = col_anno[,c("barcodes",opt$groupby)]
}

col_anno = col_anno %>% arrange(get(opt$groupby)) %>% column_to_rownames(var="barcodes")
regulonAUC_plotdata = regulonAUC_mat_out[,rownames(col_anno)]
bks <- unique(c(seq(-2.5,0, length=100),  seq(0,2.5, length=100)))
color_map = CustomCol2(1:length(unique(col_anno[[opt$groupby]])))
if (is.factor(col_anno[[opt$groupby]])) {
    names <- levels(col_anno[[opt$groupby]])
} else {
    names <- unique(col_anno[[opt$groupby]])
}
names(color_map) <- names

color_use = list()
#for (i in names(color_map)){
#    color_use[[opt$groupby]][[i]]=color_map[[i]]
#}
color_use[[opt$groupby]] <- color_map    #different with  local script

plot = pheatmap::pheatmap( regulonAUC_plotdata,
                    scale = "row",
                    cluster_cols=F,
                    cluster_rows=F,
                    show_colnames= F,
                    color=colorRampPalette(c("#406AA8","white","#D91216"))(200),
                    annotation_col = col_anno,
                    annotation_colors = color_use,
                    treeheight_col=10,
                    border_color=NA,breaks=bks,fontsize_row=6)

save_ggplots(file.path(output_dir, paste0("1.1.regulon_activity_heatmap_groupby_cells", collapse = "")),
	plot = plot,
	width = 8,
	height = 6,
	dpi = 1000,
	limitsize = F
)

# to remove outliers before mean calculation
groupby_1 =c()   # select group with only one cell
for(i in names(table(cellInfo[[opt$groupby]])) ){
	if( table(cellInfo[[opt$groupby]])[i]==1 ){
		groupby_1= i
	}
}
if(is.null(groupby_1)){
	regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo[[opt$groupby]]),function(cells){
		rowmean = apply( regulonAUC_mat_out[,cells], 1, function(x){
			mean(x[!RemoveOutlier( x, nmads = 3, type = "higher")])
		})
	} )
}else{
	cellInfo_rm1 = cellInfo[-which(cellInfo[['group']] == i),]
	cellInfo_rm1[[opt$groupby]] = droplevels(cellInfo_rm1[[opt$groupby]])
	regulonActivity_byCellType0 <- sapply(split(rownames(cellInfo_rm1), cellInfo_rm1[[opt$groupby]]),function(cells){
		rowmean = apply( regulonAUC_mat_out[,cells], 1, function(x){
			mean(x[!RemoveOutlier( x, nmads = 3, type = "higher")])
		})
	} )
	regulonAUC_mat_out_1 = regulonAUC_mat_out[,rownames(cellInfo[which(cellInfo[[opt$groupby]] == groupby_1),] )]
	regulonActivity_byCellType <- cbind(regulonActivity_byCellType0, regulonAUC_mat_out_1)
	colnames(regulonActivity_byCellType) <- c(colnames(regulonActivity_byCellType0), groupby_1)
}


# caculate the mean activity of each regulon in each group and plot heatmap
regulonActivity_byCellType_processed <- t(scale(t(regulonActivity_byCellType), center = T, scale=ifelse(length(unique(cellInfo[[opt$groupby]])) ==  2,F,T)))
regulonActivity_byCellType_processed <- na.omit(regulonActivity_byCellType_processed)
#regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[which(rowSums(regulonActivity_byCellType_processed) != 0),]
df = as.data.frame(regulonActivity_byCellType_processed) %>% tibble::rownames_to_column(var = "regulon")
write.table(df,file.path(output_dir, "1.2.centered_regulon_activity_groupby_design.xls"),
            sep = "\t", col.names =T, row.names =F, quote =F)

# regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[which(rowSums(abs(regulonActivity_byCellType_processed)) != 0),]
regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[which(rowSums(abs(regulonActivity_byCellType_processed)) > max(abs(regulonActivity_byCellType_processed))/4),]
if (dim(regulonActivity_byCellType_processed)[1]<11) {
    hig <- 4.5
} else {
    hig <- dim(regulonActivity_byCellType_processed)[1]*0.4
}
plot = pheatmap::pheatmap( regulonActivity_byCellType_processed,
                    # scale = "row",
                    cellwidth = 18,
                    cellheight = 18,
                    color=colorRampPalette(c("#406AA8","white","#D91216"))(299),
                    # annotation_col = factor(cellInfo[[opt$groupby]]),
                    angle_col = 45,
                    treeheight_col=20, 
                    treeheight_row=20,
                    border_color=NA)
save_ggplots(file.path(output_dir, paste0("1.3.regulon_activity_heatmap", collapse = "")),
	plot = plot,
	width = dim(regulonActivity_byCellType_processed)[2]/2+3,
	height = hig,
	dpi = 1000,
	limitsize = F
)

#####################################
# RSS
#####################################
seurat_ob = BinaryCount(seurat_ob, method = opt$binmethod, auc = regulonAUC, nCores = opt$ncores)
#if (!is.null(opt[["auc"]])) saveRDSMC(seurat_ob, opt$input)
rss_df = RunRSS(seurat_ob, group.by = opt$groupby)
rss_df_out = rss_df %>% subset(!grepl("extended",regulon))
write.table(rss_df_out, file.path(output_dir, "2.1.regulon_RSS_annotation.xls"),
            sep = "\t", col.names = T, row.names =F, quote =F)
gg_rss_rank = RSSRanking(rss_df, group.by = opt$groupby, top_genes = topGenes, plot_extended = non_extended )
save_ggplots(file.path(output_dir, paste0("2.2.RSS_ranking_plot", collapse = "")),
	plot = gg_rss_rank,
	width = 6,
	height = ceiling(length(unique(rss_df[,opt$groupby]))/2) * 4,
	dpi = 1000,
	limitsize = F
)

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
plot = pheatmap::pheatmap( rss_df_wide_specific,
                    cellwidth = 18,
                    cellheight = 18,
                    color=colorRampPalette(c("#406AA8","white","#D91216"))(299),
                    angle_col = 45,
                    treeheight_col=20,
                    treeheight_row=20, 
                    border_color=NA)

save_ggplots(file.path(output_dir, paste0("2.3.RSS_heatmap", collapse = "")),
	plot = plot,
	width = dim(rss_df_wide_specific)[2]/2+3,
	height = dim(rss_df_wide_specific)[1]*0.3,
	dpi = 1000,
	limitsize = F
)
