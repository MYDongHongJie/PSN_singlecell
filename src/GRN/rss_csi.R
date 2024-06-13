#!/usr/bin/env Rscript
# Processing and visualization of SCENIC results
#

#========import packages=====================================
rm(list=ls())
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scFunctions"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("ggridges"))

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i"), type = "character",
                 help = "The regulon activity scores in matrix format that from the SCENIC output file 3.4_regulonAUC.Rds." ),
    make_option( c("--RDS", "-v"), type = "character", default = "TRUE",
                 help = "the seurat object saved as R object in RDS format." ),
    make_option( c("--topGenes", "-t"), type = "integer", default = 5,
                 help = "top_genes Number of top genes to label in plot of rrs ranking." ),
    make_option( c("--celltype", "-c"), type = "character",
                 help = "the column named cell_type that assigns groupings to cells in the metadata of seurat object." ),
    make_option( c("--nclust", "-r"), type = "integer", default = 10,
                 help = "the number of regulon clusters in the csi heatmap." ),
    make_option( c("--output", "-o"), type = "character", default = "./",
                 help = "the output directory of QC results.", metavar = "outdir" )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
# parse the command line parameters
#=================================================================================
#if ( is.null(opt$RDS) ){
 #   stop("the seurat object is NOT AVAILABLE!")
#}else{
#    seurat_ob = readRDS(opt$RDS)
#    seurat_ob = UpdataSeuratObject(seurat_ob)
#}

if ( is.null(opt$topGenes) ){
    topGenes = 6
}else{
    topGenes = opt$topGenes + 1
}

if ( is.null(opt$celltype) ){
    print("NO column of cell type is available!")
}else{
    celltype = opt$celltype
}

if ( is.null(opt$nclust) ){
    nclust = 10
}else{
    nclust = opt$nclust
}

if ( is.null(opt$output) ){
    print("NO output directory specified, the current directory will be used!")
    output_dir = getwd()
}else{
    if (file.exists(opt$output)){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir)
    }
}

#=================================================================================
# Regulon analysis
#=================================================================================

if ( !is.null(opt$input) ){

    regulonAUC = readRDS(opt$input)
    metadata_sub = readRDS(opt$RDS)
#    metadata_sub = seurat_ob@meta.data

    # Binary regulons
    kmeans_thresholds = auc_thresh_kmeans(regulonAUC)
    binary_regulons = binarize_regulons(regulonAUC,kmeans_thresholds)

    joined_bin_reg <- binary_regulons %>% reduce(left_join,by="cells")
    rownames(joined_bin_reg) <- joined_bin_reg$cells
    joined_bin_reg <- joined_bin_reg[2:ncol(joined_bin_reg)]
    binary_regulons_trans <- as.matrix(t(joined_bin_reg))

    # Calculate Regulon specificity score (RSS) from binary regulon activity.
    rrs_df <- calculate_rrs(metadata_sub,
              binary_regulons = binary_regulons_trans, cell_type_column = celltype)
    rrs_df_wide = rrs_df %>% spread(cell_type,RSS)
    plot_rrs_ranking = plot_rrs_ranking(rrs_df, "RMPH", ggrepel_force = 1, ggrepel_point_padding = 0.2, top_genes = topGenes, plot_extended = FALSE)
    ggsave(file.path(output_dir, "rss_ranking.pdf"), plot = plot_rrs_ranking)
    ggsave(file.path(output_dir, "rss_ranking.png"), dpi = 1000, plot = plot_rrs_ranking)

    rrs_df_nona <- subset(rrs_df,RSS > 0)
    rrs_df_nona_plot <- ggplot(rrs_df_nona,aes(RSS,cell_type, fill = cell_type)) +
        geom_density_ridges(scale = 5, alpha = 0.75) +
        geom_vline(xintercept = 0.1) +
        theme(legend.position = "none")
    ggsave(file.path(output_dir, "rss_df_nona.pdf"), plot = rrs_df_nona_plot)
    ggsave(file.path(output_dir, "rss_df_nona.png"), dpi = 1000, plot = rrs_df_nona_plot)

    rrs_df_wide = rrs_df %>% spread(cell_type,RSS)
    rownames(rrs_df_wide) = rrs_df_wide$regulon
    rrs_df_wide = rrs_df_wide[,2:ncol(rrs_df_wide)]
    rrs_df_wide_specific = rrs_df_wide[apply(rrs_df_wide,MARGIN = 1 ,FUN =  function(x) any(x > 0.4)),]
    pdf(file.path(output_dir, "rss_heatmap.pdf"), width = 8)
    pheatmap::pheatmap(rrs_df_wide_specific,
                       treeheight_row=10, treeheight_col=10,
                       border_color=NA,
                       color = viridis(n = 10),
                       fontsize_row = 6, angle_col = 45)
    dev.off()
    png(file.path(output_dir, "rss_heatmap.png"))
    pheatmap::pheatmap(rrs_df_wide_specific,
                       treeheight_row=10, treeheight_col=10,
                       border_color=NA,
                       color = viridis(n = 10),
                       fontsize_row = 6, angle_col = 45)
    dev.off()

    # Calculate connection specificity index (CSI) for all regulons
    regulons_csi = calculate_csi(regulonAUC, calc_extended = FALSE)
    csi_test_mat = regulons_csi %>% spread(regulon_2,CSI)
    future_rownames = csi_test_mat$regulon_1
    csi_test_mat = as.matrix(csi_test_mat[,2:ncol(csi_test_mat)])
    rownames(csi_test_mat) = future_rownames
    pdf(file.path(output_dir, "csi_heatmap.pdf"), width = 8)
    pheatmap::pheatmap(csi_test_mat, show_colnames = FALSE,
                       color = viridis(n=10),
                       border_color = NA,
                       cutree_cols = 10, cutree_rows = 10,
                       fontsize_row = 4,
                       cluster_cols = TRUE, cluster_rows = TRUE,
                       treeheight_row = nclust,
                       treeheight_col = nclust,
                       clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean",
                       widt = 2000, height = 3200)
    dev.off()
    png(file.path(output_dir, "csi_heatmap.png"))
    pheatmap::pheatmap(csi_test_mat, show_colnames = FALSE,
                       color = viridis(n=10),
                       border_color = NA,
                       cutree_cols = 10, cutree_rows = 10,
                       fontsize_row = 4,
                       cluster_cols = TRUE, cluster_rows = TRUE,
                       treeheight_row = nclust,
                       treeheight_col = nclust,
                       clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean",
                       widt = 2000, height = 3200)
    dev.off()

    csi_csi_wide <- regulons_csi %>% spread(regulon_2,CSI)
    future_rownames <- csi_csi_wide$regulon_1
    csi_csi_wide <- as.matrix(csi_csi_wide[,2:ncol(csi_csi_wide)])
    rownames(csi_csi_wide) <- future_rownames
    regulons_hclust <- hclust(dist(csi_csi_wide,method = "euclidean"))
    clusters <- cutree(regulons_hclust,k= nclust)
    clusters_df <- data.frame("regulon" = names(clusters),
                              "csi_cluster" = clusters)
    #Check how many regulons are in each cluster
    clusters_df_stats <- clusters_df %>%
        group_by(csi_cluster) %>%
        mutate("regulon" = as.character(regulon)) %>%
        tally()
    clusters_df_stats_plot = ggplot(clusters_df_stats,aes(as.factor(csi_cluster),n,fill=as.factor(csi_cluster))) +
        geom_bar(color= "black",stat="identity") +
        theme(legend.position="none") +
        scale_fill_brewer(palette = "Set3") +
        labs(x = "HC clusters", y = "# Regulons")
    ggsave(file.path(output_dir, "clusters_df_stats.pdf"), plot = clusters_df_stats_plot)
    ggsave(file.path(output_dir, "clusters_df_stats.png"), dpi = 1000, plot = clusters_df_stats_plot)

    #Calculate activity scores for each csi module based on the specificity scores of all the regulons in that module.
    csi_cluster_activity_wide <- calc_csi_module_activity(clusters_df,
                                     regulonAUC,
                                     metadata_sub)
    pdf(file.path(output_dir, "csi_cluster_activity.pdf"), width = 8)
    pheatmap(csi_cluster_activity_wide,
             show_colnames = TRUE,
             color = viridis(n = 10),
             angle_col = 45,
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             clustering_distance_cols = "euclidean",
             clustering_distance_rows = "euclidean")
    dev.off()
    png(file.path(output_dir, "csi_cluster_activity.png"))
    pheatmap(csi_cluster_activity_wide,
             show_colnames = TRUE,
             color = viridis(n = 10),
             angle_col = 45,
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             clustering_distance_cols = "euclidean",
             clustering_distance_rows = "euclidean")
    dev.off()
}
