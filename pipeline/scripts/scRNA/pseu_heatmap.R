suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("monocle"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("viridis"))

opt=list()
opt$input="/public/scRNA/works/yangfan/2-7-HT2021-13068-14903/1.Fibro_clusters_2346_pseudotime_results/pseudotime_results.rds"
setwd("/public/scRNA/works/liuhongyan/Project/tamp/HT2020-14209-test/heatmap_test")
gbm_cds <- readRDS(opt$input)

genes <- as.factor(subset(gbm_cds@featureData@data, use_for_ordering == TRUE)$gene_short_name)
to_be_tested <- row.names(subset(fData(gbm_cds), gene_short_name %in% levels(genes)))



CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}


plot_genes_branched_heatmap <- function(cds_subset, branch_point = 1, branch_states = NULL,
                                        branch_labels = c("Cell fate 1", "Cell fate 2"), cluster_rows = TRUE,
                                        hclust_method = "ward.D2", num_clusters = 6, hmcols = NULL,
                                        branch_colors = c("#979797", "#F05662", "#7990C8"), add_annotation_row = NULL,
                                        add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE,
                                        scale_max = 3, scale_min = -3, norm_method = c("log", "vstExprs"),
                                        trend_formula = "~sm.ns(Pseudotime, df=3) * Branch", return_heatmap = FALSE,
                                        cores = 1, ...) {
    cds <- NA
    new_cds <- buildBranchCellDataSet(cds_subset,
        branch_states = branch_states,
        branch_point = branch_point, progenitor_method = "duplicate",
        ...
    )
    new_cds@dispFitInfo <- cds_subset@dispFitInfo
    if (is.null(branch_states)) {
        progenitor_state <- subset(pData(cds_subset), Pseudotime ==
            0)[, "State"]
        branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
    }
    col_gap_ind <- 101
    newdataA <- data.frame(
        Pseudotime = seq(0, 100, length.out = 100),
        Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1])
    )
    newdataB <- data.frame(
        Pseudotime = seq(0, 100, length.out = 100),
        Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2])
    )
    BranchAB_exprs <- genSmoothCurves(new_cds[, ],
        cores = cores,
        trend_formula = trend_formula, relative_expr = T, new_data = rbind(
            newdataA,
            newdataB
        )
    )
    BranchA_exprs <- BranchAB_exprs[, 1:100]
    BranchB_exprs <- BranchAB_exprs[, 101:200]
    common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
        setdiff(pData(new_cds)$State, branch_states), ])
    BranchP_num <- (100 - floor(max(pData(new_cds)[
        common_ancestor_cells,
        "Pseudotime"
    ])))
    BranchA_num <- floor(max(pData(new_cds)[
        common_ancestor_cells,
        "Pseudotime"
    ]))
    BranchB_num <- BranchA_num
    norm_method <- match.arg(norm_method)
    if (norm_method == "vstExprs") {
        BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
        BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
    } else if (norm_method == "log") {
        BranchA_exprs <- log10(BranchA_exprs + 1)
        BranchB_exprs <- log10(BranchB_exprs + 1)
    }
    heatmap_matrix <- cbind(
        BranchA_exprs[, (col_gap_ind - 1):1],
        BranchB_exprs
    )
    heatmap_matrix <- heatmap_matrix[!apply(
        heatmap_matrix, 1,
        sd
    ) == 0, ]
    heatmap_matrix <- Matrix::t(scale(Matrix::t(heatmap_matrix),
        center = TRUE
    ))
    heatmap_matrix <- heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
        FALSE, ]
    heatmap_matrix[is.nan(heatmap_matrix)] <- 0
    heatmap_matrix[heatmap_matrix > scale_max] <- scale_max
    heatmap_matrix[heatmap_matrix < scale_min] <- scale_min
    heatmap_matrix_ori <- heatmap_matrix
    heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[
        ,
        1
    ]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix))) / 2)
    row_dist[is.na(row_dist)] <- 1
    exp_rng <- range(heatmap_matrix)
    bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
    
    if (is.null(hmcols)) {
        hmcols <- monocle:::blue2green2red(length(bks) - 1)
    }else{
        hmcols=colorRampPalette(c("#406AA8", "white", "#D91216"))(length(bks) - 1)
    }
    ph <- pheatmap(heatmap_matrix,
        useRaster = T, cluster_cols = FALSE,
        cluster_rows = TRUE, show_rownames = F, show_colnames = F,
        clustering_distance_rows = row_dist, clustering_method = hclust_method,
        cutree_rows = num_clusters, silent = TRUE, filename = NA,
        breaks = bks, color=hmcols
    )
    annotation_row <- data.frame(Module = factor(cutree(
        ph$tree_row,
        num_clusters
    )))
    if (!is.null(add_annotation_row)) {
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
    }
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    annotation_col <- data.frame(
        row.names = c(1:ncol(heatmap_matrix)),
        `Cell Type` = c(rep(branch_labels[1], BranchA_num), rep(
            "Pre-branch",
            2 * BranchP_num
        ), rep(branch_labels[2], BranchB_num))
    )
    colnames(annotation_col) <- "Cell Type"
    if (!is.null(add_annotation_col)) {
        annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), ])$gene_short_name, 1])
    }
    names(branch_colors) <- c(
        "Pre-branch", branch_labels[1],
        branch_labels[2]
    )
    #module_colors <- scales::dscale(factor(1:num_clusters), scales::hue_pal(l = 75))
    module_colors <- CustomCol2(1:num_clusters)
    names(module_colors) <- 1:num_clusters
    annotation_colors <- list(`Cell Type` = branch_colors, Module = module_colors)
    names(annotation_colors$`Cell Type`) <- c("Pre-branch", branch_labels)
    if (use_gene_short_name == TRUE) {
        if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
            feature_label <- as.character(fData(cds_subset)[
                row.names(heatmap_matrix),
                "gene_short_name"
            ])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
            row_ann_labels <- as.character(fData(cds_subset)[
                row.names(annotation_row),
                "gene_short_name"
            ])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        } else {
            feature_label <- row.names(heatmap_matrix)
            row_ann_labels <- row.names(annotation_row)
        }
    } else {
        feature_label <- row.names(heatmap_matrix)
        row_ann_labels <- row.names(annotation_row)
    }
    row.names(heatmap_matrix) <- feature_label
    row.names(annotation_row) <- row_ann_labels
    ph_res <- pheatmap(heatmap_matrix[, ],
        useRaster = T, cluster_cols = FALSE,
        cluster_rows = TRUE, show_rownames = show_rownames, show_colnames = F,
        clustering_distance_rows = row_dist, clustering_method = hclust_method,
        cutree_rows = num_clusters, annotation_row = annotation_row,
        annotation_col = annotation_col, annotation_colors = annotation_colors,
        gaps_col = col_gap_ind, treeheight_row = 20, breaks = bks,
        fontsize = 6, color = hmcols, border_color = NA, silent = TRUE
    )
    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    if (return_heatmap) {
        return(list(
            BranchA_exprs = BranchA_exprs, BranchB_exprs = BranchB_exprs,
            heatmap_matrix = heatmap_matrix, heatmap_matrix_ori = heatmap_matrix_ori,
            ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist,
            hmcols = hmcols, annotation_colors = annotation_colors, bks = bks,
            annotation_row = annotation_row, annotation_col = annotation_col,
            ph_res = ph_res
        ))
    }
}


# change heatmap annotation from Cluster to Module
plot_pseudotime_heatmap <- function(cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2",
                                    num_clusters = 6, hmcols = NULL, add_annotation_row = NULL,
                                    add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE,
                                    norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3,
                                    trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
                                    cores = 1) {
    num_clusters <- min(num_clusters, nrow(cds_subset))
    pseudocount <- 1
    newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime),
        max(pData(cds_subset)$Pseudotime),
        length.out = 100
    ))
    m <- genSmoothCurves(cds_subset,
        cores = cores, trend_formula = trend_formula,
        relative_expr = T, new_data = newdata
    )
    m <- m[!apply(m, 1, sum) == 0, ]
    norm_method <- match.arg(norm_method)
    if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) ==
        FALSE) {
        m <- vstExprs(cds_subset, expr_matrix = m)
    } else if (norm_method == "log") {
        m <- log10(m + pseudocount)
    }
    m <- m[!apply(m, 1, sd) == 0, ]
    m <- Matrix::t(scale(Matrix::t(m), center = TRUE))
    m <- m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] <- 0
    m[m > scale_max] <- scale_max
    m[m < scale_min] <- scale_min
    heatmap_matrix <- m
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix))) / 2)
    row_dist[is.na(row_dist)] <- 1
    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- monocle:::blue2green2red(length(bks) - 1)
    } else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    ph <- pheatmap(heatmap_matrix,
        useRaster = T, cluster_cols = FALSE,
        cluster_rows = cluster_rows, show_rownames = F, show_colnames = F,
        clustering_distance_rows = row_dist, clustering_method = hclust_method,
        cutree_rows = num_clusters, silent = TRUE, filename = NA,
        breaks = bks, border_color = NA, color = hmcols
    )
    if (cluster_rows) {
        annotation_row <- data.frame(Module = factor(cutree(
            ph$tree_row, #
            num_clusters
        )))
    } else {
        annotation_row <- NULL
    }
    if (!is.null(add_annotation_row)) {
        old_colnames_length <- ncol(annotation_row)
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
        colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
    }
    if (!is.null(add_annotation_col)) {
        if (nrow(add_annotation_col) != 100) {
            stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
        }
        annotation_col <- add_annotation_col
    } else {
        annotation_col <- NA
    }
    if (use_gene_short_name == TRUE) {
        if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
            feature_label <- as.character(fData(cds_subset)[
                row.names(heatmap_matrix),
                "gene_short_name"
            ])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
            row_ann_labels <- as.character(fData(cds_subset)[
                row.names(annotation_row),
                "gene_short_name"
            ])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        } else {
            feature_label <- row.names(heatmap_matrix)
            row_ann_labels <- row.names(annotation_row)
        }
    } else {
        feature_label <- row.names(heatmap_matrix)
        if (!is.null(annotation_row)) {
              row_ann_labels <- row.names(annotation_row)
          }
    }
    row.names(heatmap_matrix) <- feature_label
    if (!is.null(annotation_row)) {
          row.names(annotation_row) <- row_ann_labels
      }
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    ph_res <- pheatmap(heatmap_matrix[, ],
        useRaster = T, cluster_cols = FALSE,
        cluster_rows = cluster_rows, show_rownames = show_rownames,
        show_colnames = F, clustering_distance_rows = row_dist,
        clustering_method = hclust_method, cutree_rows = num_clusters,
        annotation_row = annotation_row, annotation_col = annotation_col,
        treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols,
        border_color = NA, silent = TRUE, filename = NA
    )
    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    if (return_heatmap) {
        return(list(
            heatmap_matrix = heatmap_matrix, ph = ph, row_dist = row_dist,
            hmcols = hmcols, annotation_row = annotation_row, bks = bks,
            annotation_col = annotation_col, ph_res = ph_res
        ))
    }
}



opt$branchpoint=1
opt$CORES=10
clusters_num=4
showname=TRUE
branchpoint <- opt$branchpoint

gbm_cds <- gbm_cds[to_be_tested, ]

new_cds <- buildBranchCellDataSet(gbm_cds, branch_point = branchpoint, progenitor_method = "duplicate")

# cell_fate1
cell_fate1 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]), ]$State)
# cell_fate2
cell_fate2 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]), ]$State)
branch_labels <- c(paste("State", paste(sort(setdiff(cell_fate1, cell_fate2)), collapse = "-")), paste("State", paste(sort(setdiff(cell_fate2, cell_fate1)), collapse = "-")))

branch_heatmap <- plot_genes_branched_heatmap(gbm_cds, branch_point = branchpoint, cores = opt$CORES, cluster_rows = T, num_clusters = clusters_num, use_gene_short_name = T, show_rownames = showname, branch_labels = branch_labels, return_heatmap = T, hmcols="seurat")

output_dir="./"
#source("add.flag_branch.r")
source("add.flag.r")
library(grid)
## gene_name<-sample(rownames(gbm_cds),8)

gene_name<-head(rownames(gbm_cds),13)[-11]
# p1<-pheatmap(branch_heatmap)
p1 <- add.flag(branch_heatmap$ph_res,
         kept.labels = gene_name,
         repel.degree = 0)

ggsave(file.path(output_dir, "pseudotime_heatmap_branchtime_ori.pdf"), plot = branch_heatmap$ph_res)
ggsave(file.path(output_dir, "pseudotime_heatmap_branchtime_ori.png"), plot = branch_heatmap$ph_res, dpi = 1000)

# heatmap1 <- branch_heatmap$ph_res$gtable
# new.label <- heatmap1$grobs[[which(heatmap1$layout$name == "row_names")]] 
# kept.labels=gene_name
# new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")
# new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")


## branch_heatmap$ph_res
ggsave(file.path(output_dir, "pseudotime_heatmap_branchtime.pdf"), plot = p1)
ggsave(file.path(output_dir, "pseudotime_heatmap_branchtime.png"), plot = p1, dpi = 1000)

gene_clusters <- cutree(branch_heatmap$ph_res$tree_row, k = clusters_num)


### palette <- colorRampPalette(c("#406AA8", "white", "#D91216"))(n=299)
### regular ###
heatmap <- plot_pseudotime_heatmap(gbm_cds, cores = opt$CORES, cluster_rows = T, num_clusters = clusters_num, show_rownames = TRUE, return_heatmap = T)

output_dir="./"
source("add.flag.r")
library(grid)
gene_name<-sample(rownames(gbm_cds),10)
# p1<-pheatmap(branch_heatmap)
p1 <- add.flag(heatmap$ph_res,
         kept.labels = gene_name,
         repel.degree = 1)


            ggsave(file.path(output_dir, "pseudotime_heatmap.pdf"), plot = p1)
            ggsave(file.path(output_dir, "pseudotime_heatmap.png"), plot = p1, dpi = 1000)


source("add.flag.r")
heatmap1 <- heatmap$gtable
new.label <- heatmap1$grobs[[which(heatmap1$layout$name == "row_names")]] 
kept.labels=gene_name
new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")
new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")



# source("useful_R_function/add_flag.r")
library(grid)
gene_name<-sample(rownames(gbm_cds),4)
# p1<-pheatmap(branch_heatmap)
p1 <- add.flag(branch_heatmap,
         kept.labels = gene_name,
         repel.degree = 0.2)


heatmap <- branch_heatmap$gtable


add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),to = min(selected.range) - k*(min(selected.range) - min(full.range)), length.out = sum(d.select)), "npc"))
  }
  #new.y.positions <- repelled.y(new.label$y,
  #                              d.select = new.label$label != "")
  new.y.positions <-new.label$y[which(new.label$label != "")]  #基因标签位置,可能会有重叠,不知道该如何优化。。
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}


