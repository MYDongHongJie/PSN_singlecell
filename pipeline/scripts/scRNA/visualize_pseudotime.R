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
library("ggnewscale") ## 增加 2 种颜色标尺 
#===========================function definition=============================
# expressplot line funtion
plot_genes_in_pseudotime_line <- function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL,
    ncol = 1, panel_order = NULL, color_by = "State", trend_formula = paste0("~sm.ns(Pseudotime, df=3) * ", opt$groupby),
    label_by_short_name = TRUE, relative_expr = TRUE, vertical_jitter = NULL,
    horizontal_jitter = NULL){
    f_id <- NA
    Cell <- NA
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial",
        "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- cds_exprs$expression
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores = 1,
        trend_formula = trend_formula, relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame(expectation = model_expectation[x$f_id,
        x$Cell]))
    cds_exprs <- merge(cds_exprs, expectation)
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label,
            levels = panel_order)
    }
     q <- ggplot(aes(Pseudotime, expression), data = cds_exprs) +
 					  scale_y_log10() +
                      geom_point(aes_string(color = "Pseudotime"), size = I(cell_size), position = position_jitter(horizontal_jitter, vertical_jitter)) + 
                      scale_colour_viridis_c(option = "inferno") +
                      ##scale_fill_viridis_c(option = "inferno") +
                      new_scale("color") + 
 					  geom_line(aes_string(x = 'Pseudotime', y = 'expectation',color = color_by),size = 1, data = cds_exprs)+
                      expand_limits(y = min_expr)+
 					  facet_wrap(~feature_label, nrow = nrow, ncol = ncol, scales = "free_y")+
                      ylab("Expression") + xlab("Pseudotime") +
                      theme_bw() +
                      theme(panel.grid=element_blank(), panel.border=element_blank(),axis.line = element_line(size=1,colour="black"))
    q
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
    }
    ph <- pheatmap(heatmap_matrix,
        useRaster = T, cluster_cols = FALSE,
        cluster_rows = TRUE, show_rownames = F, show_colnames = F,
        clustering_distance_rows = row_dist, clustering_method = hclust_method,
        cutree_rows = num_clusters, silent = TRUE, filename = NA,
        breaks = bks, color = hmcols
    )
    annotation_row <- data.frame(Module = factor(cutree(
        ph$tree_row, #
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
    module_colors <- scales::dscale(factor(1:num_clusters), scales::hue_pal(l = 75))
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

# change width of the pseudotime evolutionary tree from 0.1 to 0.05
# adding param width = 0.05 to geom_jitter
plot_complex_cell_trajectory <- function(cds,
                                         x = 1,
                                         y = 2,
                                         root_states = NULL,
                                         color_by = "State",
                                         show_tree = TRUE,
                                         show_backbone = TRUE,
                                         backbone_color = "black",
                                         markers = NULL,
                                         show_cell_names = FALSE,
                                         cell_size = 1.5,
                                         cell_link_size = 0.75,
                                         cell_name_size = 2,
                                         show_branch_points = TRUE,
                                         ...) {
    gene_short_name <- NA
    sample_name <- NA
    data_dim_1 <- NA
    data_dim_2 <- NA

    # TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
    lib_info_with_pseudo <- pData(cds)

    if (is.null(cds@dim_reduce_type)) {
        stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
    }

    if (cds@dim_reduce_type == "ICA") {
        reduced_dim_coords <- reducedDimS(cds)
    } else if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", "SGL-tree")) {
        reduced_dim_coords <- reducedDimK(cds)
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    } else {
        stop("Error: unrecognized dimensionality reduction method.")
    }

    if (is.null(reduced_dim_coords)) {
        stop("You must first call reduceDimension() before using this function")
    }

    dp_mst <- minSpanningTree(cds)


    if (is.null(root_states)) {
        if (is.null(lib_info_with_pseudo$Pseudotime)) {
            root_cell <- row.names(lib_info_with_pseudo)[degree(dp_mst) == 1][1]
        } else {
              root_cell <- row.names(subset(lib_info_with_pseudo, Pseudotime == 0))
          }

        if (cds@dim_reduce_type != "ICA") {
              root_cell <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_cell, ]]
          }
    } else {
        candidate_root_cells <- row.names(subset(pData(cds), State %in% root_states))
        if (cds@dim_reduce_type == "ICA") {
            root_cell <- candidate_root_cells[which(degree(dp_mst, candidate_root_cells) == 1)]
        } else {
            Y_candidate_root_cells <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[candidate_root_cells, ]]
            root_cell <- Y_candidate_root_cells[which(degree(dp_mst, Y_candidate_root_cells) == 1)]
        }
    }

    tree_coords <- layout_as_tree(dp_mst, root = root_cell)

    ica_space_df <- data.frame(tree_coords)
    row.names(ica_space_df) <- colnames(reduced_dim_coords)
    colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")

    ica_space_df$sample_name <- row.names(ica_space_df)


    if (is.null(dp_mst)) {
        stop("You must first call orderCells() before using this function")
    }

    edge_list <- as.data.frame(get.edgelist(dp_mst))
    colnames(edge_list) <- c("source", "target")

    edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", by.y = "source", all = TRUE)
    edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1" = "source_prin_graph_dim_1", "prin_graph_dim_2" = "source_prin_graph_dim_2"))
    edge_df <- merge(edge_df, ica_space_df[, c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], by.x = "target", by.y = "sample_name", all = TRUE)
    edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1" = "target_prin_graph_dim_1", "prin_graph_dim_2" = "target_prin_graph_dim_2"))

    if (cds@dim_reduce_type == "ICA") {
        S_matrix <- tree_coords[, ] # colnames(cds)
    } else if (cds@dim_reduce_type %in% c("DDRTree", "SimplePPT", "SGL-tree")) {
        S_matrix <- tree_coords[closest_vertex, ]
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    }

    data_df <- data.frame(S_matrix)
    row.names(data_df) <- colnames(reducedDimS(cds))
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- row.names(data_df)
    data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "row.names")

    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in% markers)
        if (nrow(markers_fData) >= 1) {
            markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData), ])))
            colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
            markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y = "row.names")
            # print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
        data_df <- merge(data_df, markers_exprs, by.x = "sample_name", by.y = "cell_id")
        # print (head(edge_df))
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, I(cell_size))) +
            facet_wrap(~feature_label)
    } else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    }
    if (show_tree) {
        g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", yend = "target_prin_graph_dim_2"), size = cell_link_size, linetype = "solid", na.rm = TRUE, data = edge_df)
    }

    # FIXME: setting size here overrides the marker expression funtionality.
    # Don't do it!
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
        if (class(data_df[, color_by]) == "numeric") {
            g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), size = I(cell_size), na.rm = TRUE, height = 5, width = 0.03) +
                scale_color_viridis(name = paste0("log10(", color_by, ")"), ...)
        } else {
            g <- g + geom_jitter(aes_string(color = color_by), size = I(cell_size), na.rm = TRUE, height = 5, width = 0.03)
        }
    } else {
        if (class(data_df[, color_by]) == "numeric") {
            g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), size = I(cell_size), na.rm = TRUE, height = 5, width = 0.03) +
                scale_color_viridis(name = paste0("log10(", color_by, " + 0.1)"), ...)
        } else {
            g <- g + geom_jitter(aes_string(color = color_by), size = I(cell_size), na.rm = TRUE, height = 5, width = 0.03)
        }
    }

    if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
        mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
        branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[, c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
        branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, mst_branch_nodes)
        branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), ]

        g <- g + geom_point(aes_string(x = "source_prin_graph_dim_1", y = "source_prin_graph_dim_2"),
            size = 2 * cell_size, na.rm = TRUE, data = branch_point_df
        ) +
            geom_text(aes_string(x = "source_prin_graph_dim_1", y = "source_prin_graph_dim_2", label = "branch_point_idx"),
                size = 1.5 * cell_size, color = "white", na.rm = TRUE, data = branch_point_df
            )
    }
    if (show_cell_names) {
        g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
    }
    g <- g +
        # scale_color_brewer(palette="Set1") +
        theme(strip.background = element_rect(colour = "white", fill = "white")) +
        theme(panel.border = element_blank()) +
        # theme(axis.line.x = element_line(size=0.25, color="black")) +
        # theme(axis.line.y = element_line(size=0.25, color="black")) +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
        theme(panel.background = element_rect(fill = "white")) +
        theme(legend.key = element_blank()) +
        xlab("") +
        ylab("") +
        theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) +
        # guides(color = guide_legend(label.position = "top")) +
        theme(legend.key = element_blank()) +
        theme(panel.background = element_rect(fill = "white")) +
        theme(
            line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
        )
    g
}


plot_tf_auc_in_pseudotime <- function(cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL,
                                      ncol = 1, panel_order = NULL, color_by = "State", trend_formula = "~ sm.ns(Pseudotime, df=3)",
                                      label_by_short_name = TRUE, relative_expr = TRUE, vertical_jitter = NULL,
                                      horizontal_jitter = NULL) {
    f_id <- NA
    Cell <- NA
    integer_expression <- FALSE
    relative_expr <- FALSE

    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    } else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    # cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names") ## 注释掉
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- cds_exprs$expression
    } else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        } else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    } else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset,
        cores = 1,
        trend_formula = trend_formula, relative_expr = T, new_data = new_data
    )
    colnames(model_expectation) <- colnames(cds_subset)
    rownames(model_expectation) <- rownames(cds_subset) ### 增加
    expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) {
        data.frame(expectation = model_expectation[
            x$f_id,
            x$Cell
        ])
    })
    cds_exprs <- merge(cds_exprs, expectation)
    # cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    # cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr

    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label,
            levels = panel_order
        )
    }
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by),
            size = I(cell_size),
            position = position_jitter(horizontal_jitter, vertical_jitter)
        )
    } else {
        q <- q + geom_point(size = I(cell_size), position = position_jitter(
            horizontal_jitter,
            vertical_jitter
        ))
    }
    # q <- q + geom_line(aes(x = Pseudotime, y = expectation),
    #    data = cds_exprs)

    q <- q + geom_smooth(method = loess, formula = y ~ x, se = FALSE, col = "black", size = 0.5)
    q <- q + scale_y_log10() + facet_wrap(~feature_label,
        nrow = nrow,
        ncol = ncol, scales = "free_y"
    )

    # if (min_expr < 1) {
    #     q <- q + expand_limits(y = c(min_expr, 1))
    # }
    if (relative_expr) {
        q <- q + ylab("Relative Expression")
    } else {
        q <- q + ylab("AUC")
    }
    q <- q + xlab("Pseudo-time")
    q <- q + monocle:::monocle_theme_opts() ### 增加调用包
    q
}

CustomCol2 <- function(n) {
    my_palette <- c(
        "#7fc97f", "#beaed4", "#fdc086", "#386cb0", "#f0027f", "#a34e3b", "#666666", "#1b9e77", "#d95f02", "#7570b3",
        "#d01b2a", "#43acde", "#efbd25", "#492d73", "#cd4275", "#2f8400", "#d9d73b", "#aed4ff", "#ecb9e5", "#813139",
        "#743fd2", "#434b7e", "#e6908e", "#214a00", "#ef6100", "#7d9974", "#e63c66", "#cf48c7", "#ffe40a", "#a76e93",
        "#d9874a", "#adc64a", "#5466df", "#d544a1", "#54d665", "#5e99c7", "#006874", "#d2ad2c", "#b5d7a5", "#9e8442",
        "#4e1737", "#e482a7", "#6f451d", "#2ccfe4", "#ae6174", "#a666be", "#a32b2b", "#ffff99", "#3fdacb", "#bf5b17"
    )
    return(my_palette[n])
}

# ===========================command line parameters setting=============================
option_list <- list(
    make_option(c("--input", "-i"),
        type = "character",
        help = "The pseudotime result in RDS format."
    ),
    make_option(c("--genelist", "-g"),
        type = "character", default = NULL,
        help = "The genelist used to plot heatmap with header."
    ),
    make_option(c("--vismethod", "-m"),
        type = "character", default = NULL,
        help = "the visulization methods for genes.The methods can be all, heatmap, expressplot, trajectoryplot, treeplot, module."
    ),
    make_option(c("--output", "-o"),
        type = "character", default = "./",
        help = "the output directory of Clustering results."
    ),
    make_option(c("--root"),
        type = "integer", default = NULL,
        help = "[Otional] The state designate as root."
    ),
    make_option(c("--show_branch"),
        type = "logical", default = FALSE,
        help = " Whether to display cell_trajectory_color_by_Pseudotime branch, Used with the root parameter.Prepare for branch heatmap"
    ),
    make_option(c("--clusters", "-n"),
        type = "integer", default = 3,
        help = "[Heatmap] Number of clusters for the heatmap of branch genes."
    ),
    make_option(c("--branchpoint", "-b"),
        type = "integer", default = NULL,
        help = "[Heatmap] The branch point used in the branch-heatmap."
    ),
    make_option(c("--showname"),
        type = "logical", default = NULL,
        help = "[Heatmap] Whether to display the row name."
    ),
    make_option(c("--CORES", "-j"),
        type = "integer", default = 5,
        help = "[Heatmap,module] the core number used to run this script."
    ),
    make_option(c("--groupby", "-c"),
        type = "character", default = "clusters",
        help = "[expressplot, treeplot] The grouppinig variable in pseudotime_gene_express. e.g. orig.ident"
    ),
    make_option(c("--numcol", "-l"),
        type = "integer", default = 2,
        help = "[expressplot] the number of columns used when laying out the panels for each gene's expression."
    ),
    make_option(c("--pointsize", "-s"),
        type = "double", default = 0.5,
        help = "[expressplot & trajectoryplot] The point size in the plot."
    ),
    make_option(c("--module"),
        type = "character", default = NULL,
        help = "[module expressplot] the module information (for genes) from the heatmap clustering."
    ),
    make_option(c("--topn"),
        type = "integer", default = 25,
        help = "the number of top DEGs to visualizse."
    ),
    make_option(c("--toptype"),
        type = "character", default = "both",
        help = "choose from up,down,or both DEGs to visualizse."
    ),
    make_option(c("--input_sc", "-v"),
        type = "character",
        help = "The scenic result with auc data in RDS format."
    ),
    make_option(c("--sc_genelist", "-t"),
        type = "character", default = NULL,
        help = "The TF genelist used to plot with header."
    ),
    make_option(c("--anno", "-a"),
        type = "character", default = NULL,
        help = "[Heatmap] Annotation file."
    ),
    make_option(c("--nbin"),
        type = "integer", default = 50,
        help = "[bin] number of bins."
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# =================================================================================
# parse the command line parameters
# =================================================================================
if (is.null(opt$output)) {
    print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
    output_dir <- getwd()
} else {
    output_dir <- opt$output
    if (!file.exists(output_dir)) {
        dir.create(output_dir, recursive = T)
    }
}

if (is.null(opt$vismethod)) {
    print("NO vismethod is provided. Only cell_trajectory_color_by_Pseudotime will be ploted if --root is specified.")
    vismethods <- NULL
} else if (opt$vismethod == "all") {
    vismethods <- c("heatmap", "expressplot", "trajectoryplot", "treeplot", "module")
} else {
    vismethods <- unlist(strsplit(opt$vismethod, ","))
}


if (is.null(opt$showname)) {
    if (!is.null(opt$genelist)) {
        if (opt$genelist == "ordering") {
            showname <- FALSE
        } else {
            showname <- TRUE
        }
    }
} else {
    showname <- opt$showname
}

# =================================================================================
# main workflow
# =================================================================================
gbm_cds <- readRDS(opt$input)

## change the root state
if (!is.null(opt$root)) {
    gbm_cds <- orderCells(gbm_cds, root_state = opt$root) # root_state =分支编
    p <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = opt$pointsize, show_branch_points = opt$show_branch) + scale_colour_viridis_c(option = "inferno")
    if (opt$show_branch) {
        ggsave(file.path(output_dir, "cell_trajectory_color_by_Pseudotime_show_branch_points.pdf"), plot = p)
        ggsave(file.path(output_dir, "cell_trajectory_color_by_Pseudotime_show_branch_points.png"), plot = p, dpi = 1000)
    } else {
        ggsave(file.path(output_dir, "cell_trajectory_color_by_Pseudotime.pdf"), plot = p)
        ggsave(file.path(output_dir, "cell_trajectory_color_by_Pseudotime.png"), plot = p, dpi = 1000)
        saveRDS(gbm_cds, opt$input)
    }
}

## subset the gbm_cds if genelist id provided
if (!is.null(opt$genelist)) {
    if (opt$genelist == "ordering") {
        genes <- as.factor(subset(gbm_cds@featureData@data, use_for_ordering == TRUE)$gene_short_name)
        to_be_tested <- row.names(subset(fData(gbm_cds), gene_short_name %in% levels(genes)))
    } else {
        gene <- read.delim(opt$genelist, sep = "\t")
        if (dim(gene)[2] > 1) {
            up <- filter(gene, FoldChange > 1) %>%
                arrange(desc(log2FoldChange)) %>%
                top_n(opt$topn, log2FoldChange) %>%
                select(gene)
            down <- filter(gene, FoldChange < 1) %>%
                arrange(log2FoldChange) %>%
                top_n(as.numeric(paste0("-", opt$topn)), log2FoldChange) %>%
                select(gene)
            if (opt$toptype == "up") {
                gene <- up
            } else if (opt$toptype == "down") {
                gene <- down
            } else {
                gene <- rbind(up, down)
            }
            gene[, 1] <- factor(as.character(gene[, 1]))
        }
        to_be_tested <- Seurat::CaseMatch(search = levels(gene[, 1]), match = row.names(fData(gbm_cds)))
    }
    gbm_cds <- gbm_cds[to_be_tested, ]
}
# =================================================================================
# visualization
# =================================================================================
for (vismethod in vismethods) {
    if (vismethod == "bin") {
        md=pData(gbm_cds)
        if( ! opt$groupby %in% colnames(md) ){
            stop("Column ",opt$groupby," not found in data")
        }
        md$Pseudotime_bin =  floor(gbm_cds$Pseudotime*(opt$nbin-1)/max(gbm_cds$Pseudotime))+1
        xaxis="Pseudotime_bin"
        propby=opt$groupby
        field4emeta = colnames(md)
        cluster_ids <- md[, propby]
        sample_ids <- md[, xaxis]
        counts <- table(cluster_ids, sample_ids)
        df <- reshape2::melt(t(round(t(counts)/colSums(counts) * 100, 2)), 
            varnames = c(propby, xaxis), value.name = "freq")
        df[, xaxis] = factor(df[, xaxis])
        df[, propby] = factor(df[, propby])
        # plot
        p = ggplot(df) + geom_bar(aes_string(y = "freq", 
            x = xaxis, fill = propby), position = "stack", stat = "identity" ,width=1) + 
            scale_fill_manual(values = CustomCol2(1:50)) + 
            scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + 
            labs(x = paste0("Pseudotime (",opt$nbin,"bins)"), y = "Cells per bin[%]") + 
            theme_bw() + 
            theme(panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                color = NA), panel.border = element_blank(), axis.ticks.x = element_blank(), 
                axis.text = element_text(color = "black"), axis.text.x = element_blank())
        ggsave(file.path(output_dir, paste0("Proportion_of_cells_in_pseudotime_bins_by_",propby,".pdf")), plot = p)
        ggsave(file.path(output_dir, paste0("Proportion_of_cells_in_pseudotime_bins_by_",propby,".png")), plot = p)
	}
    if (vismethod == "heatmap") {
        if (is.null(opt$anno)) {
            stop("Please provide Annotation file with --anno parameter.")
        } else {
            annotation <- opt$anno
        }
        if (is.null(opt$clusters)) {
            print("The number of clusters for the heatmap,the default will be used.")
            clusters_num <- 3
        } else {
            clusters_num <- opt$clusters
        }

        if (is.null(opt$branchpoint)) {
            print("NO branchpoint is provided, the heatmap without branchtime will be plotted.")
            branchpoint <- NULL
        } else {
            branchpoint <- opt$branchpoint
        }
        if (is.null(branchpoint)) { # heatmap without branchtime
            heatmap <- plot_pseudotime_heatmap(gbm_cds, cores = opt$CORES, cluster_rows = T, num_clusters = clusters_num, show_rownames = showname, return_heatmap = T)
            ggsave(file.path(output_dir, "pseudotime_heatmap.pdf"), plot = heatmap$ph_res)
            ggsave(file.path(output_dir, "pseudotime_heatmap.png"), plot = heatmap$ph_res, dpi = 1000)
            gene_clusters <- cutree(heatmap$ph_res$tree_row, k = clusters_num)
        } else { # heatmap with branchtime
            new_cds <- buildBranchCellDataSet(gbm_cds, branch_point = branchpoint, progenitor_method = "duplicate")
            # cell_fate1
            cell_fate1 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]), ]$State)
            # cell_fate2
            cell_fate2 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]), ]$State)
            branch_labels <- c(paste("State", paste(sort(setdiff(cell_fate1, cell_fate2)), collapse = "-")), paste("State", paste(sort(setdiff(cell_fate2, cell_fate1)), collapse = "-")))
            branch_heatmap <- plot_genes_branched_heatmap(gbm_cds, branch_point = branchpoint, cores = opt$CORES, cluster_rows = T, num_clusters = clusters_num, use_gene_short_name = T, show_rownames = showname, branch_labels = branch_labels, return_heatmap = T)
            ggsave(file.path(output_dir, "pseudotime_heatmap_branchtime.pdf"), plot = branch_heatmap$ph_res)
            ggsave(file.path(output_dir, "pseudotime_heatmap_branchtime.png"), plot = branch_heatmap$ph_res, dpi = 1000)
            gene_clusters <- cutree(branch_heatmap$ph_res$tree_row, k = clusters_num)
        }
        ### annotation
        gene_clustering <- data.frame(gene_clusters)
        gene_clustering[, 1] <- as.character(gene_clustering[, 1])
        colnames(gene_clustering) <- "gene_module"
        if (opt$genelist == "ordering") {
            ord_info <- featureData(gbm_cds)@data %>%
                filter(use_for_ordering == "TRUE") %>%
                dplyr::select(gene_short_name, vst.variance.standardized, pval, qval) %>%
                tibble::column_to_rownames(var = "gene_short_name")
            gene_inter <- intersect(rownames(gene_clustering), rownames(ord_info))
            ord_info <- ord_info[gene_inter, ]
            gene_clustering <- cbind(gene_clustering, ord_info)
        }
        gene_clustering <- gene_clustering %>% tibble::rownames_to_column(var = "gene")
        anno <- read.delim(annotation, stringsAsFactors = F, sep = "\t", quote = "")
        gene_clustering <- left_join(gene_clustering, anno, by = c("gene" = "id"))
        gene_clustering[is.na(gene_clustering)] <- "--"
        write.table(gene_clustering, file.path(output_dir, "pseudotime_heatmap_gene_module_anno.xls"), sep = "\t", quote = F, row.names = F)
    }

    if (vismethod == "expressplot") {
        for (i in seq(1, length(to_be_tested), 10)) {
            j <- min(i + 9, length(to_be_tested))
            to_be_tested_sub <- to_be_tested[i:j]
            cds_subset <- gbm_cds[to_be_tested_sub, ]
            nlevel <- length(unique(cds_subset@phenoData@data[, opt$groupby]))
            if (is.null(opt$branchpoint)) {
                print("NO branchpoint is provided, the expressplot without branchtime will be plotted.")
                branchpoint <- NULL
            } else {
                branchpoint <- opt$branchpoint
            }
            if (is.null(branchpoint)) { # expressplot without branchtime
                p <- plot_genes_in_pseudotime(cds_subset, color_by = opt$groupby, cell_size = opt$pointsize, ncol = opt$numcol) + scale_colour_manual(values = CustomCol2(1:nlevel)) + guides(colour = guide_legend(override.aes = list(size = 1.5)))
                ggsave(file.path(output_dir, paste0("pseudotime_gene_express_", basename(opt$genelist), ifelse(i == 1, "", -ceiling(i / 10)), ".pdf")), height = ifelse(length(to_be_tested_sub) == 1, 7, ceiling(length(to_be_tested_sub) / opt$numcol) * 2))
                ggsave(file.path(output_dir, paste0("pseudotime_gene_express_", basename(opt$genelist), ifelse(i == 1, "", -ceiling(i / 10)), ".png")), height = ifelse(length(to_be_tested_sub) == 1, 7, ceiling(length(to_be_tested_sub) / opt$numcol) * 2), dpi = 1000)
            } else { # expressplot with branchtime
                new_cds <- buildBranchCellDataSet(cds_subset, branch_point = branchpoint, progenitor_method = "duplicate")
                # cell_fate1
                cell_fate1 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]), ]$State)
                # cell_fate2
                cell_fate2 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]), ]$State)
                branch_labels <- c(paste("State", paste(sort(setdiff(cell_fate1, cell_fate2)), collapse = "-")), paste("State", paste(sort(setdiff(cell_fate2, cell_fate1)), collapse = "-")))
                p <- plot_genes_branched_pseudotime(cds_subset, color_by = opt$groupby, branch_point = branchpoint, cell_size = opt$pointsize, ncol = opt$numcol, branch_labels = branch_labels) + scale_colour_manual(values = CustomCol2(1:nlevel)) + guides(colour = guide_legend(override.aes = list(size = 1.5)))
                ggsave(file.path(output_dir, paste0("pseudotime_gene_express_branchtime", ifelse(i == 1, "", -ceiling(i / 10)), ".pdf")), height = ceiling(length(to_be_tested_sub) / opt$numcol) * 2)
                ggsave(file.path(output_dir, paste0("pseudotime_gene_express_branchtime", ifelse(i == 1, "", -ceiling(i / 10)), ".png")), height = ceiling(length(to_be_tested_sub) / opt$numcol) * 2, dpi = 1000)
            }
        }
    }

	if ( vismethod == "expressplot_line" ){
		for( i in seq(1, length(to_be_tested), 10) ){
			j = min(i+9, length(to_be_tested))
			to_be_tested_sub = to_be_tested[i:j]
			cds_subset <- gbm_cds[to_be_tested_sub,]
			nlevel = length(unique(cds_subset@phenoData@data[,opt$groupby]))
			# expressplot without branchtime
			p = plot_genes_in_pseudotime_line(cds_subset, color_by= opt$groupby ,cell_size = opt$pointsize , ncol = opt$numcol) + scale_colour_manual( values = CustomCol2(1:nlevel) ) + guides(colour = guide_legend(override.aes = list(size=1.5))) 
			ggsave(file.path(output_dir, paste0("pseudotime_gene_express_line_",basename(opt$genelist),ifelse(i==1, "", -ceiling(i/10)), ".pdf")), height = ifelse( length(to_be_tested_sub)==1, 7, ceiling(length(to_be_tested_sub)/opt$numcol)*2 ) )
			ggsave(file.path(output_dir, paste0("pseudotime_gene_express_line_",basename(opt$genelist),ifelse(i==1, "", -ceiling(i/10)), ".png")), height = ifelse( length(to_be_tested_sub)==1, 7, ceiling(length(to_be_tested_sub)/opt$numcol)*2 ), dpi = 1000)
		}
	}

    if (vismethod == "trajectoryplot") {
        for (i in to_be_tested) {
            p <- plot_cell_trajectory(gbm_cds, markers = i, use_color_gradient = T, show_branch_points = F, show_tree = F, cell_size = opt$pointsize) + theme(legend.text = element_text(size = 10)) + scale_color_gradientn(colours = c("grey", "yellow", "red"))
            ggsave(file.path(output_dir, paste0(i, ".pdf")))
            ggsave(file.path(output_dir, paste0(i, ".png")), dpi = 1000)
        }
    }


    if (vismethod == "module") {
        # retrieve annotation colors
        if ("heatmap" %in% vismethods) {
            gene_clustering <- data.frame(gene_clusters)
            colnames(gene_clustering) <- "Gene_Module"
        } else if (!is.null(opt$module)) {
            gene_clustering <- data.frame(Module = read.table(opt$module, header = T, row.names = 1, sep = "\t"))
        } else {
            stop("no gene clustering/module provided.")
        }
        if (is.null(opt$clusters)) {
            print("The number of clusters for the heatmap,the default will be used.")
            clusters_num <- 3
        } else {
            clusters_num <- opt$clusters
        }
        gene_clustering[, 1] <- as.factor(gene_clustering[, 1])
        # colors=unlist(pheatmap:::generate_annotation_colours(gene_clustering,NULL,T),use.names=F)
        colors <- scales::dscale(factor(1:clusters_num), scales::hue_pal(l = 75))
        print("calculating for each module")
        genelist <- list()
        cds_exprs <- list()
        p <- list()
        if (is.null(opt$branchpoint)) {
            print("NO branchpoint is provided, the heatmap without branchtime will be plotted.")
            branchpoint <- NULL
        } else {
            branchpoint <- opt$branchpoint
        }
        if (is.null(branchpoint)) { # without branch
            for (i in unique(gene_clustering[, 1])) {
                genelist[[i]] <- as.factor(rownames(subset(gene_clustering, Gene_Module == i)))
                cds_subset <- gbm_cds[row.names(subset(fData(gbm_cds), gene_short_name %in% levels(genelist[[i]]))), ]
                cds_exprs[[i]] <- exprs(cds_subset)
                cds_exprs[[i]] <- Matrix::t(Matrix::t(cds_exprs[[i]]) / sizeFactors(cds_subset))
                cds_exprs[[i]] <- reshape2::melt(round(as.matrix(cds_exprs[[i]])))
                min_expr <- cds_subset@lowerDetectionLimit
                colnames(cds_exprs[[i]]) <- c("f_id", "Cell", "expression")
                cds_exprs[[i]] <- aggregate(. ~ Cell, data = cds_exprs[[i]], mean) %>% select(f_id, Cell, expression)
                cds_exprs[[i]]$f_id <- "gene_set"
                cds_pData <- pData(cds_subset)
                cds_fData <- fData(cds_subset)
                cds_exprs[[i]] <- merge(cds_exprs[[i]], cds_pData, by.x = "Cell", by.y = "row.names")
                cds_exprs[[i]]$adjusted_expression <- cds_exprs[[i]]$expression
                cds_exprs[[i]]$feature_label <- cds_exprs[[i]]$f_id
                cds_exprs[[i]]$f_id <- as.character(cds_exprs[[i]]$f_id)
                cds_exprs[[i]]$feature_label <- factor(cds_exprs[[i]]$feature_label)
                new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
                model_expectation <- genSmoothCurves(cds_subset, cores = opt$CORES, relative_expr = T, new_data = new_data)
                colnames(model_expectation) <- colnames(cds_subset)
                model_expectation <- t(colMeans(model_expectation))
                rownames(model_expectation) <- "gene_set"
                expectation <- ddply(cds_exprs[[i]], .(f_id, Cell), function(x) data.frame(expectation = model_expectation[x$f_id, x$Cell]))
                cds_exprs[[i]] <- merge(cds_exprs[[i]], expectation)
                cds_exprs[[i]]$expression[cds_exprs[[i]]$expression < min_expr] <- min_expr
                cds_exprs[[i]]$expectation[cds_exprs[[i]]$expectation < min_expr] <- min_expr
                saveRDS(cds_exprs, "cds_exprs.rds") # genSmoothCurves() takes time. saving the expectation result just in case
                p[[i]] <- ggplot() +
                    geom_line(aes(x = Pseudotime, y = expectation), size = 1, colour = colors[as.integer(i)], data = cds_exprs[[i]]) +
                    expand_limits(y = c(min_expr, 1)) +
                    ylab("Relative Expression") +
                    xlab("Pseudo-time") +
                    scale_y_log10() +
                    theme_bw() +
                    theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(size = 1, colour = "black"))
                ggsave(file.path(output_dir, paste0("expressplot_module", i, ".pdf")), plot = p[[i]], width = 5, height = 3)
                ggsave(file.path(output_dir, paste0("expressplot_module", i, ".png")), plot = p[[i]], width = 5, height = 3)
            }
        } else { ## branchtime
            new_cds <- buildBranchCellDataSet(gbm_cds, branch_point = branchpoint, progenitor_method = "duplicate")
            cell_fate1 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]), ]$State)
            cell_fate2 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]), ]$State)
            branch_labels <- c(paste("State", paste(sort(setdiff(cell_fate1, cell_fate2)), collapse = "-")), paste("State", paste(sort(setdiff(cell_fate2, cell_fate1)), collapse = "-")))
            new_cds <- buildBranchCellDataSet(gbm_cds, branch_point = branchpoint, branch_labels = branch_labels, progenitor_method = "duplicate")
            for (i in unique(gene_clustering[, 1])) {
                genelist[[i]] <- as.factor(rownames(subset(gene_clustering, Gene_Module == i)))
                cds_subset <- new_cds[row.names(subset(fData(new_cds), gene_short_name %in% levels(genelist[[i]]))), ]
                Branch <- NA
                CM <- exprs(cds_subset)
                CM <- Matrix::t(Matrix::t(CM) / sizeFactors(cds_subset))
                cds_exprs[[i]] <- reshape2::melt(round(as.matrix(CM)))
                min_expr <- cds_subset@lowerDetectionLimit
                colnames(cds_exprs[[i]]) <- c("f_id", "Cell", "expression")
                cds_exprs[[i]] <- aggregate(. ~ Cell, data = cds_exprs[[i]], mean) %>% select(f_id, Cell, expression)
                cds_exprs[[i]]$f_id <- "gene_set"
                cds_pData <- pData(cds_subset)
                cds_fData <- fData(cds_subset)
                cds_exprs[[i]] <- merge(cds_exprs[[i]], cds_pData, by.x = "Cell", by.y = "row.names")
                cds_exprs[[i]]$adjusted_expression <- round(cds_exprs[[i]]$expression)
                cds_exprs[[i]]$feature_label <- cds_exprs[[i]]$f_id
                cds_exprs[[i]]$feature_label <- as.factor(cds_exprs[[i]]$feature_label)
                cds_exprs[[i]]$Branch <- as.factor(cds_exprs[[i]]$Branch)
                new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)
                full_model_expectation <- genSmoothCurves(cds_subset,
                    cores = opt$CORES,
                    trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch", relative_expr = T, new_data = new_data
                )
                colnames(full_model_expectation) <- colnames(cds_subset)
                full_model_expectation2 <- t(colMeans(full_model_expectation))
                rownames(full_model_expectation2) <- "gene_set"
                cds_exprs[[i]]$full_model_expectation <- apply(cds_exprs[[i]], 1, function(x) full_model_expectation2[x[2], x[1]])
                cds_exprs[[i]]$feature_label <- cds_exprs[[i]]$f_id
                cds_exprs[[i]]$feature_label <- factor(cds_exprs[[i]]$feature_label)
                cds_exprs[[i]]$expression[is.na(cds_exprs[[i]]$expression)] <- min_expr
                cds_exprs[[i]]$expression[cds_exprs[[i]]$expression < min_expr] <- min_expr
                cds_exprs[[i]]$full_model_expectation[is.na(cds_exprs[[i]]$full_model_expectation)] <- min_expr
                cds_exprs[[i]]$full_model_expectation[cds_exprs[[i]]$full_model_expectation < min_expr] <- min_expr
                cds_exprs[[i]]$State <- as.factor(cds_exprs[[i]]$State)
                cds_exprs[[i]]$Branch <- as.factor(cds_exprs[[i]]$Branch)
                saveRDS(cds_exprs, "cds_exprs_branched.rds")
                p[[i]] <- ggplot(aes(Pseudotime, expression), data = cds_exprs[[i]]) +
                    scale_y_log10() +
                    expand_limits(y = min_expr) +
                    geom_line(aes_string(
                        x = "Pseudotime", y = "full_model_expectation",
                        linetype = "Branch"
                    ), size = 1, color = colors[as.integer(i)], data = cds_exprs[[i]]) +
                    ylab("Expression") +
                    xlab("Pseudotime (stretched)") +
                    theme_bw() +
                    theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(size = 1, colour = "black"))
                ggsave(file.path(output_dir, paste0("expressplot_module", i, "_branchtime.pdf")), plot = p[[i]], width = 5, height = 3)
                ggsave(file.path(output_dir, paste0("expressplot_module", i, "_branchtime.png")), plot = p[[i]], width = 5, height = 3)
            }
        }
    }

    if (vismethod == "treeplot") {
        nlevel <- length(unique(gbm_cds@phenoData@data[, opt$groupby]))
        p <- plot_complex_cell_trajectory(gbm_cds, color_by = opt$groupby, show_branch_points = T, cell_size = opt$pointsize, cell_link_size = 0.3) + scale_colour_manual(values = CustomCol2(1:nlevel)) + guides(colour = guide_legend(override.aes = list(size = 3)))
        ggsave(file.path(output_dir, paste0("pseudotime_treeplot_colorby_", opt$groupby, ".pdf")), plot = p)
        ggsave(file.path(output_dir, paste0("pseudotime_treeplot_colorby_", opt$groupby, ".png")), plot = p, dpi = 1000)
        p_facet <- plot_complex_cell_trajectory(gbm_cds, color_by = opt$groupby, show_branch_points = T, cell_size = opt$pointsize, cell_link_size = 0.3) + scale_colour_manual(values = CustomCol2(1:nlevel)) + guides(colour = guide_legend(override.aes = list(size = 3))) + facet_wrap(eval(expr = parse(text = paste0("~ ", opt$groupby))))
        ggsave(file.path(output_dir, paste0("pseudotime_treeplot_colorby_", opt$groupby, "_splited.pdf")), plot = p_facet, width = 10, height = 10)
        ggsave(file.path(output_dir, paste0("pseudotime_treeplot_colorby_", opt$groupby, "_splited.png")), plot = p_facet, dpi = 1000, width = 10, height = 10)
    }
}


if (!is.null(opt$input_sc) && !is.null(opt$sc_genelist)) {
    auc <- readRDSMC(opt$input_sc)

    ## keep barcode number same
    inter <- intersect(colnames(gbm_cds), colnames(auc))
    auc <- auc[, inter]
    gbm_cds <- gbm_cds[, inter]

    ### Replace the expression matrix in pseudotime RDS with AUC matrix
    e.scenic <- rlang::env(exprs = auc@assays$data$AUC)
    assayData(gbm_cds) <- e.scenic

    ### genelist tf as gene "FOXL1 (16g)"
    gene <- read.delim(opt$sc_genelist, sep = "\t")

    to_be_tested <- as.vector(gene[, 1])
    for (i in seq(1, length(to_be_tested), 10)) {
        j <- min(i + 9, length(to_be_tested))
        to_be_tested_sub <- to_be_tested[i:j]
        cds_subset <- gbm_cds[to_be_tested_sub, ]
        nlevel <- length(unique(cds_subset@phenoData@data[, opt$groupby]))
        p <- plot_tf_auc_in_pseudotime(cds_subset, color_by = opt$groupby, cell_size = opt$pointsize, ncol = opt$numcol, min_expr = 0.0001) + scale_colour_manual(values = CustomCol2(1:nlevel)) + guides(colour = guide_legend(override.aes = list(size = 1.5)))

        ggsave(file.path(output_dir, paste0("pseudotime_TF_AUC_", basename(opt$sc_genelist), ifelse(i == 1, "", -ceiling(i / 10)), ".pdf")), height = ceiling(length(to_be_tested_sub) / opt$numcol) * 2)
        ggsave(file.path(output_dir, paste0("pseudotime_TF_AUC_", basename(opt$sc_genelist), ifelse(i == 1, "", -ceiling(i / 10)), ".png")), height = ceiling(length(to_be_tested_sub) / opt$numcol) * 2, dpi = 1000)
    }
}