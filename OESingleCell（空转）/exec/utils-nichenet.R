#为了不加载seurat包而调用函数，因此在函数内部增加seurat::xxxx调用
get_expressed_genes_updata <- function (ident, seurat_obj, pct = 0.1, assay_oi = NULL){
  if (!"RNA" %in% names(seurat_obj@assays)) {
      if ("Spatial" %in% names(seurat_obj@assays)) {
          if (class(seurat_obj@assays$Spatial@data) != "matrix" &
              class(seurat_obj@assays$Spatial@data) != "dgCMatrix") {
              warning("Spatial Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$Spatial@data' for default or 'seurat_obj@assays$SCT@data' for when the single-cell transform pipeline was applied")
          }
          if (sum(dim(seurat_obj@assays$Spatial@data)) == 0) {
              stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$Spatial@data'")
          }
      }
  }else {
      if (class(seurat_obj@assays$RNA@data) != "matrix" & class(seurat_obj@assays$RNA@data) != "dgCMatrix") {
          warning("Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data or seurat_obj@assays$SCT@data for when the single-cell transform pipeline was applied")
      }
      if ("integrated" %in% names(seurat_obj@assays)) {
          if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$integrated@data)) ==0)
              stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data")
      }
      else if ("SCT" %in% names(seurat_obj@assays)) {
          if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$SCT@data)) ==0) {
              stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$SCT@data' for data corrected via SCT")
          }
      }
      else {
          if (sum(dim(seurat_obj@assays$RNA@data)) == 0) {
              stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data'")
          }
      }
  }
  if (sum(ident %in% unique(Seurat::Idents(seurat_obj))) != length(ident)) {
      stop("One or more provided cell clusters is not part of the 'Idents' of your Seurat object")
  }
  if (!is.null(assay_oi)) {
      if (!assay_oi %in% Seurat::Assays(seurat_obj)) {
          stop("assay_oi should be an assay of your Seurat object")
      }
  }
  cells_oi <- Seurat::Idents(seurat_obj) %>%
    .[Seurat::Idents(seurat_obj) %in% ident] %>%
    names()
  if (!is.null(assay_oi)) {
      cells_oi_in_matrix <- intersect(colnames(seurat_obj[[assay_oi]]@data),
                                      cells_oi)
      exprs_mat <- seurat_obj[[assay_oi]]@data %>% .[, cells_oi_in_matrix]
  }else {
      if ("integrated" %in% names(seurat_obj@assays)) {
          warning("Seurat object is result from the Seurat integration workflow. The expressed genes are now defined based on the integrated slot. You can change this via the assay_oi parameter of the get_expressed_genes() functions. Recommended assays: RNA or SCT")
          cells_oi_in_matrix <- intersect(colnames(seurat_obj@assays$integrated@data),
                                          cells_oi)
          if (length(cells_oi_in_matrix) != length(cells_oi))
            stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$integrated@data). Please check that the expression matrix contains cells in columns and genes in rows.")
          exprs_mat <- seurat_obj@assays$integrated@data %>%
            .[, cells_oi_in_matrix]
      }
      else if ("SCT" %in% names(seurat_obj@assays) & !"Spatial" %in%
          names(seurat_obj@assays)) {
          warning("Seurat object is result from the Seurat single-cell transform workflow. The expressed genes are defined based on the SCT slot. You can change this via the assay_oi parameter of the get_expressed_genes() functions. Recommended assays: RNA or SCT")
          cells_oi_in_matrix <- intersect(colnames(seurat_obj@assays$SCT@data),
                                          cells_oi)
          if (length(cells_oi_in_matrix) != length(cells_oi))
            stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$SCT@data). Please check that the expression matrix contains cells in columns and genes in rows.")
          exprs_mat <- seurat_obj@assays$SCT@data %>% .[, cells_oi_in_matrix]
      }
      else if ("Spatial" %in% names(seurat_obj@assays) & !"SCT" %in%
          names(seurat_obj@assays)) {
          warning("Seurat object is result from the Seurat spatial object. The expressed genes are defined based on the Spatial slot. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! ;-) )")
          cells_oi_in_matrix <- intersect(colnames(seurat_obj@assays$Spatial@data),
                                          cells_oi)
          if (length(cells_oi_in_matrix) != length(cells_oi))
            stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial@data). Please check that the expression matrix contains cells in columns and genes in rows.")
          exprs_mat <- seurat_obj@assays$Spatial@data %>% .[,
            cells_oi_in_matrix]
      }
      else if ("Spatial" %in% names(seurat_obj@assays) & "SCT" %in%
          names(seurat_obj@assays)) {
          warning("Seurat object is result from the Seurat spatial object, followed by the SCT workflow. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! The expressed genes are defined based on the SCT slot, but this can be changed via the assay_oi parameter.")
          cells_oi_in_matrix <- intersect(colnames(seurat_obj@assays$SCT@data),
                                          cells_oi)
          if (length(cells_oi_in_matrix) != length(cells_oi))
            stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial@data). Please check that the expression matrix contains cells in columns and genes in rows.")
          exprs_mat <- seurat_obj@assays$SCT@data %>% .[, cells_oi_in_matrix]
      }
      else {
          if (sum(cells_oi %in% colnames(seurat_obj@assays$RNA@data)) ==
            0)
            stop("None of the cells are in colnames of 'seurat_obj@assays$RNA@data'. The expression matrix should contain cells in columns and genes in rows.")
          cells_oi_in_matrix <- intersect(colnames(seurat_obj@assays$RNA@data),
                                          cells_oi)
          if (length(cells_oi_in_matrix) != length(cells_oi))
            stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$RNA@data). Please check that the expression matrix contains cells in columns and genes in rows.")
          exprs_mat <- seurat_obj@assays$RNA@data %>% .[, cells_oi_in_matrix]
      }
  }
    n_cells_oi_in_matrix <- length(cells_oi_in_matrix)

    genes <- exprs_mat %>%
      apply(1, function(x) {
          sum(x > 0) / n_cells_oi_in_matrix
      }) %>%
      .[. >= pct] %>%
      names()

  return(genes)
}


sub_nichenet <- subparsers$add_parser(
    "nichenet",
    # description = docstring,
    # formatter_class = 'argparse.RawTextHelpFormatter',
    formatter_class = "lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
    argument_default = "True",
    help = "The goal of NicheNet is to study intercellular communication from a computational perspective."
)
sub_nichenet$add_argument(
    "-s",
    "--species",
    type = "character",
    default = "human",
    choices=c("human","mouse"),
    help = "the species for the count matrix:human, mouse.[default: %(default)s]"
)
sub_nichenet$add_argument(
    "-l",
    "--celltype",
    type = "character",
    help = "The cell type annotation column name to use in seurat metadata."
)
sub_nichenet$add_argument(
    "--receiver",
    type = "character",
    help = "The cell type annotation column name to use in seurat metadata."
)
sub_nichenet$add_argument(
    "--sender",
    type = "character",
    help = "The cell type annotation column name to use in seurat metadata."
)
sub_nichenet$add_argument(
    "-g",
    "--genelist",
    type = "character",
    help = "[OPTIONAL]The cell type annotation column name to use in seurat metadata."
)
sub_nichenet$add_argument(
    "-d",
    "--diff",
    type = "character",
    help = "[OPTIONAL] colname of  matadata for diff,variable1,variable2."
)
sub_nichenet$add_argument(
    "-n",
    "--topn",
    type = "integer",
    help = "the top activity ligand."
)
sub_nichenet$add_argument(
    "-m",
    "--database",
    # type = "charactor",
    default = "/data/database/NicheNet",
    help = "nichenet database  directory, include:ligand_target_matrix.rds,lr_network.rds,weighted_networks.rds.[default: %(default)s]"
)
sub_nichenet$add_argument(
    "--pct",
    type = "double",
    default = 0.10,
    help = "Consider genes expressed if they are expressed in at least a specific fraction of cells of a cluster.[default: %(default)s]"
)
sub_nichenet$add_argument(
    "--ligand_visiual_cutoff",
    type = "double",
    default = 0.25,
    help =paste0("cutoff for prepare_ligand_target_visualization:Quantile cutoff on the ligand-target scores of the ",
    "input weighted ligand-target network. Scores under this cutoff will be set to .[default: %(default)s]")
)

args <- commandArgs(TRUE)
if ("nichenet" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "nichenet") {
        # ==============================================================================================================
        futile.logger::flog.info("step1：导入seurat对象")
        suppressMessages(data_ob <- OESingleCell::ReadX(
            input = opt$input,
            informat = opt$informat,
            assays = assays[1],
            data.use = dataslots,
            verbose = F
        ))
        # ==============================================================================================================
        futile.logger::flog.info("step2: 导入nichenet数据库并根据情况进行转换")
        ligand_target_matrix <- readRDS(glue::glue("{opt$database}/ligand_target_matrix.rds"))
        lr_network <- readRDS(glue::glue("{opt$database}/lr_network.rds"))
        weighted_networks <- readRDS(glue::glue("{opt$database}/weighted_networks.rds"))
        weighted_networks_lr <- weighted_networks$lr_sig %>%
            dplyr::inner_join(
                lr_network %>% dplyr::distinct(from, to),
                by = c("from", "to")
            )

        ## =====
        futile.logger::flog.info("nichenet数据库，默认物种为人，如果为小鼠，需要进行转换")
        if (opt$species == "mouse") {
            lr_network <- lr_network %>%
                dplyr::mutate(
                    from = convert_human_to_mouse_symbols(from),
                    to = convert_human_to_mouse_symbols(to)
                ) %>%
                tidyr::drop_na()
            colnames(ligand_target_matrix) <- ligand_target_matrix %>%
                colnames() %>%
                convert_human_to_mouse_symbols()
            rownames(ligand_target_matrix) <- ligand_target_matrix %>%
                rownames() %>%
                convert_human_to_mouse_symbols()
            ligand_target_matrix <- nichenetr::ligand_target_matrix %>%
                .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
            weighted_networks_lr <- weighted_networks_lr %>%
                dplyr::mutate(
                    from = convert_human_to_mouse_symbols(from),
                    to = convert_human_to_mouse_symbols(to)
                ) %>%
                tidyr::drop_na()
        }
        # ==============================================================================================================
        futile.logger::flog.info("step3: 定义细胞群：[sender/niche] cell population 和 [receiver/target] cell population")
        celltype <- opt$celltype
        Seurat::Idents(data_ob) <- data_ob[[celltype]]

        if(!is.null(opt$receiver)){
              receiver<- unlist(strsplit(opt$receiver, ",", perl = T))
        }else{
              receiver<- names(table(data_ob[[celltype]]))
        }

        if(!is.null(opt$sender)){
              sender_celltypes<- unlist(strsplit(opt$sender, ",", perl = T))
        }else{
              sender_celltypes<- names(table(data_ob[[celltype]]))
        }

        # ==============================================================================================================
        futile.logger::flog.info("step4:定义背景基因和感兴趣的靶基因集合")
        futile.logger::flog.info("4.1: 定义背景基因")

        expressed_genes_receiver <- receiver%>%
            unique() %>%
            lapply(get_expressed_genes_updata, data_ob, opt$pct, assay_oi = assays) %>%
            unlist() %>%
            unique()

        background_expressed_genes <- nichenetr::get_expressed_genes(
            receiver,
            data_ob,
            pct = opt$pct,
            assay_oi = assays
        ) %>%
            .[. %in% rownames(ligand_target_matrix)]

        expressed_genes_sender <- sender_celltypes %>%
            unique() %>%
            lapply(get_expressed_genes_updata, data_ob, opt$pct, assay_oi = assays) %>%
            unlist() %>%
            unique()
        # ==============================================================================================================
        futile.logger::flog.info("4.2:定义感兴趣的靶基因集合，可以自行输入基因列表，也可以输入参数采用Seurat::FindMarkers进行鉴定，选择前top200个基因， 默认筛选标准为: min.pct = 0.10 & p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25 ")
        if (!is.null(opt$genelist)) {
            geneset_oi <- read.table(opt$genelist, head = T, sep = "\t")
            if (dim(geneset_oi)[2] > 1) {
                geneset_oi <- dplyr::filter(
                    geneset_oi,
                    pvalue <= 0.05 & log2FoldChange >= 0.25
                ) %>%
                    dplyr::pull(gene) %>%
                    unique()
            } else {
                geneset_oi <- geneset_oi[, 1] %>% unique()
            }
            geneset_oi <- intersect(rownames(data_ob), geneset_oi)
        } else if (!is.null(opt$diff)) {
              seurat_obj_receiver <-  subset(data_ob, idents = receiver)
              condition <- unlist(strsplit(opt$diff, ",", perl = T))
              seurat_obj_receiver <- Seurat::SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[[condition[1]]])
              condition_oi <-  condition[2]
              condition_reference <-  condition[3]
              DE_table_receiver <-  Seurat::FindMarkers(object  = data_ob,
                                                        ident.1 = condition_oi,
                                                        ident.2 = condition_reference,
                                                        min.pct = 0.10 ) %>%
                                    tibble::rownames_to_column("gene")
              readr::write_tsv(DE_table_receiver  ,file=glue::glue("{output_dir}/input-no-filtered-all-diff-gene-for-{condition_oi}-vs-{condition_reference}.tsv") )
              geneset_oi <-  DE_table_receiver %>% dplyr::filter(p_val_adj <= 0.05 & avg_log2FC >= 0.25) %>% head(200)%>% dplyr::pull(gene)
              readr::write_tsv(DE_table_receiver%>%dplyr::filter(p_val_adj <= 0.05 & avg_log2FC >= 0.25) %>% head(200) ,file=glue::glue("{output_dir}/input-target-gene-for-{condition_oi}-vs-{condition_reference}.tsv"))
              geneset_oi <-  geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
        } else {
            futile.logger::flog.info("input gene list or setting condition")
        }

        # ==============================================================================================================
        futile.logger::flog.info("step5: 确定潜在的高活性配体")
        ligands <- lr_network %>%
            dplyr::pull(from) %>%
            unique()
        receptors <- lr_network %>%
            dplyr::pull(to) %>%
            unique()
        expressed_ligands <- intersect(ligands, expressed_genes_sender)
        expressed_receptors <- intersect(receptors, expressed_genes_receiver)
        potential_ligands <- lr_network %>%
            dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
            dplyr::pull(from) %>%
            unique()

        # ===================================================================================================================
        futile.logger::flog.info("step6: Perform NicheNet’s ligand activity analysis on the gene set of interest")
        ligand_activities <- nichenetr::predict_ligand_activities(
            geneset = geneset_oi,
            background_expressed_genes = background_expressed_genes,
            ligand_target_matrix = ligand_target_matrix,
            potential_ligands = potential_ligands
        ) %>%
            dplyr::arrange(-pearson) %>%
            dplyr::mutate(rank = rank(dplyr::desc(pearson)))

        if (!is.null(opt$topn)) {
            topn <- opt$topn
            best_upstream_ligands <- ligand_activities %>%
                dplyr::top_n(topn, pearson) %>%
                dplyr::arrange(-pearson) %>%
                dplyr::pull(test_ligand) %>%
                unique()
        } else {
            best_upstream_ligands <- ligand_activities %>%
                dplyr::arrange(-pearson) %>%
                dplyr::pull(test_ligand) %>%
                unique()
        }
        ## 导出活性配体文件
        write.table(
            data.frame(ligand_activities),
            file = paste0(output_dir, "/ligand_activities.txt"),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F
        )

        # ==============================================================================================================
        futile.logger::flog.info("step7 :Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity anal")
        active_ligand_target_links_df <- best_upstream_ligands %>%
            lapply(
                nichenetr::get_weighted_ligand_target_links,
                geneset = geneset_oi,
                ligand_target_matrix = ligand_target_matrix,
                n = 200
            ) %>%
            dplyr::bind_rows() %>%
            tidyr::drop_na()

        active_ligand_target_links <- nichenetr::prepare_ligand_target_visualization(
            ligand_target_df = active_ligand_target_links_df,
            ligand_target_matrix = ligand_target_matrix,
            cutoff = opt$ligand_visiual_cutoff
        )

        order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>%
            rev() %>%
            make.names()
        order_targets <- active_ligand_target_links_df$target %>%
            unique() %>%
            intersect(rownames(active_ligand_target_links)) %>%
            make.names()
        rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names()
        colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names()

        vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()

        write.table(
            vis_ligand_target,
            file = paste0(output_dir, "/vis_ligand_target.txt"),
            col.names = T,
            row.names = T,
            sep = "\t",
            quote = F
        )

        # ===================================================================================================================
        futile.logger::flog.info("step8: Receptors of top-ranked ligands")
        lr_network_top <- lr_network %>%
            dplyr::filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
            dplyr::distinct(from, to)
        best_upstream_receptors <- lr_network_top %>%
            dplyr::pull(to) %>%
            unique()
        lr_network_top_df_large <- weighted_networks_lr %>%
            dplyr::filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
        lr_network_top_df <- lr_network_top_df_large %>%
            tidyr::spread("from", "weight", fill = 0)
        lr_network_top_matrix <- lr_network_top_df %>%
            dplyr::select(-to) %>%
            as.matrix() %>%
            magrittr::set_rownames(lr_network_top_df$to)

        dist_receptors <- dist(lr_network_top_matrix, method = "binary")
        hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
        order_receptors <- hclust_receptors$labels[hclust_receptors$order]
        dist_ligands <- dist(lr_network_top_matrix %>% t(), method = "binary")
        hclust_ligands <- hclust(dist_ligands, method = "ward.D2")

        order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]
        order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix))
        order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

        vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]
        rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
        colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()


        write.table(
            vis_ligand_receptor_network,
            file = paste0(output_dir, "/vis_ligand_receptor_network.txt"),
            col.names = T,
            row.names = T,
            sep = "\t",
            quote = F
        )
        # ===============================================================================================================
        futile.logger::flog.info("step9:combined plot of ligand")
        # ligand activity heatmap
        ligand_pearson_matrix <- ligand_activities %>%
            dplyr::select(pearson) %>%
            as.matrix() %>%
            magrittr::set_rownames(ligand_activities$test_ligand)
        rownames(ligand_pearson_matrix) <- rownames(ligand_pearson_matrix) %>% make.names()
        colnames(ligand_pearson_matrix) <- colnames(ligand_pearson_matrix) %>% make.names()
        vis_ligand_pearson <- ligand_pearson_matrix[order_ligands, ] %>%
            as.matrix(ncol = 1) %>%
            magrittr::set_colnames("Pearson")

       ## 图一： 配体活性打分热图
        p_ligand_pearson <- vis_ligand_pearson %>%
            nichenetr::make_heatmap_ggplot(
                glue::glue("{paste0(sender_celltypes,collapse = ' ')} Prioritized ligands"),
                "Ligand activity",
                color = "darkorange",
                legend_position = "top",
                x_axis_position = "top",
                legend_title = "Pearson correlation coefficient\ntarget gene prediction ability"
            ) +
            ggplot2::theme(legend.text = element_text(size = 9),
                           axis.ticks = ggplot2::element_blank(),
                           axis.title= element_text(size=13),
                           axis.text =   element_text(size =10,colour ="black"),
                           legend.key.width=unit(0.8,"cm") ,
                           legend.direction='horizontal',
                           legend.title = element_text(size=10,hjust =0))

        # ligand expression Seurat dotplotsi
        order_ligands_adapted <- order_ligands
        # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
        order_ligands_adapted[order_ligands_adapted == "H2.M3"] <- "H2-M3"
        # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
        order_ligands_adapted[order_ligands_adapted == "H2.T23"] <- "H2-T23"
        # ===================================================================================================================
        futile.logger::flog.info("step10: rotated_dotplot of ligands")
        T_F <- colnames(data_ob)[data_ob[[celltype]][, 1] %in% sender_celltypes]
        ## 图二
        rotated_dotplot <- Seurat::DotPlot(data_ob %>%
                                                subset(cells = T_F),
                                            dot.min = 0,
                                            dot.scale =9,
                                            features = order_ligands_adapted,
                                            cols = "RdYlBu") +
                                            ggplot2::coord_flip() +
                            #ggplot2::theme_minimal()+
                            ggplot2::theme(
                                legend.key.width=unit(0.8,"cm") ,
                                legend.direction='horizontal',
                                legend.text = ggplot2::element_text(size = 9),
                                legend.title = ggplot2::element_text(size = 10),
                                axis.text.x  = ggplot2:: element_text(size =10,  angle = 90,hjust = 0),
                                axis.title= element_text(size=13),
                                axis.text.y = ggplot2::element_text(size = 10),
                                axis.ticks = ggplot2::element_blank()) +
                            Seurat::RotatedAxis()+ # flip of coordinates necessary because we want to show ligands in the rows
                            ggplot2::theme( axis.text.x = element_text(size = 10,  angle = 90,hjust = 0))+
                            ylab("Expression in Sender") +
                            xlab("") +
                            scale_y_discrete(position = "right")+
                            guides(color = guide_colorbar(title = 'Average Expression',order = 1), fill = guide_legend(order = 0))

        ## 图三：配体-靶基因热图
        p_ligand_target_network <- vis_ligand_target %>%
            nichenetr::make_heatmap_ggplot(
              glue::glue("Prioritized ligands in {paste0(sender_celltypes,collapse = ' ')}"),
              glue::glue("Predicted target genes in {paste0(receiver,collapse = ' ')}"),
                color = "purple",
                legend_position = "top",
                x_axis_position = "top",
                legend_title = "Regulatory potential"
            ) +
            ggplot2::theme(legend.text = element_text(size = 9),
                           axis.title= element_text(size=13),
                           axis.ticks = ggplot2::element_blank(),
                           axis.text =   element_text(size =10,colour ="black"),
                           legend.key.width=unit(0.8,"cm") ,
                           legend.direction='horizontal',
                           legend.title = element_text(size=10,hjust =0))+
            ggplot2::scale_fill_gradient2(
                low = "whitesmoke",
                high = "purple",
                breaks = c(0, 0.0045, 0.0090)
            )
       #=============================================
        figures_without_legend <- cowplot::plot_grid(
            ##图一
            p_ligand_pearson +
                ggplot2::theme(legend.position = "none") +
                ggplot2::theme(axis.title.x = ggplot2::element_text()),

            ##图二
            rotated_dotplot + ggplot2::theme(legend.position = "none"),

            ##图三
            p_ligand_target_network +
                ggplot2::theme(legend.position = "none") +
                ggplot2::ylab(""),
            align = " hv",
            nrow = 1,
            rel_widths = c(1.5,
                           1+length(sender_celltypes) ,
                           ncol(vis_ligand_target )*0.4)
        )
        #=============================================
        legends <- cowplot::plot_grid(
            ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
            ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
            ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
            nrow = 1,
            #align = "h",
            rel_widths = c(1.3, 1, 1)
        )

        combined_plot <- cowplot::plot_grid(figures_without_legend,
                                            legends,
                                            rel_heights = c(nrow(vis_ligand_target ),2),
                                            nrow = 2,
                                            align = "hv")
        OESingleCell::save_ggplots(
            glue::glue("{output_dir}/combined_plot_target"),
            plot = combined_plot,
            width =1.5+1+length(sender_celltypes) +ncol(vis_ligand_target )*0.4,
            height = nrow(vis_ligand_target ) * 0.5 + 1,
            to.png = FALSE
        )

        # ===================================================================================================================
        futile.logger::flog.info("step11: 拼图  配体活性+配体在sender中表达热图+配体-受体作用热图")
        vis_ligand_pearson_recp <- ligand_pearson_matrix[colnames(vis_ligand_receptor_network), ] %>%
            as.matrix(ncol = 1) %>%
            magrittr::set_colnames("Pearson")
        vis_ligand_pearson_recp <- as.matrix(vis_ligand_pearson_recp[order(vis_ligand_pearson_recp[, 1]), ])
        colnames(vis_ligand_pearson_recp) <- "Pearson"
        ## 图一 配体活性热图
        p_ligand_pearson_recp <- vis_ligand_pearson_recp %>%
            nichenetr::make_heatmap_ggplot(
                glue::glue("{paste0(sender_celltypes,collapse = ' ')} Prioritized ligands"),
                "Ligand activity",
                color = "darkorange",
                legend_position = "top",
                x_axis_position = "top",
                legend_title = "Pearson correlation coefficient\nreceptor gene prediction ability"
            ) +
            ggplot2::theme(legend.text = element_text(size = 9),
                           axis.ticks = ggplot2::element_blank(),
                           axis.title= element_text(size=13),
                           axis.text =   element_text(size =10,colour ="black"),
                           legend.key.width=unit(0.8,"cm") ,
                           legend.direction='horizontal',
                           legend.title = element_text(size=10,hjust =0))

        futile.logger::flog.info("rotated_dotplot_recp of ligands")
        ## 图二 配体表达dotplot图
        rotated_dotplot_recp <- Seurat::DotPlot(
            data_ob %>% subset(cells = T_F),
            features = rownames(vis_ligand_pearson_recp),
            cols = "RdYlBu"
        ) +
            ggplot2::coord_flip() +
            Seurat::RotatedAxis()+
            #ggplot2::theme_minimal()+
            ggplot2::theme(
                legend.key.width=unit(0.8,"cm") ,
                legend.direction='horizontal',
                legend.text = ggplot2::element_text(size = 9),
                legend.title = ggplot2::element_text(size = 10),
                axis.text.x  = ggplot2:: element_text(size =10, angle = 90,hjust = 0),
                axis.title= element_text(size=13),
                axis.text.y = ggplot2::element_text(size = 10,colour ="black"),
                axis.ticks = ggplot2::element_blank()) +
            Seurat::RotatedAxis()+ # flip of coordinates necessary because we want to show ligands in the rows when
            ggplot2::theme( axis.text.x = element_text(size = 10, colour ="black", angle = 90,hjust = 0))+
            ylab("Expression in Sender") +
            xlab("") +
            scale_y_discrete(position = "right")+
            guides(color = guide_colorbar(title = 'Average Expression',order = 1), fill = guide_legend(order = 0))

        ## 图三 配体-受体热图
        vis_ligand_receptor_network_recp <- vis_ligand_receptor_network[, rownames(vis_ligand_pearson_recp)]
        p_ligand_receptor_network_recp <- vis_ligand_receptor_network_recp %>%
            t() %>%
            nichenetr::make_heatmap_ggplot(
                 glue::glue("Prioritized ligands in {paste0(sender_celltypes,collapse = ' ')}"),
                  glue::glue("Receptors in {paste0(receiver,collapse = ' ')}"),
                 color = "mediumvioletred",
                 x_axis_position = "top",
                 legend_title = "Prior interaction potential"
            )+
            ggplot2::theme(legend.text = element_text(size = 9),
                           axis.ticks = ggplot2::element_blank(),
                           axis.title= element_text(size=13),
                           axis.text =   element_text(size =10,colour ="black"),
                           legend.key.width=unit(0.8,"cm") ,
                           legend.direction='horizontal',
                           legend.title = element_text(size=10,hjust =0))

        ## 拼图
        figures_without_legend_recp <- cowplot::plot_grid(
            p_ligand_pearson_recp +ggplot2::theme(legend.position = "none"),
            rotated_dotplot_recp +ggplot2::theme(legend.position = "none"),
            p_ligand_receptor_network_recp +ggplot2::theme(legend.position = "none")+ ggplot2::ylab(""),
            align = "hv",
            nrow = 1,
            rel_widths = c(1.5,
                           1+length(sender_celltypes) ,
                           dim(vis_ligand_receptor_network_recp)[1]*0.4)
        )
        legends_recp <- cowplot::plot_grid(
            ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_recp)),
            ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot_recp)),
            ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor_network_recp)),
            #align = "h",
            nrow = 1,
            rel_widths = c(1.3, 1, 1)
        )

        combined_plot_recp <- cowplot::plot_grid(
            figures_without_legend_recp,
            legends_recp,
            rel_heights = c( dim(vis_ligand_receptor_network_recp)[2], 2),
            nrow = 2,
            align = "hv"
        )
        ###导出图片
        OESingleCell::save_ggplots(
            glue::glue("{output_dir}/combined_plot_receptor"),
            plot = combined_plot_recp,
            width = 1.5+1+length(sender_celltypes) +dim(vis_ligand_receptor_network_recp)[1]*0.4,
            height = dim(vis_ligand_receptor_network_recp)[2] * 0.5 + 1,
            to.png = FALSE
        )
        if(!file.exists(file.path(output_dir, "NicheNet细胞通讯分析说明.doc"))){
          file.copy("/home/xfyang/document/documents/NicheNet细胞通讯分析说明.docx",
          file.path(output_dir, "NicheNet细胞通讯分析说明.docx"))}
        ## output session information===============================================================================
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}
