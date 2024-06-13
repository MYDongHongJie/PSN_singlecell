#=======================================================================================================================
mia_plot <- function(mia_score_file,
                     output_dir,
                     annotation_legend_side = "top",
                     nrow = 2,
                     #width = 6,
                     #height = 4,
                     querryprefix,
                     refprefix,
                     reorder = "yes") {
  ##
  M <- read.table(mia_score_file, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  ##
  if (reorder == "yes") {
    colnames(M) <- colnames(M) %>%
                   as.data.frame() %>%
                   tidyr::separate(".", into = c("ref", "c2"), sep = "[(]") %>%
                   dplyr::select("ref") %>%
                   dplyr::mutate(ref = stringr::str_replace(ref, querryprefix, "")) %>%
                   unlist(use.names = FALSE)
    rownames(M) <- rownames(M) %>%
                   as.data.frame() %>%
                   tidyr::separate(".", into = c("ref", "c2"), sep = "[(]") %>%
      dplyr::select("ref") %>%
      dplyr::mutate(ref = stringr::str_replace(ref, refprefix, "")) %>%
      unlist(use.names = FALSE)
    M <- M %>%
      tibble::rownames_to_column("ref") %>%
      dplyr::arrange(as.numeric(ref)) %>%
      tibble::column_to_rownames("ref")
    M <- M[, order(as.numeric(colnames(M)))]
    write.table(x = M, file = paste0(output_dir, "/", "MIA_enrichment_for_plot.xls"), sep = "\t", col.names = NA)
  }
  width <- dim(M)[[2]] * 0.9
  height <- dim(M)[[1]] * 0.7
  #==============================================================================
  score_range <- M %>%
    dplyr::summarise_all(max) %>%
    quantile %>%
    dplyr::select(`75%`) %>%
    dplyr::mutate(`75%` = `75%` / 10) %>%
    ceiling() * 10
  col_fun <- circlize::colorRamp2(c(-(score_range$`75%`), 0, score_range$`75%`), c("navyblue", "white", "red"))
  st_clusters_name <- colnames(M)
  clusters.colors <- structure(names = st_clusters_name,
                               suppressMessages(OESingleCell::SelectColors(1:length(st_clusters_name),palette="ditto")))
  colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(clusters = factor(st_clusters_name, levels = st_clusters_name)),
                                              col = list(clusters = clusters.colors),
                                              #show_annotation_name = T,
                                              #annotation_name_side ="left",
                                              #simple_anno_size_adjust=T,
                                              annotation_legend_param = list(clusters = list(
                                                column_order = factor(st_clusters_name,
                                                                      levels = st_clusters_name), #固定legend标签顺序
                                                nrow = nrow, ##设置标签排序组合
                                                direction = "horizontal" ## 方向标签方向，水平或者垂直
                                              )
                                              ))
  pdf(NULL)
  p<-ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(ComplexHeatmap::Heatmap(as.matrix(M),
                                                                                           col = col_fun,
                                                                                           cluster_rows = F,
                                                                                           cluster_columns = F,
                                                                                           show_column_names = F,
                                                                                           column_title_side = "top",
                                                                                           show_row_names = T,
                                                                                           row_names_gp = grid::gpar(fontsize = 10),
                                                                                           row_names_rot = 0,
                                                                                           row_names_side = "left",
                                                                                           row_dend_width = unit(0.8, "mm"),
                                                                                           heatmap_legend_param = list(
                                                                                             legend_gp = grid::gpar(fontsize = 15),
                                                                                             title = "     -log10(pvalue)",
                                                                                             #labels_rot = 0,
                                                                                             #title_position = "",
                                                                                             at = c(-(score_range$`75%`), 0, score_range$`75%`),
                                                                                             labels = c(glue::glue("-{score_range$`75%`} \n depletion"), 0,
                                                                                                        glue::glue("{score_range$`75%`} \n enrichment")),
                                                                                             direction = "horizontal",
                                                                                             legend_width = unit(5, "cm")
                                                                                           ),
                                                                                           top_annotation = colAnn),
       merge_legend = FALSE,
       heatmap_legend_side = "bottom",
       annotation_legend_side = annotation_legend_side
    )))
  OESingleCell::save_ggplots(plot =p,
                             filename=glue::glue("{output_dir}/celltype_heatmap"),
                             width = width,
                             height = height)
  #=====================================================================================================================

  data_long <- M %>%
    tibble::rownames_to_column('rowname') %>%
    tidyr::gather(key = 'key', value = 'value', -rowname) %>%
    dplyr::group_by(key) %>%
    dplyr::slice_max(value) %>%
    tidyr::separate(key, c("key", "tmp1"), "\\.") %>%
    tidyr::separate(rowname, c("rowname", "tmp2"), "\\(") %>%
    dplyr::select(c("rowname", "key", "value"))
  colnames(data_long) <- c("source", "target", "value")
  data_long$target <- paste0(data_long$target, " ")

  #制作nodes
  nodes <- data.frame(name = c(as.character(data_long$source),
                               as.character(data_long$target)) %>% unique()
  )
  data_long$IDsource <- match(data_long$source, nodes$name) - 1
  data_long$IDtarget <- match(data_long$target, nodes$name) - 1
  #ColourScal='d3.scaleOrdinal() .range(oeplot::ditto_pal)'

  network <- networkD3::sankeyNetwork(Links = as.data.frame(data_long),
                                      Nodes = nodes,
                                      Source = "IDsource",
                                      Target = "IDtarget",
                                      Value = "value",
                                      NodeID = "name",
                                      sinksRight = FALSE,
                                      nodeWidth = 40,
                                      fontSize = 15,
                                      nodePadding = 10)
  networkD3::saveNetwork(network, paste0(output_dir, "/", "celltype_sankeyNetwork.html"), selfcontained = TRUE)
  webshot::webshot(paste0(output_dir, "/", "celltype_sankeyNetwork.html"),
                   paste0(output_dir, "/", "celltype_sankeyNetwork.png"))
  webshot::webshot(paste0(output_dir, "/", "celltype_sankeyNetwork.html"),
                   paste0(output_dir, "/", "celltype_sankeyNetwork.pdf"))

}

# ==================MIA results visulazitio=============================================================================

sub_mia_plot = subparsers$add_parser("mia", help = "MIA results visulazition.")
sub_mia_plot$add_argument("--file", type = "character",
                          help = "[REQUIRED] MIA output file: celltype_enrichment.xls ")
sub_mia_plot$add_argument("--reorder", type = "character", default = "yes",
                          help = "[OPTIONAL]reorder the columns.[default: %(default)s]")
sub_mia_plot$add_argument("--queryprefix", type = "character", default = "st_",
                          help = "[OPTIONAL]the prefix set for output refby .[default: %(default)s] ")
sub_mia_plot$add_argument("--refprefix", type = "character", default = "sc_",
                          help = "[OPTIONAL]the prefix set for output refby .[default: %(default)s] ")
args <- commandArgs(TRUE)
if ("mia" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "mia") {
    ##==================================================================================================================
    futile.logger::flog.info("running MIA analysis results visualize:")
    mia_plot(opt$file,
             output_dir,
             reorder = opt$reorder,
             querryprefix = opt$queryprefix,
             refprefix = opt$refprefix)
    if (!file.exists(file.path(output_dir, "MIA分析说明文档.docx"))) {
        file.copy("/public/dev_scRNA/oesinglecell3_test/test_util_mia/test_2023_01/MIA分析说明文档.docx", output_dir)
        message("MIA分析说明文档.docx 完成拷贝")
    }else {
        message("MIA分析说明文档.docx 已完成拷贝, 不再重复拷贝")
    }
    ##step: output session information==================================================================================
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
  }
}