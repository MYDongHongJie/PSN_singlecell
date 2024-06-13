#=======================================================================================================================
docstring <- " example1:\\n\\n\\
sctool  -i st.rds -f rds --assay SCT -o ./ celltrek  --refobject sc.rds  --refassay RNA --groupby new_celltype "
sub_celltrek <- subparsers$add_parser(
  "celltrek",
  description = docstring,
  formatter_class = "argparse.RawTextHelpFormatter",
  # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
  argument_default = "True",
  help = "Using celltrek to transfer scRNA-seq data to spatial"
)
sub_celltrek$add_argument(
  "--refobject",
  type = "character",
  help = "the customized reference expression matrix in seurat format."
)
sub_celltrek$add_argument(
  "--refassay",
  type = "character",
  default = "RNA",
  help = "the default assay for reference data to use for the this run.[default: %(default)s]"
)
sub_celltrek$add_argument(
  "--groupby",
  type = "character",
  help = "choose which metadata info of ref object to run celltrek."
)
sub_celltrek$add_argument(
  "--reduction",
  type = "character",
  default = "pca",
  help = "Dimension reduction method, usually 'pca'.[default: %(default)s]."
)
sub_celltrek$add_argument(
  "--top_spot",
  type = "integer",
  default = 10,
  help = "maximum number of spots that one cell can be charted.[default: %(default)s]."
)
sub_celltrek$add_argument(
  "--spot_n",
  type = "integer",
  default = 10,
  help = "maximum number of cells that one spot can contain.[default: %(default)s]."
)
sub_celltrek$add_argument(
  "--spointsize",
  type = "double",
  default = 0.8,
  help = "[OPTIONAL]the point size in the spatialplot.[default: %(default)s]."
)
sub_celltrek$add_argument(
  "--palette",
  type = "character",
  default = "ditto",
  help = paste0(
    "the discrete color schema mapped to the cell annotations specified by --groupby.[default: %(default)s]",
    " Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,",
    " tableau20:20,Buen:17,UKBB:18,TF1:17,paired:12"
  )
)
sub_celltrek$add_argument(
  "--seed",
  type = "integer",
  default = 20,
  help = "set random number to generate repeatable results.[default: %(default)s]."
)
sub_celltrek$add_argument(
    "--crop",
    default = 'FALSE',
    help = "whether to crop in spatialplot, data from cytassist project should be 'TRUE'.[default: %(default)s]"
)
#=======================================================================================================================
args <- commandArgs(TRUE)
if ("celltrek" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "celltrek") {
    futile.logger::flog.info("step1: load query and reference object data")
    #load query object data=============================================================================================
    suppressMessages(data_ob <- OESingleCell::ReadX(
      input = opt$input,
      informat = opt$informat,
      assays = assays[1],
      data.use = dataslots,
      verbose = F
    ))
    #load reference object data=========================================================================================
    suppressMessages(ref_ob <- OESingleCell::ReadX(
      input = opt$refobject,
      informat = "rds",
      assays = opt$refassay,
      data.use = "data,counts",
      verbose = F
    ))

    #get subset of cells if necessary===================================================================================
    if (!is.null(opt$predicate)) {
       futile.logger::flog.info(glue::glue("step1.5: get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
       #factors <- stringr::str_split(opt$predicate, ",")[[1]]
       df <- slot(data_ob, "meta.data")
       desired_cells <- subset(df, eval(parse(text = opt$predicate)))
       #desired_cells <- subset(df, df[[factors[1]]] %in% factors[2:length(factors)])
       data_ob <- subset(data_ob, cells = rownames(desired_cells))
    }

    #run for each image=================================================================================================
    sampleid <- unique(data_ob@meta.data$sampleid)
    for (i in sampleid) {
      output_celltrek <- glue::glue("{output_dir}/celltrek_results/celltrek_{i}.rds")
      if (!file.exists(output_celltrek)) {
        #data preprocessing===============================================================================================
        futile.logger::flog.info(glue::glue("step2: data preprocessing"))
        sub_data_ob <- subset(data_ob, subset = sampleid == i)
        if (!is.null(Seurat::Images(sub_data_ob))) {
          unuse_images <- Seurat::Images(sub_data_ob)[!Seurat::Images(sub_data_ob) %in% (sub_data_ob@meta.data$sampleid %>% unique())]
          if (length(unuse_images) > 0) {
            sub_data_ob@images[unuse_images] <- NULL
          }
        }
        sub_data_ob <- SeuratObject::RenameCells(sub_data_ob, new.names = make.names(SeuratObject::Cells(sub_data_ob)))
        ref_ob <- SeuratObject::RenameCells(ref_ob, new.names = make.names(SeuratObject::Cells(ref_ob)))

        if (!"orig.ident" %in% colnames(ref_ob@meta.data)) {
          metadata <- ref_ob@meta.data %>% dplyr::mutate(orig.ident = rownames(ref_ob@meta.data))
          ref_ob <- SeuratObject::AddMetaData(ref_ob, metadata$orig.ident, "orig.ident")
        }
        SeuratObject::Idents(ref_ob) <- ref_ob[[opt$groupby]]

        #data integratation===============================================================================================
        futile.logger::flog.info(glue::glue("step3: scRNA-seq and ST data integratation"))
        traint <- CellTrek::traint(st_data = sub_data_ob,
                                   sc_data = ref_ob,
                                   sc_assay = opt$refassay,
                                   cell_names = opt$groupby)

        dimplot <- Seurat::DimPlot(traint, group.by = "type")
        OESingleCell::save_ggplots(filename = glue::glue("{output_dir}/celltrek_results/1.intergreated_results/{i}/{i}_dimplot"),
                     plot = dimplot,
                     limitsize = FALSE,
                     height = 6,
                     dpi = 600)

        #run celltrek=====================================================================================================
        futile.logger::flog.info(glue::glue("step4: run celltrek"))
        celltrek <- CellTrek::celltrek(st_sc_int = traint, #Seurat traint object.  st和sc的整合数据
                                       int_assay = 'traint', #Integration assay ('traint'). 整合数据的assay
                                       sc_data = ref_ob, #SC data, optional. 单细胞数据集
                                       sc_assay = opt$refassay, #SC assay. 应使用标准化后的assay,
                                       reduction = opt$reduction, #Dimension reduction method, usually 'pca'. 降维方法，一般使用pca
                                       intp = T, #If True, do interpolation. 是否对单细胞数据进行插值
                                       intp_pnt = 10000, #Number of interpolation points. 插值数目
                                       intp_lin = F, #If TRUE, use linear interpolation. 是否线性插值
                                       nPCs = 30, #Number of PCs used for CellTrek. 用于celltrek的维度数
                                       ntree = 1000, #Number of trees in random forest. 随机森林的决策树数目
                                       dist_thresh = 0.4, #Distance threshold.
                                       top_spot = opt$top_spot, #Maximum number of spots that one cell can be charted. 单个细胞可以映射的最大spots数目
                                       spot_n = opt$spot_n, #Maximum number of cells that one spot can contain. 单个spot可以包含的最大细胞细胞数
                                       repel_r = 5, #Repelling radius
                                       repel_iter = 10, #Repelling iterations
                                       keep_model = T #是否保留随机森林原始模型信息
        )$celltrek
        saveRDS(celltrek, glue::glue("{output_dir}/celltrek_results/celltrek_{i}.rds"))
      }else {
        celltrek <- readRDS(output_celltrek)
      }

      #plot for celltrek results========================================================================================
      futile.logger::flog.info(glue::glue("step5: visualization of celltrek results"))
      crop <- as.logical(opt$crop)
      metadata_sampleid <- celltrek@meta.data %>% dplyr::mutate(sampleid = i)
      celltrek <- SeuratObject::AddMetaData(celltrek, metadata_sampleid$sampleid, "sampleid")
      celltrek[[opt$groupby]] <- factor(celltrek[[opt$groupby]][,1], levels = unique(celltrek[[opt$groupby]][,1]))
      colors <- OESingleCell::SelectColors(unique(celltrek[[opt$groupby]][,1]), palette = opt$palette)
      p0 <- OESingleCell::SpatialPlot(celltrek,
                        pt.size.factor = opt$spointsize,
                        assay = opt$refassay,
                        group.by = opt$groupby,
                        ncol = 1,
                        images = SeuratObject::Images(celltrek),
                        image.alpha = 0.7,
                        theme_void = TRUE,
                        combine = FALSE,
                        crop = crop,
                        cols = opt$palette,
                        min.cutoff = "q10",
                        max.cutoff = "q95")
      OESingleCell::save_ggplots(filename = glue::glue("{output_dir}/celltrek_results/1.intergreated_results/{i}/{i}_all_spatialplot"),
                   plot = p0[[1]],
                   limitsize = FALSE,
                   height = 6,
                   dpi = 600)

      metadata_newgroup <- celltrek@meta.data %>% dplyr::mutate(groupby.new = celltrek[[opt$groupby]][,1])
      celltrek <- SeuratObject::AddMetaData(celltrek, metadata_newgroup$groupby.new, "groupby.new")
      for (ct in unique(celltrek$groupby.new)) {
        p <- OESingleCell::SpatialPlot(subset(celltrek, groupby.new == ct),
                         pt.size.factor = opt$spointsize,
                         assay = opt$refassay,
                         group.by = "groupby.new",
                         ncol = 1,
                         images = SeuratObject::Images(celltrek),
                         image.alpha = 0.7,
                         theme_void = TRUE,
                         combine = FALSE,
                         crop = crop,
                         cols = c(colors[[ct]], "black"),
                         min.cutoff = "q10",
                         max.cutoff = "q95")
        OESingleCell::save_ggplots(filename = glue::glue("{output_dir}/celltrek_results/1.intergreated_results/{i}/{i}_{ct}_spatialplot"),
                     plot = p[[1]],
                     limitsize = FALSE,
                     height = 6,
                     dpi = 600)
      }

      #run scoloc function==============================================================================================
      futile.logger::flog.info(glue::glue("step6: scoloc analysis"))
      graph_KL <- CellTrek::scoloc(celltrek,
                                   col_cell = opt$groupby,
                                   use_method = 'KL',
                                   eps=1e-50)
      scoloc_result <- graph_KL$mst_cons %>% tibble::rownames_to_column(var = "id")
      if (!dir.exists(glue::glue("{output_dir}/celltrek_results/2.celltrek_scoloc/{i}"))) {
        dir.create(glue::glue("{output_dir}/celltrek_results/2.celltrek_scoloc/{i}"), recursive = T)
      }
      write.table(scoloc_result,
                  glue::glue("{output_dir}/celltrek_results/2.celltrek_scoloc/{i}/{i}_scoloc_results.xls"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)

      #generate igraph object===========================================================================================
      futile.logger::flog.info(glue::glue("step7: prepare for igraph object"))
      #作图前处理，获取细胞类型信息
      rownames(graph_KL$mst_cons) <- stringi::stri_replace_all_fixed(rownames(graph_KL$mst_cons),
                                                                     c(".","/","-","_"),
                                                                     " ",
                                                                     vectorize_all = FALSE)
      colnames(graph_KL$mst_cons) <- stringi::stri_replace_all_fixed(colnames(graph_KL$mst_cons),
                                                                     c(".","/","-","_"),
                                                                     " ",
                                                                     vectorize_all = FALSE)
      graph_KL_mst_cons <- graph_KL$mst_cons
      cell_class <- celltrek@meta.data %>%
        dplyr::select(id = opt$groupby) %>%
        unique
      cell_class$id <- stringi::stri_replace_all_fixed(cell_class$id, c(".","/","-","_"), " ", vectorize_all = FALSE)
      celltrek_count <- data.frame(freq = table(celltrek[[opt$groupby]]))
      colnames(celltrek_count) <- c("freq.Var1", "freq.Freq")
      celltrek_count$freq.Var1 <- stringi::stri_replace_all_fixed(celltrek_count$freq.Var1, c(".","/","-","_"), " ", vectorize_all = FALSE)
      cell_class_new <- merge(cell_class, celltrek_count, by.x = "id", by.y = "freq.Var1")

      #定义网络图节点
      mst_cons_node <- data.frame(id = rownames(graph_KL_mst_cons), label = rownames(graph_KL_mst_cons))
      col_colmn <- OESingleCell::SelectColors(mst_cons_node$id, palette = opt$palette)
      mst_cons_node$color <- col_colmn
      size_df <- data.frame(id = cell_class_new$id,
                            value = as.numeric(cell_class_new[, "freq.Freq"]) / sum(as.numeric(cell_class_new[, "freq.Freq"])))
      mst_cons_node <- dplyr::left_join(mst_cons_node, size_df, by = "id") %>% data.frame

      #定义网络图边界
      graph_KL_mst_cons <- data.frame(id = rownames(graph_KL_mst_cons), graph_KL_mst_cons, check.names = F)
      mst_cons_edge <- reshape2::melt(graph_KL_mst_cons) %>%
        na.omit() %>%
        magrittr::set_colnames(c('from', 'to', 'value'))
      mst_cons_edge <- mst_cons_edge[mst_cons_edge$value > 0,]
      #构建网络图对象
      graph <- igraph::graph_from_data_frame(mst_cons_edge, directed = TRUE, vertices = mst_cons_node)

      #plot for igraph object===========================================================================================
      futile.logger::flog.info(glue::glue("step8: make plot for igraph object"))
      set.seed(opt$seed) #生成随机数，使图的布局可重复
      l <- igraph::layout.fruchterman.reingold(graph) #设置图的布局方式为弹簧式发散的布局
      igraph::V(graph)$size <- 50 * mst_cons_node$value + 5 #degree(graph)*2+30  #节点大小与该细胞类型的细胞总数相关
      igraph::V(graph)$color <- OESingleCell::SelectColors(unique(celltrek[[opt$groupby]][,1]), palette = opt$palette) #根据类型设置颜色,按照类型分组
      igraph::V(graph)$label.color <- 'black' #设置节点标记的颜色
      igraph::V(graph)$frame.color <- NA
      igraph::V(graph)$label.degree <- -pi / 2
      igraph::V(graph)$label.dist <- -1.2
      igraph::E(graph)$width <- 5 * igraph::E(graph)$value #根据频次列设置边宽度
      igraph::E(graph)$label <- NA #根据频次列设置边标签
      igraph::E(graph)$arrow.width <- 0
      igraph::E(graph)$arrow.size <- 0 #设置箭头大小
      pdf(file = glue::glue("{output_dir}/celltrek_results/2.celltrek_scoloc/{i}/{i}_scoloc_network_plot.pdf"))
      plot(graph, layout = l)
      dev.off()
      png(file = glue::glue("{output_dir}/celltrek_results/2.celltrek_scoloc/{i}/{i}_scoloc_network_plot.png"))
      plot(graph, layout = l)
      dev.off()
    }
    if (!file.exists(file.path(output_dir, "CellTrek-单细胞转录组联合空间转录组分析说明文档.docx"))) {
        file.copy("/public/dev_scRNA/oesinglecell3_test/document/CellTrek-单细胞转录组联合空间转录组分析说明文档.docx",
                  file.path(output_dir, "CellTrek-单细胞转录组联合空间转录组分析说明文档.docx")) }
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}
