# Title     : spatialDWLS空间转录组细胞类型鉴定
# Objective : TODO
# Created by: xiufeng.yang
# Created on: 2022/11/24

#' This function is to find ref marker for scRNAseq
#'
#' @param ref_ob  参考数据集seurat对象
#' @param ref_assay  对应assay
#' @param ref_markers 细胞类型对应的marker信息
#' @param ref_celltype 细胞类型对应列
#' @param test_mode 筛选方法选择
#' @param avg_log2FC marker基因筛选阈值
#' @param pvalue 筛选pvlaue设置
#' @param min_pct1  最小阈值
#' @param min_pct2  最大阈值
#' @param pct_fold  pct1/pct2最小值设置
#' @param n_topmarkers marker基因筛选
#' @return This function returns the result of spatialDWLS deconvolution
find_ref_marker <- function(ref_ob,
                            ref_assay = "RNA",
                            ref_marker = "auto",
                            ref_celltype = "celltype",
                            test_mode = "wilcox",
                            avg_log2FC = 1,
                            pvalue = 0.05,
                            min_pct1 = 0.25,
                            max_pct2 = 0.5,
                            pct_fold = 1,
                            n_topmarkers = 200) {
  #=====================================================================================================================
  futile.logger::flog.info(glue::glue("convert  reference scRNA-seq object to giotto_obj_scRNA "))
  ncelltypes <- length(unique(ref_ob[[ref_celltype, drop = TRUE]]))
  futile.logger::flog.info(glue::glue("Found {ncelltypes} cell types in the reference."))
  instrs <- Giotto::createGiottoInstructions(python_path = stringr::str_split(Sys.which("python"), "'")[[1]])
  sc_ob <- Giotto::createGiottoObject(instructions = instrs,
                                      raw_exprs = OESingleCell::GetAssayData(ref_ob, assay = ref_assay, slot = "counts"),
                                      norm_expr = OESingleCell::GetAssayData(ref_ob, assay = ref_assay, slot = "data"),
                                      cell_metadata = ref_ob[[ref_celltype]])
  #=====================================================================================================================
  futile.logger::flog.info("FindAllMarkers for reference scRNA-seq object data per celltype,This step use DEG function implemented in seurat")
  if (ref_marker != "auto") {
    refmarkers <- readr::read_tsv(ref_marker)
  } else {
    futile.logger::flog.info("利用findAllMarker来进行筛选marker基因")
    ref_ob <- Seurat::SetIdent(ref_ob, value = ref_celltype)
    refmarkers <- Seurat::FindAllMarkers(ref_ob,
                                         assay = ref_assay,
                                         slot = "data",
                                         test.use = test_mode,
                                         only.pos = T,
                                         logfc.threshold = 0,
                                         min.pct = 0)
    refmarkers <- refmarkers %>%
      dplyr::mutate(gene_diff = round(refmarkers$pct.1 / refmarkers$pct.2, 3)) %>%
      dplyr::select(gene, tidyselect::everything())
    write.table(refmarkers,
                file = glue::glue("{output_dir}/refmarkers_for_{ref_celltype}.xls"),
                col.names = T,
                row.names = F,
                sep = "\t",
                quote = F)
  }
  # Use top 100 markers as marker genes (or as much as there is)
  top_markers <- lapply(unique(refmarkers$cluster) %>% as.character(), function(celltype) {
    refmarkers %>%
      dplyr::filter(cluster == celltype,
                    avg_log2FC >= avg_log2FC &
                      p_val < pvalue &
                      pct.1 > min_pct1 &
                      pct.2 < max_pct2 &
                      gene_diff > pct_fold) %>%
      dplyr::arrange(p_val, dplyr::desc(avg_log2FC), dplyr::desc(gene_diff)) %>%
      dplyr::filter(gene_diff > pct_fold) %>%
      dplyr::top_n(n_topmarkers, gene_diff) %>%
      .$gene })
  names(top_markers) <- unique(refmarkers$cluster) %>% as.character()
  gene_list <- unlist(top_markers, use.names = FALSE)
  #Calculate median expression value of signature genes in each cell type
  norm_exp <- 2^(sc_ob@norm_expr) - 1
  id <- unique(refmarkers$cluster)
  ExprSubset <- norm_exp[gene_list,]
  Sig_exp <- NULL
  for (i in unique(id)) {
    Sig_exp <- cbind(Sig_exp, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
  }
  colnames(Sig_exp) <- unique(id)
  return(Sig_exp)
}

#' This function is to deconvolute ST data using the spatialDWLS method
#'
#' @param Sig_exp  参考数据集的细胞类型marker的list列表
#' @param spatial_ob  this is a seurat object
#' @param spatial_assay "Spatial"
#' @param slice  样本切片信息
#' @return This function returns the result of spatialDWLS deconvolution
#'
spatialDWLS_run <- function(Sig_exp,
                            spatial_ob,
                            spatial_assay = "Spatial",
                            n_cell = 20,
                            slice,
                            output_dir) {
  if (file.exists(glue::glue("{output_dir}/{slice}_Gitto_SpatilDWLS.rds"))) {
    st_data <- readRDS(glue::glue("{output_dir}/{slice}_Gitto_SpatilDWLS.rds"))
  }else {
    instrs <- Giotto::createGiottoInstructions(python_path = stringr::str_split(Sys.which("python"), "'")[[1]])
    #=====================================================================================================================
    futile.logger::flog.info("3):创建Giotto 对象")
    spatial_count <- OESingleCell::GetAssayData(spatial_ob, assay = spatial_assay, slot = "counts")
    spatial_locs <- Seurat::GetTissueCoordinates(spatial_ob, image = slice)
    #spatial_count  <- spatial_ob@assays$Spatial@counts %>% as.data.frame()
    #spatial_locs <-  Seurat::GetTissueCoordinates(spatial_ob,image = slice)
    #data pre-processing
    st_data <- Giotto::createGiottoObject(raw_exprs = spatial_count,
                                          spatial_locs = spatial_locs,
                                          instructions = instrs)

    futile.logger::flog.info("4):对空间数据进行过滤，均一化和聚类")
    st_data <- Giotto::filterGiotto(gobject = st_data,
                                    expression_threshold = 0,
                                    gene_det_in_min_cells = 0,
                                    min_det_genes_per_cell = 0,
                                    expression_values = c('raw'),
                                    verbose = T)
    st_data <- Giotto::normalizeGiotto(gobject = st_data)
    st_data <- Giotto::calculateHVG(gobject = st_data)
    gene_metadata <- Giotto::fDataDT(st_data)
    st_featgenes <- gene_metadata[which(gene_metadata$hvg == "yes"),]$gene_ID
    st_data <- Giotto::runPCA(gobject = st_data, genes_to_use = st_featgenes, scale_unit = F) ##主成分分析
    st_data <- Giotto::createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 10) ## 近邻网络构建
    st_data <- Giotto::doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000) ## leiden聚类分析

    futile.logger::flog.info("5):输入st_data及Sig_exp, 运行spatialDWLS")
    st_data <- Giotto::runDWLSDeconv(gobject = st_data, #gitto对象
                                     expression_values = c("normalized"), ## 采用的data为标准化后数据
                                     logbase = 2, ##log标准化的方法
                                     cluster_column = "leiden_clus",
                                     sign_matrix = Sig_exp,
                                     n_cell = n_cell, ##每个spots包含细胞数目
                                     cutoff = 2, ##
                                     name = NULL, ##结果存储对象名称
                                     return_gobject = TRUE)
    saveRDS(st_data, file = glue::glue("{output_dir}/{slice}_Gitto_SpatilDWLS.rds"))
  }
  decon_mtrx <- st_data@spatial_enrichment$DWLS %>%
    tibble::column_to_rownames("cell_ID")

  return(decon_mtrx)
}

#######################################################################################################################
docstring <- "example1:\\n\\
sctool  -i  st.rds   -f rds  -o results  -d rds   --assay  Spatial  --dataslot counts spatialDWLS  --refexp sc.rds --refassay RNA --refcelltype  celltype --refmarker refmarkers.xls  \\n\\
example2:\\n\\
sctool  -i  st.rds   -f rds  -o results  -d rds   --assay  Spatial  --dataslot counts spatialDWLS  --refexp sc.rds --refassay RNA --refcelltype  celltype --refmarker auto   "
sub_spatialDWLS <- subparsers$add_parser("spatialDWLS",
                                         description = docstring,
                                         formatter_class = "argparse.RawTextHelpFormatter",
                                         # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                         argument_default = "True",
                                         help = "Using CARD to refer Spaital transcriptome cell type")
sub_spatialDWLS$add_argument("--refexp", type = "character", default = NULL,
                             help = paste0("the customized reference expression matrix in seurat format, which may come from the",
                                           " quantification results of microarray, Bulk RNA-seq and scRNA-seq sequencing."))
sub_spatialDWLS$add_argument("--refassay", type = "character", default = NULL,
                             help = paste0("the ref assay in data object to use. When it comes to multimodal assay, this is the assay used to",
                                           "initialize the object, all the other assays will merged into it. "))
sub_spatialDWLS$add_argument("--refcelltype", type = "character", help = "the cell type annotation column id in the reference expression matrix.")
sub_spatialDWLS$add_argument("--refmarker", default = "auto", help = "the  refmarker  file  for reference data to use for the this run")

sub_spatialDWLS$add_argument("--min_pct1", type = "double", default = 0.25,
                             help = "the minimium ratio of cells expressing one specific gene in a cluster.[default: %(default)s] ")
sub_spatialDWLS$add_argument("--max_pct2", type = "double", default = 0.5,
                             help = "the maximiium ratio of cells expressing one specific gene in all other clusters.[default: %(default)s] ")
sub_spatialDWLS$add_argument("--pct_fold", type = "double", default = 1,
                             help = "the minimiu fold of pct1 for gene in a specific cluster against pct2 for all other cluster.[default: %(default)s] ")
sub_spatialDWLS$add_argument("--avg_log2FC", type = "double", default = 0.58,
                             help = "The average log2FC of the gene UMI count in its cluster against the all other clusters.[default: %(default)s] ")
sub_spatialDWLS$add_argument("--pvalue", type = "double", default = 0.05,
                             help = "the P-value of the gene differential expression.[default: %(default)s] ")
sub_spatialDWLS$add_argument("--test", type = "character", default = "wilcox",
                             help = "test methods used to cell markers. Options are: wilcox,presto,venice, t, bimod,poisson, negbinom, MAST, DESeq2, DESeq2LRT,limma, edgeR.[default: %(default)s] ")
sub_spatialDWLS$add_argument("--n_topmarkers", type = "character", default = "200",
                             help = "top n markers setting for reference object.[default: %(default)s] ")

sub_spatialDWLS$add_argument("--n_cell", default = "20",
                             help = "setting  n cells for each spot.[default: %(default)s] ")
# === Subcmd: spatialDWLS, deconvolute the cell type composition for spatial transcriptomics  using scRNA-seq data ========
args <- commandArgs(TRUE)
if ("spatialDWLS" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "spatialDWLS") {
    set.seed(100)
    # ==================================================================================================================
    futile.logger::flog.info("step1:导入空间转录组seurat对象：{opt$input}")
    suppressMessages(data_ob <- OESingleCell::ReadX(
      input = opt$input,
      informat = opt$informat,
      assays = assays[1],
      data.use = dataslots,
      verbose = F
    ))
    Seurat::DefaultAssay(data_ob) <- assays[1]
    cellmeta <- OESingleCell::colData(data_ob)
    if (!is.null(Seurat::Images(data_ob))) {
      unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique)]
      if (length(unuse_images) > 0) { data_ob@images[unuse_images] <- NULL }
    }
    # ==================================================================================================================
    futile.logger::flog.info("step2:导入参考单细胞转录组seurat对象：{opt$refexp}, 并依据提供参数，筛选细胞类型特异性marker")
    suppressMessages(ref_ob <- OESingleCell::ReadX(
      input = opt$refexp,
      informat = "rds",
      assays = opt$refassay,
      data.use = "counts,data",
      verbose = F
    ))

    # the input is a seurat object which may contain more than one sample
    SeuratObject::DefaultAssay(ref_ob) <- opt$refassay
    ref_ob <- SeuratObject::SetIdent(ref_ob, value = opt$refcelltype)
    Sig_exp <- find_ref_marker(ref_ob = ref_ob,
                               ref_assay = opt$refassay,
                               ref_marker = opt$refmarker,
                               ref_celltype = opt$refcelltype,
                               test_mode = opt$test,
                               avg_log2FC = opt$avg_log2FC,
                               pvalue = opt$pvalue,
                               min_pct1 = opt$min_pct1,
                               max_pct2 = opt$max_pct2,
                               pct_fold = opt$pct_fold,
                               n_topmarkers = opt$n_topmarkers)

    # ==================================================================================================================
    futile.logger::flog.info("step3:分别针对每个空间样本切片进行spatialDWLS分析")
    images <- Seurat::Images(data_ob)
    celltype_out <- list()
    # parallel::mclapply(images, function(slice) {存在bug,暂不可行
    for (slice in images) {
      data_ob_sub <- base::subset(data_ob, subset = sampleid == as.character(slice))
      if (!is.null(Seurat::Images(data_ob_sub))) {
        unuse_images <- Seurat::Images(data_ob_sub)[!Seurat::Images(data_ob_sub) %in% slice]
        if (length(unuse_images) > 0) { data_ob_sub@images[unuse_images] <- NULL }
      }
      celltype_out[[slice]] <- spatialDWLS_run(Sig_exp = Sig_exp,
                                               spatial_ob = data_ob_sub,
                                               spatial_assay = "Spatial",
                                               n_cell = as.numeric(opt$n_cell),
                                               slice = slice,
                                               output_dir = output_dir)
    }
    #}, mc.cores = 4)
    # ==================================================================================================================
    futile.logger::flog.info("step4: 将spatialDWLS输出的各个样本进行汇总，并存入misc@spatialDWLS_results对象中")
    out_celltype <- dplyr::bind_rows(celltype_out) %>%
      data.frame() %>%
      tibble::rownames_to_column("barcodes")
    write.table(
      out_celltype,
      file = glue::glue("{output_dir}/all-celltype.xls"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    data_ob <- subset(data_ob, cells = out_celltype$barcodes)
    ## save results to seurat object's mis part
    data_ob@misc$spatialDWLS_results <- out_celltype
    maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
    metadata <- data_ob@misc$spatialDWLS_results %>%
      tibble::column_to_rownames("barcodes") %>%
      dplyr::mutate(
        top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
        top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
      )
    data_ob <- Seurat::AddMetaData(data_ob, metadata = metadata)
    if (tolower(opt$outformat) == "h5seurat" & as.logical(opt$update)) {
      SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob)
    } else {
      OESingleCell::SaveX(
        data_ob,
        output = opt$output,
        update = FALSE,
        outformat = opt$outformat,
        prefix = opt$prefix
      )
    }
    #==============================================================================================================
    ## save session information
    futile.logger::flog.info("step5：保存日志文件")
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}
