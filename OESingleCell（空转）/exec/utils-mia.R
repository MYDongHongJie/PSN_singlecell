#' filter genes specified by user or according to the cell number where one gene can be detected
#' @param object the data object: SingleCellExperiment or Seurat
#' @param slot the matrix slot to use, default to 'counts'.
#' @param min.value the minimum counts across each cell
#' @param min.cells the minimum cell number one gene can be detected above zero.
#' @param filter.genes a vector of gene names to be removed directly.
#'
#' @importFrom Matrix rowSums
#' @importFrom Seurat GetAssayData CaseMatch
#'
#' @export
FilterGenes <- function(
  object,
  assay = "RNA",
  slot = "counts",
  min.value = 1,
  min.cells = 0,
  filter.genes = NULL
) {
  genes.use <- rownames(object)
  Seurat::DefaultAssay(object = object) <- assay
  if (min.cells > 0) {
    num.cells <- Matrix::rowSums(Seurat::GetAssayData(object, assay = assay, slot = slot) >= min.value)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    counts <- Seurat::GetAssayData(object, assay = assay, slot = slot)[genes.use,]
    object[[assay]] <- Seurat::CreateAssayObject(counts = counts)
    # object = subset( object, features = genes.use)
    # object = SetAssayData(object, new.data = GetAssayData(object, slot ="data")[genes.use,])
    # object@data <- object@data[genes.use, ] # Seurat V2.x
  }
  if (!is.null(filter.genes)) {
    filter.genes <- Seurat::CaseMatch(search = filter.genes, match = rownames(object))
    genes.use <- setdiff(genes.use, filter.genes) #keep genes not in filter.genes
    object <- object[genes.use,]
  }
  # object <- LogSeuratCommand(object)
  return(object)
}

marker_gene_matrix <- function(
  seurat_object_file,
  format,
  assay,
  marker_file,
  filter = "yes",
  groupby = "cluster",
  test.use = "t",
  logfc.threshfold = 1,
  min_pct=0.5,
  max_pct=0.5,
  prefix="ref",
  out_dir,
  pct_fold,
  pval=0.05,
  recorrect_umi = TRUE,
  predicate = NULL,
  sub_obj = NULL
) {
  suppressMessages(seurat_ob<- OESingleCell::ReadX(input =seurat_object_file,
                                                   informat = format,
                                                   assays = assay,
                                                   data.use = "data",
                                                   verbose = F))
  ##  get raw marker genes
  if(is.null(marker_file)){
      if (groupby %in% colnames(seurat_ob@meta.data)) {
          seurat_ob <- Seurat::SetIdent(seurat_ob, value = groupby)
      }
      if (!is.null(seurat_ob@assays$SCT)) {
          if ((length(seurat_ob@assays$SCT@SCTModel.list) > 1) & is.null(seurat_ob@misc$PreSCT)) {
              ## https://rdrr.io/cran/Seurat/man/PrepSCTFindMarkers.html
              futile.logger::flog.info("Prepare object to run differential expression on SCT assay ")
              seurat_ob <- Seurat::PrepSCTFindMarkers(seurat_ob, assay = "SCT", verbose = TRUE)
              seurat_ob@misc$PreSCT <- TRUE
          }
      }
      if (!is.null(predicate)) {
          futile.logger::flog.info(glue::glue("step1.5: get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
          #factors <- stringr::str_split(opt$predicate, ",")[[1]]
          df <- slot(seurat_ob, "meta.data")
          desired_cells <- subset(df, eval(parse(text = predicate)))
          print(head(desired_cells))
          seurat_ob <- subset(seurat_ob, cells = rownames(desired_cells))
          seurat_ob@meta.data$sampleid <- factor(seurat_ob@meta.data$sampleid, levels = unique(seurat_ob@meta.data$sampleid))
      }
      if (!is.null(Seurat::Images(seurat_ob))) {
          unuse_images <- Seurat::Images(seurat_ob)[!Seurat::Images(seurat_ob) %in% (seurat_ob@meta.data$sampleid %>% unique())]
          if (length(unuse_images) > 0) {
              seurat_ob@images[unuse_images] <- NULL
          }
      }
      if (!is.null(sub_obj)) {
          futile.logger::flog.info("Get subseted object after PreSCT")
          suppressMessages(seurat_sub <- OESingleCell::ReadX(
                  input = sub_obj,
                  informat = format,
                  assays = assay,
                  verbose = FALSE
          ))
          seurat_ob <- subset(seurat_ob, cells = rownames(seurat_sub@meta.data))
          if (!is.null(Seurat::Images(seurat_ob))) {
              unuse_images <- Seurat::Images(seurat_ob)[!Seurat::Images(seurat_ob) %in% (seurat_ob@meta.data$sampleid %>% unique())]
              if (length(unuse_images) > 0) {
                  seurat_ob@images[unuse_images] <- NULL
              }
          }
          seurat_ob@meta.data <- seurat_sub@meta.data
          if (!groupby == "clusters") {
              SeuratObject::Idents(seurat_ob) <- seurat_ob[[groupby]]
          }
      }
      recorrect_umi <- as.logical(recorrect_umi)
      markers_raw <- Seurat::FindAllMarkers(object = seurat_ob,
                                            slot = "data",
                                            assay = assay,
                                            only.pos = T,
                                            test.use = test.use,
                                            recorrect_umi = recorrect_umi,
                                            logfc.threshold = 0,
                                            min.pct = 0)
      markers_raw <- markers_raw %>%
              # dplyr::rename_with(~groupby, cluster) %>%
              dplyr::mutate(gene_diff = round(markers_raw$pct.1 / markers_raw$pct.2, 3)) %>%
              dplyr::select(gene, tidyselect::everything())
      write.table(markers_raw,
                  file = paste0(out_dir, "/", prefix, test.use, "_all_markers_for_each_cluster.xls"),
                  col.names = T,
                  row.names = F,
                  sep = "\t",
                  quote = F)
  }else {
      markers_raw <- readr::read_tsv(file = marker_file, col_types = readr::cols(cluster = readr::col_character()))
      if (!"cluster" %in% colnames(markers_raw)) {
          markers_raw <- markers_raw %>% dplyr::mutate(cluster = groupby)
      }
      if (!"avg_log2FC" %in% colnames(markers_raw)) {
          markers_raw <- markers_raw %>% dplyr::mutate(avg_log2FC = avg_logFC) #适配seurat3产出的marker.list
      }
  }

    genes <- rownames(Seurat::GetAssayData(FilterGenes(seurat_ob, assay = assay, min.cells = 1, filter.genes = NULL),
                                           assay = assay))
    if (filter == "yes") {
        markers <- markers_raw %>%
                dplyr::group_by(cluster) %>%
                dplyr::filter(avg_log2FC >= logfc.threshfold) %>%
                dplyr::filter(p_val < pval) %>%
                dplyr::filter(pct.1 > min_pct & pct.2 < max_pct) %>%
                dplyr::arrange(p_val, plyr::desc(avg_log2FC), plyr::desc(gene_diff)) %>%
                dplyr::filter(gene_diff > pct_fold)
    }else if (filter == "no") {
        markers <- markers_raw
    }
    match <- Seurat::CaseMatch(search = as.vector(markers$gene), match = rownames(seurat_ob))
    markers <- markers %>% dplyr::filter(gene %in% names(match))
    write.table(x = markers,
                file = paste0(out_dir, "/", prefix, "marker.xls"),
                col.names = T,
                row.names = F,
                sep = "\t",
                quote = F)
    cluster_info <- markers %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarize(n = dplyr::n()) %>%
            dplyr::mutate(cluster_new = paste0(prefix, cluster, "(", n, " genes)"))
    markers <- markers %>%
            dplyr::left_join(cluster_info, by = "cluster") %>%
            dplyr::mutate(cluster = cluster_new)
    groups <- unique(markers$cluster)
    seurat_cluster <- tibble::tibble(gene = rep(genes, length(groups)),
                                     cluster = rep(groups, each = length(genes)))
    groupby_markers <- markers %>%
            dplyr:::select(gene, cluster) %>%
            unique() %>%
            dplyr::mutate(logic = 1) %>%
            dplyr::right_join(y = seurat_cluster, by = c("gene" = "gene", "cluster" = "cluster")) %>%
            tidyr::spread(key = "cluster", value = "logic", fill = 0) %>%
            dplyr::arrange(gene) %>%
            tibble::column_to_rownames(var = "gene") %>%
            as.matrix()
    return(groupby_markers)
}

#=======================================================================================================================
run_mia <- function(ref_marker, query_marker, output_dir) {
  M <- matrix(ncol = ncol(ref_marker), nrow = ncol(query_marker))
  M_enr <- matrix(ncol = ncol(ref_marker), nrow = ncol(query_marker))
  M_dep <- matrix(ncol = ncol(ref_marker), nrow = ncol(query_marker))
  geneList <- rownames(query_marker)
  b <- length(geneList)
  for (region in 1:ncol(query_marker)) {
    #region=1
    G <- rownames(query_marker[which(query_marker[, region] == 1),])
    for (celltype in 1:ncol(ref_marker)) {
      #celltype=2
      C <- rownames(ref_marker[which(ref_marker[, celltype] == 1),])
      a <- length(intersect(G, C))
      c <- length(G)
      d <- length(C)
      M_enr[region, celltype] <- round(-1 * log10(phyper(c - a - 1, c, b - c, b - d)), 2)
      M_dep[region, celltype] <- round(-1 * log10(1 - (phyper(c - a - 1, c, b - c, b - d))), 2)

      if (M_enr[region, celltype] < M_dep[region, celltype]) {
        M[region, celltype] <- -1 * M_dep[region, celltype]
      }
      else {
        M[region, celltype] <- M_enr[region, celltype]
      }
    }
  }
  M <- t(M)
  rownames(M) <- colnames(ref_marker)
  colnames(M) <- colnames(query_marker)
  M[is.infinite(M)] <- max(M[M != max(M)])
  print(M)
  write.table(x = M, file = paste0(output_dir,"/","celltype_enrichment.xls"), sep = "\t", col.names = NA)
}

# ==================Multimodal intersection analysis (MIA) subcommond  =================================================
sub_mia <- subparsers$add_parser("mia", help = "celltype mapping by Multimodal intersection analysis (MIA) ")

# sub_mia$add_argument("--queryexp", type = "character",
#             help = "[REQUIRED]the expression matrix in seurat format" )
# sub_mia$add_argument( "--queryassay", type = "character", default = "SCT",
#             help = "[REQUIRED]the default assay for queryexp.[default: %(default)s] ")

sub_mia$add_argument("--queryby", type = "character", default="clusters",
                     help = "[REQUIRED]the query column id in the reference expression matrix.[default: %(default)s] ")
sub_mia$add_argument("--querymarker", type = "character",
                     help = "[REQUIRED]the marker gene  for queryby(eg: cluster)" )
sub_mia$add_argument("--queryfilter", type = "character", default = "yes",
                     help = "[REQUIRED]filter the marker gene or not ,yes or no [default: %(default)s]" )
sub_mia$add_argument("--queryprefix", type = "character", default = "st_",
                     help = "[REQUIRED]the prefix set for output queryby .[default: %(default)s] ")

sub_mia$add_argument("--refexp", type = "character",
                     help = "[REQUIRED]the customized reference expression matrix in seurat format. ")
sub_mia$add_argument("--refassay", type = "character", default = "RNA",
                     help = "[OPTIONAL]the default assay for refexp.[default: %(default)s] ")
sub_mia$add_argument("--reformat", type = "character", default="rds",
                     help = "[REQUIRED]the customized reference expression matrix format.[default: %(default)s] ")
sub_mia$add_argument("--refby", type = "character", default="cell_type",
                     help = "[REQUIRED]the reference column id in the reference expression matrix.[default: %(default)s] " )
sub_mia$add_argument("--refmarker", type = "character",
                     help = "[REQUIRED] the marker gene  for refby (eg: cluster)")
sub_mia$add_argument("--refilter", type = "character", default = "yes",
                     help = "[OPTIONAL] filter the marker gene or not ,yes or no [default: %(default)s]" )
sub_mia$add_argument("--refprefix", type = "character", default = "ref_",
                     help = "[REQUIRED]the prefix set for output refby .[default: %(default)s] ")


sub_mia$add_argument("--min_pct1", "-t", type = "double", default = 0.5,
                     help ="the minimium ratio of cells expressing one specific gene in a cluster.[default: %(default)s] ")
sub_mia$add_argument("--max_pct2", "-T", type = "double", default = 0.5,
                     help ="the maximiium ratio of cells expressing one specific gene in all other clusters.[default: %(default)s] ")
sub_mia$add_argument("--pct_fold", "-c", type = "double", default = 2,
                     help ="the minimiu fold of pct1 for gene in a specific cluster against pct2 for all other cluster.[default: %(default)s] ")
sub_mia$add_argument("--avg_log2FC", "-k", type = "double", default = 1,
                     help = "The average log2FC of the gene UMI count in its cluster against the all other clusters.[default: %(default)s] ")
sub_mia$add_argument("--pvalue", "-p", type = "double", default = 0.05,
                     help = "the P-value of the gene differential expression.[default: %(default)s] ")
sub_mia$add_argument("--test", "-e", type = "character", default = "wilcox",
                     help = "test methods used to cell markers. Options are: wilcox,presto,venice, t, bimod,poisson, negbinom, MAST, DESeq2, DESeq2LRT,limma, edgeR.[default: %(default)s] ")
sub_mia$add_argument(
    "--sub_obj",
    default = NULL,
    help = "if the inputted rds had been subclustered or subsetted, it should be given by this parameter.[default: %(default)s]"
)
sub_mia$add_argument(
    "--recorrect_umi",
    default = TRUE,
    help = "whether to recorrect umi when run FindAllMarkers. If the object had run PrepSCTFindMarkers and subset, it should be 'FALSE'.[default: %(default)s]"
)
# 'wilcox' : Wilcoxon rank sum test (default).
# 'presto' : performace improved Wilcoxon rank sum test for big data.
# 'venice' : package from bioTuning with function VeniceAllMarkers imtating the FindAllMarkers in Seurat but much faster.")
# 't' : Student\'s t-test.
# 'bimod' : Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013).
# 'poisson' : Likelihood ratio test assuming an underlying poisson distribution. Use only for UMI-based datasets.
# 'negbinom' : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets.
# 'MAST' : GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015).
# 'DESeq2' : DE based on a model using the negative binomial distribution (Love et al, Genome Biology, 2014).
# ============== Subcmd: findallmarkers Find All markers for each cell group ================================================
args <- commandArgs(TRUE)
if ("mia" %in% args) {
  #parser parameter
  opt<- intial_setting()
  if (opt$sub_name == "mia") {
    futile.logger::flog.info("step1: refenence markers tiding:")#=========================================================
    ref_genes <-marker_gene_matrix(seurat_object_file=opt$refexp,
                                   format = opt$reformat,
                                   assay = opt$refassay,
                                   marker_file=opt$refmarker,
                                   filter = opt$refilter,
                                   groupby= opt$refby,
                                   test.use= opt$test,
                                   logfc.threshfold=opt$avg_log2FC,
                                   min_pct=opt$min_pct1,
                                   max_pct=opt$max_pct2,
                                   pct_fold=opt$pct_fold,
                                   prefix = opt$refprefix,
                                   pval=opt$pvalue,
                                   out_dir=output_dir)
    futile.logger::flog.info("step2: query markers tiding:")#=============================================================
    query_genes <-marker_gene_matrix(seurat_object_file=opt$input,
                                     format = opt$informat,
                                     assay = opt$assay,
                                     filter = opt$queryfilter,
                                     marker_file=opt$querymarker,
                                     groupby= opt$queryby,
                                     prefix = opt$queryprefix,
                                     test.use= opt$test,
                                     logfc.threshfold=opt$avg_log2FC,
                                     min_pct=opt$min_pct1,
                                     max_pct=opt$max_pct2,
                                     pct_fold=opt$pct_fold,
                                     pval=opt$pvalue,
                                     out_dir=output_dir,
                                     recorrect_umi = opt$recorrect_umi,
                                     predicate = opt$predicate,
                                     sub_obj = opt$sub_obj)
    futile.logger::flog.info("step3: running MIA analysis:")#=============================================================
    run_mia(ref_genes, query_genes, output_dir)
    if (!file.exists(file.path(output_dir, "MIA分析说明文档.docx"))) {
        file.copy("/public/dev_scRNA/oesinglecell3_test/document/MIA分析说明文档.docx",
                  file.path(output_dir, "MIA分析说明文档.docx")) }
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
  }
}



