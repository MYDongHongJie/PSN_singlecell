docstring <- "example1:\\n\\n\\
sctool -i seurat.rds -f rds -o ./ gsva -m metadata.csv -g path/to/gmt -c 1000 -k Poisson -s 2 -S 10000 -x TRUE -d sampleid:A:B -p 0.05 -r 10 \\n\\n\\
example2: \\n\\n\\
sctool -i seurat.rds -f rds -o ./ gsva -m metadata.csv -g path/to/gmt -c 1000 -k Poisson -s 2 -S 10000 -x TRUE \\n\\n\\
example3: \\n\\n\\
sctool -i seurat.rds -f rds -o ./ gsva -e enrich_score_matrix.xls -d sampleid:A:B -p 0.05 -r 10"
sub_gsva <- subparsers$add_parser("gsva", description = docstring,
                                  #formatter_class ="lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                  formatter_class = "argparse.RawTextHelpFormatter", argument_default = "True",
                                  help = "Gene Set Variation Analysis.")
sub_gsva$add_argument("-m", "--metadata", type = "character", default = NULL,
                      help = "the sample metadata which must include sample id in this assay design.[default: %(default)s].")
sub_gsva$add_argument("-M", "--method", type = "character", default = "gsva",
                      help = paste0("Method to employ in the estimation of gene-set enrichment scores per sample. By default this is set ",
                                    "to gsva (Hanzelmann et al, 2013) and other options are ssgsea (Barbie et al, 2009), zscore (Lee et al, 2008)",
                                    "or plage (Tomfohr et al, 2005)."))
sub_gsva$add_argument("-g", "--gmt", type = "character", default = NULL,
                      help = "the gene sets in gmt format from the MsigDB/KEGG/GO database etc.")
sub_gsva$add_argument("-n", "--splitby", type = "character", default = NULL,
                      help = paste0("[OPTIONAL]Parallelized processing by spliting the cell size using one groupping factors in the metadata",
                                    "in case of too many cell for calculation"))
sub_gsva$add_argument("-c", "--chunkby", type = "integer", default = NULL,
                      help = paste0("[OPTIONAL]Parallelized processing by even chunked cell numbers for improving performace in case of",
                                    "too many cell for calculation"))
sub_gsva$add_argument("-k", "--kcdf", type = "character", default = "Gaussian",
                      help = paste0("Character string denoting the kernel to use during the non-parametric estimation of the cumulative ",
                                    "distribution function of expression levels across samples when method='gsva'. ",
                                    "The option can be Guassian, Possion and none. By default, kcdf='Gaussian' which is suitable when",
                                    " input expression values are continuous, such as microarray fluorescent units in logarithmic scale",
                                    "RNA-seq log- CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as ",
                                    "those derived from RNA-seq experiments, then this argument should be set to kcdf='Poisson'. This",
                                    "argument supersedes arguments rnaseq and kernel, which are deprecated and will be removed in the next",
                                    "release.[default: %(default)s]."))
sub_gsva$add_argument("-a", "--abs_rank", type = "logical", default = F,
                      help = paste0("Flag used only when mx_diff=TRUE. When abs_ranking=FALSE(default) a modified Kuiper statistic is used",
                                    "to calculate enrichment scores, taking the magnitude difference between the largest positive and negative",
                                    "random walk deviations. When abs.ranking=TRUE the original Kuiper statistic that sums the largest",
                                    "positive and negative random walk deviations, is used. In this latter case, gene sets with genes ",
                                    "enriched on either extreme (high or low) will be regarded as'highly' activated.[default: %(default)s]."))
sub_gsva$add_argument("-s", "--min_sz", type = "integer", default = 2,
                      help = "Minimum size of the resulting gene sets.[default: %(default)s].")
sub_gsva$add_argument("-S", "--max_sz", type = "integer", default = 100000,
                      help = "Maximum size of the resulting gene sets.[default: %(default)s].")
sub_gsva$add_argument("-x", "--mx_diff", type = "logical", default = T,
                      help = paste0("Offers two approaches to calculate the enrichment (ES) from the KS random walk statistic.mx_diff=FALSE:",
                                    "ES is calculated as the maximum distance of the random walk from 0. mx_diff=TRUE (default):",
                                    " ES is calculated as the magnitude difference between the largest positive and negative random walk ",
                                    "deviations.[default: %(default)s]."))
sub_gsva$add_argument("-t", "--tau", type = "double", default = 1,
                      help = paste0("Exponent defining the weight of the tail in the random walk performed by both the gsva ",
                                    "(Hanzelmann et al., 2013) and the ssgsea (Barbie et al., 2009) methods. By default, this tau=1 when ",
                                    "method='gsva' and tau=0.25 when method='ssgsea' just as specified by Barbie et al. (2009) where this",
                                    " parameter is called alpha.[default: %(default)s]."))
sub_gsva$add_argument("-q", "--WHICH_GROUP", type = "character", default = NULL,
                      help = paste0("[Optional]select the groupping column in metadata used to do subsettig of cell if necessary.",
                                    "[default: %(default)s]."))
sub_gsva$add_argument("-u", "--WHICH_CELLS", type = "character", default = NULL,
                      help = paste0("[Optional] the level id list in selected cell groupping used to celltyping.For some cell  clusters",
                                    " with high hertergensity this may be useful for sub-celltyping combined with the option -l/--LEVEL",
                                    " with single mode. If not specified with cell clusters's ID, all cells will be used.[default: %(default)s]."))
sub_gsva$add_argument("-d", "--contrast", type = "character", default = "NULL",
                      help = paste0("[Required]levels of a factor used to compare with for final differenetial results.The format is ",
                                    "Factor:interesting_level:reference_level."))
sub_gsva$add_argument("-F", "--fdr", type = "double", default = NULL,
                      help = "[Optional]the adjust P value threshold.")
sub_gsva$add_argument("-p", "--pval", type = "double", default = 0.05,
                      help = "[Optional]the P value threshold.[default: %(default)s].")
sub_gsva$add_argument("-r", "--topn", type = "integer", default = 10,
                      help = "the maximium number of ranked pathway terms on the top for each cluster.[default: %(default)s].")
sub_gsva$add_argument("-e", "--enrichmatrix", type = "character", default = NULL,
                      help = "the GSVA enrichment score matrix.")
args <- commandArgs(TRUE)
#=================================================================================
#function definition
#=================================================================================
`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`
## Define function to get commare pairs ================================================================================
get_compare_pairs <- function(contrasts, all_levels) {
  if (contrasts[2] == "all" & contrasts[3] != "all") {
    case_levels <- all_levels[-which(all_levels == contrasts[3])] #delete the reference level
    all_comparisions <- paste(contrasts[1], case_levels, contrasts[3], sep = ":")
  }else if (contrasts[2] != "all" & contrasts[3] == "all") {
    ref_levels <- all_levels[-which(all_levels == contrasts[2])] #delete the interested level
    all_comparisions <- paste(contrasts[1], contrasts[2], ref_levels, sep = ":")
  }else if (contrasts[2] == "all" & contrasts[3] == "all") {
    all_comparisions <- lapply(all_levels, function(x)
      paste(contrasts[1],
            x,
            paste0(all_levels[-which(all_levels == x)], collapse = ","),
            sep = ":"))
    all_comparisions <- unlist(all_comparisions)
  }else {
    all_comparisions <- opt$contrast
  }
  return(all_comparisions)
}

## find diff marker ====================================================================================================
run_findmarker <- function(seurat_ob, contrastsx) {
  Seurat::Idents(seurat_ob) <- contrastsx[1]
  pathway_diffexp <- Seurat::FindMarkers(seurat_ob,
                                         ident.1 = as.character(unlist(strsplit(contrastsx[2], ",", perl = T))),
                                         ident.2 = as.character(unlist(strsplit(contrastsx[3], ",", perl = T))),
                                         test.use = "limma",
                                         min.pct = 0,
                                         group.by = contrastsx[1],
                                         logfc.threshold = -Inf,
                                         assay = "GSVA",
                                         slot = "data")
  pathway_diffexp <- pathway_diffexp %>% tibble::rownames_to_column(var = "geneset")
  return(pathway_diffexp)
}

## get diffexp result ==================================================================================================
get_diffexp_result <- function(seurat_ob, contrasts, all_comparisions) {
  if (contrasts[2] == "all" & contrasts[3] == "all") {
    gsva_results <- c()
    for (contrastx in all_comparisions) {
      contrastsx <- unlist(strsplit(contrastx, ':', perl = T))
      pathway_diffexp <- run_findmarker(seurat_ob, contrastsx)
      pathway_diffexp$cluster <- as.character(contrastsx[2])
      gsva_results <- rbind(gsva_results, pathway_diffexp)
    }
    gsva_results <- gsva_results %>%
      dplyr::select(geneset, logFC, AveExpr, t, p_val, adj.P.Val, B, cluster)
    colnames(gsva_results) <- c("geneset", "logFC", "avgExp", "t", "pval", "FDR", "B", "Interested_group")
    write.table(gsva_results,
                file.path(output_dir, "GSVA_results.xls"),
                quote = F,
                sep = "\t",
                col.names = T,
                row.names = F)
    ####
    if (contrasts[1] == "clusters") {
      gsva_results$Interested_group <- factor(gsva_results$Interested_group,
                                              levels = sort(unique(as.numeric(seurat_ob@meta.data$clusters))))
    }else {
      gsva_results$Interested_group <- as.factor(gsva_results$Interested_group)
    }
    if (is.null(FDR_cutoff)) {
      plot_term <- gsva_results %>%
        subset(pval < pval_cutoff & t > 0) %>%
        dplyr::group_by(Interested_group) %>%
        dplyr::top_n(opt$topn, t) %>%
        dplyr::arrange(Interested_group)
    }else {
      plot_term <- gsva_results %>%
        subset(FDR < FDR_cutoff & t > 0) %>%
        dplyr::group_by(Interested_group) %>%
        dplyr::top_n(opt$topn, t) %>%
        dplyr::arrange(Interested_group)
    }
    write.table(plot_term, file.path(output_dir, paste0("GSVA_top", opt$topn, "_results.xls")),
                quote = F, sep = "\t", col.names = T, row.names = F)
    ###=========================================================================================================================
    plot_data <- tidyr::spread(as.data.frame(gsva_results[, c("geneset", "Interested_group", "t")]),
                               Interested_group,
                               t) %>%
      tibble::column_to_rownames(var = "geneset")
    plot_data <- plot_data[as.vector(unique(plot_term$geneset)),]
    plot_data<- plot_data %>% t %>%  scale(center = TRUE, scale = TRUE)  %>%t
    pdf(NULL)
    p<-ComplexHeatmap::Heatmap(as.matrix(plot_data), # 输入数据为矩阵
              name = 'heatmap for GSVA', # 热图图例名称
              col = grDevices::colorRampPalette(c("blue", "white", "red"))(n = 299), #配色方案
              border = FALSE, # 显示边框
              cluster_rows = F, #按照行进行聚类
              cluster_columns = F, #按照列进行聚类
              #column_km = 2, # 划分列聚类
              #row_km = 2, # 划分行聚类
              show_column_names = T, #显示列名
              show_row_names = T, #显示行名与否
              column_names_side = "bottom", , #列名位置 top bottom
              row_names_rot = 0, ##行名旋转角度
              column_names_rot = 90, ##列名旋转角度
              row_names_side = "right", #行名位置 left right
              column_names_gp = grid::gpar(fontsize = 8), ##列字体大小
              row_names_gp = grid::gpar(fontsize = 8), ##行字体大小
              #alpha=0.5, #设置图片颜色透明度
              # width = ncol(gene_data)*unit(5, "mm"), # 格子的宽度
              # height = nrow(gene_data)*unit(2, "mm"),# 格子的高度
              top_annotation = NULL, # 添加左侧注释信息
              right_annotation = NULL, # 添加右侧注释信息
              left_annotation = NULL, # 添加左侧注释信息
              bottom_annotation = NULL, # 添加下方注释信息
              heatmap_legend_param = list(title = "",
                                          title_position = "topcenter", # 标题相对图例的位置 topcenter, topleft, leftcenter, lefttop.
                                          #at=c(0,10), #图例范围
                                          labels_gp = grid::gpar(fontsize = 10), ##图例字体
                                          legend_direction = "vertical", ##图例放置
                                          legend_height = unit(5, "cm") #图例长度
              ))
	  OESingleCell::save_ggplots(file.path(output_dir, paste0("top", opt$topn, "_gsva_term.pdf")),
                                 plot= ggplotify::as.ggplot(p),
                                 width=dim(plot_data)[2]*0.1+max(string::str_length(colnames(plot_data)))*1,
                                 height=dim(plot_data)[1]*0.1+max(string::str_length(rownames(plot_data)))*0.05)

    write.table(as.data.frame(plot_data) %>% tibble::rownames_to_column(var = "geneset"),
                file.path(output_dir, paste0("top", opt$topn, "_term_t_value.xls")),
                quote = F,
                col.names = T,
                row.names = F,
                sep = "\t")
  }else {
    for (contrastx in all_comparisions) {
      contrastsx <- unlist(strsplit(contrastx, ':', perl = T))
      pathway_diffexp <- run_findmarker(seurat_ob, contrastsx)
      pathway_diffexp <- pathway_diffexp %>% dplyr::select(geneset, logFC, AveExpr, t, p_val, adj.P.Val, B)
      colnames(pathway_diffexp) <- c("geneset", "logFC", "avgExp", "t", "pval", "FDR", "B")
      write.table(pathway_diffexp,
                  file.path(output_dir, paste0("diffexp_genesets_GSVA_score4", contrastsx[1], "_", contrastsx[2], "-vs-", contrastsx[3], ".xls")),
                  quote = F, sep = '\t', row.names = F)
      ####
      pathway_diffexp$just <- ifelse(pathway_diffexp$t < 0, 0, 1)
      pathway_diffexp$is.just <- pathway_diffexp$just == 1
      if (is.null(FDR_cutoff)) {
        pathway_diffexp$is.sig <- pathway_diffexp$pval < pval_cutoff
      }else {
        pathway_diffexp$is.sig <- pathway_diffexp$FDR < FDR_cutoff
      }
      pathway2vis <- pathway_diffexp %>%
        dplyr::filter(pval < pval_cutoff) %>%
        dplyr::group_by(is.just) %>%
        dplyr::arrange(abs(t)) %>%
        dplyr::top_n(opt$topn, abs(t))
      pp <- ggplot(pathway2vis, aes(reorder(geneset, t), t)) +
        geom_col(aes(fill = is.just)) +
        scale_fill_manual(values = c("#6CC570", "#2A5078")) +
        coord_flip() +
        labs(x = "Pathway", y = "t value of GSVA Score") +
        theme_minimal() +
        geom_text(aes(x = geneset, y = 0, label = geneset), hjust = pathway2vis$just, size = 3.5) +
        theme(axis.text.y = element_blank()) +
        theme(panel.grid = element_blank())
      pp <- pp + labs(fill = paste0("t value > 0"))
      ggsave(file.path(output_dir, paste0("diffexp_genesets_GSVA_score4", contrastsx[1], "_", contrastsx[2], "-vs-", contrastsx[3], "_barplot.pdf")),
             plot = pp, width = 12, height = 8, bg = "white")
      ggsave(file.path(output_dir, paste0("diffexp_genesets_GSVA_score4", contrastsx[1], "_", contrastsx[2], "-vs-", contrastsx[3], "_barplot.png")),
             plot = pp, width = 12, height = 8, bg = "white")
    }
  }
}

#=================================================================================
#input file and parameter
#=================================================================================
if ("gsva" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "gsva") {
    if (!is.null(opt$metadata)) { assay_metadata <- read.csv(opt$metadata, sep = ",", header = T, row.names = 1) }
    if (is.null(opt$mx_dff)) {
      if (opt$method == "gsva") { mx.diff <- T } else { mx.diff <- F }
    }else { mx.diff <- as.logical(toupper(opt$mx_diff)) }
    if (is.null(opt$tau)) {
      if (opt$method == "gsva") { tau <- 1 } else { tau <- 0.25 }
    }else { tau <- as.numeric(opt$tau) }
    if (!is.null(opt$gmt)) { gene_sets <- GSEABase::getGmt(opt$gmt) }
    if (!is.null(opt$enrichmatrix)) {
      gsva_scores <- read.table(opt$enrichmatrix,
                                sep = "\t",
                                header = T,
                                row.names = 1,
                                check.names = F,
                                quote = '"') }
    if (is.null(opt$fdr)) {
      FDR_cutoff <- NULL
      if (is.null(opt$pval)) {
        stop("Warning:NO FDR or P value is specified!")
      }else { pval_cutoff <- opt$pval }
    }else { FDR_cutoff <- opt$fdr }
    ##==============================================================================================================
    if ((opt$informat == "rds" || opt$format == "h5seurat")) {
      seurat_ob <- OESingleCell::ReadX(input = opt$input,
                                       informat = opt$informat,
                                       assays = opt$assay,
                                       data.use = dataslots,
                                       verbose = F)
      ###
      if (!is.null(opt$metadata)) {
        sampleidx <- gsub("(_|-)[ATGC]{16,}.*", "", cellnames, perl = T)
        additional_cell_meta <- vector()
        for (colidx in colnames(assay_metadata)) {
          additional_cell_meta <- cbind(additional_cell_meta, as.vector(assay_metadata[sampleidx, colidx])) }
        colnames(additional_cell_meta) <- colnames(assay_metadata)
        rownames(additional_cell_meta) <- cellnames
        additional_cell_meta <- as.data.frame(additional_cell_meta)
        seurat_ob <- SeuratObject::AddMetaData(seurat_ob, additional_cell_meta)
        assay_metadata <- seurat_ob@meta.data %>%
          tibble::rownames_to_column(var = "cellbarcode") %>%
          dplyr::select(cellbarcode, dplyr::everything())
        write.table(assay_metadata,
                    file = file.path(output_dir, "metadata4each_cell.xls", quote = F, col.names = T, row.names = F, sep = ","))
      }
      ###
      if (is.null(seurat_ob@meta.data$clusters)) {
        seurat_ob <- SeuratObject::StashIdent(seurat_ob, save.name = "clusters")
      }else {
        Seurat::Idents(object = seurat_ob) <- seurat_ob[["clusters"]] }
      if (!is.null(opt$WHICH_CELLS)) {
        cluster_list <- unlist(strsplit(opt$WHICH_CELLS, ",", perl = T))
        seurat_ob <- SubsetData(seurat_ob, subset.name = opt$WHICH_GROUP, accept.value = cluster_list)
      }
      cellnames <- Seurat::Cells(seurat_ob)
      cell_num <- length(cellnames)
      ####
    }
    #===============================================================================================================
    ### get gsva_score for diff formart input file
    ######
    if (is.null(opt$enrichmatrix)) {
      if (opt$informat == "raw") {
        if (!"assay_metadata" %in% ls()) { stop("Please provide metadata file for raw matrix using -m") }
        exprs_matrix <- read.table(opt$input, sep = "\t", header = T, row.names = 1)
        exprs_matrix <- as.matrix(exprs_matrix)
        eset <- Biobase::ExpressionSet(assayData = expers_matrix, phenoData = Biobase:AnnotatedDataFrame(assay_metadata))
        gsva_scores <- GSVA::gsva(exprs(eset),
                                  gene_sets,
                                  method = opt$method,
                                  kcdf = opt$kcdf,
                                  min.sz = opt$min_sz,
                                  max.sz = opt$max_sz,
                                  mx.diff = mx.diff)
      }

      ######
      if (opt$informat == "rds" || opt$format == "h5seurat") {
        if (!is.null(opt$chunkby)) {
          if (cell_num < opt$chunkby) { chunkby <- cell_num }else { chunkby <- opt$chunkby }
          intervals <- seq(1, cell_num, by = chunkby)
          doParallel::registerDoParallel(cores = length(intervals))
          gsva_enrichment <- foreach::foreach(x = intervals) %dopar% GSVA::gsva(
                     Seurat::GetAssayData(seurat_ob,
                                  assay = opt$assay,
                                  slot = "counts")[, cellnames[seq(x, min(x + chunkby - 1, cell_num))]] %>% as.matrix,
                      gene_sets,
                      method = opt$method,
                      kcdf = opt$kcdf,
                      min.sz = opt$min_sz,
                      max.sz = opt$max_sz,
                      mx.diff = mx.diff)
          doParallel::stopImplicitCluster()
          gsva_scores <- as.data.frame(gsva_enrichment[[1]], check.names = F)
          for (indx in 2:length(gsva_enrichment)) {
            gsva_scores <- transform(merge(gsva_scores, as.data.frame(gsva_enrichment[[indx]]), by = 0, all = F),
                                     row.names = Row.names, Row.names = NULL) }
          colnames(gsva_scores) <- gsub("^X", "", colnames(gsva_scores))
        }else if (!is.null(opt$splitby)) {
          seurat_ob_list <- Seurat::SplitObject(seurat_ob, split.by = opt$splitby)
          doParallel::registerDoParallel(cores = length(seurat_ob_list))
          gsva_enrichment <- foreach::foreach(seurat_subset = seurat_ob_list) %dopar% GSVA::gsva(
            as.matrix(Seurat::GetAssayData(seurat_ob, assay = opt$assay, slot = "counts")),
            gene_sets,
            method = opt$method,
            kcdf = opt$kcdf,
            min.sz = opt$min_sz,
            max.sz = opt$max_sz,
            mx.diff = mx.diff)
          doParallel::stopImplicitCluster()
          gsva_scores <- as.data.frame(gsva_enrichment[[1]], check.names = F)
          for (indx in 2:length(gsva_enrichment)) {
            gsva_scores <- transform(merge(gsva_scores, as.data.frame(gsva_enrichment[[indx]]),
                                           by = 0, all = F),
                                     row.names = Row.names,
                                     Row.names = NULL)
          }
          colnames(gsva_scores) <- gsub("^X", "", colnames(gsva_scores))
        }else {
          gsva_scores <- GSVA::gsva(as.matrix(as.matrix(Seurat::GetAssayData(seurat_ob, assay = opt$assay, slot = "counts"))),
                                    gene_sets,
                                    method = opt$method,
                                    kcdf = opt$kcdf,
                                    min.sz = opt$min_sz,
                                    max.sz = opt$max_sz,
                                    mx.diff = mx.diff)
        }
      }

      ####### save  gsva_scores results
      gsva_scores <- as.data.frame(gsva_scores)
      gsva_scores$geneset <- rownames(gsva_scores)
      gsva_scores <- gsva_scores %>% dplyr::select(geneset, dplyr::everything())
      gsva_scores_path <- file.path(output_dir, "GSVA_enrichment_results.xls")
      write.table(gsva_scores, file = gsva_scores_path, col.names = T, row.names = F, sep = "\t")
    }
    #=====================================================================================================================
    #### diff analy
    if (!is.null(opt$contrast)) {
      seurat_ob[['GSVA']] <- Seurat::CreateAssayObject(counts = gsva_scores)
      seurat_ob <- Seurat::ScaleData(seurat_ob, assay = "GSVA")
      assay_metadata <- seurat_ob@meta.data
      contrasts <- unlist(strsplit(opt$contrast, ":", perl = T))
      all_levels <- as.vector(unique(assay_metadata[, contrasts[1]]))
      all_comparisions <- get_compare_pairs(contrasts, all_levels)
      get_diffexp_result(seurat_ob, contrasts, all_comparisions)
    }
    ## save session information
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}