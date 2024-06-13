sub_barcode_rank <- subparsers$add_parser("barcode_rank_gene_umi", help = "绘制barcode拐点图")
sub_barcode_rank$add_argument("-m", "--starsolo", type = "character", help = "starsolo  output dir ,eg:result/2.STAR/XHR-bm-ctr-74/XHR-bm-ctr-74_Solo.out/GeneFull/")
sub_barcode_rank$add_argument("-s", "--prefix", type = "character", help = "sample name id,eg:XHR-bm-ctr-74")
sub_barcode_rank$add_argument("-f", "--filtered", type = "character", help = "filtered or re-filtered")
args <- commandArgs(TRUE)
if ("barcode_rank_gene_umi" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "barcode_rank_gene_umi") {
    futile.logger::flog.info("step1:导入starsolo结果") #==================================================================
    counts_mtx <- Seurat::Read10X(glue::glue("{opt$starsolo}/raw/"))
    #sce<- DropletUtils::read10xCounts(glue::glue("{opt$starsolo}/raw/"), col.names = T)
    data_ob <- Seurat::CreateSeuratObject(counts = counts_mtx, assay = "RNA")
    #data_ob<-Seurat::as.Seurat(sce, counts = "counts",data = "logcounts")
    rm(counts)
    data_ob <- subset(data_ob, subset = nCount_RNA > 0)
    filtered_barcode <- data.table::fread(glue::glue("{opt$starsolo}/{opt$filtered}/barcodes.tsv.gz"), header = F)

    futile.logger::flog.info("step2:绘制Rank_Barcode-vs-Gene_countͼ") #===================================================
    barcode_gene <- Seurat::FetchData(data_ob, vars = "nFeature_RNA") %>% dplyr::arrange(nFeature_RNA)
    barcode_gene <- barcode_gene %>%
      as.data.frame %>%
      tibble::rownames_to_column("barcode") %>%
      dplyr::mutate(kept = ifelse(barcode %in% filtered_barcode$V1, "True", "False")) %>%
      dplyr::mutate(kept = factor(kept, levels = c("True", "False")))
    barcode_gene$rank <- rank(-barcode_gene$nFeature_RNA)
    plot1 <- ggplot(barcode_gene, aes(x = rank, y = nFeature_RNA)) +
      geom_point(shape = 1, aes(colour = kept)) +
      scale_x_log10(labels = scales::number) +
      scale_y_log10(labels = scales::number) +
      scale_colour_manual(values = c("#3558AE", "grey")) +
      ggtitle(glue::glue("{opt$prefix}_Gene")) +
      ylab("Total Counts") +
      xlab("Rank Barcode") +
      theme_bw(base_size = 20) +
      theme(panel.grid = element_blank(), legend.position = c(.12, .12), plot.title = element_text(hjust = 0.5)) +
      theme(legend.background = element_rect(color = "black", linetype = "solid", size = 0.1))

    OESingleCell::save_ggplots(plot1, filename = glue::glue("{output_dir}/Barcode_vs_Gene_Counts_for_{opt$prefix}"), width = 8, height = 8)

    futile.logger::flog.info("step3:绘制Rank_Barcode_vs_UMI_Countsͼ") #=============================================
    barcode_gene <- Seurat::FetchData(data_ob, vars = "nCount_RNA") %>% dplyr::arrange(nCount_RNA)
    barcode_gene <- barcode_gene %>%
      as.data.frame %>%
      tibble::rownames_to_column("barcode") %>%
      dplyr::mutate(kept = ifelse(barcode %in% filtered_barcode$V1, "True", "False")) %>%
      dplyr::mutate(kept = factor(kept, levels = c("True", "False")))
    barcode_gene$rank <- rank(-barcode_gene$nCount_RNA)
    plot2 <- ggplot(barcode_gene, aes(x = rank, y = nCount_RNA)) +
      geom_point(shape = 1, aes(colour = kept)) +
      scale_x_log10(labels = scales::number) +
      scale_y_log10(labels = scales::number) +
      scale_colour_manual(values = c("#3558AE", "grey")) +
      ggtitle(glue::glue("{opt$prefix}_UMI")) +
      ylab("Total Counts") +
      xlab("Rank Barcode") +
      theme_bw(base_size = 20) +
      theme(panel.grid = element_blank(), legend.position = c(.12, .12), plot.title = element_text(hjust = 0.5)) +
      theme(legend.background = element_rect(color = "black", linetype = "solid", size = 0.1))
    OESingleCell::save_ggplots(plot2, filename = glue::glue("{output_dir}/Barcode_vs_UMI_Counts_for_{opt$prefix}"), width = 8, height = 8)

    futile.logger::flog.info("step4：saving session information ") #=============================================
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  } }
