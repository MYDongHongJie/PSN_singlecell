sub_qcmetrics = subparsers$add_parser("qcmetrics", help = "integration of multiple visualizion for st-piplint.")
sub_qcmetrics$add_argument("-l","--level", type = "character", default = "cell",
        help = "the analysis object level: cell or spot.[default: %(default)s]")
sub_qcmetrics$add_argument("--metadata", "-m",
  type = "character", default = NULL,
  help = "[Optional]the additional sample metadata in case of assay design changes. [default: %(default)s]")
sub_qcmetrics$add_argument("--colorby", "-c",
  type = "character", default = "sampleid",
  help = "[Otional]visualize cells' metadata by coloring cells in different color according to cell grouping levels. [default: %(default)s]")
sub_qcmetrics$add_argument("--prefix", "-d",
  type = "character", default = "beforeQC",
  help = "[Optional]the prefix of the output file. [default: %(default)s]")
sub_qcmetrics$add_argument("--metrics", "-v",
  type = "character", default = "nGene,nUMI,percent.mito,percent.ribo,batch,cellcycle",
  help = "The QC metrics list of source of hetergenesity to visualize. [default: %(default)s]")
sub_qcmetrics$add_argument("--perplexity", "-P",
  type = "integer", default = 30,
  help = "[Optional]The value of the perplexity used for tSNE.[default: %(default)s]")
sub_qcmetrics$add_argument("--components", "-n",
  type = "integer", default = 20,
  help = "[Optional]the appropriate number of statistically significant components to use for clustering,which derived from the JackStraw result. [default: %(default)s]")
sub_qcmetrics$add_argument("--ptsize",
  type = "double", default = 1,
  help = "[OPTIONAL]setting the dot point size on the violin plot. Set this to 0 when you prefer no points. [default: %(default)s]")
sub_qcmetrics$add_argument("--alpha",
  type = "double", default = 0.6,
  help = "[OPTIONAL]setting the dot opacity on the violin plot. [default: %(default)s]")
sub_qcmetrics$add_argument("--pointsize", "-s",
  type = "double", default = 0.5,
  help = "[OPTIONAL]the point size in the plot. [default: %(default)s]")
sub_qcmetrics$add_argument("--colshema",
  type = "character", default = "red,grey",
  help = "the color schema for groupping coloring. [default: %(default)s]")
sub_qcmetrics$add_argument( "--mito", "-M" , type = "character", default = NULL, help = "mitochondrial gene list.")
args <- commandArgs(TRUE)
if ( "qcmetrics" %in% args ){
  opt<-intial_setting()
  prefix <- opt$prefix
  groupfactor <- opt$colorby
  QC_metrics <- unlist(strsplit(opt$metrics, ",", perl = T))
  pt.size <- as.numeric(opt$ptsize)
  assay2use <- Hmisc::capitalize(assays)
  alpha2use <- opt$alpha
# read in the 10X data

  seurat_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslot, # "data","scale.data"
                                 reductions = opt$reduct,
                                 graphs = FALSE,
                                 images = opt$image, # no graph object needed here
                                 verbose = FALSE)

  # set the default assay as the specified
  Seurat::DefaultAssay(seurat_ob) <- assay2use

# update the metedata in the singlecell_ob@meta.data with new additional sample metadata
  if (!is.null(opt$metadata)) {
    additional_metadata <- read.csv(opt$metadata, sep = ",", header = T)
    rownames(additional_metadata) <- additional_metadata$sampleid
    cellnames <- Cells(seurat_ob)
    sampleidx <- gsub("(-|_)[ATGC]{16,}.*", "", cellnames, perl = T) # the index order is the same as the row index of the assay metadata
    # integrate the additional metadata from the assay design
    additional_cell_meta <- vector()
    for (colidx in colnames(additional_metadata)) {
      additional_cell_meta <- cbind(additional_cell_meta, as.vector(additional_metadata[sampleidx, colidx]))
    }
    colnames(additional_cell_meta) <- colnames(additional_metadata)
    rownames(additional_cell_meta) <- cellnames
    additional_cell_meta <- as.data.frame(additional_cell_meta)
    seurat_ob <- Seurat::AddMetaData(seurat_ob, metadata = additional_cell_meta)
  }

  # setting the maxumium chunck mermory usage much bigger in case of big data
  options(future.globals.maxSize = Inf)
  future::plan("multiprocess", workers = opt$ncores) # parallization start from here
  if (!length(Seurat::VariableFeatures(seurat_ob))) { # if calculated, ignore it.
    seurat_ob <- Seurat::NormalizeData(
      object = seurat_ob,
      normalization.method = "LogNormalize", scale.factor = 10000
    )
    seurat_ob <- Seurat::FindVariableFeatures(
      object = seurat_ob,
      mean.function = "FastExpMean",
      dispersion.function = "FastLogVMR", do.plot = F
    )
    seurat_ob <- Seurat::ScaleData(object = seurat_ob, verbose = F) # takes some time
  }

  # check the components to use
  suppressWarnings({
    if (is.null(opt$components)) {
      if (is.null(SeuratObject::Misc(seurat_ob, "components_num"))) {
        print("NO previous components calculation is AVAILABLE, 20 will be used as default.")
        SeuratObject::Misc(seurat_ob, "components_num") <- 20
      } else {
        components_num <- SeuratObject::Misc(seurat_ob, "components_num")
      }
    } else {
      SeuratObject::Misc(seurat_ob, "components_num") <- opt$components
      components_num <- opt$components
    }
  })

  # check the value of perplexity
  suppressWarnings({
    if (is.null(opt$perplexity)) {
      if (is.null(SeuratObject::Misc(seurat_ob, "perplexity"))) {
        print("NO previous perplexity is AVAILABLE, 30 will be used as default.")
        SeuratObject::Misc(seurat_ob, "perplexity") <- 30
      } else {
        perplexity <- SeuratObject::Misc(seurat_ob, "perplexity")
      }
    } else {
      perplexity <- opt$perplexity
      SeuratObject::Misc(seurat_ob, "perplexity") <- opt$perplexity
    }
  })

  if (!"pca" %in% names(Seurat::Key(seurat_ob))) { # if calcuated, ignore it
    seurat_ob <- Seurat::RunPCA(object = seurat_ob, features = Seurat::VariableFeatures(seurat_ob), verbose = F)
  }
  if (!"umap" %in% names(Seurat::Key(seurat_ob))) { # if calcuated, ignore it
    seurat_ob <- Seurat::RunUMAP(
      object = seurat_ob, dims = 1:components_num,
      reduction = "pca",
      # assay = assay2use,
      perplexity = perplexity, force.recalc = T, check_duplicates = F,
      do.fast = T, dim.embed = 2
    )
  }
  # The number of genes and UMIs are automatically calculated for every object by Seurat.
  # For non-UMI data, nUMI represents the sum of the non-normalized values within a cell.
  # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
  # ,since this represents non-transformed and non-log-normalized counts.
  # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.

  if ("nGene" %in% QC_metrics) {
    print(unique(seurat_ob$sampleid))
    nfeatures <- paste0("nFeature_", assay2use)
    vln4nGene <- Seurat::VlnPlot(
      object = seurat_ob, features = nfeatures, group.by = groupfactor, ncol = 1,
      pt.size = pt.size
    ) + ggplot2::ggtitle("nGene") + ggplot2::labs(xlab = "") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + Seurat::NoLegend()
    if (assay2use == "Spatial") {
    gsp <- OESingleCell::SpatialPlot(seurat_ob,
                          features = nfeatures,
                          min.cutoff = "q0",
                          max.cutoff = "q90",
                          alpha = c(0.1, 1),
                          pt.size.factor = opt$pointsize,
                          pt.alpha=FALSE,
                          stroke=0.2,combine=FALSE
                          )

      SpatialColors <- colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
      nfeatures_range<-seurat_ob@meta.data %>% dplyr::select(c(nfeatures)) %>% range()
      gsp <- lapply(gsp, function(x) x
              + ggplot2::scale_fill_gradientn(colors= SpatialColors(n = 100),
                                     limits = c(nfeatures_range[1],nfeatures_range[2])))
      gsb <- do.call(ggpubr::ggarrange,
                     c(gsp, list(ncol = length(gsp),
                                 nrow = 1,
                                 common.legend = TRUE,
                                 legend = "right",
                                 align = "none")))
      wit_rot <- length(gsp)
      plot <- cowplot::plot_grid(vln4nGene, gsb,rel_widths= c(1, wit_rot))
      width <-  3 * length(gsp) + length(vln4nGene)
      height <- 4
    }
    else {
      plot <- vln4nUMI
      width<-8
      heigh<-6
    }
    OESingleCell::save_ggplots(file.path(output_dir,
                       paste0(prefix, "_total_genes4each_", opt$level, "_on_violin_plot", collapse = "")),
             plot = plot,
             width = width,
             height = height,
             dpi=300,bg='white')


  }
  if ("nUMI" %in% QC_metrics) {
    nCount <- paste0("nCount_", assay2use)
    vln4nUMI <- Seurat::VlnPlot(
      object = seurat_ob, features = nCount, group.by = groupfactor, ncol = 1,
      pt = pt.size
    ) + ggplot2::ggtitle("nUMI") + ggplot2::labs(xlab = "") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + Seurat::NoLegend()
    if (assay2use == "Spatial") {
      gsp <- OESingleCell::SpatialPlot(seurat_ob,
                features =  nCount,
                min.cutoff = "q0",
                max.cutoff = "q90",
                alpha = c(0.1, 1),
                pt.size.factor = opt$pointsize,
                pt.alpha=TRUE,
                stroke=0.2,combine=FALSE
                )

      SpatialColors <- colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
      nCount_range<-seurat_ob@meta.data %>% dplyr::select(c(nCount)) %>% range()
      gsp <- lapply(gsp, function(x) x
              + ggplot2::scale_fill_gradientn(colors= SpatialColors(n = 100),
                                     limits = c(nCount_range[1],nCount_range[2])))
      gsb <- do.call(ggpubr::ggarrange,
                     c(gsp, list(ncol = length(gsp),
                                 nrow = 1,
                                 common.legend = TRUE,
                                 legend = "right",
                                 align = "none")))
      plot <- cowplot::plot_grid(vln4nUMI, gsb,rel_widths= c(1, length(gsp)))
      width <-  3 * length(gsp) + length(vln4nGene)
      height <- 4
    }
    else {
      plot <- vln4nUMI
      width<-8
      heigh<-6
    }
    OESingleCell::save_ggplots(file.path(output_dir,
                     paste0(prefix, "_total_UMIs4each_", opt$level, "_on_violin_plot")),
           plot = plot,
           width = width,
           height = height,
           dpi=300,bg='white')
  }
      if ("batch" %in% QC_metrics) {
      # to visulaize the metadata of each cell on the PCA plot

      umap4batch <- Seurat::DimPlot(object = seurat_ob,
                            dims = c(1, 2),
                            reduction = "umap",
                            label = F,
                            pt.size = opt$pointsize,
                            group.by = "batchid") +
                    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      OESingleCell::save_ggplots(file.path(output_dir,
                     paste0(prefix, "_visualize_batch_effect_on_umap_plot", collapse = "")),
                   plot = umap4batch,
                   dpi=300,bg='white')
    }

    # We calculate the percentage of mitochondrial genes here and store it in percent.mito using AddMetaData.
    if ("percent.mito" %in% QC_metrics) {
      if ( !is.null(opt$mito) ){
                                mito.genes <- utils::read.delim(opt$mito)
                                raw.counts = Seurat::GetAssayData(seurat_ob, slot = "counts")
                                percent.mito <- Matrix::colSums(raw.counts[mito.genes[,1], ])/Matrix::colSums(raw.counts)
                        }else{
                                mito.genes <- grep(pattern = "^(MT|mt)(-|_)", x = rownames(seurat_ob), value = T,perl=T)
                                if(length(mito.genes)==0){
                                  stop("Can not fetch mito genes and please input mito gene list")
                                }else{
                                  raw.counts = Seurat::GetAssayData(seurat_ob, slot = "counts")
                                  percent.mito <- Matrix::colSums(raw.counts[mito.genes, ])/Matrix::colSums(raw.counts)

                                }
                        }
      seurat_ob <- Seurat::AddMetaData(object = seurat_ob, metadata = percent.mito, col.name = "percent.mito")
      vln4mito <- Seurat::VlnPlot(
        object = seurat_ob, features = "percent.mito",
        group.by = groupfactor, ncol = 1, pt.size = pt.size
      ) + ggplot2::labs(xlab = "") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + Seurat::NoLegend()
      if (assay2use == "Spatial") {

        gsp <- OESingleCell::SpatialPlot(seurat_ob,
                                  features = "percent.mito",
                                  min.cutoff = "q0",
                                  max.cutoff = "q90",
                                  alpha = c(0.1, 1),
                                  pt.size.factor = opt$pointsize,
                                  pt.alpha=TRUE,
                                  stroke=0.2,combine=FALSE)

        SpatialColors <- colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
        percent.mito_range<-seurat_ob@meta.data %>% dplyr::select(c( "percent.mito")) %>% range()
        gsp <- lapply(gsp, function(x) x
                + ggplot2::scale_fill_gradientn(colors= SpatialColors(n = 100),
                                       limits = c(percent.mito_range[1],percent.mito_range[2])))
        gsb <- do.call(ggpubr::ggarrange,
                       c(gsp,
                         list(ncol = length(gsp),
                              nrow = 1,
                              common.legend = TRUE,
                              legend = "right",
                              align = "hv")))
        plot <- cowplot::plot_grid(vln4mito, gsb,rel_widths= c(1, length(gsp)))
        width <-  3 * length(gsp) + length(vln4nGene)
        height <- 4
      }
      else {
        plot <- vln4mito
        width<-8
        heigh<-6
      }
      OESingleCell::save_ggplots(file.path(output_dir,
                 paste0(prefix, "_mitochondroin_transcript_ratio_in_each_", opt$level, "_violin_plot")),
                   plot = plot,
                   width = width,
                   height = height,
                   dpi=300,bg='white')


    }
    # We calculate the percentage of ribosome genes here and store it in percent.ribo using AddMetaData.
    if ("percent.ribo" %in% QC_metrics) {
        # Calculate percent ribosomal genes.
        ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", rownames(seurat_ob), value = TRUE, ignore.case = T)
        percent.ribo <- Seurat::PercentageFeatureSet(seurat_ob, features = ribo.genes)
        seurat_ob <- Seurat::AddMetaData(object = seurat_ob, metadata = percent.ribo, col.name = "percent.ribo")
        vln4ribo <- Seurat::VlnPlot(object = seurat_ob,
                            features = "percent.ribo",
                            group.by = groupfactor,
                            ncol = 1,
                            pt.size = pt.size) + labs(xlab = "") + theme(plot.title = element_text(hjust = 0.5))
        gg_ribo_pca <- Seurat::FeaturePlot(seurat_ob,
          features = "percent.ribo", reduction = "pca",
          cols = c("grey", "red"), pt.size = pt.size
        )
        gg_ribo_umap <- Seurat::FeaturePlot(seurat_ob,
          features = "percent.ribo", reduction = "umap",
          cols = c("grey", "red"), pt.size = pt.size
        )
        if (assay2use == "Spatial") {

          gsp <- OESingleCell::SpatialPlot(seurat_ob,
                                          features =  "percent.ribo",
                                          min.cutoff = "q0",
                                          max.cutoff = "q90",
                                          alpha = c(0.1, 1),
                                          pt.size.factor = opt$pointsize,
                                          pt.alpha=TRUE,
                                          stroke=0.2,combine=FALSE)
          gsb <- do.call(ggpubr::ggarrange,
                         c(gsp,
                           list(ncol = length(gsp),
                                nrow = 1,
                                common.legend = TRUE,
                                legend = "right",
                                align = "hv")))
          wit_rot <- length(gsp)
          plot1 <- cowplot::plot_grid(vln4ribo, gsb,rel_widths= c(1, wit_rot ))
          plot2 <- cowplot::plot_grid(gg_ribo_pca, gsb,rel_widths= c(1,wit_rot ))
          plot3 <- cowplot::plot_grid(gg_ribo_umap, gsb,rel_widths= c(1,wit_rot ))
          width <-  3 * length(gsp) + length(vln4ribo)
          height <- 4
        }
        else {
          plot1 <- vln4ribo
          plot2 <- gg_ribo_pca
          plot3 <- gg_ribo_umap
          width<-8
          heigh<-6
        }
        OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_ribosome_transcript_ratio_in_each_", opt$level, "_violin_plot")),
                                  plot=plot1,
                                  width=width,
                                  height=height,
                                  dpi=300,bg='white')
        OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_ribosome_transcript_ratio_in_each_", opt$level, "_PCA_plot")),
                                  plot=plot1,
                                  width=width,
                                  height=height,
                                  dpi=300,bg='white')
        OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_ribosome_transcript_ratio_in_each_", opt$level, "_umap_plot")),
                                  plot=plot1,
                                  width=width,
                                  height=height,
                                  dpi=300,bg='white')

      }
  if ("top50" %in% QC_metrics) {
      top50_gene <- apply(seurat_ob@raw.data, 2, function(x) sum(x[order(x, decreasing = T)][1:50]) / sum(x))
      seurat_ob <- Seurat::AddMetaData(object = seurat_ob, metadata = top50_gene, col.name = "top50")
      vln4top50 <- Seurat::VlnPlot(
                  object = seurat_ob, features = "top50",
                  group.by = groupfactor, ncol = 1, pt.size = pt.size) +
        ggplot2::labs(xlab = "") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      gg_top50_pca <- Seurat::FeaturePlot(seurat_ob, features = "top50", reduction = "pca", cols = c("grey", "red"), pt.size = pt.size)
      gg_top50_umap <- Seurat::FeaturePlot(seurat_ob, features = "top50", reduction = "umap", cols = c("grey", "red"), pt.size = pt.size)
      if (assay2use == "Spatial") {
        gsp <- OESingleCell::SpatialPlot(seurat_ob,
                                  features =  "top50",
                                  min.cutoff = "q0",
                                  max.cutoff = "q90",
                                  alpha = c(0.1, 1),
                                  pt.size.factor = opt$pointsize,
                                  pt.alpha=TRUE,
                                  stroke=0.2,combine=FALSE)


        gsb <- do.call(ggpubr::ggarrange,
                       c(gsp,
                         list(nrow = length(gsp),
                              ncol = 1,
                              common.legend = TRUE,
                              legend = "right",
                              align = "hv")))
        wit_rot <- length(gsp)
        plot1 <- cowplot::plot_grid(vln4top50, gsb,rel_widths= c(1, wit_rot))
        plot2 <- cowplot::plot_grid(gg_top50_pca, gsb,rel_widths= c(1, wit_rot))
        plot3 <- cowplot::plot_grid(gg_top50_umap, gsb,rel_widths= c(1, wit_rot))
        width <-  3 * length(gsp) + 5 * (length(gsp)/2)
        height <- 6
      }
      else {
        plot1 <- vln4top50
        plot2 <- gg_top50_pca
        plot3 <- gg_top50_umap
        width<-8
        heigh<-6
      }
      OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_transcript4top50_gene_ratio_in_each_", opt$level, "_violin_plot")),
                                  plot=plot1,
                                  width=width,
                                  height=height,
                                  dpi=300,bg='white')
      OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_transcript4top50_gene_ratio_in_each_", opt$level, "_PCA_plot")),
                                  plot=plot1,
                                  width=width,
                                  height=height,
                                  dpi=300,bg='white')
      OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_transcript4top50_gene_ratio_in_each_", opt$level, "_umap_plot")),
                                  plot=plot1,
                                  width=width,
                                  height=height,
                                  dpi=300,bg='white')

  }
      if ("UMI.per.gene" %in% QC_metrics) {
      tmp <- seurat_ob@meta.data$nCount_RNA / seurat_ob@meta.data$nFeature_RNA
      names(tmp) <- row.names(seurat_ob@meta.data)
      seurat_ob <- Seurat::AddMetaData(seurat_ob, metadata = tmp, "UMI.per.Gene")
      vln4UPG <- Seurat::VlnPlot(object = seurat_ob,
                         features = "UMI.per.Gene",
                         group.by = groupfactor,
                         ncol = 1,
                         pt.size = pt.size) + labs(xlab = "") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + Seurat::NoLegend()
      gg_umi4gene_pca <- Seurat::FeaturePlot(seurat_ob,
                                     features = "UMI.per.Gene",
                                     reduction = "pca",
                                     cols = c("grey", "red"),
                                     pt.size = pt.size)
      gg_umi4gene_umap <- Seurat::FeaturePlot(seurat_ob,
                                      features = "UMI.per.Gene",
                                      reduction = "umap",
                                      cols = c("grey", "red"),
                                      pt.size = pt.size)
      if (assay2use == "Spatial") {

      gsp <- OESingleCell::SpatialPlot(seurat_ob,
                                  features =  "UMI.per.Gene",
                                  min.cutoff = "q0",
                                  max.cutoff = "q90",
                                  alpha = c(0.1, 1),
                                  pt.size.factor = opt$pointsize,
                                  pt.alpha=TRUE,
                                  stroke=0.2,combine=FALSE)

      gsb <- do.call(ggpubr::ggarrange,
                     c(gsp, list(nrow = 1,
                                 ncol = 1,
                                 common.legend = TRUE,
                                 legend = "right",
                                 align = "hv")))
        plot1 <- cowplot::plot_grid(vln4UPG, gsb,rel_widths= c(1, length(gsp)))
        plot2 <- cowplot::plot_grid(gg_umi4gene_pca, gsb,rel_widths= c(1, length(gsp)))
        plot3 <- cowplot::plot_grid(gg_umi4gene_umap, gsb,rel_widths= c(1, length(gsp)))
        width <-  3 * length(gsp) + 5 * (length(gsp)/1)
        height <- 6
      }
      else {
        plot1 <- vln4UPG
        plot2 <- gg_umi4gene_pca
        plot3 <- gg_umi4gene_umap
        width<-8
        heigh<-6
      }


      OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_mean_UMIs4each_gene_in_each_", opt$level, "_violin_plot")),
                                plot=plot1,
                                width=width,
                                height=height,
                                dpi=300,bg='white')
      OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_mean_UMIs4each_gene_in_each_", opt$level, "_PCA_plot")),
                                plot=plot1,
                                width=width,
                                height=height,
                                dpi=300,bg='white')
      OESingleCell::save_ggplots(file.path(output_dir,paste0(prefix, "_mean_UMIs4each_gene_in_each_", opt$level, "_umap_plot")),
                                plot=plot1,
                                width=width,
                                height=height,
                                dpi=300,bg='white')
    }














}


# scVis  \
# -i ./result/sctools/spatial.rds \
# -f rds \
# -o ./ \
# --image TRUE \
# --assay Spatial \
# qcmetrics
# -v "nGene,nUMI,percent.mito,batch" \
# --ptsize 0.1  \
# --alpha 0.6 \
# --pointsize 1.5 \
# -c sampleid \
# -M ./MT_genelist.txt \
# -l spot


