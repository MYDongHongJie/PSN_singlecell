sub_neighbor <- subparsers$add_parser("neighbor",  help = "neighbor analysis")
sub_neighbor$add_argument("--center", type = "character",
                          help = "the center celltype or tissue area.")
sub_neighbor$add_argument("--celltype", type = "character",
                          help = "the selected colname of matedata.")
sub_neighbor$add_argument("--palette", type = "character",default="blindless",
                          help = "color palette for celltype")
sub_neighbor$add_argument("--layer_num", type = "integer", default = 2,
                          help = "the layer number of neighor.[default: %(default)s]")
sub_neighbor$add_argument("--plotmarker", type = "character", default = FALSE,
                          help = "plot the DGE of different distance class[default: %(default)s]")
sub_neighbor$add_argument("--glm_nb", type = "character", default = FALSE,
                          help = "Negative two-term regression.[default: %(default)s]")
sub_neighbor$add_argument("--specie", type = "character", default = "human",
                          help = "for GO reference.[default: %(default)s]")
sub_neighbor$add_argument(
    "--crop",
    default = 'FALSE',
    help = "whether to crop in spatialplot, data from cytassist project should be 'TRUE'.[default: %(default)s]"
)
##必须注意：celltype_col中的细胞类型必须可以在names(meta。data)中匹配到。 GO富集暂时只有人和小鼠
args <- commandArgs(TRUE)
if ("neighbor" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "neighbor") {
    #===================================================================================================================
    futile.logger::flog.info("step1:read the specified assay and data slot in data object into memory")
    suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                    informat = opt$informat,
                                                    assays = assays,
                                                    data.use = dataslot,
                                                    verbose = F))
   # data_ob[[opt$celltype]]<-as.factor(data_ob[[opt$celltype]])
    
    #print(data_ob@meta.data[opt$celltype])
    if ('slice' %in% names(data_ob@meta.data) == FALSE) { data_ob$slice <- data_ob$sampleid }
    ##根据所有细胞类型，合理设置各个细胞类型配色方案
    colors_all <- OESingleCell::SelectColors(unique(data_ob[[opt$celltype]][,1]), 
                                             palette = opt$palette)
    print(colors_all)
                                    
    ##remove unused sample image in data_ob@image
    if (!is.null(opt$predicate)) {
      futile.logger::flog.info(glue::glue("get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
      df <- OESingleCell::colData(data_ob)
      desired_cells <- subset(df, eval(parse(text = opt$predicate)))
      data_ob <- subset(data_ob,cells=rownames(desired_cells))
      if (!is.null(Seurat::Images(data_ob))) {
        unuse_images<-Seurat::Images(data_ob)[ ! Seurat::Images(data_ob)  %in% (data_ob@meta.data$sampleid%>%unique)]
        if(length(unuse_images)>0){ data_ob@images[unuse_images]<-NULL}
        }
    }
    data_ob$sampleid <- factor(data_ob$sampleid, levels = unique(data_ob$sampleid)) 
    ##根据实际细胞类型，合理设置各个细胞类型配色方案

   
    #select_celltype<- unique(data_ob[[opt$celltype]])
    #print(select_celltype) 
    colors<-colors_all#[select_celltype]

    if(length(Seurat::Images(data_ob))>1){
      futile.logger::flog.info("报错：该分析一次只针对一个样本，请根据predicate参数来提取单个样本")
      quit()
    }
    #===================================================================================================================
    futile.logger::flog.info("step2: 根据参考组织计算欧式距离并绘图")
    center <- opt$center
    locat <- data_ob@images[[1]]@coordinates[, c("row", "col")]
    distest <- dist(locat, p = 2)
    center_barcode <- rownames(data_ob@meta.data)[grep(center, unlist(data_ob[[ opt$celltype]]))]
    distest2 <- as.data.frame(as.matrix(distest)) %>% .[rownames(data_ob@meta.data), center_barcode]
    distest2$dist <- apply(distest2, 1, min)
    center_distance <- paste0(center, "_distance")
    data_ob[[center_distance]] <- distest2$dist
    plot_data <- data.frame(celltype = factor(data_ob@meta.data[,  opt$celltype],
                                              levels = unique(data_ob@meta.data[,  opt$celltype])),
                            dist = distest2$dist)
    ##空间欧氏距离图
    crop <- as.logical(opt$crop)
    p1_slice <- OESingleCell::SpatialPlot(seurat_ob = data_ob,
                                          features = paste0(center, "_distance"),
                                          min.cutoff = 0.1,
                                          combine = FALSE,
                                          ncol = 2,
                                          cols=rev(OESingleCell::SelectColors(palette = "spectral",is.discrete = FALSE)),
                                          alpha =1,
                                          pt.size.factor = 1,
                                          image.alpha = 0.01,
                                          images = Seurat::Images(data_ob)[1],
                                          crop = crop
    )
    ##密度曲线图
    p1_density <- ggplot2::ggplot(plot_data) +
                  ggplot2::stat_density(ggplot2::aes(x = dist, colour = celltype),
                                        geom = "line",
                                        position = "identity",
                                        size = 0.2) +
                  ggplot2::scale_color_manual(values = colors) +
                  ggplot2::theme_set(ggplot2::theme_bw()) +
                  ggplot2::theme(panel.grid = ggplot2::element_blank()) +
                  ggplot2::xlab(paste0("distance to ", center, " of ", Seurat::Images(data_ob)[1])) +
                  ggplot2::scale_y_continuous(limits = c(0, 0.2))
    ##merge +密度曲线图
    layout <- glue::glue("{strrep('A',1)}{strrep('#',1)}\n{strrep('B',2)}")
    p1_slice_density <- patchwork::wrap_plots(A = p1_slice[[1]],
                                              B = p1_density,
                                              widths=c(1+stringr::str_length(paste0(center, "_distance"))/8,2),
                                              design = layout)

    OESingleCell::save_ggplots(glue::glue("{output_dir}/2.distance_analysis/{center}_distance_and_density"),
                               plot =  p1_slice_density ,
                               width = 10,
                               height = 8 ,
                               dpi = 1000)
    #===================================================================================================================
    futile.logger::flog.info("step3: 近邻点确认及绘图")
    #num取外围层数，推荐一层，两层
    get_loc <- function(loc, barcode, num) {
      row <- loc[barcode,]$row
      col <- loc[barcode,]$col
      col_list <- list()
      for (i in 1:(num + 1)) {
        part_col <- loc[loc$row %in% c((row - (i - 1)):(row + (i - 1))) & loc$col %in% c((col - 2 * num + (i - 1)):(col + 2 * num - (i - 1))),]
        col_list[[i]] <- rownames(part_col)
      }
      col_list <- unique(unlist(col_list))
      #     part_col <- loc[loc$row %in% c((row-1*num):(row+1*num))&loc$col %in% c((col-2*num):(col+2*num)),]
      #     return(rownames(part_col))
      #     part_col$spot_type <- "inter_spot"
      #     part_col[barcode,"spot_type"] <- "intra_spot"
      #     part_col=tibble::rownames_to_column(part_col,"barcode")
    }
    barcodes <- rownames(data_ob@meta.data[data_ob[[opt$celltype]][,1] == center,])
    allbarcode <- lapply(barcodes, get_loc, loc = locat, num = opt$layer_num)
    allbarcode <- unique(unlist(allbarcode))
    spottype <- data.frame(celltype = data_ob[[opt$celltype]])
    spottype$spottype <- NA
    spottype[rownames(spottype) %in% allbarcode & spottype[, 1] == center, 2] <- "intra_spot" #内点
    spottype[rownames(spottype) %in% allbarcode & spottype[, 1] != center, 2] <- "inter_spot" #间点
    spottype[is.na(spottype[, 2]), 2] <- "distal_spot" #外点
    center_spot <- paste0(center, "_spot_type")
    data_ob[[center_spot]] <- factor(spottype[, 2], levels = c("intra_spot", "inter_spot", "distal_spot"))
 
    p2_celltype <- OESingleCell::SpatialPlot(seurat_ob = data_ob,
                                             group.by = opt$celltype,
                                             cols = colors,
                                             combine = FALSE,
                                             pt.size.factor = 1,
                                             image.alpha = 0.01,
                                             crop = crop
                                             )
    p2_spottype <- OESingleCell::SpatialPlot(seurat_ob = data_ob,
                                             group.by = center_spot,
                                             cols = c("#FF5733", "#ED9481", "grey"),
                                             # min.cutoff = 0.1,
                                             combine = FALSE,
                                             ncol = 2,
                                             alpha = c(1.1, 1.5),
                                             pt.size.factor = 1,
                                             image.alpha = 0.01,
                                             images = Seurat::Images(data_ob)[1],
                                             crop = crop
    )
    p2_celltype_spottype <- patchwork::wrap_plots(A = p2_celltype[[1]],
                                                  B = p2_spottype[[1]],
                                                  widths=c(4+(data_ob@meta.data[[opt$celltype]] %>%
                                                              unique %>%
                                                              stringr::str_length()%>% max%>%+5) /10,
                                                           4+stringr::str_length(center_spot)/10),
                                                  design = "AB")
    # p2_celltype_spottype <- patchwork::wrap_plots(p2_celltype[[1]],
    #                                           p2_spottype[[1]],
    #                                           ncol=2,
    #                                           widths=c(1+stringi::stri_length("celltype")/8,
    #                                                   1+stringi::stri_length(center_spot)/8) )
    OESingleCell::save_ggplots(glue::glue("{output_dir}/1.co-region/spottype_and_celltype"),
                               plot = p2_celltype_spottype,
                               width = 8 +stringr::str_length("celltype")/10+stringr::str_length(center_spot)/10,
                               height = 4 ,
                               dpi = 1000)
    #===================================================================================================================
    futile.logger::flog.info("step4: 分面boxplot图-celltype")
    celltype <- c(names(table(data_ob[[opt$celltype]])), paste0(center, "_spot_type")) %>% data_ob[[.]]
    my_comparisons <- list(c("intra_spot", "inter_spot"), c("intra_spot", "distal_spot"), c("inter_spot", "distal_spot"))
    box_list <- list()
    for (i in 1:(dim(celltype)[2] - 1)) {
      p <- ggpubr::ggboxplot(celltype[, c(i, dim(celltype)[2])],
                             y = names(celltype)[i],
                             x = center_spot,
                             fill  = center_spot,
                             palette =c("#FF5733", "#ED9481", "grey"),
                             add = "jitter", shape = NA) +
           ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format")+
           ggplot2::theme(text = element_text(size = 10))+
           ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
           ggplot2::guides(fill = FALSE)
      box_list[[i]] <- p
    }
    p3 <- patchwork::wrap_plots(box_list)
    OESingleCell::save_ggplots(glue::glue("{output_dir}/1.co-region/spot_type_of_different_celltype"),
                               plot = p3,
                               width = 15,
                               height = 15,
                               dpi = 300)
    #====================================================================================================================
    futile.logger::flog.info("step5: 分面boxplot图-region")
 
    spots_type <-unique(data_ob@meta.data[,center_spot]) #names(table(data_ob[[center_spot]]))
  
    spot_celltype <- list()
    c_name <- unique(data_ob@meta.data[,opt$celltype])[!unique(data_ob@meta.data[,opt$celltype]) %in% c(opt$center)] %>%as.character
    print(c_name)
    for (i in 1:length(spots_type)) {
      data <- data_ob@meta.data %>%
              tibble::rownames_to_column("barcodes")%>%
              dplyr::filter(!!rlang::sym(center_spot)== spots_type[i]) %>% 
              dplyr::select(c("barcodes",c_name))  %>%
              tibble::column_to_rownames("barcodes") %>%
              reshape2::melt(, variable.name = 'celltype', value.name = 'ratio')
      #colors<-colors[c_name]
      p <- ggpubr::ggboxplot(data,
                             y = "ratio",
                             x = "celltype",
                             fill = "celltype",
                             palette = colors[c_name],
                             add = "jitter",
                             shape = NA,
                             title = spots_type[i]) +
           ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, size = 10, hjust = 1))+
           ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
           ggplot2::guides(fill = FALSE)
      spot_celltype[[i]] <- p
    }
    p4 <- patchwork::wrap_plots(spot_celltype)
    OESingleCell::save_ggplots(glue::glue("{output_dir}/1.co-region/celltype_of_different_spottype"),
                               plot = p4,
                               width = 9,
                               height = 5,
                               dpi = 300)
    ##==================================================================================================================
    if (opt$plotmarker == TRUE) {
      ##对于中心点，间点，远端点finderallmarker，并注释
      suppressPackageStartupMessages(library('clusterProfiler'))
      suppressPackageStartupMessages(library('ComplexHeatmap'))
      suppressPackageStartupMessages(library('circlize'))
      reflib <- switch(opt$specie,
                       "mouse" = { suppressPackageStartupMessages(library('org.Mm.eg.db'))
                         lib <- "org.Mm.eg.db" },
                       "human" = { suppressPackageStartupMessages(library('org.Hs.eg.db'))
                         lib <- "org.Hs.eg.db" })
      distance <- data_ob[[center_distance]][data_ob[[center_distance]][, 1] != 0, 1]
      center_distance_class <- paste0(center_distance, "_class")
      data_ob[[center_distance_class]] <- "start"
      data_ob[[center_distance_class]][data_ob[[center_distance]][, 1] <= quantile(distance, c(2 / 3)) & data_ob[[center_distance]][, 1] > quantile(distance, c(1 / 3)), 1] <- "middle"
      data_ob[[center_distance_class]][data_ob[[center_distance]][, 1] <= max(distance) & data_ob[[center_distance]][, 1] > quantile(distance, c(2 / 3)), 1] <- "end"
      data_ob[[center_distance_class]][data_ob[[center_distance]][, 1] == 0, 1] <- "reference"
      data_ob[[center_distance_class]] <- factor(data_ob[[center_distance_class]][, 1], levels = c("reference", "start", "middle", "end"))
      colors4 <- c('#33a02c', "#ff7f00", "#b6baf2", "#6a3d9a")
      names(colors4) <- c("reference", "start", "middle", "end")
      print(names(colors4))
      distance_class <- OESingleCell::SpatialPlot(seurat_ob = data_ob,
                                                  group.by = center_distance_class,
                                                  cols = colors4,
                                                  min.cutoff = 0.1,
                                                  combine = FALSE,
                                                  ncol = 2,
                                                  alpha = c(0.5, 2),
                                                  pt.size.factor = 1.5,
                                                  image.alpha = 0.01,
                                                  images = Seurat::Images(data_ob)[1],
                                                  crop = crop
      )
      OESingleCell::save_ggplots(glue::glue("{output_dir}/2.distance_analysis/{center}_distance_class_plot"),
                                 plot = distance_class[[1]],
                                 width = 7,
                                 height = 5,
                                 dpi = 1000)
      Seurat::DefaultAssay(data_ob) <- assays
      Seurat::Idents(data_ob) <- data_ob[[center_distance_class]][, 1]
      marker <- Seurat::FindAllMarkers(data_ob,
                                       assay = assays,
                                       slot = "data",
                                       only.pos = T,
                                       logfc.threshold = 0.25,
                                       min.pct = 0.25)
      # markers_raw <- Seurat::FindAllMarkers(object = data_ob,
      #                                 slot = "data",
      #                                 assay = assays,
      #                                 only.pos = T,
      #                                 test.use = "wilcox",
      #                                 logfc.threshold = 0,
      #                                 min.pct = 0 )
      # markers_raw <- markers_raw %>%
      #                # dplyr::rename(cluster=groupby) %>%
      #                dplyr::mutate( gene_diff = round(markers_raw$pct.1 / markers_raw$pct.2, 3)) %>%
      #                dplyr::select( gene, everything() )
      # print(head(markers_raw))
      # marker <-  markers_raw  %>%
      #            dplyr::group_by (cluster) %>%
      #            dplyr::filter(avg_log2FC >= 0.25) %>%
      #            dplyr::filter( p_val< 0.05 ) %>%
      #            dplyr::filter( pct.1 > 0.25  & pct.2 < 0.5) %>%
      #            dplyr::arrange(p_val,desc(avg_log2FC),desc(gene_diff)) %>%
      #            dplyr::filter(gene_diff > 1 )
      write.table(marker ,
            file = glue::glue("{output_dir}/2.distance_analysis/markers_for_4_defined_region.xls"),
            col.names =T,
            row.names = F,
            sep = "\t",
            quote=F)
      ego_list <- list()
      ego_net <- list()

      for (i in unique(marker$cluster)) {
        gene <- marker[marker$cluster == i, "gene"]
        ego_list[[i]] <- clusterProfiler::enrichGO(gene = gene,
                                                   OrgDb = lib,
                                                   ont = "ALL",
                                                   keyType = "SYMBOL",
                                                   pAdjustMethod = "BH",
                                                   pvalueCutoff = 1,
                                                   qvalueCutoff = 1)
        write.table(ego_list[[i]] ,
                    file=glue::glue("{output_dir}/2.distance_analysis/{i}_markergene_go_output.xls"),
                    sep = "\t", quote = F)
        if (dim(data.frame(ego_list[[i]]))[1] > 0) {
          ego_net[[i]] <- clusterProfiler::cnetplot(ego_list[[i]],
                                                    showCategory = 5,
                                                    categorySize = "pvalue",
                                                    circular = TRUE,
                                                    colorEdge = TRUE) }
      }
      ego_list <- ego_list[lapply(ego_list, function(x) { dim(x)[1] }) > 0]
      data_list <- list()
      for (i in names(ego_list)) {
        data <- data.frame(ego_list[[i]])
        data$spot_class <- i
        data_list[[i]] <- data
      }
      data_list <- do.call(rbind, data_list)
      write.table(data_list, file=glue::glue("{output_dir}/2.distance_analysis/markergene_go_output.xls"), sep = "\t", quote = F)
      go_data_gene_df <- list()
      for (i in unique(marker$cluster)) {
        go_data <- data.frame((ego_list[[i]]))
        go_data_gene <- go_data[go_data$Description %in% ego_net[[i]]$data$name[1:5], "geneID"] %>%
                        lapply(., strsplit, split = "/") %>%
                        unlist %>%
                        unique
        if (length(go_data_gene) > 0) {
          go_data_gene_df[[i]] <- data.frame(Seurat::GetAssayData(data_ob[go_data_gene,], slot = "scale.data"), genetype = i)
        }
      }
      names(go_data_gene_df) <- NULL
      ##整理好的用于画热图主体
      go_data_gene_df <- do.call(rbind, go_data_gene_df)
      gene_ann <- data.frame(gene_annot = go_data_gene_df[, dim(go_data_gene_df)[2]])
      rownames(gene_ann) <- rownames(go_data_gene_df)
      print(dim(go_data_gene_df)[2])
      go_data_gene_df <- go_data_gene_df[, -dim(go_data_gene_df)[2]] %>% .[, order(data_ob[[center_distance]][, 1])]
      gene_ann[, 1] <- factor(gene_ann[, 1], levels = unique(gene_ann[, 1]))
      ##注释
      ha_top <- HeatmapAnnotation(distance = data_ob[[center_distance]][, 1][order(data_ob[[center_distance]][, 1])])
      ha_bottom <- HeatmapAnnotation(distance_class = data_ob[[center_distance_class]][, 1][order(data_ob[[center_distance]][, 1])], col = list(distance_class = colors4))
      ha_row <- rowAnnotation(class = gene_ann[, 1], col = list(class = colors4[names(ego_list)]))

      heatmap1 <- ComplexHeatmap::Heatmap(go_data_gene_df,
                                          top_annotation = ha_top,
                                          left_annotation = ha_row,
                                          bottom_annotation = ha_bottom,
                                          cluster_rows = FALSE,
                                          cluster_columns = FALSE,
                                          show_column_names = F,
                                          show_row_names = F,
                                          row_split = gene_ann[, 1],
                                          row_title_rot = 0)
      heatmap1 <- ggplotify::as.ggplot(heatmap1)
      OESingleCell::save_ggplots(glue::glue("{output_dir}/2.distance_analysis/marker_heatmap"),
                                 plot = heatmap1,
                                 width = 8,
                                 height = 8,
                                 dpi = 200, bg = 'white')


      get_map <- function(x) {
        test_list <- list()
        for (i in 1:dim(x)[1]) {
          test_list[[i]] <- unlist(strsplit(x[i, 2], split = "/"))

        }
        names(test_list) <- x[, 1]
        test_data <- data.frame(to = unlist(test_list)) %>% mutate(from = gsub("[0-9]", "", rownames(.)))
        return(test_data)

      }

      colors_light <- c('#d1c4ff', "#ff785b", "#F7D4CC", "#cccccc")
      names(colors_light) <- names(colors4)

      map_list <- list()
      grid.col_list <- list()
      for (i in names(ego_list)) {
        go_data <- data.frame((ego_list[[i]]))
        test <- go_data[go_data$Description %in% ego_net[[i]]$data$name[1:5], c("Description", "geneID")]
        mapdata <- get_map(test)
        map_list[[i]] <- mapdata
        test_col_go <- colorRampPalette(c(colors_light[i], colors4[i]))(5)
        names(test_col_go) <- unique(mapdata$from)
        test_col_go <- test_col_go[mapdata$from]
        test_col_gene <- test_col_go
        names(test_col_gene) <- mapdata$to
        grid.col_list[[i]] <- c(test_col_go, test_col_gene)
      }


      # par(mfrow=c(2,2))
      pdf(glue::glue("{output_dir}/2.distance_analysis/marker_go.pdf"), w = 8, h = 8)
      # heatmap1
      for (i in names(ego_list)) {
        circlize::chordDiagram(map_list[[i]], annotationTrack = "grid", grid.col = grid.col_list[[i]],
                               preAllocateTracks = list(track.height = 0.075))
        # we go back to the first track and customize sector labels
        circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
          circlize::circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
        }, bg.border = NA)
        title(paste0("top 5 GO annotation of ", i))
        circos.clear()

      }
      dev.off()
    }

    ##==================================================================================================================
    if (opt$glm_nb == TRUE) {
      suppressPackageStartupMessages(library('clusterProfiler'))
      reflib <- switch(opt$specie,
                       "mouse" = { suppressPackageStartupMessages(library('org.Mm.eg.db'))
                         lib <- "org.Mm.eg.db" },
                       "human" = { suppressPackageStartupMessages(library('org.Hs.eg.db'))
                         lib <- "org.Hs.eg.db" })
      suppressPackageStartupMessages(library('circlize'))
      suppressPackageStartupMessages(library('ComplexHeatmap'))
      output_dir2 <- file.path(output_dir, "2.distance_analysis/glmresult")
      if (!file.exists(output_dir2)) { dir.create(output_dir2, recursive = T) }
      output_dir2 <- normalizePath(output_dir2)
      data <- Seurat::GetAssayData(data_ob)
      future::plan("multiprocess", workers = 20)
      rate <- apply(data, 1, function(x) { length(which(x != 0)) })
      rate2 <- (rate / dim(data)[2]) * 100
      data_f1 <- data[rate2 > 10,]
      data_f1_cv <- apply(data_f1, 1, function(x) { 100 * sd(x) / mean(x) })
      data_f2 <- data.frame(t(data.frame(data_f1[order(data_f1_cv, decreasing = TRUE)[1:5000],])))
      rownames(data_f2) <- gsub("\\.", "-", rownames(data_f2))
      data_f2$dis <- data_ob[[center_distance]][rownames(data_f2),]

      get_glmnb <- function(x) {
        data <- data_f2[, c(x, "dis")]
        names(data) <- c("gene", "dis")
        test <- MASS::glm.nb(gene ~ dis, data = data)
        result <- data.frame(Estimate = summary(test)$coefficients[2, 1],
                             Pr = summary(test)$coefficients[2, 4])
        rownames(result) <- x
        return(result)
      }

      data_glmnb <- future.apply::future_lapply(names(data_f2), get_glmnb)
      save(data_glmnb, file = file.path(output_dir2, "data_glmnb.RData")) ##给你加一个保险
      data_glmnb <- do.call(rbind, data_glmnb)
      write.table(data_glmnb,
                  file = file.path(output_dir2, paste0(center_distance, "_glm_result.txt")),
                  sep = "\t",
                  quote = F)
      ####只选取p值<0.01的基因
      data_glmnb_filter_gene <- rownames(data_glmnb[data_glmnb[, 2] < 0.01,])
      Seurat::DefaultAssay(data_ob) <- assays
      data_glmnb_filter <- Seurat::GetAssayData(data_ob[data_glmnb_filter_gene,], slot = "data") %>% data.frame()
      future::plan("multiprocess", workers = opt$ncores)
      dis <- factoextra::get_dist(data_glmnb_filter, method = "pearson")
      hc <- hclust(dis, method = "ward.D2")
      p <- factoextra::fviz_dend(hc, k = 8, cex = 0.7, rect = TRUE)
      OESingleCell::save_ggplots(file.path(output_dir2, "hclust_of_distancegenes"),
                                 plot = p,
                                 width = 5,
                                 height = 5,
                                 dpi = 200)
      save(dis, hc, p, file = file.path(output_dir2, "dis_hc_p.RData")) ##再加一个保险
      ##提取聚簇标签
      cluster_im <- ggplot2::ggplot_build(p)$data[[2]]
      group <- unique(cluster_im$group)%>% sort()
      ##对聚簇富集分析
      ego_list <- list()

      for (i in group) {
        gene <- cluster_im[cluster_im$group == i, "label"]
        ego_list[[i]] <- clusterProfiler::enrichGO(gene = gene,
                                                   OrgDb = lib,
                                                   ont = "ALL",
                                                   keyType = "SYMBOL",
                                                   pAdjustMethod = "BH",
                                                   pvalueCutoff = 1,
                                                   qvalueCutoff = 1)
        write.table(ego_list[[i]] , file=glue::glue("{output_dir2}/enrichment_{i}_markergene_go.xls"), sep = "\t", quote = F)
      }
      names(ego_list) <- group
      ego_list <- ego_list[unlist(lapply(ego_list, function(x) { dim(x)[1] > 0 }))]
      ego_list_output <- lapply(ego_list, data.frame)
      for (i in names(ego_list_output)) {

        ego_list_output[[i]]$cluster <- i
      }
      ego_list_output %>%
        do.call(rbind, .) %>%
        write.table(., file = file.path(output_dir2, "go_for_distancegene.txt"), sep = "\t", quote = FALSE)

      get_map <- function(x) {
        test_list <- list()
        for (i in 1:dim(x)[1]) {
          test_list[[i]] <- unlist(strsplit(x[i, 2], split = "/"))

        }
        names(test_list) <- x[, 1]
        test_data <- data.frame(to = unlist(test_list)) %>% mutate(from = gsub("[0-9]", "", rownames(.)))
        return(test_data)

      }

      top_goerich <- function(x, name) {
        if (any(x@result$p.adjust < 0.01)) {
          data <- data.frame(x) %>%
                  dplyr::filter(p.adjust < 0.01) %>%
                  dplyr::top_n(5, Count) %>%
                  dplyr::select(c(Description, geneID))
          data <- get_map(data)
          data$cluster <- name
          return(data)
        }else {
          return(NULL)
        }

      }
      geo_data <- list()
      for (i in names(ego_list)) {
        i <- as.character(i)
        data <- top_goerich(x = ego_list[[i]], name = i)
        geo_data[[i]] <- data

      }
      go_data_gene_df <- do.call(rbind, geo_data)
      gene_ann <- data.frame(gene = go_data_gene_df$to, gene_annot = go_data_gene_df[, dim(go_data_gene_df)[2]])
      gene_ann[, 2] <- factor(gene_ann[, 2], levels = unique(gene_ann[, 2]))
      dis_order <- data_ob[[center_distance]][, 1][order(data_ob[[center_distance]][, 1])]
      go_data_gene_df <- Seurat::GetAssayData(data_ob[go_data_gene_df$to,], slot = "scale.data")
      go_data_gene_df <- go_data_gene_df[, order(data_ob[[center_distance]][, 1])]
      go_data_gene_df <- go_data_gene_df[gene_ann[, 1],] %>% data.frame()
      ###Loess smoothing
      #https://github.com/IStevant/XX-XY-mouse-gonad-scRNA-seq/blob/b2068562edd17331b85021f02b4b2483860ed55a/scripts/analysis_functions.R
      ##平滑热图参考 https://github.com/IStevant/XX-XY-mouse-gonad-scRNA-seq/blob/master/scripts/XX_analysis_dm.R
      smooth_gene_exp <- function(data = data, pseudotime = pseudotime, span = 0.75) {
        smooth_data <- data
        for (gene in 1:nrow(data)) {
          gene_exp <- t(data[gene,])
          smooth <- loess(formula = gene_exp ~ pseudotime, span = span)
          smooth_data[gene,] <- predict(smooth, newdata = pseudotime)
        }
        return(smooth_data)
      }

      go_data_gene_df <- smooth_gene_exp(go_data_gene_df, dis_order, span = 0.4) %>%
                         t() %>%
                         scale() %>%
                         t() %>%
                         data.frame()
      cold <- colorRampPalette(c('#f7fcf0', '#41b6c4', '#253494', '#081d58', '#081d58'))
      warm <- colorRampPalette(c('#ffffb2', '#fecc5c', '#e31a1c', '#800026', '#800026'))
      mypalette <- c(rev(cold(21)), warm(20))
      heatmapcol <- OESingleCell::SelectColors(gene_ann[, 2], palette = "customecol2")
      ha_row <- ComplexHeatmap::rowAnnotation(class = gene_ann[, 2], col = list(class = heatmapcol))
      ha_top <- ComplexHeatmap::HeatmapAnnotation(distance = dis_order)
      heatmap <- ComplexHeatmap::Heatmap(go_data_gene_df,
                                         name = 'gene exp',
                                         top_annotation = ha_top,
                                         left_annotation = ha_row,
                                         cluster_rows = FALSE,
                                         cluster_columns = FALSE,
                                         show_column_names = F,
                                         show_row_names = F,
                                         row_split = gene_ann[, 2],
                                         row_title_rot = 0,
                                         col = circlize::colorRamp2(c(-4, -1, 0, 1, 4),
                                                                    c(mypalette[15],
                                                                      mypalette[19],
                                                                      mypalette[21],
                                                                      mypalette[22],
                                                                      mypalette[41]))) %>% ggplotify::as.ggplot()

      OESingleCell::save_ggplots(file.path(output_dir2, "distance_gene_heatmap"),
                                 plot = heatmap,
                                 width = 5,
                                 height = 5,
                                 dpi = 100)
      ###这也是保险
      save(gene_ann, go_data_gene_df, dis_order, heatmapcol, file = file.path(output_dir2, "fordis_heatmap.RData"))
      ##淡紫色#ddd2f9
      #geo_data 用于go circlize图
      #heatmapcol 簇颜色
      col_dis_list <- list()
      for (i in sort(names(geo_data))) {
        col_t <- colorRampPalette(c("#ddd2f9", heatmapcol[i]))(length(unique(geo_data[[i]]$from)))
        names(col_t) <- unique(geo_data[[i]]$from)
        num <- match(geo_data[[i]]$from, names(col_t))
        geo_data[[i]]$col <- col_t[num]
        col_gene <- geo_data[[i]]$col
        names(col_gene) <- geo_data[[i]]$to
        col_dis_list[[i]] <- c(col_t, col_gene)
      }
      pdf(file.path(output_dir2, "top_5_distance_gene_GO.pdf"), h = 10, w = 8)
      for (i in sort(names(geo_data))) {
        circlize::chordDiagram(geo_data[[i]], annotationTrack = "grid", grid.col = col_dis_list[[i]],
                               preAllocateTracks = list(track.height = 0.075))
        # we go back to the first track and customize sector labels
        circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
          circlize::circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
        }, bg.border = NA)
        title(glue::glue("top 5 GO annotation of cluster {i}"))
        circos.clear()

      }
      dev.off()

      ###切片集合
      go_top5 <- do.call(rbind, geo_data)
      data <- go_top5 %>%
        split(., .$cluster) %>%
        lapply(., function(x) { x <- split(x, x$from) })
      for (i in 1:length(data)) {
        names(data[[i]]) <- gsub(" ", "_", names(data[[i]])) %>%
          gsub("-", "_", .) %>%
          gsub(",", "", .)
      }
      for (j in 1:length(data)) {
        for (i in names(data[[j]])) {
          geneset <- data[[j]][[i]]$to
          data_ob <- Seurat::AddModuleScore(data_ob, features = list(geneset), name = i)
        }
      }

      spplot_list <- function(go_gene, data_ob) {
        name_t <- paste0(names(go_gene), "1")
        p_list <- list()
        for (i in name_t) {
          p_test <- OESingleCell::SpatialPlot(seurat_ob = data_ob, 
                                              features = i, 
                                              cols = "horizon",
                                              min.cutoff = 0.1,
                                              max.cutoff = 1,
                                              combine = FALSE,
                                              ncol = 2,
                                              alpha = c(0.1, 0.5),
                                              pt.size.factor = 0.8,
                                              image.alpha = 0.01,
                                              images = Seurat::Images(data_ob)[1],
                                              crop = crop
          )
          titstr <- unlist(strsplit(i, split = "_"))
          if (length(titstr) > 3) {
            titstr <- paste(paste(titstr[1:3], collapse = " "),
                            paste(titstr[4:length(titstr)], collapse = " "),
                            sep = " \n ")
          }else {
            titstr <- paste(titstr, collapse = " ")
          }

          p_list[[i]] <- p_test[[1]] +
            ggplot2::labs(title = titstr) +
            ggplot2::guides(fill = ggplot2::guide_colorbar(title = ""))

        }
        return(p_list)
      }

      image_list <- lapply(data, spplot_list, data_ob)
      for (i in names(image_list)) {
        plot_top5 <- do.call(ggpubr::ggarrange, c(image_list[[i]],
                                                  list(ncol = 3,
                                                       nrow = ceiling(length(image_list[[i]]) / 3)),
                                                  common.legend = TRUE,
                                                  legend = "right",
                                                  align = "hv"))
        OESingleCell::save_ggplots(glue::glue("{output_dir2}/cluster_{i}_top5_go"),
                                   plot = plot_top5,
                                   width = 3 * 3,
                                   height = 3 * ceiling(length(image_list[[i]]) / 3),
                                   dpi = 200)
      }
      ###go circle图
    }

    ##==================================================================================================================
    ### save seurat object
    OESingleCell::SaveX(data_ob,
                        output = output_dir,
                        update = FALSE,
                        outformat = opt$outformat,
                        prefix = opt$prefix)
    if (!file.exists(file.path(output_dir, "基于层数及欧式距离定义的空间近邻分析说明文档.docx"))) {
        file.copy("/public/dev_scRNA/oesinglecell3_test/document/基于层数及欧式距离定义的空间近邻分析说明文档.docx",
                  file.path(output_dir, "基于层数及欧式距离定义的空间近邻分析说明文档.docx")) }
    ## save session informations
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
}