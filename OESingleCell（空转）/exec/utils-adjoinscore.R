sub_adjoin <- subparsers$add_parser("adjoin", help = "Identify and quantify cluster-cluster or incluster neighbours.")
sub_adjoin$add_argument("--sample",type="character",default=NULL,
                         help = "select sample to analyze.Default is all.[default: %(default)s].")
sub_adjoin$add_argument("--col_name",type="character",default = "clusters",
                         help = "[OPTIONAL]the colname of metadata to plot.[default: %(default)s]")
sub_adjoin$add_argument("--imageinf",type="character",
                         help = "[OPTIONAL]the image information used to created STobject.[default: %(default)s]")
sub_adjoin$add_argument("--inclust",type="character",default = "TRUE",
                         help = "Visium neighbourhood analysis of cluster community formation or aggregation .[default: %(default)s]")
sub_adjoin$add_argument("--betclust",type="character",default = "TRUE",
                         help = "Identify and quantify cluster-cluster neighbours.[default: %(default)s]")
args <- commandArgs(TRUE)
if ( "adjoin" %in% args ){
  opt<-intial_setting()
  library(ggplot2)
  library(STutility)
  col_name <- opt$col_name
  object <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = opt$assay,
                                 verbose = F,
                                 images=opt$image)
  if(!is.null(opt$sample)){
      sample <- unlist(strsplit(sample,split=","))
      difsample <- setdiff(names(object@images),sample)
      object <- subset(object,sampleid %in% sample)
      for(i in difsample){
          object@images[[i]] <- NULL
      }
  }
  print('stp1')
  output_dir_homotypic <- file.path(output_dir, 'homotypic')
  dir.create(output_dir_homotypic, recursive = T)

  col <- OESingleCell::SelectColors(unique(object[[col_name]][,1]),'ditto')
  object$celltype_col <- col[object[[col_name]][,1]]#确定颜色


  rename_sample <- seq(unique(object$sampleid))
  names(rename_sample) <- unique(object$sampleid)
  object$sampleid_num <- rename_sample[object$sampleid]

  inf_data <- read.table(opt$imageinf,sep=",",quote = "")

  imgs <- inf_data[grep(".png",inf_data[,2]),] %>%
          tibble::rownames_to_column() %>%
          tibble::column_to_rownames("V1") %>%
          .[names(rename_sample),'V2']

  imdata <- list()
  for (i in names(rename_sample)){
      imdata[[i]] <- data.frame(x=object@images[[i]]@coordinates$col,
                                y=object@images[[i]]@coordinates$row,
                                adj_x=object@images[[i]]@coordinates$col,
                                adj_y=object@images[[i]]@coordinates$row,
                                pixel_x=object@images[[i]]@coordinates$imagecol*object@images[[1]]@scale.factors$hires ,
                                pixel_y=object@images[[i]]@coordinates$imagerow*object@images[[1]]@scale.factors$hires ,
                                sample=rename_sample[i]%>% as.character ,
                                original_x=object@images[[i]]@coordinates$imagecol*object@images[[1]]@scale.factors$hires,
                                original_y=object@images[[i]]@coordinates$imagerow*object@images[[1]]@scale.factors$hires)
  #sampleid = names(object@images[i]) )
  rownames(imdata[[i]]) <- rownames(object@images[[i]]@coordinates)

  }
  imdata <- dplyr::bind_rows( imdata)
  print('stp2')
  object@tools$Staffli <- STutility::CreateStaffliObject(imgs = imgs,
                                             meta.data = imdata,
                                             platforms = rep('Visium',length(rename_sample)))

  object@meta.data$id <- seq(1:dim(object)[2])
  object@meta.data$labels <- "Default"
  print('stp3')
  ##add jsonfile information

  jsonfile <- inf_data[grep(".json",inf_data[,2]),] %>%
          tibble::rownames_to_column() %>%
          tibble::column_to_rownames("V1") %>%
          .[names(rename_sample),'V2']
  names(jsonfile) <- unname(rename_sample)

  scaleVisium <- NULL
  scaleVisium <- scaleVisium %||% sapply(jsonfile, function(f) {jsonlite::read_json(f)$tissue_hires_scalef})
  spot.diamater.list <- sapply(jsonfile, function(f) {jsonlite::read_json(f)$spot_diameter_fullres})*scaleVisium
  ranges.list <- list()
  for(i in object@tools$Staffli@samplenames){
      ranges.list[[i]] <- sapply(imdata[imdata$sample==i, c("pixel_x", "pixel_y")], range)
  }

  dims.list <- lapply(object@tools$Staffli@samplenames, function(i) {
        rs <- ranges.list[[i]]
        data.frame(min_x = rs[1, 1], max_x = rs[2, 1], min_y = rs[1, 2], max_y = rs[2, 2], spot_diameter = spot.diamater.list[i])
      })
  object@tools$Staffli@dims <- dims.list

  object <- STutility::LoadImages(object, time.resolve = T, verbose = T, xdim = 100)
  if(opt$inclust=='TRUE'){
    ##多样本簇内构建网络
    spatnet_init <- STutility::GetSpatNet(object)

    ##注释barcodes的细胞类型，样本来源
    spatnet <- do.call(rbind, lapply(seq_along(spatnet_init), function(i) {
      spnet <- spatnet_init[[i]]
      spnet$cluster_from <- object[[]][spnet$from, col_name]
      spnet$cluster_to <- object[[]][spnet$to, col_name]
      spnet$sample <- paste0(i)
      return(spnet)
    }))
    # spatnet$cluster_from <- as.factor(as.numeric(spatnet$cluster_from))
    # spatnet$cluster_to <- as.factor(as.numeric(spatnet$cluster_to))

    colors_clusters <- object@meta.data[, c("clusters", "sampleid_num","celltype_col",col_name)] %>%
      dplyr::distinct() %>%
      dplyr::arrange(clusters)#根据clusters排序


    ####Q:为什么横向在视觉效果上看起来比斜对面连线短，却没有连线（即from自己to自己）因为横向点与点是挨着的，没有连线说明即使在短距内，也没贴邻
    #####依然存在其他细胞类型点
    num <- Seurat::Images(object) %>% length()
    for(cluster_view in unique(object[[col_name]][,1])){
      spatnet.subset <- subset(spatnet, cluster_from %in% cluster_view)
      spatnet.subset_conly <- subset(spatnet.subset, cluster_to %in% cluster_view)
      if(dim(spatnet.subset_conly)[1]>0){
        degree <- table(spatnet.subset_conly$from) %>% data.frame %>% tibble::column_to_rownames('Var1')
        spatnet.subset$degree <- degree[spatnet.subset$from,1]
        spatnet.subset$degree[is.na(spatnet.subset$degree)] <- 0.5
        spatnet.subset$degree <- as.character(spatnet.subset$degree)
        cell_color <- colors_clusters[colors_clusters[[col_name]]==cluster_view,"celltype_col"] %>% unique()
        cell_color <- colorRampPalette(c("white",cell_color))(7)
        names(cell_color) <- as.character(c(0.5,1:6))

        turn_sample <- data.frame(object$sampleid_num,object$sampleid)[,c(1,2)] %>%
              unique %>%
              tibble::rownames_to_column('barcode') %>%
              tibble::column_to_rownames('object.sampleid_num')
        spatnet.subset$sampleid <- turn_sample[spatnet.subset$sample,2]
        spatnet.subset_conly$sampleid <- turn_sample[spatnet.subset_conly$sample,2]

        p_net <- ggplot2::ggplot() +
          ggplot2::geom_segment(data = spatnet.subset_conly, ggplot2::aes(x = start_x, xend = end_x, y = -start_y, yend = -end_y), size=0.1) +
          ggplot2::geom_point(data = spatnet.subset, ggplot2::aes(start_x, -start_y,color = degree), size = .2) +
          ggplot2::labs(color=cluster_view) +
          ggplot2::scale_color_manual(values =cell_color) +
          ggplot2::facet_wrap(~sampleid)+
         theme(panel.background = element_rect(fill = "grey80",
                                      colour = "grey80",size = 0.5))+
         theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
          w=2*num
          h=2*num
          print(cluster_view)

        OESingleCell::save_ggplots(file.path(output_dir_homotypic , paste0(cluster_view,"obs_without_nourmalized.pdf")),
                                  plot=p_net,
                                  width=w,
                                  height=h,
                                  dpi=200,bg='white')
      }
    }
    #dev.off()
    n_samples <- unique(object$sampleid_num)

    k_avg_list <- list()
    k_avg_df <- data.frame(row.names = paste0("S", n_samples))

    length(unique(c(spatnet.subset$from)))

    ##计算观测平均值

    for(cluster_view in unique(object[[col_name]][,1])){##细胞类型
      for(s in n_samples){
        spatnet.subset <- subset(spatnet, cluster_from %in% cluster_view)
        spatnet.subset <- subset(spatnet.subset, sample == s)
        spatnet.subset_conly <- subset(spatnet.subset, cluster_to %in% cluster_view)

        N <- length(unique(c(spatnet.subset$from))) #点
        ki <- spatnet.subset_conly %>%
          dplyr::group_by(from) %>%
          dplyr::count(name = "ki")
        L <- sum(ki$ki)/2 ####边
        k_avg <- (2*L)/N #每一个
        k_avg_list[[paste0("S", s)]] <- k_avg
      }
      k_avg_df_add <- as.data.frame(unlist(k_avg_list))
      colnames(k_avg_df_add) <- paste0("kavg_cluster_", cluster_view)
      k_avg_df <- cbind(k_avg_df, k_avg_df_add)
    }

    k_avg_df$sample <- rownames(k_avg_df)

    write.table(k_avg_df,file=file.path(output_dir_homotypic ,'observed_k_avg_dif.txt'),sep="\t",quote = F)


    ##构建预测值
    n_perm <- 50 ###打乱标签50次
    avgk_df_perm_list <- list()
    for(i in seq(1, n_perm)){
      avgk_df_perm_list[[i]] <- OESingleCell::RandomClusteringSpatNet(object,
                                                        column.cluster.id = col_name,
                                                        column.sample.id = "sampleid_num",
                                                        se.SpatNet = spatnet_init,
                                                        random.seed = i
                                                        )
    }

    avgk_df_perm_avg <- data.frame(row.names = row.names(avgk_df_perm_list[[1]]))
    avgk_df_perm_sd <- data.frame(row.names = row.names(avgk_df_perm_list[[1]]))
    avgk_df_perm_me <- data.frame(row.names = row.names(avgk_df_perm_list[[1]]))

    n_cols <- dim(avgk_df_perm_list[[1]])[2]
    for(c in seq(1:n_cols)){
      cluster_avgks <- do.call(rbind, lapply(avgk_df_perm_list, `[[`, c))##筛选出每一次随机中循环到的簇
      avgk_df_perm_avg$c_avg <- colMeans(cluster_avgks)
      avgk_df_perm_sd$c_sd <-  matrixStats::colSds(cluster_avgks)
      avgk_df_perm_me$c_me <- qnorm(.95) * ( avgk_df_perm_avg$c_avg / sqrt(length(avgk_df_perm_list)) )  # margin of error at 90% CI
      colnames(avgk_df_perm_avg)[colnames(avgk_df_perm_avg)=="c_avg"] <-
        colnames(avgk_df_perm_sd)[colnames(avgk_df_perm_sd)=="c_sd"] <-
        colnames(avgk_df_perm_me)[colnames(avgk_df_perm_me)=="c_me"] <-
        colnames(avgk_df_perm_list[[1]])[c]
    }
    ##以上，c_avg为每一个簇的50次随机的平均，sd为每个簇中50此随机的标准差。以此类推

    se_stats <- object@meta.data %>%
      dplyr::group_by(sampleid, sampleid_num,clusters ,!!!syms(c(col_name))) %>%
      dplyr::count(name = "spots")

    ###实际值减去模拟值
    k_avg_perm_diff <- k_avg_df[rownames(avgk_df_perm_avg), colnames(avgk_df_perm_avg)] - avgk_df_perm_avg
    k_avg_perm_zscore <- (k_avg_df[rownames(avgk_df_perm_avg), colnames(avgk_df_perm_avg)] - avgk_df_perm_avg) / avgk_df_perm_sd

    write.table(k_avg_perm_diff,file=file.path(output_dir_homotypic ,'k_avg_perm_diff.txt'),sep="\t",quote = F)
    write.table(k_avg_perm_zscore,file=file.path(output_dir_homotypic ,'k_avg_perm_zscore.txt'),sep="\t",quote = F)

    k_avg_perm_zscore

    summary_df_kavg_diff <- as.data.frame(t(do.call(cbind, lapply(k_avg_perm_diff, summary))))
    summary_df_kavg_diff$cluster <- as.character(sub(pattern = "kavg_cluster_", "", rownames(summary_df_kavg_diff))) %>% factor(.,levels=.)

    ##细胞类型在不同样本中，实际与模拟差值最大是哪个？
    p1 <- ggplot(summary_df_kavg_diff, aes(x=reorder(cluster, Max.), y=Max.)) +
      geom_col(fill=col[summary_df_kavg_diff$cluster[order(summary_df_kavg_diff$Max.,decreasing = TRUE)]])+
      labs(x="") +
      coord_flip() +
      theme_classic()

    ##不同细胞类型在几个样本中的平均（实际与模拟差值）
    p2 <- ggplot(summary_df_kavg_diff, aes(x=reorder(cluster, Mean), y=Mean)) +
      geom_col(fill=col[summary_df_kavg_diff$cluster[order(summary_df_kavg_diff$Mean,decreasing = TRUE)]])+
      labs(x="") +
      coord_flip() +
      theme_classic()

    p=p1 + p2 + patchwork::plot_annotation(title = '<k> obs-exp difference') & theme(plot.title = element_text(hjust=0.5))
    OESingleCell::save_ggplots(file.path(output_dir_homotypic ,'obs_exp_difference'),
                                plot=p,
                                width=10,
                                height=6,
                                dpi=200,bg='white')


    df_kavg_diff <- as.data.frame(t(k_avg_perm_diff))
    df_kavg_diff$cluster <- as.character(sub(pattern = "kavg_cluster_", "", rownames(df_kavg_diff)))
    df_kavg_diff_long <- reshape2::melt(df_kavg_diff, id.vars = c("cluster"))

    color <- col %>% .[df_kavg_diff_long$cluster]

    ##不同细胞类型在不同样本中的（实际与观测值之差）的箱式图
    #实际-平均（50次模拟）
    p3 <- ggplot(df_kavg_diff_long, aes(x=reorder(cluster, value),fill=cluster, y=value)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      geom_point(color=color, size=1) +scale_fill_manual(values=c(color))+
      labs(x="", y="<k>", title="<k> obs-exp diff.") +
      coord_flip() +
      theme_classic()

    p4 <- ggplot(df_kavg_diff_long, aes(x=reorder(cluster, value), y=value)) +
      geom_col(data = summary_df_kavg_diff, aes(x=reorder(cluster, Mean), y=Mean),fill="white",color = "black", width = .8) +
      geom_hline(yintercept = 0) +
      # geom_boxplot() +
      geom_point(color=color, size=1) +
      labs(x="", y="<k>", title="<k> obs-exp diff.") +
      coord_flip() +
      theme_classic()

    cowplot::plot_grid(p4, p3,rel_widths= c(1, 2.5)) %>%
    OESingleCell::save_ggplots(file.path(output_dir_homotypic ,'obs_exp_diff'),
                                    plot=.,
                                    width=8,
                                    height=6,
                                    dpi=200,bg='white')
    summary_df_kavg_zscore <- as.data.frame(t(do.call(cbind, lapply(k_avg_perm_zscore, summary))))
    summary_df_kavg_zscore$cluster <- as.character(sub(pattern = "kavg_cluster_", "", rownames(summary_df_kavg_zscore))) %>% factor(.,levels=.)
    colcor <- col[as.character(summary_df_kavg_zscore$cluster[order(summary_df_kavg_zscore$Max.,decreasing = T)])]
    ##zscore后的值

    p5 <- ggplot(summary_df_kavg_zscore[order(summary_df_kavg_zscore$Max.,decreasing = T), ], aes(x=reorder(cluster, Max.), y=Max.)) +
      geom_col(fill=col[summary_df_kavg_zscore[order(summary_df_kavg_zscore$Max.,decreasing = T),'cluster']]) +
      labs(x="") +
      coord_flip() +
      theme_classic()
    p6 <- ggplot(summary_df_kavg_zscore[order(summary_df_kavg_zscore$Mean,decreasing = T), ],aes(x=reorder(cluster, Mean), y=Mean)) +
      geom_col(fill=col[summary_df_kavg_zscore[order(summary_df_kavg_zscore$Mean,decreasing = T),'cluster']]) +
      labs(x="") +
      coord_flip() +
      theme_classic()
    p5_6 <- p5+p6 + patchwork::plot_annotation(title = '<k> obs-exp z-score') & theme(plot.title = element_text(hjust=0.5))
    OESingleCell::save_ggplots(file.path(output_dir_homotypic ,'obs_exp_difference_zscore'),
                                    plot=p5_6 ,
                                    width=9,
                                    height=5,
                                    dpi=200,bg='white')
    df_kavg_zscore <- as.data.frame(t(k_avg_perm_zscore))
    df_kavg_zscore$cluster <- as.character(sub(pattern = "kavg_cluster_", "", rownames(df_kavg_zscore)))
    df_kavg_zscore_long <- reshape2::melt(df_kavg_zscore, id.vars = c("cluster"))
    p7 <- ggplot2::ggplot(df_kavg_zscore_long, ggplot2::aes(x=reorder(cluster, value), fill=cluster,y=value)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_boxplot() +
      ggplot2::geom_point(color=color, size=.5) +ggplot2::scale_fill_manual(values=c(col[df_kavg_zscore_long$cluster]))+
      ggplot2::labs(x="", y="z-score", title="<k> obs-exp diff.") +
      ggplot2::coord_flip() +
      ggplot2::theme_classic()
    OESingleCell::save_ggplots(file.path(output_dir_homotypic ,'obs_exp_diff_zscore'),
                                    plot=p7,
                                    width=8,
                                    height=6,dpi=200,bg='white')

  }

###簇间

  if(opt$betclust=='TRUE'){
    output_dir_heterotype <- file.path(output_dir, 'heterotype')
    dir.create(output_dir_heterotype, recursive = T)
    c_include <- unique(object[[col_name]][,1])

    ##构建随机模拟数据50次，每次：获取每一种细胞类型的接邻数量
    n_perm <- 50
    perm_adj_mat_list <- list()
    for(i in 1:n_perm){
        ##模拟标签（50次）为了可复原结果，以i设置seed
      object <- OESingleCell::RandomiseClusteringIDs(se.object = object,
                                                     column.cluster.id = col_name,
                                                     column.sample.id = 'sampleid_num',
                                                     random.seed = i)
      object <- Seurat::SetIdent(object, value = "clusters_perm")
       ##STutility::RegionNeighbours所产生的接邻标签列将以nbs_开头。为了防止后续被混淆，已有nbs_匹配列将被清除
      for(column_rm in grep(pattern = "nbs_", colnames(object[[]]), value = T)){
        object[[column_rm]] <- NULL
      }
      for(c in c_include){
        object <- STutility::RegionNeighbours(object, id = c, keep.within.id = T, verbose = TRUE)##在Ident中去选取id(id=c)的簇，根据此簇筛选出近邻
      }
      perm_nbs_df <- object[[]][, c("clusters_perm", grep(pattern = "nbs_", colnames(object[[]]), value = T))]
      perm_adj_mat <- OESingleCell::CreateAdjMatrix(nbs.df = perm_nbs_df,
                                                    cluster.include = c_include,
                                                    column.clusters.id = "clusters_perm")
      perm_adj_mat_list[[i]] <- perm_adj_mat##循环50次，length(perm_adj_mat_list) 为50
    }
    n_cols <- dim(perm_adj_mat_list[[1]])[1]
    perm_adj_mat_avg <- perm_adj_mat_sd <- matrix(0L, nrow = n_cols, ncol = n_cols)
    for(i in 1:n_cols){
      for(j in 1:n_cols){
        list_ij <- c()
        for(list_n in 1:length(perm_adj_mat_list)){
          list_ij <- c(list_ij, perm_adj_mat_list[[list_n]][i,j])
        }
        perm_adj_mat_avg[i,j] <- mean(list_ij)###50个矩阵中同一位置上的mean
        perm_adj_mat_sd[i,j] <- sd(list_ij)###50个矩阵中同一位置上的sd
      }
    }

    #' Check if random data is normally dist.
    perm_adj_mat_df <- data.frame(matrix(unlist(perm_adj_mat_list), nrow=length(perm_adj_mat_list), byrow=T))  # dims: n_permX(18*18) = 50x324
    perm_adj_mat_df$perm_n <- rownames(perm_adj_mat_df)
    perm_adj_mat_df_long <- reshape2::melt(perm_adj_mat_df, id.vars = c("perm_n"))

    #' OBSERVED VALUES: Run RegionNeighbours()##真实数据构建
    Seurat::Idents(object) <- object[[col_name]][,1]
    for(column_rm in grep(pattern = "nbs", colnames(object[[]]), value = T)){
      object[[column_rm]] <- NULL
    }
    for(c in c_include){
      object <- STutility::RegionNeighbours(object, id = c, keep.within.id = T, verbose = TRUE)
    }

    #' Create df with nbs info
    nbs_df <- object[[]][, c(col_name, grep(pattern = "nbs", colnames(object[[]]), value = T))]

    plot_coln <- grep(pattern = "nbs", colnames(object[[]]), value = T)

    sample_num <- length(Seurat::Images(object))
        for(i in plot_coln){
        object[[i]][,1] <- factor(object[[i]][,1],levels <- unique(object[[i]][,1]))
        cell=strsplit(i,split = "nbs_") %>% unlist %>% .[2]
        col_i <- c(col[cell],'#4D3D87','gray80')
        names(col_i) <- c(cell,i,NA)

        p_sp <- OESingleCell::SpatialPlot(object,group.by=i,
                              combine=TRUE, pt.size.factor=1,alpha = c(0.1, 0.5),
                              col=col_i)
         OESingleCell::save_ggplots(file.path(output_dir_heterotype,paste0(i,"_adjmap")),
                         plot = p_sp,
                         width = 5*sample_num,
                         height = 5,
                         dpi=100)

    }
    sample_num <- length(Seurat::Images(object))
        for(i in plot_coln){
        object[[i]][,1] <- factor(object[[i]][,1],levels <- unique(object[[i]][,1]))
        cell=strsplit(i,split = "nbs_") %>% unlist %>% .[2]
        col_i <- c(col[cell],'#4D3D87','gray80')
        names(col_i) <- c(cell,i,NA)

        p_sp <- OESingleCell::SpatialPlot(object,group.by=i,
                              combine=TRUE, pt.size.factor=1.2,alpha = c(0.1, 0.5),
                              col=col_i)
         OESingleCell::save_ggplots(file.path(output_dir_heterotype,paste0(i,"_adjmap")),
                         plot = p_sp,
                         width = 5*sample_num,
                         height = 5,
                         dpi=100)

    }
    # c_count <- nbs_df %>%
    #   dplyr::group_by(celltype) %>%
    #   dplyr::count()
    c_count <- data.frame(table(nbs_df[,col_name]))##此种细胞类型
    names(c_count) <- c("cluster", "n")
    rownames(c_count) <- c_count$cluster
    # colnames(c_count) <- c("cluster", "n")
    c_count$id <- paste0("C", c_count$cluster)

    nbs_adjmat <- OESingleCell::CreateAdjMatrix(nbs.df = nbs_df,
                                  cluster.include = c_include,
                                  column.clusters.id =col_name
                                  )##实际观测值的adjmat


    #' Normalize method 1
    ##如：A与B接邻，则A与B的接邻标准化方式为：（A,B接邻数）/(A的接邻总和+B的接邻总和)
    nbs_adjmat_norm <- nbs_adjmat
    for (i in seq(1:dim(nbs_adjmat_norm)[1])) {
      for (j in seq(1:dim(nbs_adjmat_norm)[1])) {
        e_sum <- sum(nbs_adjmat_norm[i,i], nbs_adjmat_norm[j,j])##i的所有的接邻点和j的所有的接邻点之和
        if (i!=j) {
          nbs_adjmat_norm[i,j] <- round(nbs_adjmat_norm[i,j]/e_sum, 4)*100##保留4位小数并*100
        }
      }
    }

    #' Normalize method 1
    ##如：A与B接邻，则A与B的接邻标准化方式为：（A,B接邻数）/(A的接邻总和+B的接邻总和)
    nbs_adjmat_norm <- nbs_adjmat
    for (i in seq(1:dim(nbs_adjmat_norm)[1])) {
      for (j in seq(1:dim(nbs_adjmat_norm)[1])) {
        e_sum <- sum(nbs_adjmat_norm[i,i], nbs_adjmat_norm[j,j])##i的所有的接邻点和j的所有的接邻点之和
        if (i!=j) {
          nbs_adjmat_norm[i,j] <- round(nbs_adjmat_norm[i,j]/e_sum, 4)*100##保留4位小数并*100
        }
      }
    }

    #' "Normalize"/Standardise method 2: Using permuted values（（obs-exp）/exp）这是文章方法 （推荐）
    nbs_adjmat_permscore <- round( ((nbs_adjmat - perm_adj_mat_avg) / perm_adj_mat_sd) , digits = 3)#保留三位有效小数
    diag(nbs_adjmat_permscore) <- diag(nbs_adjmat)##对角元素依然为总接邻数量
    rownames(nbs_adjmat) <- rownames(nbs_adjmat_norm) <- rownames(nbs_adjmat_permscore) <- c_include
    colnames(nbs_adjmat) <- colnames(nbs_adjmat_norm) <- colnames(nbs_adjmat_permscore) <- c_include
    rownames(perm_adj_mat_avg) <- rownames(perm_adj_mat_sd) <- colnames(perm_adj_mat_avg) <- colnames(perm_adj_mat_sd) <- c_include
    write.table(nbs_adjmat,file=file.path(output_dir_heterotype,'observed_nbs_adjmat.txt'),sep="\t",quote = F)
    write.table(nbs_adjmat_norm,file=file.path(output_dir_heterotype,'nbs_adjmat_norm_method1.txt'),sep="\t",quote = F)
    write.table(nbs_adjmat_permscore,file=file.path(output_dir_heterotype,'nbs_adjmat_permscore_method2.txt'),sep="\t",quote = F)

    nbs_adjmat_numeric <- as.numeric(nbs_adjmat_norm)[! (as.numeric(nbs_adjmat_norm) %in% diag(nbs_adjmat_norm))]
    p_adjmat <- pheatmap::pheatmap(nbs_adjmat_norm,
                  cluster_rows = F, cluster_cols = F,
                  breaks = seq(min(nbs_adjmat_numeric),max(nbs_adjmat_numeric),3),
                  color = colorRampPalette(c("#4d2664", "white", "#376431"))(length(seq(min(nbs_adjmat_numeric),max(nbs_adjmat_numeric),3))),
                  border_color = "white",
                  na_col = "grey80")
    OESingleCell::save_ggplots(file.path(output_dir_heterotype,'nbs_adjmat_norm'),
                            plot=p_adjmat,
                            width=5,
                            height=5,
                            dpi=200,bg='white')
        ##每个样本的实际连通度
    nbs_adjmat_list_subject <- OESingleCell::CreateAdjMatrixPerSubject(nbs.df = nbs_df,
                                                         se.metadata = object@meta.data,
                                                         column.subject.id = "sampleid_num",
                                                         cluster.include = c_include,
                                                         column.clusters.id = col_name
                                                         )
    ##针对每一张切片（样本）
    #' PER SUBJECT, EXPECTED VALUES
    n_perm <- 50
    perm_adj_mat_list_subject <- list()
    for(i in 1:n_perm){
      object <- OESingleCell::RandomiseClusteringIDs(se.object = object, column.cluster.id = col_name, column.sample.id = "sampleid_num", random.seed = i)
      object <- Seurat::SetIdent(object, value = "clusters_perm")
      for(column_rm in grep(pattern = "nbs_", colnames(object[[]]), value = T)){
        object[[column_rm]] <- NULL
      }
      for(c in c_include){
        object <- STutility::RegionNeighbours(object, id = c, keep.within.id = T, verbose = TRUE)
      }
      perm_nbs_df <- object[[]][, c("clusters_perm", grep(pattern = "nbs_", colnames(object[[]]), value = T))]
      perm_adj_mat_list_subject_id <- OESingleCell::CreateAdjMatrixPerSubject(nbs.df = perm_nbs_df,
                                                               se.metadata = object@meta.data,
                                                               cluster.include = c_include,
                                                               column.clusters.id = "clusters_perm",
                                                               column.subject.id = "sampleid_num")
      perm_adj_mat_list_subject[[i]] <- perm_adj_mat_list_subject_id
    }
    n_cols <- dim(perm_adj_mat_list_subject[[1]][[1]])[1]
    perm_adj_mat_avg_list_subject <- perm_adj_mat_sd_list_subject <- list()
    perm_adj_mat_avg_subject <- perm_adj_mat_sd_subject <- matrix(0L, nrow = n_cols, ncol = n_cols)


    for (subject_id in names(perm_adj_mat_list_subject[[1]])) {
      perm_adj_mat_list_subject_id <- list()
      for (i in 1:length(perm_adj_mat_list_subject)) {
        perm_adj_mat_list_subject_id[[i]] <- perm_adj_mat_list_subject[[i]][[subject_id]]
      }
      for (i in 1:n_cols) {
        for (j in 1:n_cols) {
          list_ij <- c()
          for (list_n in 1:length(perm_adj_mat_list_subject_id)) {
            list_ij <- c(list_ij, perm_adj_mat_list_subject_id[[list_n]][i,j])
          }
          perm_adj_mat_avg_subject[i,j] <- mean(list_ij)
          perm_adj_mat_sd_subject[i,j] <- sd(list_ij)
        }
        perm_adj_mat_avg_list_subject[[subject_id]] <- perm_adj_mat_avg_subject
        perm_adj_mat_sd_list_subject[[subject_id]] <- perm_adj_mat_sd_subject
      }
    }
    #' PER SUBJECT, (OBS - EXP_AVG)/EXP_SD
    nbs_adjmat_permscore_subject <- list()
    for (subject_id in names(nbs_adjmat_list_subject)) {
      nbs_adjmat_permscore_subject[[subject_id]] <- round( (nbs_adjmat_list_subject[[subject_id]] - perm_adj_mat_avg_list_subject[[subject_id]]) / perm_adj_mat_sd_list_subject[[subject_id]], digits = 3)
      diag(nbs_adjmat_permscore_subject[[subject_id]]) <- diag(nbs_adjmat_list_subject[[subject_id]])

      #' Set row/column names
      rownames(nbs_adjmat_permscore_subject[[subject_id]]) <- colnames(nbs_adjmat_permscore_subject[[subject_id]]) <- c_include
    }
    write.table(do.call(rbind,nbs_adjmat_permscore_subject),file=file.path(output_dir_heterotype,'nbs_adjmat_permscore_subject.txt'),sep="\t",quote=F)

    color_low2 <- "#62376e"
    nbs_adjmat_df <- nbs_adjmat_norm
    g <- igraph::graph.adjacency(nbs_adjmat_df, mode = "undirected", weighted = TRUE, diag = F)
    e_sum <- diag(nbs_adjmat_df)
    df <- data.frame(name = names(e_sum),
                     size = e_sum,
                     norm_size = e_sum/max(e_sum)*40,
                     group=names(e_sum),
                     color = col,
                     c_size = c_count[names(e_sum), 'n'],
                     stringsAsFactors = F)
    links <- data.frame(
      source=igraph::as_edgelist(g, names = TRUE)[,1],
      target=igraph::as_edgelist(g, names = TRUE)[,2])

    g_df_norm <- data.frame(links,
                            weight = igraph::E(g)$weight,
                            weight_minmax = OESingleCell::minmax_norm(abs(igraph::E(g)$weight)))

    g2 <- igraph::graph_from_data_frame(d = links, vertices = df, directed = F)

    g2 <- igraph::set_edge_attr(g2, "weight", value=igraph::E(g)$weight)
    igraph::E(g2)$width <- igraph::E(g2)$weight*1.5

    pdf(file.path(output_dir_heterotype,'nbs_adjmat_df_net.pdf'),w=15,h=15)
    plot(g2,
    #      layout=layout_in_circle,
         vertex.label.family = "Helvetica", vertex.label.color = "black", vertex.label.cex = 1.5,
         vertex.size = igraph::V(g2)$norm_size,
         vertex.color = grDevices::adjustcolor(igraph::V(g2)$color, alpha.f=.9),
         vertex.frame.color="white",
         edge.color = grDevices::adjustcolor("grey70", alpha.f = .5),
         edge.curved=0.1)
    dev.off()


    g3 <- g2
    igraph::E(g3)$width <- igraph::E(g2)$weight
    for(clusters_plot in names(e_sum)){
      e_pairs <- c()
      for(c in names(e_sum)[!names(e_sum)==clusters_plot]){
        e_pairs <- c(e_pairs, clusters_plot, c)
      }
      ecol <- rep(NA, igraph::ecount(g3))
      ecol[ igraph::get.edge.ids(g, vp = e_pairs)] <- "gray50"

      node_ids <- names( igraph::V(g3))
      node_sizes <- data.frame(node = node_ids,
                               size =  igraph::V(g3)$norm_size,
                               stringsAsFactors = F)
      LS <-  igraph::layout_as_star(g3, center = clusters_plot, order = order(node_sizes$size, decreasing = T))

    #   fname <- paste0("nbs_analysis.graph_adipocytes_", cluster_plot, ".", ANALYSIS)
    #   cairo_ps(filename = file.path(DIR_FIG, paste0(fname, ".eps")), width = 5.5, height = 5.5)
     pdf(file.path(output_dir_heterotype,paste0(clusters_plot,'_nbs_adjmat_sub_net.pdf')),w=15,h=15)
        plot(g3,
           layout=LS,
           vertex.label.family = "Helvetica",
           vertex.label.color = "black",
           vertex.label.cex = 1.5,
           vertex.size = igraph::V(g3)$norm_size*1.25,
           vertex.color = grDevices::adjustcolor(igraph::V(g3)$color, alpha.f=1),
           vertex.frame.color="white",
           edge.color = grDevices::adjustcolor(ecol, alpha.f = .5),
           edge.curved=0)
       dev.off()
    }

    #' PLOT HEATMAP
    hm_df <- nbs_adjmat_permscore##校正后
    diag(hm_df) <- NA
    hm_df[lower.tri(hm_df)] <- NA
    hm_df <- hm_df[,rev(colnames(hm_df))]

    breaksList <- seq(-6, 6, by = .5)
    p1 <- pheatmap::pheatmap(hm_df,
                  cluster_rows = F, cluster_cols = F,
                  breaks = breaksList,
                  color = colorRampPalette(c("#4d2664", "white", "#376431"))(length(breaksList)-1),
                  border_color = "white",
                  na_col = "grey80")
    OESingleCell::save_ggplots(file.path(output_dir_heterotype,'nbs_adjmat_permscore'),
                            plot=p1,
                            width=5,
                            height=5,
                            dpi=200,bg='white')
    thesample <- names(rename_sample)
    names(thesample) <- as.character(1:length(rename_sample))

     ##按样本切分后
    for (subject_id in names(nbs_adjmat_permscore_subject)) {
      hm_df <- nbs_adjmat_permscore_subject[[subject_id]]

      diag(hm_df) <- NA
      hm_df[lower.tri(hm_df)] <- NA
      hm_df[is.infinite(hm_df)] <- NA
      hm_df <- hm_df[,rev(colnames(hm_df))]
      breaksList <- seq(-6, 6, by = .5)
       psub <- pheatmap::pheatmap(hm_df,
               cluster_rows = F,
               cluster_cols = F,
               breaks = breaksList,
               color = colorRampPalette(c("#4d2664", "white", "#376431"))(length(breaksList)-1),
               border_color = "white",
               na_col = "grey80",
               main  = thesample[subject_id])
       OESingleCell::save_ggplots(file.path(output_dir_heterotype,paste0(thesample[subject_id],'nbs_adjmat_permscore')),
                        plot=psub,
                        width=5,
                        height=5,
                        dpi=200,bg='white')
    }




  }







}


