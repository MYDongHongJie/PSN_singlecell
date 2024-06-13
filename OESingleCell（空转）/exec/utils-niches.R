sub_niches <- subparsers$add_parser("niches", help = "do niches in scRNA or Spatial.")
sub_niches$add_argument("--gene_int",type="character",default='/public/dev_scRNA/yfang/niche/niches_human_omnipath_interactions.RData',
                         help = "the interaction files")
sub_niches$add_argument("--method",type="character",default = 'NeighborhoodToCell',
                         help = "the method for interaction.It can be CellToSystem ,SystemToCell,CellToCellSpatial,CellToNeighborhood,NeighborhoodToCell.[default: %(default)s]")
sub_niches$add_argument("--latent_time_file",type="character",default = NULL,
                         help = "the latent_time_file .[default: %(default)s]")
sub_niches$add_argument("--plotmarker",type="logical",default = TRUE,
                         help = "if to plot marker .[default: %(default)s]")
sub_niches$add_argument("--reduction",type="character",default = 'SCT_umap',
                         help = "the reduction for plotmarker  .[default: %(default)s]")

args <- commandArgs(TRUE)
if ( "niches"  %in% args ){
  opt<-intial_setting()
  if (opt$sub_name == "niches") {
    library(Seurat)
    methodlist <- c("CellToCell",##for single cell
    "CellToSystem" ,#for single cell,and sending cells are centre cells
    "SystemToCell",
    "CellToCellSpatial" ,
    "CellToNeighborhood",
    "NeighborhoodToCell")
    method <- methodlist %in% opt$method
    object <- readRDS(opt$input)
    objectcolnames <- gsub("-","_",colnames(object))
    object <- SeuratObject::RenameCells(object, new.names = objectcolnames)

    Seurat::Idents(object) <- object$clusters

    objectlist <- OESingleCell::SplitObject(object,'sampleid')

    if(opt$image==TRUE){
    objectlist <- lapply(objectlist,function(data){
        sample=unique(data$sampleid)
        data@meta.data$x=data@images[[sample]]@coordinates$row
        data@meta.data$y=data@images[[sample]]@coordinates$col
        data <- SeuratWrappers::RunALRA(data)
        return(data)

    })}

    intersection=load(opt$gene_int)
    get_niche <- function(object,LRdata,method){
      NICHES_output <- NICHES::RunNICHES(object = object,
                             LR.database = 'custom',
                             custom_LR_database=LRdata,
                             species = "human",
                             assay = "alra",
                             position.x = 'x',
                             position.y = 'y',
                             rad.set = 2, # Geometry dependent
                             min.cells.per.ident = 0,
                             min.cells.per.gene = NULL,
                             #meta.data.to.map = 'clusters',
                             CellToCell =method[1],CellToSystem =method[2],SystemToCell = method[3],
                             CellToCellSpatial = method[4],CellToNeighborhood =method[5],NeighborhoodToCell =method[6])
      #对RUNniche结果后的L_R互作对聚类
      niche <- NICHES_output[[methodlist[method]]]
      Seurat::Idents(niche) <- niche[['ReceivingType']]
      niche <- Seurat::ScaleData(niche)
      niche <- Seurat::FindVariableFeatures(niche,selection.method = "disp")
      niche <- Seurat::RunPCA(niche)

      niche <- Seurat::RunUMAP(niche,dims = 1:10)
      return(niche)

  }


    niche_list <- lapply(objectlist,function(x){
        obj <- x
        niche_output <- get_niche(object=obj,LRdata=niches_intersect,method=method)
        return(niche_output)
    })

    for(i in niche_list){
      #对重聚类后的结果仍以原聚簇信息（颜色）显示
        plot <-  Seurat::DimPlot(i,reduction = 'umap',pt.size = 0.5,shuffle = T, label = T) +ggplot2::ggtitle('Cellular Microenvironment')
        sample <- i[['orig.ident']][1,]
          OESingleCell::save_ggplots(
                            file.path(output_dir,paste(sample,'_dimplot')) ,
                            plot = plot,
                            dpi = 300 ,
                            limitsize = F,
                            width = 7,
                            height = 7)


    }
    output_dir_marker <- paste0(output_dir,"/marker")
    if (! file.exists( output_dir_marker)) { dir.create(output_dir_marker, recursive = T) }
    output_dir_marker<-normalizePath(output_dir_marker)
    top_marker_list <- lapply(niche_list,function(x){
        global_DEGs <- Seurat::FindAllMarkers(
                object = x,
                only.pos = T,
                test.use = 'wilcox',
                logfc.threshold = 0,
                min.pct = 0.25
            )

        global_DEGs <- global_DEGs  %>%
            dplyr::mutate(gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>%
            dplyr::select(gene, dplyr::everything())

        topn_markers <- global_DEGs %>%
                        dplyr::group_by(cluster) %>%
                        dplyr::filter(avg_log2FC >= 1 &
                            p_val < 0.05 &
                            pct.1 > 0.1 &
                            pct.2 < 2 &
                            gene_diff > 2) %>%
                        dplyr::arrange(p_val, dplyr::desc(avg_log2FC), dplyr::desc(gene_diff)) %>%
                        dplyr::filter(gene_diff > 2) %>%
                        dplyr::top_n(10, gene_diff)
         readr::write_tsv(
                global_DEGs,
                file = glue::glue("{output_dir_marker}/{x[['orig.ident']][1,1]}_all_markers_for_each_cluster_anno.xls"),
            )
         readr::write_tsv(
                topn_markers,
                file = glue::glue("{output_dir_marker}/{x[['orig.ident']][1,1]}_top10_markers_for_each_cluster_anno.xls"),)
        x <- x %>% Seurat::ScaleData(assay = "SCT", features = unique(topn_markers$gene))
        heatmap <- Seurat::DoHeatmap(x,features = unique(topn_markers$gene))+ ggplot2::scale_fill_gradientn(colors = c("grey","white", "blue"))
         l=length(unique(x$ReceivingType))
        OESingleCell::save_ggplots(
                            file.path(output_dir_marker,paste(x[['orig.ident']][1,1],'_top10_markers_heatmap')) ,
                            plot =  heatmap,
                            dpi = 300 ,
                            limitsize = F,
                            width = l,
                            height = 1.5*l)
        return(topn_markers)


    })



    niche_list2 <- list()
    for (x in 1:length(niche_list)){
        niches.data <- Seurat::GetAssayData(object =  niche_list[[x]][[methodlist[method]]], slot = 'data')
        colnames(niches.data) <- niche_list[[x]][['ReceivingCell']]$ReceivingCell
        sampcol <- paste0(niche_list[[x]][['orig.ident']][1,1],'_',methodlist[method])
        print(sampcol)# <- strsplit(sampcol,split="—")[[1]][2]
        print(objectlist[[x]][['sampleid']][1,1])
        objectlist[[x]][[methodlist[method]]] <- Seurat:: CreateAssayObject(data = niches.data )
        DefaultAssay(objectlist[[x]]) <- methodlist[method]
        objectlist[[x]] <- Seurat::ScaleData(objectlist[[x]])
        niches.data <- Seurat::ScaleData(niches.data)
        niches.data <- Seurat::CreateSeuratObject(counts = niches.data)
        niche_list2[[objectlist[[x]][['sampleid']][1,1]]] <- niches.data

    }
    mergedata <- merge( x = niche_list2[[1]], y = niche_list2[2:length(x = niche_list2)] )
    mergedata <- Seurat::GetAssayData(mergedata)

    object[[methodlist[method]]] <- CreateAssayObject(data = mergedata )

    DefaultAssay(object) <- methodlist[method]

    saveRDS(object,file=file.path(output_dir,'object_with_nightcell.rds'))

    for (i in 1:length(objectlist)){
        sample <- objectlist[[i]][['sampleid']][1,1]
        allimages <- Seurat::Images(objectlist[[i]])
        othersample <- allimages[!c(allimages %in% sample)]
        for (j in othersample){
            objectlist[[i]][[j]] <- NULL
        }
        saveRDS(objectlist[[i]],file=file.path(output_dir,paste0(sample,'_object_with_interaction.rds')))
    }
#调用marker作图
    if(opt$plotmarker==TRUE){
      for (i in 1:length(objectlist)){
        sample <- objectlist[[i]][['sampleid']][1,1]
        rdsfile <- file.path(output_dir,paste0(sample,'_object_with_interaction.rds'))
        markerfile <- file.path(output_dir_marker,
                                paste0("Neighborhood—",sample,'_top10_markers_for_each_cluster_anno.xls'))
        plot_output <- file.path(output_dir_marker,paste0(sample,'_markerplot'))

        if (! file.exists(plot_output)) {  dir.create(plot_output, recursive = T) }
        plot_output<-normalizePath(plot_output)

        plotassay <- methodlist[method]
        print('plot markergenes')
        cmd <- glue::glue('scVis -i {rdsfile} -o {plot_output} --assay {plotassay} --image TRUE markervis -l {markerfile} ',
                            '-n 10 -c avg_log2FC -m vlnplot,featureplot -s 0 --reduct {opt$reduction} -p 1.2 -a 1 --plotlegend none')
        futile.logger::flog.info(glue::glue("running:{cmd}"))
        system(cmd)



      }

    }

    smooth_gene_exp <- function(data = data, pseudotime = pseudotime, span = 0.75) {
        smooth_data <- data
        for (gene in 1:nrow(data)) {
            gene_exp <- t(data[gene, ])
            smooth <- loess(formula = gene_exp ~ pseudotime, span = span)
            smooth_data[gene, ] <- predict(smooth, newdata = pseudotime)
        }
        return(smooth_data)
    }
    output_dir_latenttime <- paste0(output_dir,'/cor_with_latenttime')
    if (! file.exists(output_dir_latenttime)) { dir.create(output_dir_latenttime, recursive = T) }
    output_dir_latenttime<-normalizePath(output_dir_latenttime)
    if(!is.null(opt$latent_time_file)){
        latent_time <- read.table(opt$latent_time_file,sep="\t",quote="",header=1,row.names=1)
        rownames(latent_time ) <- gsub("-","_",rownames(latent_time ))
        for(i in 1:length(objectlist)){
            nichedata <- data.frame(Seurat::GetAssayData(objectlist[[i]],slot = "data", assay =  methodlist[method]))
            nichedata <- nichedata[rowSums(nichedata)>0,]
            interbarcodes <- intersect(colnames(nichedata),rownames(latent_time))
            latent_time_sub  <- latent_time[interbarcodes,]
            latent_time_sub <- latent_time_sub[order(latent_time_sub$latent_time),]
            nichedata <- nichedata[,rownames(latent_time_sub)]
            clusterstype <- data.frame(objectlist[[i]]$clusters)[rownames(latent_time_sub),]
            names(clusterstype ) <- rownames(latent_time_sub)
            L_R_cor <- sapply(1:dim(nichedata)[1],function(x){
            data <- as.numeric(nichedata[x,])
            result <- cor.test(data,latent_time_sub$latent_time,method='pearson')
            result <- c(p=result$p.value,cor=result$estimate)
            })
            nichedata_filter=nichedata[L_R_cor[1,]<0.05&abs(L_R_cor[2,])>0.4,]
            L_R_cor_filter <- L_R_cor[,L_R_cor[1,]<0.05&abs(L_R_cor[2,])>0.4]
            L_R_cor_filter_output <- data.frame(t(L_R_cor_filter))
            rownames(L_R_cor_filter_output) <- rownames(nichedata_filter)
            write.table(L_R_cor_filter_output,
                        file=file.path(output_dir_latenttime,paste0(objectlist[[i]][['sampleid']][1,1],'_cor_wiht_latenttime.txt')),
                        sep="\t",quote=F)
            nichedata_heatmap <- smooth_gene_exp(nichedata_filter, pseudotime = latent_time_sub$latent_time, span = 0.4) %>%
                            t() %>%
                            scale() %>%
                            t() %>%
                            as.matrix()
            min_time <- min(latent_time_sub$latent_time)
            max_tim <- max(latent_time_sub$latent_time)
            med <- median(latent_time_sub$latent_time)
            col_dis='#899923'
            col_type <- OESingleCell::SelectColors(levels(clusterstype))
            ha =  ComplexHeatmap::HeatmapAnnotation(
                Latent_time= latent_time_sub$latent_time,
                ReceivingType = unname(clusterstype[rownames(latent_time_sub)]),
                col = list(Latent_time = circlize::colorRamp2(c(min_time,med,max_tim),c("white","#FF7A92",'#F0242E')),
                           ReceivingType = col_type[unname(clusterstype[rownames(latent_time_sub)])])

            )
            ha_row <- ComplexHeatmap::rowAnnotation(cor = L_R_cor_filter[2,])


            heatmap1 <- ComplexHeatmap::Heatmap(
                    nichedata_heatmap,
                    top_annotation = ha,
                    #left_annotation = ha_row,
                    cluster_rows = T,
                    cluster_columns = FALSE,
                    show_column_names = F,
                    show_row_names = T,
                    #row_split = row_ann$group,
                    row_title_rot = 0,
                    row_names_gp = grid::gpar(fontsize = 6)
                )%>% ggplotify::as.ggplot()
            file_heatmap <- file.path(output_dir_latenttime,paste0(objectlist[[i]][['sampleid']][1,1],'_cor_wiht_latenttimef'))
            h=ifelse(dim(nichedata_heatmap)[1]*0.09>50,50,dim(nichedata_heatmap)[1]*0.09)
            OESingleCell::save_ggplots(
                        file_heatmap,
                        plot = heatmap1,
                        width = length(col_type),
                        height = h,
                        dpi = 300
                    )
        }

    }















  }

}
