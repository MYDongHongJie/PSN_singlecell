#.libPaths("/PERSONALBIO/work/singlecell/s00/software/R.library/library")
#.libPaths(c("/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/sc/lib/R/library","/PERSONALBIO/work/singlecell/s00/software/R.library/library"))
options(bitmapType='cairo')

source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/color/color.R")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/plot.r")

Theme <- function(){
    thememe <- theme(axis.title = element_text(face = 'bold',
                                      size = "16",color = "black"),
            # legend.position = 'right',
            axis.text.x = element_text(color = "black",face = 'bold',angle = 90,
                                       size = 10,  hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(face = 'bold',size =10,color = "black"),
            legend.text = element_text(face = 'bold',color = "black",size = 10),
            legend.title = element_text(face = 'bold',color = "black",size = 10),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position="right",
            strip.text.x = element_text( face="bold",size = 10),
            strip.text.y = element_text( face="bold",size = 10),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            strip.background = element_rect(colour = "white", fill = "grey"),
            plot.title = element_text(face = "bold",color = "black",lineheight=.8,  hjust=0.5, size =5))
    return(thememe)
}

Topn_heatmap_marker<-function(immune.combined,collapseby = opt$collapseby,Topn = opt$Topn,seurat_exp_cluster_dir=seurat_exp_cluster_dir,markers=markers,colors =colors ){
        suppressPackageStartupMessages(library(ComplexHeatmap))
        cluster_topn_markers=markers %>% group_by(cluster)%>% top_n(n = Topn, wt = avg_log2FC)
        out_dir = seurat_exp_cluster_dir
        write.table(cluster_topn_markers, file.path(out_dir,glue::glue("cluster_top{Topn}_markers_avg_heatmap.xls")),quote = F, row.names = F, sep = "\t")

        markers2vis = cluster_topn_markers$gene
        count = as.matrix(GetAssayData(immune.combined, slot = "data")[markers2vis,])
        metadata = immune.combined@meta.data
        metadata$id = rownames(metadata)
        collapsed_count = vector()
        if ( !collapseby %in% colnames(metadata) ){
            stop("NO specified column found!")
        }
        collapsed_group = metadata %>% group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))
        for ( cells in collapsed_group$cellid ){
            samplex = unlist(strsplit(cells, ",", perl =T))
            collapsed_count= cbind(collapsed_count,rowMeans( count[,samplex,drop=F] ))
        }
        collapsed_count = as.matrix( collapsed_count )
        collapsed_group = as.data.frame(collapsed_group)
        colnames(collapsed_count) = as.matrix(collapsed_group[,1])
        data = tibble::rownames_to_column(as.data.frame(collapsed_count),"GeneID")
        #out_dir = file.path(seurat_exp_cluster_dir,glue::glue('Top{Topn}_Marker_heatmap'))
        out_dir = seurat_exp_cluster_dir
        if (!file.exists(out_dir)){
                dir.create(out_dir,recursive = TRUE)
        }
        #write.table(data, file.path(out_dir,glue::glue("Top{Topn}_avgExpression_marker_heatmap_data.xls")),quote = F, row.names = F, sep = "\t")
        if (dim(collapsed_count)[2]>2) {
            df <- t(scale(t(collapsed_count)))
                scale_data = tibble::rownames_to_column(as.data.frame(df),"GeneID")
                #write.table(scale_data, file.path(out_dir,glue::glue("Top{Topn}_avgExpression_marker_heatmap_scaledata.xls")),quote = F, row.names = F, sep = "\t")
        }else{
                df <- collapsed_count
        }
        #增加Annobar
        if(is.factor(metadata[,collapseby]) ){
            Anno_df = levels(metadata[,collapseby])
        }else{
            Anno_df = as.vector(colnames(collapsed_count))
        }
        Anno_df = factor(Anno_df,levels=Anno_df)
        col <- setNames(colors[1:length(Anno_df)], Anno_df)
                ncol = ceiling(length(unique(metadata[[collapseby]]))/10)
        Ha = HeatmapAnnotation(CellType=Anno_df,col=list(CellType =col ),simple_anno_size = unit(0.3, "cm"),show_annotation_name=F,annotation_legend_param =list(ncol=ncol))
        palette <-colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))(299)

        plot2 =Ha %v% Heatmap(df,
                row_names_gp = gpar(fontsize = 4),
                column_names_gp = gpar(fontsize = 8),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names= FALSE,
                col = palette,
                heatmap_legend_param = list(title = "Expression-level",title_position ="leftcenter-rot"),
                column_names_rot = 45,
                name=" ")
        ggheatmap <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot2, padding = grid::unit(c(2, 3, 2, 3), "mm"))))
        ggsave(plot=ggheatmap,filename=file.path(out_dir,glue::glue("cluster_top{Topn}_markers_avg_heatmap.pdf")),width=3.8+ncol*0.2,height=4)
        ggsave(plot=ggheatmap,filename=file.path(out_dir,glue::glue("cluster_top{Topn}_markers_avg_heatmap.png")),width=3.8+ncol*0.2,height=4,dpi=1000)
}


Seurat.Plot <- function(immune.combined,colors=colors,seurat_exp_cluster_dir=seurat_exp_cluster_dir,markers=markers,defaultassay="RNA",topn=30){
		cluster_overviwe_dir = file.path(seurat_exp_cluster_dir,"cluster_overview")
    if(!dir.exists(cluster_overviwe_dir)){
        dir.create(cluster_overviwe_dir,recursive = T)
    }
    tsne_dir = file.path(seurat_exp_cluster_dir,"tsne")
    if(!dir.exists(tsne_dir)){
        dir.create(tsne_dir,recursive = T)
    }
    umap_dir = file.path(seurat_exp_cluster_dir,"umap")
    if(!dir.exists(umap_dir)){
        dir.create(umap_dir,recursive = T)
    }
       
    #展示每个cluster的marker基因
    DefaultAssay(immune.combined) <- defaultassay
    nsample <- length(unique(immune.combined$sample))
    ngroup <- length(unique(immune.combined$group))
    ncluster <- length(unique(immune.combined$seurat_clusters))
    if("celltype" %in% colnames(immune.combined@meta.data)){
        ncluster <- length(unique(immune.combined$celltype))
        collapseby <- "celltype"
        #按样本/组展示细胞占比。柱状图。横坐标是样本/组，纵坐标是细胞占比
        temp <- data.frame(table(immune.combined$celltype,immune.combined$sample))
    }else{
        collapseby <- "seurat_clusters"
        ncluster <- length(unique(immune.combined$seurat_clusters))
        temp <- data.frame(table(immune.combined$seurat_clusters,immune.combined$sample))
    }

    #按样本/组展示细胞占比。柱状图。横坐标是样本/组，纵坐标是细胞占比
	colnames(temp) <- c("Cluster","Sample","number")
	p <- ggplot(temp,aes(x=Sample,fill=Cluster,y=number))+
	  geom_bar(stat = "identity",position ='fill',colour="black")+
	  labs(x='',y='')+
	  guides(fill = guide_legend(title = '',reverse=TRUE))+
	  scale_fill_manual(values = colors)+
	  theme_bw()+ Theme()
	ggsave(p,file=paste(cluster_overviwe_dir,"cluster_sample_cellcounts.png",sep="/"),width = 8,height = 7)
	ggsave(p,file=paste(cluster_overviwe_dir,"cluster_sample_cellcounts.pdf",sep="/"),width = 8,height = 7)
    if("celltype" %in% colnames(immune.combined@meta.data)){
        temp <- data.frame(table(immune.combined$celltype,immune.combined$group))
    }else{
        temp <- data.frame(table(immune.combined$seurat_clusters,immune.combined$group))
    }
    colnames(temp) <- c("Cluster","Group","number")
    p <- ggplot(temp,aes(x=Group,fill=Cluster,y=number))+
      geom_bar(stat = "identity",position ='fill',colour="black")+
      labs(x='',y='')+
      guides(fill = guide_legend(title = '',reverse=TRUE))+
      scale_fill_manual(values = colors)+
      theme_bw()+ Theme()
    ggsave(p,file=paste(cluster_overviwe_dir,"cluster_group_cellcounts.png",sep="/"),width = 8,height = 7)
    ggsave(p,file=paste(cluster_overviwe_dir,"cluster_group_cellcounts.pdf",sep="/"),width = 8,height = 7)
    

    #分样本展示
    p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample",cols=colors)
    ggsave(p1,file=paste(umap_dir,"cluster_umap.sample.pdf",sep="/"),width = 4.8,height = 4,limitsize = FALSE)
    ggsave(p1,file=paste(umap_dir,"cluster_umap.sample.png",sep="/"),width = 4.8,height = 4,limitsize = FALSE)
    p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "sample",cols=colors)
    ggsave(p1,file=paste(tsne_dir,"cluster_tsne.sample.pdf",sep="/"),width = 4.8,height = 4,limitsize = FALSE)
    ggsave(p1,file=paste(tsne_dir,"cluster_tsne.sample.png",sep="/"),width = 4.8,height = 4,limitsize = FALSE)

    #分组展示
    p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "group",cols=colors)
    ggsave(p1,file=paste(umap_dir,"cluster_umap.group.pdf",sep="/"),width = 4.8,height = 4,limitsize = FALSE)
    ggsave(p1,file=paste(umap_dir,"cluster_umap.group.png",sep="/"),width = 4.8,height = 4,limitsize = FALSE)
    p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "group",cols=colors)
    ggsave(p1,file=paste(tsne_dir,"cluster_tsne.group.pdf",sep="/"),width = 4.8,height = 4,limitsize = FALSE)
    ggsave(p1,file=paste(tsne_dir,"cluster_tsne.group.png",sep="/"),width = 4.8,height = 4,limitsize = FALSE)

    #分cluster展示
    p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,cols=colors) +NoLegend()
    ggsave(paste(umap_dir,"cluster_umap.label.pdf",sep="/"),width = 4,height = 4,limitsize = FALSE)
    ggsave(paste(umap_dir,"cluster_umap.label.png",sep="/"),width = 4,height = 4,limitsize = FALSE)
    p2 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE,cols=colors) +NoLegend()
    ggsave(paste(tsne_dir,"cluster_tsne.label.pdf",sep="/"),width = 4,height = 4,limitsize = FALSE)
    ggsave(paste(tsne_dir,"cluster_tsne.label.png",sep="/"),width = 4,height = 4,limitsize = FALSE)
    p2 <- DimPlot(immune.combined, reduction = "umap", label = FALSE,cols=colors)
    ggsave(paste(umap_dir,"cluster_umap.pdf",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6,limitsize = FALSE)
    ggsave(paste(umap_dir,"cluster_umap.png",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6,limitsize = FALSE)
    p2 <- DimPlot(immune.combined, reduction = "tsne", label = FALSE,cols=colors)
    ggsave(paste(tsne_dir,"cluster_tsne.pdf",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6,limitsize = FALSE)
    ggsave(paste(tsne_dir,"cluster_tsne.png",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6,limitsize = FALSE)


    if(nsample==1){swidth <- 5}else{swidth <- 9}
    sheight <- ceiling(nsample/2)*4
    if(ngroup==1){gwidth <- 5}else{gwidth <- 9}
    gheight <- ceiling(ngroup/2)*4
    #分样本展示 按样本分页;分组展示，按组分页
    DimPlot(immune.combined, reduction = "umap", split.by = "sample",ncol=2,cols=colors)
    ggsave(paste(umap_dir,"cluster_umap_splitsample.pdf",sep="/"),width = swidth,height = sheight,limitsize = FALSE)
    ggsave(paste(umap_dir,"cluster_umap_splitsample.png",sep="/"),width = swidth,height = sheight,limitsize = FALSE)
    DimPlot(immune.combined, reduction = "tsne", split.by = "sample",ncol=2,cols=colors)
    ggsave(paste(tsne_dir,"cluster_tsne_splitsample.pdf",sep="/"),width = swidth,height = sheight,limitsize = FALSE)
    ggsave(paste(tsne_dir,"cluster_tsne_splitsample.png",sep="/"),width = swidth,height = sheight,limitsize = FALSE)

    DimPlot(immune.combined, reduction = "umap", split.by = "group",ncol=2,cols=colors)
    ggsave(paste(umap_dir,"cluster_umap_splitgroup.pdf",sep="/"),width = gwidth,height = gheight,limitsize = FALSE)
    ggsave(paste(umap_dir,"cluster_umap_splitgroup.png",sep="/"),width = gwidth,height = gheight,limitsize = FALSE)
    DimPlot(immune.combined, reduction = "tsne", split.by = "group",ncol=2,cols=colors)
    ggsave(paste(tsne_dir,"cluster_tsne_splitgroup.pdf",sep="/"),width = gwidth,height = gheight,limitsize = FALSE)
    ggsave(paste(tsne_dir,"cluster_tsne_splitgroup.png",sep="/"),width = gwidth,height = gheight,limitsize = FALSE)

    #marke基因展示
    Marker_dir = file.path(seurat_exp_cluster_dir,"../Marker")
    if(!file.exists(Marker_dir)){dir.create(Marker_dir)}
   
    all_top10_markers=markers %>% filter(avg_log2FC != Inf) %>%group_by(cluster)  %>%  top_n(n = 1, wt = avg_log2FC) %>% as.data.frame() %>% dplyr::distinct(.,gene,.keep_all = T) %>% head(n = 10)
    height= ceiling(length(all_top10_markers$gene)/5)*4

    cluster_top10_markers=markers %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)
    tmp <- subset(immune.combined,downsample=1000)
		tmp = ScaleData(tmp,features= unique(cluster_top10_markers$gene))
    DoHeatmap(tmp, features = unique(cluster_top10_markers$gene),group.colors = colors) + NoLegend() +theme(axis.text.y = element_text(size = 4)) 

    ggsave(paste(Marker_dir,"cluster_top10_markers_heatmap.pdf",sep="/"),width = 14,height = 8)
    ggsave(paste(Marker_dir,"cluster_top10_markers_heatmap.png",sep="/"),width = 14,height = 8)
    Topn_heatmap_marker(immune.combined=immune.combined,markers=markers,collapseby=collapseby,Topn=topn,seurat_exp_cluster_dir=Marker_dir,colors=colors)
    #png(paste(seurat_exp_cluster_dir,"cluster_top10_markers_heatmap.png",sep="/"),width =1200,height = 800)
    #DoHeatmap(tmp, features = unique(cluster_top10_markers$gene)) + NoLegend() +theme(axis.text.y = element_text(size = 4))
    #dev.off()
}
