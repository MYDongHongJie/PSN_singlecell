
run_cluster<-function(immune.combined,seurat_exp_cluster_dir,type,idents,colors,avg_log2FC='0.25'){
    if(!dir.exists(seurat_exp_cluster_dir)){
        dir.create(seurat_exp_cluster_dir,recursive = T)
    }
    DefaultAssay(immune.combined) <- "RNA"
    Idents(immune.combined)<-idents
    #immune.combined <- NormalizeData(immune.combined)
    #immune.combined <- ScaleData(immune.combined, feature=rownames(immune.combined),verbose = FALSE)

    markers <- FindAllMarkers(immune.combined, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10)
		
    all_top10_markers=markers %>%  top_n(n = 50, wt = avg_log2FC) %>%dplyr::distinct(.,gene,.keep_all = T) %>% top_n(n = 10, wt = avg_log2FC)
		
    height= ceiling(length(all_top10_markers$gene)/5)*3
    for( clust_num in  unique(Idents(immune.combined))){
        cluster_dir=file.path(seurat_exp_cluster_dir,"Each_celltype_marker",paste("cluster",clust_num,sep="_"))
        if(!file.exists(cluster_dir)){
        dir.create(cluster_dir,recursive = T)
        }
				enrichment_dir = file.path(seurat_exp_cluster_dir,"../06_Enrichment",paste("cluster",clust_num,sep="_"))
				index=paste0("c",clust_num)
        upcluster_dir_enrich=paste(enrichment_dir,"up",sep="/")
        downcluster_dir_enrich=paste(enrichment_dir,"down",sep="/")
        allcluster_dir_enrich=paste(enrichment_dir,"all",sep="/")
        if(!file.exists(upcluster_dir_enrich)){dir.create(upcluster_dir_enrich,recursive =TRUE)}
        if(!file.exists(downcluster_dir_enrich)){dir.create(downcluster_dir_enrich,recursive =TRUE)}
        if(!file.exists(allcluster_dir_enrich)){dir.create(allcluster_dir_enrich,recursive =TRUE)}

        cluster_markers=subset(markers,cluster==clust_num)
        rownames(cluster_markers)<-cluster_markers$gene
        if(nrow(cluster_markers)>1){
          #genelist=cluster_markers$gene
            up =subset(cluster_markers,p_val < 0.05 & avg_log2FC > 0.25)
            down =subset(cluster_markers,p_val < 0.05 & avg_log2FC < -0.25)
            all  =subset(cluster_markers,p_val < 0.05 & abs(avg_log2FC) > 0.25)
            upgenelist=up$gene
            downgenelist=down$gene
            allgenelist = all$gene
            try(enrichment(species=type,outDir=upcluster_dir_enrich,geneList=upgenelist))
            try(enrichment(species=type,outDir=downcluster_dir_enrich,geneList=downgenelist))
            try(enrichment(species=type,outDir=allcluster_dir_enrich,geneList=allgenelist))
            write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
            top10_markers=cluster_markers %>%  top_n(n = 10, wt = avg_log2FC)
            VlnPlot(immune.combined, features = top10_markers$gene,pt.size = 0 ,ncol=5,cols=colors)
            ggsave(paste(cluster_dir,paste0(index,"_top10_vlnplot.pdf"),sep="/"),width =20,height = height,bg='#ffffff')
            ggsave(paste(cluster_dir,paste0(index,"_top10_vlnplot.png"),sep="/"),width =20,height = height,bg='#ffffff')

            FeaturePlot(immune.combined, features = top10_markers$gene, min.cutoff = "q9",ncol=5,order=T,cols=c("lightgrey", "red"))
            ggsave(paste(cluster_dir,paste0(index,"_top10_umap.pdf"),sep="/"),width = 20,height = height,bg='#ffffff')
            ggsave(paste(cluster_dir,paste0(index,"_top10_umap.png"),sep="/"),width = 20,height = height,bg='#ffffff')
            #DotPlot(immune.combined, features = rev(unique(top10_markers$gene)))
            #ggsave(paste(cluster_dir,"top10_dotplot.pdf",sep="/"),width = 10,height = height,bg='#ffffff')
            #ggsave(paste(cluster_dir,"top10_dotplot.png",sep="/"),width = 10,height = height,bg='#ffffff')
            gc(TRUE)
        }
    }
  write.table(markers,paste(seurat_exp_cluster_dir,"allmarkers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
  return(markers)
}

Topmarker_plot <-function(immune.combined,markers,seurat_exp_cluster_dir,type,idents,colors,avg_log2FC='0.25'){
	 #marke基因展示
	  marker_number<-table(markers$cluster)%>% reshape2::melt()
    marker_number$Cluster <- as.factor(marker_number$Var1)
		number_of_bar <- nrow(marker_number)
		marker_number$id <- seq(1, number_of_bar)
		angle <- 90 - 360 * (marker_number$id-0.5) /number_of_bar

		marker_number$hjust <- ifelse( angle < -90, 1, 0)
		marker_number$angle <- ifelse(angle < -90, angle+180, angle)

		p <- ggplot(marker_number, aes(x=as.factor(id), y=value, fill=Cluster)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
		geom_bar(stat="identity", alpha=0.5)+scale_fill_manual(values=colors) +
		ylim(-(min(marker_number$value)-10),NA) +
		theme_minimal() +
		theme(
			#legend.position = "none",
			axis.text = element_blank(),
			axis.title = element_blank(),
			panel.grid = element_blank(),
			plot.margin = unit(rep(0.1,4), "cm") 
		)  +
		coord_polar() + 
		geom_text(data=marker_number, aes(x=id, y=value+10, label=as.vector(value), hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3, angle= marker_number$angle, inherit.aes = FALSE ) 

    ggsave(plot = p,paste(seurat_exp_cluster_dir,"marker_number.pdf",sep="/"),width = 8,height = 8,bg='white')
    ggsave(plot = p,paste(seurat_exp_cluster_dir,"marker_number.png",sep="/"),width =8,height = 8,dpi=1000,bg='white')



    #top基因展示
    all_top10_markers=markers %>% filter(avg_log2FC != Inf) %>%group_by(cluster)  %>%  top_n(n = 1, wt = avg_log2FC) %>% as.data.frame() %>% dplyr::distinct(.,gene,.keep_all = T) %>% head(n = 10)
    height= ceiling(length(all_top10_markers$gene)/5)*4
    VlnPlot(immune.combined, features = all_top10_markers$gene,pt.size = 0.1 ,ncol=5,cols=colors)
    ggsave(paste(seurat_exp_cluster_dir,"cluster_top1_vlnplot.pdf",sep="/"),width = 25,height = height,bg='#ffffff')
    ggsave(paste(seurat_exp_cluster_dir,"cluster_top1_vlnplot.png",sep="/"),width = 25,height = height,bg='#ffffff')
		if (length(unique(markers$cluster))<31){
			eachclusters_top5_markers=markers %>% filter(avg_log2FC != Inf) %>%group_by(cluster)  %>%  top_n(n = 5, wt = avg_log2FC) %>% as.data.frame() %>% dplyr::distinct(.,gene,.keep_all = T)
			data_ob =immune.combined
		}else {
			 eachclusters_top5_markers=markers %>% filter( cluster <= 20, avg_log2FC != Inf) %>%group_by(cluster)  %>%  top_n(n = 5, wt = avg_log2FC) %>% as.data.frame() %>% dplyr::distinct(.,gene,.keep_all = T)
			 data_ob = subset(immune.combined,seurat_clusters %in% c(0:20))
		}
		 
    #DotPlot(immune.combined, features = unique(all_top10_markers$gene))
		if(length(eachclusters_top5_markers$gene) > 9){
				direction ="vertical"
				box = "vertical"
		}else{
				direction ="vertical"
				box = "horizontal"
		}

		ggdots = Seurat::DotPlot(object = data_ob,dot.scale=4,
                    features = eachclusters_top5_markers$gene ) +							
										Seurat::RotatedAxis() +
                    guides(color = guide_colorbar(title = "Exp avg"), 
                            size = guide_legend(title = "Pct Exp %")) +		
                    #ggplot2::coord_flip() + 
                    ggplot2::scale_colour_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100) ) + 
                    theme( legend.position = "right", # 设置图例位置在底部
                            legend.direction = direction, # 设置图例排列方向为水平
                            legend.box = box, # 设置图例框的方向为水平
                            legend.box.just = "center", # 设置图例框居中显示
                            legend.spacing.x = unit(0.3, "cm"), # 设置图例条目之间的水平间距
														axis.text.x = element_text(size=9)
                        )

    ggsave(plot = ggdots,paste(seurat_exp_cluster_dir,"cluster_top5_dotplot.pdf",sep="/"),width = length(eachclusters_top5_markers$gene)*0.16+1.2,height = 0.3*length(unique(data_ob@meta.data$seurat_clusters))-1,limitsize = FALSE)
    ggsave(plot = ggdots,paste(seurat_exp_cluster_dir,"cluster_top5_dotplot.png",sep="/"),width = length(eachclusters_top5_markers$gene)*0.16+1.2,height = 0.3*length(unique(data_ob@meta.data$seurat_clusters))-1,limitsize = FALSE,bg='#ffffff')

    #FeaturePlot(immune.combined, features = all_top10_markers$gene, min.cutoff = "q9",ncol=5,order=T,cols=c("lightgrey", "red"))
    FeaturePlot(immune.combined, features = all_top10_markers$gene,ncol=5,order=T,cols=c("lightgrey", "red"))
    ggsave(paste(seurat_exp_cluster_dir,"cluster_top1_umap.pdf",sep="/"),width = 18,height = height,bg='#ffffff')
    ggsave(paste(seurat_exp_cluster_dir,"cluster_top1_umap.png",sep="/"),width = 18,height = height,bg='#ffffff')



    cluster_top10_markers=markers %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)
    #tmp <- subset(immune.combined,downsample=1000)
		DefaultAssay(immune.combined) = "RNA" ##TODO
		immune.combined <- ScaleData(immune.combined, feature=cluster_top10_markers$gene,verbose = FALSE)## 更改为只针对marker基因进行scale,以减小内存压力
    #DoHeatmap(tmp, features = unique(cluster_top10_markers$gene),group.colors = colors) + NoLegend() +theme(axis.text.y = element_text(size = 4)) 
		ggheat = Seurat::DoHeatmap( object = immune.combined,
                        features = cluster_top10_markers$gene,
                        group.colors = colors,
                        group.by = 'seurat_clusters', group.bar = T, label = F, draw.lines = T) +theme(axis.text.y = element_text(size = 4)) 
  
		ncluster = length(unique(cluster_top10_markers$cluster))
    ggsave(paste(seurat_exp_cluster_dir,"cluster_top10_markers_heatmap.pdf",sep="/"),width = ceiling(ncluster/20)*8,height = 14)
    ggsave(paste(seurat_exp_cluster_dir,"cluster_top10_markers_heatmap.png",sep="/"),width = ceiling(ncluster/20)*8,height = 14)
    
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
