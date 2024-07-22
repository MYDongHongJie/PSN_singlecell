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




Seurat.Plot <- function(immune.combined,colors=colors,seurat_exp_cluster_dir=seurat_exp_cluster_dir,markers=markers,defaultassay="RNA",topn=30,avg_log2FC){
    DefaultAssay(immune.combined) <- defaultassay
    nsample <- length(unique(immune.combined$sample))
    ngroup <- length(unique(immune.combined$group))
    ncluster <- length(unique(immune.combined$seurat_clusters))
		cluster_overviwe = file.path(seurat_exp_cluster_dir,'2.cluster_overview')
		if(!file.exists(cluster_overviwe)){dir.create(cluster_overviwe)}
		tsne_res = file.path(seurat_exp_cluster_dir,'3.tsne')
		if(!file.exists(tsne_res)){dir.create(tsne_res)}
		umap_res = file.path(seurat_exp_cluster_dir,'4.umap')
		if(!file.exists(umap_res)){dir.create(umap_res)}

    if("celltype" %in% colnames(immune.combined@meta.data)){
        ncluster <- length(unique(immune.combined$celltype))
        #collapseby <- "celltype"
        p1 = plot.cluster.std(data = immune.combined,clusters =  "celltype", xlab = "Cluster number", log =FALSE, group = "group",legend.title = "Group",widths = c(2,2),colors=colors)
        ggsave(paste(cluster_overviwe,"cluster_number.pdf",sep="/"),width = 10,height = 6)
        ggsave(paste(cluster_overviwe,"cluster_number.png",sep="/"),width = 10,height = 6)

        #按样本/组展示细胞占比。柱状图。横坐标是样本/组，纵坐标是细胞占比
        temp <- data.frame(table(immune.combined$celltype,immune.combined$sample))
    }else{
        #collapseby <- "seurat_clusters"
        ncluster <- length(unique(immune.combined$seurat_clusters))
        p1 = plot.cluster.std(data = immune.combined,clusters =  "seurat_clusters", xlab = "Cluster number", log =FALSE, group = "group",legend.title = "Group",widths = c(2,2),colors=colors)
        ggsave(paste(cluster_overviwe,"cluster_number.pdf",sep="/"),width = 6,height = 6)
        ggsave(paste(cluster_overviwe,"cluster_number.png",sep="/"),width = 6,height = 6)
        temp <- data.frame(table(immune.combined$seurat_clusters,immune.combined$sample))
    }

    #按样本/组展示细胞占比。柱状图。横坐标是样本/组，纵坐标是细胞占比
	colnames(temp) <- c("Cluster","Sample","number")
	p <- ggplot(temp,aes(x=Sample,fill=Cluster,y=number))+
	  geom_bar(stat = "identity",position ='fill')+
	  labs(x='',y='')+
		scale_y_continuous(expand = c(0,0))+
	  guides(fill = guide_legend(title = '',reverse=TRUE))+
	  scale_fill_manual(values = colors)+
	  theme_bw() + 
		theme(panel.grid.major = element_blank(), 
		strip.background = element_rect(fill = NA,color = NA), 
		panel.border = element_blank(), 
		axis.ticks.x = element_blank(), 
    axis.text = element_text(color = "black"), 
		axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)
		)+theme_cowplot()
	ggsave(p,file=paste(cluster_overviwe,"cluster_sample_cellcounts.png",sep="/"),width = 8,height = 7,bg='white')
	ggsave(p,file=paste(cluster_overviwe,"cluster_sample_cellcounts.pdf",sep="/"),width = 8,height = 7,bg='white')
    if("celltype" %in% colnames(immune.combined@meta.data)){
        temp <- data.frame(table(immune.combined$celltype,immune.combined$group))
    }else{
        temp <- data.frame(table(immune.combined$seurat_clusters,immune.combined$group))
    }
    colnames(temp) <- c("Cluster","Group","number")
    p <- ggplot(temp,aes(x=Group,fill=Cluster,y=number))+
      geom_bar(stat = "identity",position ='fill')+
      labs(x='',y='')+
			scale_y_continuous(expand = c(0,0))+
      guides(fill = guide_legend(title = '',reverse=TRUE))+
      scale_fill_manual(values = colors)+ 
			theme(panel.grid.major = element_blank(), 
			strip.background = element_rect(fill = NA,color = NA), 
			panel.border = element_blank(), 
			axis.ticks.x = element_blank(), 
			axis.text = element_text(color = "black"), 
			axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)
			)+theme_cowplot()
    ggsave(p,file=paste(cluster_overviwe,"cluster_group_cellcounts.png",sep="/"),width = 8,height = 7,bg='white')
    ggsave(p,file=paste(cluster_overviwe,"cluster_group_cellcounts.pdf",sep="/"),width = 8,height = 7,bg='white')
    
    #分样本展示
    p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample",cols=colors)
    ggsave(p1,file=paste(umap_res,"cluster_umap.sample.pdf",sep="/"),width = 4.8,height = 4)
    ggsave(p1,file=paste(umap_res,"cluster_umap.sample.png",sep="/"),width = 4.8,height = 4)
    p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "sample",cols=colors)
    ggsave(p1,file=paste(tsne_res,"cluster_tsne.sample.pdf",sep="/"),width = 4.8,height = 4)
    ggsave(p1,file=paste(tsne_res,"cluster_tsne.sample.png",sep="/"),width = 4.8,height = 4)

    #分组展示
    p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "group",cols=colors)
    ggsave(p1,file=paste(umap_res,"cluster_umap.group.pdf",sep="/"),width = 4.8,height = 4)
    ggsave(p1,file=paste(umap_res,"cluster_umap.group.png",sep="/"),width = 4.8,height = 4)
    p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "group",cols=colors)
    ggsave(p1,file=paste(tsne_res,"cluster_tsne.group.pdf",sep="/"),width = 4.8,height = 4)
    ggsave(p1,file=paste(tsne_res,"cluster_tsne.group.png",sep="/"),width = 4.8,height = 4)

    #分cluster展示
    p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,cols=colors) +NoLegend()
    ggsave(paste(umap_res,"cluster_umap.label.pdf",sep="/"),width = 4,height = 4)
    ggsave(paste(umap_res,"cluster_umap.label.png",sep="/"),width = 4,height = 4)
    p2 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE,cols=colors) +NoLegend()
    ggsave(paste(tsne_res,"cluster_tsne.label.pdf",sep="/"),width = 4,height = 4)
    ggsave(paste(tsne_res,"cluster_tsne.label.png",sep="/"),width = 4,height = 4)
    p2 <- DimPlot(immune.combined, reduction = "umap", label = FALSE,cols=colors)
    ggsave(paste(umap_res,"cluster_umap.pdf",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6)
    ggsave(paste(umap_res,"cluster_umap.png",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6)
    p2 <- DimPlot(immune.combined, reduction = "tsne", label = FALSE,cols=colors)
    ggsave(paste(tsne_res,"cluster_tsne.pdf",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6)
    ggsave(paste(tsne_res,"cluster_tsne.png",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6)


    if(nsample==1){swidth <- 4}else{swidth <- 8}
    sheight <- ceiling(nsample/2)*4
    if(ngroup==1){gwidth <- 4}else{gwidth <- 8}
    gheight <- ceiling(ngroup/2)*4
    #分样本展示 按样本分页;分组展示，按组分页
    DimPlot(immune.combined, reduction = "umap", split.by = "sample",ncol=2,cols=colors,label=TRUE)+NoLegend()
    ggsave(paste(umap_res,"cluster_umap_splitsample.pdf",sep="/"),width = swidth,height = sheight)
    ggsave(paste(umap_res,"cluster_umap_splitsample.png",sep="/"),width = swidth,height = sheight)
    DimPlot(immune.combined, reduction = "tsne", split.by = "sample",ncol=2,cols=colors,label=TRUE)+NoLegend()
    ggsave(paste(tsne_res,"cluster_tsne_splitsample.pdf",sep="/"),width = swidth,height = sheight)
    ggsave(paste(tsne_res,"cluster_tsne_splitsample.png",sep="/"),width = swidth,height = sheight)

    DimPlot(immune.combined, reduction = "umap", split.by = "group",ncol=2,cols=colors,label=TRUE)+NoLegend()
    ggsave(paste(umap_res,"cluster_umap_splitgroup.pdf",sep="/"),width = gwidth,height = gheight)
    ggsave(paste(umap_res,"cluster_umap_splitgroup.png",sep="/"),width = gwidth,height = gheight)
    DimPlot(immune.combined, reduction = "tsne", split.by = "group",ncol=2,cols=colors,label=TRUE)+NoLegend()
    ggsave(paste(tsne_res,"cluster_tsne_splitgroup.pdf",sep="/"),width = gwidth,height = gheight)
    ggsave(paste(tsne_res,"cluster_tsne_splitgroup.png",sep="/"),width = gwidth,height = gheight)

### 结果展示
    #暂定放在process文件夹
		preprocess = file.path(seurat_exp_cluster_dir,'1.preprocess')
		if(!file.exists(preprocess)){
    	dir.create(preprocess)
    }
		scFeature=c('nCount_RNA','nFeature_RNA',"percent.mt")
		
		for (feature in scFeature){
			if (feature %in% colnames(immune.combined@meta.data)){
				nfeatureplot = FeaturePlot(immune.combined, reduction = "umap", label = F,features=feature)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
    		ggsave(plot=nfeatureplot,filename = file.path(preprocess,paste0(feature,"_umap.pdf")),width = 7,height = 7)
				ggsave(plot=nfeatureplot,filename = file.path(preprocess,paste0(feature,"_umap.png")),width = 7,height = 7,dpi=1000)
			}
		}

    #补充一张cellcyle的图
		if ('cells_cell_cycle.pdf' %in% dir(preprocess)){
			p<-dittoBarPlot(immune.combined,"Phase", group.by = "seurat_clusters",color.panel=colors)
      ggsave(p,filename=paste(preprocess,"cell_cycle_per_cluster_barplot.pdf",sep="/"))
      ggsave(p,filename=paste(preprocess,"cell_cycle_per_cluster_barplot.png",sep="/"))
			p2 <- DimPlot(immune.combined,group.by="Phase",cols=colors,reduction='umap',label=F)+theme(axis.text.x = element_text(size=14,face="bold"),
		axis.text.y = element_text(size=14,face="bold"),
		axis.title.x = element_text(size=16,face="bold"),
		axis.title.y = element_text(size=16,face="bold"),
		aspect.ratio = 1/1)
			ggsave(p2,filename=paste(preprocess,"cell_cycle_per_cluster_umap.png",sep="/"),dpi=1000,bg = 'white')
			ggsave(p2,filename=paste(preprocess,"cell_cycle_per_cluster_umap.pdf",sep="/"),bg = 'white')
		}
	
}
