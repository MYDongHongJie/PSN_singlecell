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

Seurat.Plot <- function(immune.combined,colors=colors,seurat_exp_cluster_dir=seurat_exp_cluster_dir,markers=markers,defaultassay="RNA"){
    DefaultAssay(immune.combined) <- defaultassay
    nsample <- length(unique(immune.combined$sample))
    ngroup <- length(unique(immune.combined$group))
    ncluster <- length(unique(immune.combined$seurat_clusters))
    if("celltype" %in% colnames(immune.combined@meta.data)){
        ncluster <- length(unique(immune.combined$celltype))
        p1 = plot.cluster.std(data = immune.combined,clusters =  "celltype", xlab = "Cluster number", log =FALSE, group = "group",legend.title = "Group",widths = c(2,2),colors=colors)
        ggsave(paste(seurat_exp_cluster_dir,"cluster_number.pdf",sep="/"),width = 10,height = 6,limitsize = FALSE)
        ggsave(paste(seurat_exp_cluster_dir,"cluster_number.png",sep="/"),width = 10,height = 6,limitsize = FALSE)

        #按样本/组展示细胞占比。柱状图。横坐标是样本/组，纵坐标是细胞占比
        temp <- data.frame(table(immune.combined$celltype,immune.combined$sample))
    }else{
        ncluster <- length(unique(immune.combined$seurat_clusters))
        p1 = plot.cluster.std(data = immune.combined,clusters =  "seurat_clusters", xlab = "Cluster number", log =FALSE, group = "group",legend.title = "Group",widths = c(2,2),colors=colors)
        ggsave(paste(seurat_exp_cluster_dir,"cluster_number.pdf",sep="/"),width = 6,height = 6,limitsize = FALSE)
        ggsave(paste(seurat_exp_cluster_dir,"cluster_number.png",sep="/"),width = 6,height = 6,limitsize = FALSE)
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
	ggsave(p,file=paste(seurat_exp_cluster_dir,"cluster_sample_cellcounts.png",sep="/"),width = 8,height = 7)
	ggsave(p,file=paste(seurat_exp_cluster_dir,"cluster_sample_cellcounts.pdf",sep="/"),width = 8,height = 7)
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
    ggsave(p,file=paste(seurat_exp_cluster_dir,"cluster_group_cellcounts.png",sep="/"),width = 8,height = 7)
    ggsave(p,file=paste(seurat_exp_cluster_dir,"cluster_group_cellcounts.pdf",sep="/"),width = 8,height = 7)
    




    #分样本展示
    p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample",cols=colors)
    ggsave(p1,file=paste(seurat_exp_cluster_dir,"cluster_umap.sample.pdf",sep="/"),width = 6.8,height = 6,limitsize = FALSE)
    ggsave(p1,file=paste(seurat_exp_cluster_dir,"cluster_umap.sample.png",sep="/"),width = 6.8,height = 6,limitsize = FALSE)
    p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "sample",cols=colors)
    ggsave(p1,file=paste(seurat_exp_cluster_dir,"cluster_tsne.sample.pdf",sep="/"),width = 6.8,height = 6,limitsize = FALSE)
    ggsave(p1,file=paste(seurat_exp_cluster_dir,"cluster_tsne.sample.png",sep="/"),width = 6.8,height = 6,limitsize = FALSE)

    #分组展示
    p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "group",cols=colors)
    ggsave(p1,file=paste(seurat_exp_cluster_dir,"cluster_umap.group.pdf",sep="/"),width = 6.8,height = 6,limitsize = FALSE)
    ggsave(p1,file=paste(seurat_exp_cluster_dir,"cluster_umap.group.png",sep="/"),width = 6.8,height = 6,limitsize = FALSE)
    p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "group",cols=colors)
    ggsave(p1,file=paste(seurat_exp_cluster_dir,"cluster_tsne.group.pdf",sep="/"),width = 6.8,height = 6,limitsize = FALSE)
    ggsave(p1,file=paste(seurat_exp_cluster_dir,"cluster_tsne.group.png",sep="/"),width = 6.8,height = 6,limitsize = FALSE)

    #分cluster展示
    p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,cols=colors) +NoLegend()
    ggsave(paste(seurat_exp_cluster_dir,"cluster_umap.label.pdf",sep="/"),width = 6,height = 6,limitsize = FALSE)
    ggsave(paste(seurat_exp_cluster_dir,"cluster_umap.label.png",sep="/"),width = 6,height = 6,limitsize = FALSE)
   
    p2 <- DimPlot(immune.combined, reduction = "umap", label = FALSE,cols=colors)
    ggsave(paste(seurat_exp_cluster_dir,"cluster_umap.pdf",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6,limitsize = FALSE)
    ggsave(paste(seurat_exp_cluster_dir,"cluster_umap.png",sep="/"),width = 6+ceiling(ncluster/13)*0.5,height = 6,limitsize = FALSE)
    

    sout_width=ceiling(sqrt(nsample))+1
    sheight <- nsample/sout_width*6
    #swidth <-ceiling(sout_width/2)*5 + ceiling(nsample/13)*0.5
    swidth <- sout_width*5 + ceiling(nsample/13)*0.5
    gout_width=ceiling(sqrt(ngroup))
    gheight <- ngroup/gout_width*6
    #gwidth <-ceiling(gout_width/2)*5 + ceiling(ngroup/13)*0.5
    gwidth <-gout_width*5 + ceiling(ngroup/13)*0.5
    #分样本展示 按样本分页;分组展示，按组分页
    DimPlot(immune.combined, reduction = "umap", split.by = "sample",ncol=sout_width,cols=colors)
    ggsave(paste(seurat_exp_cluster_dir,"cluster_umap_splitsample.pdf",sep="/"),width = swidth,height = sheight,limitsize = FALSE)
    ggsave(paste(seurat_exp_cluster_dir,"cluster_umap_splitsample.png",sep="/"),width = swidth,height = sheight,limitsize = FALSE)
    
    DimPlot(immune.combined, reduction = "umap", split.by = "group",ncol=gout_width,cols=colors,label.size = 2)
    
    ggsave(paste(seurat_exp_cluster_dir,"cluster_umap_splitgroup.pdf",sep="/"),width = gwidth,height = gheight,limitsize = FALSE)
    ggsave(paste(seurat_exp_cluster_dir,"cluster_umap_splitgroup.png",sep="/"),width = gwidth,height = gheight,limitsize = FALSE)
    

    #png(paste(seurat_exp_cluster_dir,"cluster_top10_markers_heatmap.png",sep="/"),width =1200,height = 800)
    #DoHeatmap(tmp, features = unique(cluster_top10_markers$gene)) + NoLegend() +theme(axis.text.y = element_text(size = 4))
    #dev.off()
}
