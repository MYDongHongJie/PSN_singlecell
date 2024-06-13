library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

seurat_ob<-readRDS("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20220811/OE2018H2930S/v3_object.clustering_resolution0.3.rds")

p=DimPlot(seurat_ob,group.by = "sampleid",reduction = "tsne") + scale_color_manual(values=c("#e5ddb2","#d48571","#5e8ac0"))
ggsave("groupby_sampleid.png",plot=p,width =7.4)
ggsave("groupby_sampleid.pdf",plot=p,width =7.4)

plot_data = seurat_ob@meta.data[,c("sampleid","clusters")] %>% group_by(.dots = c("sampleid","clusters")) %>% dplyr::summarize(cell_number = n()) %>% mutate(freq = round((cell_number / sum(cell_number)) * 100,2))
plot_data$Freq_group = ""
for (i in 1:24){
    cluster_sub = as.character(plot_data[i,"clusters"])
    total = sum(subset(plot_data,clusters == cluster_sub)$cell_number)
    plot_data$Freq_group[i] = (plot_data$cell_number[i]/total)
}
plot_data$Freq_group=as.numeric(plot_data$Freq_group)
plot_data$clusters=factor(plot_data$clusters, levels=c(7:0))
plot_data$sampleid=gsub("62CAFs","CAF-62",plot_data$sampleid)
plot_data$sampleid=gsub("71CAFs","CAF-71",plot_data$sampleid)
plot_data$sampleid=gsub("95CAFs","CAF-95",plot_data$sampleid)
dot = ggplot(plot_data) + geom_point(pch = 16,aes(x = sampleid, y = clusters, color = sampleid,size = Freq_group)) + 
      scale_color_manual(values=c("CAF-62"="#e5ddb2","CAF-71"="#d48571","CAF-95"="#5e8ac0")) + scale_size(range = c(2, 10),limits = c(0,1),breaks = c(0,0.2,1)) + 
      theme_linedraw() + theme(panel.grid.major = element_blank(),panel.border = element_blank(), 
      axis.title.y = element_blank(), axis.ticks.y = element_blank(),axis.title.x=element_blank(), axis.ticks.x = element_blank(), 
      axis.text.x = element_text(size=12, angle=45, hjust=1),axis.text.y = element_text(size=14),legend.text = element_text(size=13) )+
      guides(color=guide_legend(title=NULL,override.aes = list(size = 7)), size = guide_legend(title=NULL), shape = NULL)

ggsave("dotplot.png",plot=dot,width = 3.5, height = 5)
ggsave("dotplot.pdf",plot=dot,width = 3.5, height = 5)

p=p+theme(legend.position= "none" )

aa=plot_grid(p, dot, scale = c(1,0.8),rel_widths = c(2, 1.2))
ggsave("dimplot_dotplot.png",plot=aa,width = 10)
ggsave("dimplot_dotplot.pdf",plot=aa,width = 10)

########################################################

seurat_seurat<-readRDS("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/HT2021-14371_20210812/new_clusters_CAF.rds")
seurat_seurat@meta.data[which(seurat_seurat@meta.data$new_clusters %in% c(1,5,10)),"new_celltype"] = "iCAF"
seurat_seurat@meta.data[which(seurat_seurat@meta.data$new_clusters %in% c(2,3,6,7,8)),"new_celltype"] = "myCAF"
seurat_seurat@meta.data[which(seurat_seurat@meta.data$new_clusters %in% c(4)),"new_celltype"] = "unCAF"
seurat_seurat@meta.data[which(seurat_seurat@meta.data$new_clusters %in% c(9)),"new_celltype"] = "unkown"

seurat_seurat@meta.data$new_celltype=factor(seurat_seurat@meta.data$new_celltype, levels=c("iCAF", "myCAF","unCAF","unkown"))
## visualization 
toid="new_celltype"
nlevel = length(unique(seurat_seurat@meta.data[,toid]))
opt=list();opt$reduct="tsne"
output_dir="HT2021-14371_20210812"
dir.create(output_dir,recursive=T)
if (dim(seurat_seurat)[2] < 500){
  pointsize=1.5
}else{
  pointsize=0.5
}
ggtsne = DimPlot(object = seurat_seurat, reduction = opt$reduct,pt.size = pointsize,group.by=toid) +
    ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_colour_manual( values = c(rgb(245,174,89,maxColorValue = 255),rgb(137,201,111,maxColorValue = 255),c(rgb(184,101,185,maxColorValue = 255),"#5e8ac0")))
ggplot2::ggsave(file.path(output_dir,paste0(toid,".pdf")), ggtsne, width = max(nchar( as.vector(unique(seurat_seurat@meta.data[,toid]))))/15+7)
ggplot2::ggsave(file.path(output_dir,paste0(toid,".png")), ggtsne, width = max(nchar( as.vector(unique(seurat_seurat@meta.data[,toid]))))/15+7, dpi=1000)

simplified_meta = seurat_seurat@meta.data %>%
                               dplyr::rename( "Barcode" = "orig.ident") %>%
                               dplyr::select( Barcode, sampleid, clusters,group,!!toid)
write.table(simplified_meta, quote = F,sep =",",row.names = F,
             file.path(output_dir,paste0(toid,".metadata.csv",collapse = ""))) 
saveRDS(seurat_seurat,file.path(output_dir,paste0(toid,Sys.Date(),".rds")))


feat = FeaturePlot(seurat_seurat,features = "TRPA1",cols = c("grey","red"), reduction= opt$reduct,ncol = 2, pt.size = 0.7, order= T, max.cutoff = 3 ) +
       theme( plot.title = element_text(hjust = 0.5)) 
ggplot2::ggsave(file.path(output_dir,paste0("TRPA1_FeaturePlot.pdf")), feat, width = max(nchar( as.vector(unique(seurat_seurat@meta.data[,toid]))))/15+7)
ggplot2::ggsave(file.path(output_dir,paste0("TRPA1_FeaturePlot.png")), feat, width = max(nchar( as.vector(unique(seurat_seurat@meta.data[,toid]))))/15+7, dpi=1000)
########################################################
library(pheatmap)
library(RColorBrewer)
heatmap1<-read.delim("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/HT2021-14371_20210812/heatmap1.xls",sep="\t",row.names=1)
heatmap2<-read.delim("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/HT2021-14371_20210812/heatmap2.xls",sep="\t",row.names=1)

heatmap_data = as.matrix(t(heatmap1))
heatmap_data = as.matrix(t(heatmap2))

bk <- c(seq(-1,0,by=0.01),seq(0.01,2,by=0.01))

ht1 =  pheatmap(heatmap_data, cellheight=30,cellwidth=30,
                colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk)),
                cluster_rows = F, border=FALSE, scale = "column", legend_breaks = seq(-1,2,1),
                breaks = bk, treeheight_col = 10, angle_col = "45")
ggsave("HT2021-14371_20210812/heatmap1.png", ht1, height=3 ,width = 10)
ggsave("HT2021-14371_20210812/heatmap1.pdf", ht1, height=3 ,width = 10)

ggsave("HT2021-14371_20210812/heatmap2.png", ht1, height=4 ,width = 12)
ggsave("HT2021-14371_20210812/heatmap2.pdf", ht1, height=4 ,width = 12)

########################################################
genelist = read.delim("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/HT2021-14371_20210812/genelist.txt")
genes = as.vector(genelist$Gene)
feat = FeaturePlot(seurat_seurat, features = genes, cols = c("grey","red"), reduction= opt$reduct,ncol = 4, pt.size = 0.6, order= T ) +
       theme( plot.title = element_text(hjust = 0.5)) 
ggplot2::ggsave(file.path(output_dir,paste0("FeaturePlot.pdf")), feat,width = 24, height = 38 )
ggplot2::ggsave(file.path(output_dir,paste0("FeaturePlot.png")), feat,width = 24, height = 38 )

########################################################
library(monocle)

gbm_cds126 <- readRDS("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/pseudo_plot/2.figure2-pseudo_plot/2.cluster126/pseudotime_results.rds")
gbm_cds035 <- readRDS("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/pseudo_plot/2.figure2-pseudo_plot/1.cluster035/pseudotime_results.rds")
output_dir = "pseudo_plot"

gbm_cds126 <- orderCells(gbm_cds126, root_state = 2) # root_state =分支编
p <- plot_cell_trajectory(gbm_cds126, color_by = "Pseudotime", show_branch_points = F) + scale_colour_viridis_c(option = "inferno")
ggsave(file.path(output_dir, "visualize_126/cell_trajectory_color_by_Pseudotime_126.pdf"), plot = p)
ggsave(file.path(output_dir, "visualize_126/cell_trajectory_color_by_Pseudotime_126.png"), plot = p, dpi = 1000)
saveRDS(gbm_cds126, "/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/pseudo_plot/2.figure2-pseudo_plot/2.cluster126/pseudotime_results.rds")

gbm_cds035 <- orderCells(gbm_cds035, root_state = 5) # root_state =分支编
p <- plot_cell_trajectory(gbm_cds035, color_by = "Pseudotime", show_branch_points = F) + scale_colour_viridis_c(option = "inferno")
ggsave(file.path(output_dir, "visualize_035/cell_trajectory_color_by_Pseudotime_035.pdf"), plot = p)
ggsave(file.path(output_dir, "visualize_035/cell_trajectory_color_by_Pseudotime_035.png"), plot = p, dpi = 1000)
saveRDS(gbm_cds035, "/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/pseudo_plot/2.figure2-pseudo_plot/1.cluster035/pseudotime_results.rds")

genes035 = read.table("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/pseudo_plot/genelist035.txt",sep="\t",header=T)
genes126 = read.table("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2020-13912/houxu-20221018/pseudo_plot/genelist126.txt",sep="\t",header=T)
genes035$limit = as.vector(genes035$limit)
genes126$limit = as.vector(genes126$limit)
for (i in as.vector(genes126$genes)) {
    # limits = as.numeric(unlist(strsplit(genes126[which(genes126$genes == i),"limit"],",")))
    p <- plot_cell_trajectory(gbm_cds126, markers = i, use_color_gradient = T, show_branch_points = F, show_tree = F) + 
         theme(legend.text = element_text(size = 10)) + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100) )#, limits=limits[c(1,3)], breaks = limits
    ggsave(file.path(output_dir, paste0("visualize_126/",i, ".pdf")))
    ggsave(file.path(output_dir, paste0("visualize_126/",i, ".png")), dpi = 1000)

}

for (i in as.vector(genes035$genes)) {
    p <- plot_cell_trajectory(gbm_cds035, markers = i, use_color_gradient = T, show_branch_points = F, show_tree = F) + 
         theme(legend.text = element_text(size = 10)) + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100) )#,limits=c(0.5,2),breaks=c(1,2,0.5)
    ggsave(file.path(output_dir, paste0("visualize_035/",i, ".pdf")))
    ggsave(file.path(output_dir, paste0("visualize_035/",i, ".png")), dpi = 1000)
}

for (i in as.vector(genes126$genes)) {
    limits = as.numeric(unlist(strsplit(genes126[which(genes126$genes == i),"limit"],",")))
    p <- plot_cell_trajectory_cutoff(gbm_cds126, markers = i, use_color_gradient = T, show_branch_points = F, show_tree = F,min.cutoff = limits[1], max.cutoff = limits[3]) + 
         theme(legend.text = element_text(size = 10)) + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100) )
    ggsave(file.path(output_dir, paste0("visualize_126/",i, "_limit.pdf")))
    ggsave(file.path(output_dir, paste0("visualize_126/",i, "_limit.png")), dpi = 1000)
}

for (i in as.vector(genes035$genes)) {
    limits = as.numeric(unlist(strsplit(genes035[which(genes035$genes == i),"limit"],",")))
    p <- plot_cell_trajectory_cutoff(gbm_cds035, markers = i, use_color_gradient = T, show_branch_points = F, show_tree = F, min.cutoff = limits[1], max.cutoff = limits[3]) + 
         theme(legend.text = element_text(size = 10)) + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100) )
    ggsave(file.path(output_dir, paste0("visualize_035/",i, "_limit.pdf")))
    ggsave(file.path(output_dir, paste0("visualize_035/",i, "_limit.png")), dpi = 1000)
}
