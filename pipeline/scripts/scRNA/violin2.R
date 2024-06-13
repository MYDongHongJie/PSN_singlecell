# module purge && module load OESingleCell/2.0.0

suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("Seurat")) 
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("future")) ## 可以学习下

## ,text= element_text(family="Arial") ## png 可以保存， pdf 无法保存。
setwd("/public/scRNA/works/liuhongyan/Project/scRNA/done_report/HT2020-11162-3-human/Further_analysis_2022-03-09/1.picture/1.QC")

rds="/public/scRNA/works/liuhongyan/Project/scRNA/done_report/HT2020-11162-3-human/SA20210208/singlecell_object.clustering_resolution0.1.rds"
seurat_ob = readRDSMC( rds, cores = availableCores() )

# 一、质控部分（4张小图）
# 1.每个样本的基因数-小提琴图
# 2.每个样本的UMIs-小提琴图 
# 3.UMIs与线粒体基因-散点图
# 4.UMIs和基因的关系-散点图

## 样本颜色编号：#1F7784, #FF7F0E, #2CA02C, #D52122
### gene
vln4nGene = VlnPlot(object = seurat_ob, features = "nFeature_RNA",x.lab.rot = T, alpha = 0.6, cols = c("#1F7784", "#FF7F0E", "#2CA02C", "#D52122") ,
                            group.by = "sampleid",ncol = 1, pt.size = 0)+ggtitle("") + ylab("Number of genes") +
                            theme( legend.position="none", axis.title.x= element_blank(),axis.text.y=element_text(size=20), axis.text=element_text(size=22), title=element_text(size=24),plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"), axis.ticks.length =unit(0.4, "lines"),axis.text.x=element_text(angle=0,hjust=0.5) ) 
output_dir="plot"
dir.create(output_dir)
ggsave(file.path(output_dir,paste0("afterQC_total_genes4each_cell_on_violin_plot.png", collapse = "")),
        dpi=1000, plot = vln4nGene)
ggsave(file.path(output_dir,paste0("afterQC_total_genes4each_cell_on_violin_plot.pdf", collapse = "")),
        plot = vln4nGene)


## umi 
assay2use="RNA"
nCount = paste0( "nCount_", assay2use)
vln4nUMI = VlnPlot(object = seurat_ob, features = nCount ,alpha = 0.6, cols = c("#1F7784", "#FF7F0E", "#2CA02C", "#D52122") ,
                    group.by = "sampleid", ncol = 1,
                    x.lab.rot = T,pt.size = 0) +ggtitle("") +ylab("Number of UMIs") +
                    theme( legend.position="none", axis.title.x= element_blank(),axis.text.y=element_text(size=20), axis.text=element_text(size=22), title=element_text(size=24), plot.margin = unit(c(0, 0.5, 0, 0.5), "cm") , axis.ticks.length =unit(0.4, "lines"),axis.text.x=element_text(angle=0,hjust=0.5) )


ggsave(file.path(output_dir,paste0("afterQC_total_UMIs4each_cell_on_violin_plot.png")),
        dpi=1000, plot = vln4nUMI)
ggsave(file.path(output_dir,paste0("afterQC_total_UMIs4each_cell_on_violin_plot.pdf")), plot = vln4nUMI)


## 线粒体比例和UMI 散点图

meta=seurat_ob@meta.data

## plot.margin = unit(c(0, 0.2, 0, 0), "cm")   unit中的四个值代表top, right, bottom, left
# p <- ggplot(meta, aes(x=nCount_RNA,y=percent.mito)) + geom_point(aes(color=sampleid))  +
#         scale_colour_manual(values=c("#1F7784", "#FF7F0E", "#2CA02C", "#D52122") ) + theme_classic() +
#         theme( legend.position="top", legend.title=element_blank() , legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=24),plot.margin = unit(c(0, 1, 0, 1), "cm") , axis.ticks.length =unit(0.4, "lines") ) +
#         xlab("Number of UMIs")+ ylab("% of mitochondrial genes") 

# ## theme_cowplot()
# p <- ggplot(meta, aes(x=nCount_RNA,y=percent.mito)) + geom_point(aes(color=sampleid))  +
#         scale_colour_manual(values=c("#1F7784", "#FF7F0E", "#2CA02C", "#D52122") ) + theme_cowplot() +
#         theme( legend.position="top", legend.title=element_blank() , legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=24),plot.margin = unit(c(0, 1, 0, 1), "cm") , axis.ticks.length =unit(0.4, "lines") ) +
#         xlab("Number of UMIs")+ ylab("% of mitochondrial genes")
p <- FeatureScatter(object = seurat_ob, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size=1, group.by="sampleid", cols=c("#1F7784", "#FF7F0E", "#2CA02C", "#D52122") ) + ggtitle("") + 
        theme( legend.title=element_blank() , legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=24),plot.margin = unit(c(0, 1, 0, 1), "cm") , axis.ticks.length =unit(0.4, "lines"),legend.position = "top", legend.direction = "horizontal",legend.justification=c(0.5,1) ) +
        xlab("Number of UMIs")+ ylab("% of mitochondrial genes")


ggsave(file.path(output_dir,paste0("UMIs_percent.mito.png")),
        dpi=1000, plot = p)
ggsave(file.path(output_dir,paste0("UMIs_percent.mito.pdf")), plot = p)

## UMIs和基因的关系-散点图
# p <- ggplot(meta, aes(x=nCount_RNA,y=nFeature_RNA)) + geom_point(aes(color=sampleid))  +
#         scale_colour_manual(values=c("#1F7784", "#FF7F0E", "#2CA02C", "#D52122") ) + theme_classic() +
#         theme( legend.position="top", legend.title=element_blank() , legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=24),plot.margin = unit(c(0, 1, 0, 1), "cm") , axis.ticks.length =unit(0.4, "lines")) +
#         xlab("Number of UMIs")+ ylab("Number of genes") 


p <- FeatureScatter(object = seurat_ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=1, group.by="sampleid", cols=c("#1F7784", "#FF7F0E", "#2CA02C", "#D52122") ) + ggtitle("") + 
        theme( legend.title=element_blank() , legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=24),plot.margin = unit(c(0, 1, 0, 1), "cm") , axis.ticks.length =unit(0.4, "lines"),legend.position = "top", legend.direction = "horizontal",legend.justification=c(0.5,1) ) +
        xlab("Number of UMIs")+ ylab("Number of genes") 

ggsave(file.path(output_dir,paste0("UMIs_genes.png")),
        dpi=1000, plot = p)
ggsave(file.path(output_dir,paste0("UMIs_genes.pdf")), plot = p)

