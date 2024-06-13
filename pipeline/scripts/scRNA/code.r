.libPaths("/home/ziqingzhen/R/x86_64-conda_cos6-linux-gnu-library/3.6")
library(Seurat)
library(readr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(paletteer)
#BiocManager::install("paletteer")

pbmc = read_rds("Fibroblast_new_celltypes_20220308.rds")
pbmc@meta.data$new_sampleid = ""
pbmc@meta.data[which(pbmc@meta.data$sampleid%in%"P1_M"),"new_sampleid"] = "P1-Adj"
pbmc@meta.data[which(pbmc@meta.data$sampleid%in%"P1_SCC"),"new_sampleid"] = "P1-Ca"
pbmc@meta.data[which(pbmc@meta.data$sampleid%in%"P2_M"),"new_sampleid"] = "P2-Adj"
pbmc@meta.data[which(pbmc@meta.data$sampleid%in%"P2_SCC"),"new_sampleid"] = "P2-Ca"
pbmc@meta.data[which(pbmc@meta.data$sampleid%in%"P3_M"),"new_sampleid"] = "P3-Adj"
pbmc@meta.data[which(pbmc@meta.data$sampleid%in%"P3_SCC"),"new_sampleid"] = "P3-Ca"
pbmc@meta.data[which(pbmc@meta.data$sampleid%in%"P4_M"),"new_sampleid"] = "P4-Adj"
pbmc@meta.data[which(pbmc@meta.data$sampleid%in%"P4_SCC"),"new_sampleid"] = "P4-Ca"

table(pbmc@meta.data$new_sampleid)

pbmc@meta.data$new_group = ""
pbmc@meta.data[which(pbmc@meta.data$group%in%"SCC"),"new_group"] = "Cancer"
pbmc@meta.data[which(pbmc@meta.data$group%in%"Mucosa"),"new_group"] = "Adjacent"
table(pbmc@meta.data$new_group)


source("./custom_seurat_functions.R")
############### 即可开始绘图

p = plot.clusters.group(data = pbmc,clusters =  "new_celltypes", xlab = "Cluster number", log = FALSE, group = "new_sampleid",group2 = "new_group",legend.title = "Patient",widths = c(5,3,4),color = 2)

pdf("part7.pdf",width =12,height = 3 )
p
dev.off()
ggsave("part7.png",plot = p,width =12,height = 3 ,dpi=1000)
