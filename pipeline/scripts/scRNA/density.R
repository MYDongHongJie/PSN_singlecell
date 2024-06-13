#####Z_CQ
data_ob = OESingleCell::ReadX(input = '/public/scRNA_works/works/lumeiyun/project/scRNA/ZOE2023010050-b1-human/houxu20230208/diff/seurat.h5seurat', informat = 'h5seurat', verbose = F)

#读up
c2 <- read.delim("/public/scRNA_works/works/lumeiyun/project/scRNA/ZOE2023010050-b1-human/houxu20230208/density_FC1.5/Z_CQ/addmodulescore_up/geneset_visualization/addmodeulescore_up.xls",sep=",",header=T)

data_ob@meta.data$up = ""

### 给所有细胞进行细胞类型修改
for (i in unique(c2$up) ){
    subset_ob= subset(c2, up==i)
    data_ob@meta.data[which(data_ob@meta.data$rawbc %in% as.vector(subset_ob$Barcode)) ,"up"] = i
}
data_ob@meta.data$up = as.numeric(data_ob@meta.data$up)
BiocManager::install("Nebulosa")
library("Nebulosa")
library("Seurat")
up = plot_density(data_ob, "up")
pdf("up.pdf")
print(up)
dev.off

#down
c3 <- read.delim("/public/scRNA_works/works/lumeiyun/project/scRNA/ZOE2023010050-b1-human/houxu20230208/density_FC1.5/Z_CQ/addmodulescore_down/geneset_visualization/addmodeulescore_down.xls",sep="\t",header=T)

data_ob@meta.data$down = ""

### 给所有细胞进行细胞类型修改
for (i in unique(c3$down) ){
    subset_ob= subset(c3, down==i)
    data_ob@meta.data[which(data_ob@meta.data$rawbc %in% as.vector(subset_ob$Barcode)) ,"down"] = i
}
data_ob@meta.data$down = as.numeric(data_ob@meta.data$down)
BiocManager::install("Nebulosa")
library("Nebulosa")
library("Seurat")
down = plot_density(data_ob, "down")
pdf("down.pdf")
print(down)
dev.off

#total
c4 <- read.delim("/public/scRNA_works/works/lumeiyun/project/scRNA/ZOE2023010050-b1-human/houxu20230208/density_FC1.5/Z_CQ/addmodulescore_total/geneset_visualization/addmodeulescore_total.xls",sep="\t",header=T)

data_ob@meta.data$total = ""

### 给所有细胞进行细胞类型修改
for (i in unique(c4$total) ){
    subset_ob= subset(c4, total==i)
    data_ob@meta.data[which(data_ob@meta.data$rawbc %in% as.vector(subset_ob$Barcode)) ,"total"] = i
}
data_ob@meta.data$total = as.numeric(data_ob@meta.data$total)
BiocManager::install("Nebulosa")
library("Nebulosa")
library("Seurat")
total = plot_density(data_ob, "total")
pdf("total.pdf")
print(total)
dev.off()
