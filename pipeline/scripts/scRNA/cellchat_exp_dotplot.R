#一、根据附件挑选出来的受配体细胞和受配体关系对，绘制目标关系对的dotplot（两个表格分别绘图，共2张图）
#要求：映射表达量情况
#  路径：DZOE2023020735_周智航老师_20230823_Fib\2.Cellchat\By_new_celltype
#工作路径
/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230906_Plot

一、根据附件挑选出来的受配体细胞和受配体关系对，绘制目标关系对的dotplot（两个表格分别绘图，共2张图），如下：
      路径：DZOE2023020735_周智航老师_20230823_Fib\2.Cellchat\By_new_celltype

#绘图表格
/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230823_Fib/2.Cellchat/Cellchat_sub_new_celltype/communication.xls

#绘图数据
/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230823_Fib/2.Cellchat/Cellchat_sub_new_celltype/cellchat_results.rds

cellchat_list = readRDS("/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230823_Fib/2.Cellchat/Cellchat_sub_new_celltype/cellchat_list.rds")

cellchat = readRDS("/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230906_Plot/1.Cellchat_plot/cellchat_results.rds")

.libPaths("/public/scRNA_works/works/guokaiqi/software/R/x86_64-conda-linux-gnu-library/4.0")
source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/CellChat_functions.R")

suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
    library(patchwork)
    library(ggalluvial)
    library(NMF)
    library(ggplot2)
    library(optparse)
    library(ComplexHeatmap)
    library(future)
    library(tidyverse)
    library(corrplot)
    library(circlize)
    library(pheatmap)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(Cairo)
})

####################绘图代码Fib_significant_interactions_bubble_plot#######################
pairLRsig <- cellchat@LR$LRsig
  group <- cellchat@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
nLR <- nrow(pairLRsig)
numCluster <- nlevels(group)



data = as.matrix(cellchat@data.signaling)
  data.use <- data/max(data)
  nC <- ncol(data.use)

data.use.avg <- aggregate(t(data.use), list(group), FUN =triMean)

data.use.avg <- t(data.use.avg[,-1])
colnames(data.use.avg) <- levels(group)

complex_input <- cellchat@DB$complex
cofactor_input <- cellchat@DB$cofactor
dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "A")
dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "I")
dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor

rownames(dataLavg) = geneL
rownames(dataRavg) = geneR

colnames(dataLavg) <- levels(group)
colnames(dataRavg) <- levels(group)

data = read.delim("/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230906_Plot/1.Cellchat_plot/Fib.txt")


dataLavg = as.data.frame(dataLavg)
dataLavg$geneID = geneL
dataLavg =dataLavg %>% select("geneID","myCAF","matrix_CAF","iCAF","EMT_like_CAF","Cancer_pre")

dataRavg = as.data.frame(dataRavg)
dataRavg$geneID = geneR
dataRavg = dataRavg%>% select("geneID","Cancer_pre")

Ligand_cell = c("myCAF","matrix_CAF","iCAF","EMT_like_CAF")
Receptor_cell = c("Cancer_pre")
Ligand = data$Ligand
Receptor = data$Receptor

selected_ligand_expr = data.frame()
k = 1
for (i in Ligand_cell){
for (j in Ligand){
print(i)
print(j)
selected_ligand_expr[k,1] = i
selected_ligand_expr[k,2] = j
expr  <- dataLavg %>% select(geneID,i) %>% filter(geneID == j)  %>% select(-geneID)
selected_ligand_expr[k,3] = expr[1,1]
k = k+1
print(k)
}
}
colnames(selected_ligand_expr)[1]  ="source"
colnames(selected_ligand_expr)[2]  ="ligand"
colnames(selected_ligand_expr)[3]  ="expr"


selected_receptor_expr = data.frame()
k = 1
for (i in Receptor_cell){
for (j in Receptor){
print(i)
print(j)
selected_receptor_expr[k,1] = i
selected_receptor_expr[k,2] = j
expr  <- dataRavg %>% select(geneID,i) %>% filter(geneID == j)  %>% select(-geneID)
selected_receptor_expr[k,3] = expr[1,1]
k = k+1
print(k)
}
}
colnames(selected_receptor_expr)[1]  ="target"
colnames(selected_receptor_expr)[2]  ="receptor"
colnames(selected_receptor_expr)[3]  ="expr"

data$ligand_receptor = paste(data$Ligand, data$Receptor, sep = "_")
head(cellchat@LR$LRsig)
test = subset(cellchat@LR$LRsig, interaction_name %in% data$ligand_receptor)
test = test %>% select(interaction_name)

p1 <- netVisual_bubble(cellchat, remove.isolate = FALSE, return.data = T,sources.use = c(1,2,3,4),targets.use = 11, pairLR.use = test, font.size = 12)
communication = p1$communication

communication$source_ligand = paste(communication$source, communication$ligand)
selected_ligand_expr$source_ligand = paste(selected_ligand_expr$source, selected_ligand_expr$ligand)
for (i in communication$source_ligand){
communication[which(communication$source_ligand ==i),"source_ligand_expr"] = selected_ligand_expr[which(selected_ligand_expr$source_ligand ==i),"expr"]
}

communication$target_receptor = paste(communication$target, communication$receptor)
selected_receptor_expr$target_receptor= paste(selected_receptor_expr$target,selected_receptor_expr$receptor)
for (i in unique(communication$target_receptor)){
communication[which(communication$target_receptor ==i),"target_receptor_expr"] =unique(selected_receptor_expr[which(selected_receptor_expr$target_receptor ==i),"expr"])
}

communication$exprs = communication$source_ligand_expr + communication$target_receptor_expr
communication$LOG2exprs = log2(communication$exprs/2+0.0001) 


df = communication

source("/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230906_Plot/1.Cellchat_plot/visualization.R")

p2 =  netVisual_bubble_plot(cellchat,communication_df = communication, remove.isolate = FALSE, return.data = T,sources.use = c(1,2,3,4),targets.use = 11, pairLR.use = test, font.size = 12)
output_dir = "/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230906_Plot/1.Cellchat_plot"
ggsave(file.path(output_dir, "Fib.pdf"),plot = p2$gg.obj,height = 7, width = 7, limitsize = FALSE)

communication$pval <- gsub("1", "p > 0.05", communication$pval)
communication$pval <- gsub("2", "0.01 < p < 0.05", communication$pval)
communication$pval <- gsub("3", "p < 0.01", communication$pval)

write.table(communication, file.path(output_dir, "Fib_cell_comm_annotation.xls"), sep = "\t", col.names = T, row.names = F, quote = F)

########另一个图
pairLRsig <- cellchat@LR$LRsig
  group <- cellchat@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
nLR <- nrow(pairLRsig)
numCluster <- nlevels(group)


data = as.matrix(cellchat@data.signaling)
  data.use <- data/max(data)
  nC <- ncol(data.use)


data.use.avg <- aggregate(t(data.use), list(group), FUN =triMean)


data.use.avg <- t(data.use.avg[,-1])
colnames(data.use.avg) <- levels(group)


complex_input <- cellchat@DB$complex
cofactor_input <- cellchat@DB$cofactor
dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "A")
dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "I")
dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor


rownames(dataLavg) = geneL
rownames(dataRavg) = geneR


colnames(dataLavg) <- levels(group)
colnames(dataRavg) <- levels(group)


data = read.delim("/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230906_Plot/1.Cellchat_plot/PMC.txt")

dataLavg = as.data.frame(dataLavg)
dataLavg$geneID = geneL
dataLavg =dataLavg %>% select("geneID","C2_LIPF","C3_REG4","C4_GAST","C5_MSLN","C6_XIST")


dataRavg = as.data.frame(dataRavg)
dataRavg$geneID = geneR
dataRavg = dataRavg%>% select("geneID","Cancer_pre")


Ligand_cell = c("C2_LIPF","C3_REG4","C4_GAST","C5_MSLN","C6_XIST")
Receptor_cell = c("Cancer_pre")
Ligand = data$Ligand
Receptor = data$Receptor


selected_ligand_expr = data.frame()
k = 1
for (i in Ligand_cell){
for (j in Ligand){
print(i)
print(j)
selected_ligand_expr[k,1] = i
selected_ligand_expr[k,2] = j
expr  <- dataLavg %>% select(geneID,i) %>% filter(geneID == j)  %>% select(-geneID)
selected_ligand_expr[k,3] = expr[1,1]
k = k+1
print(k)
}
}
colnames(selected_ligand_expr)[1]  ="source"
colnames(selected_ligand_expr)[2]  ="ligand"
colnames(selected_ligand_expr)[3]  ="expr"



selected_receptor_expr = data.frame()
k = 1
for (i in Receptor_cell){
for (j in Receptor){
print(i)
print(j)
selected_receptor_expr[k,1] = i
selected_receptor_expr[k,2] = j
expr  <- dataRavg %>% select(geneID,i) %>% filter(geneID == j)  %>% select(-geneID)
selected_receptor_expr[k,3] = expr[1,1]
k = k+1
print(k)
}
}
colnames(selected_receptor_expr)[1]  ="target"
colnames(selected_receptor_expr)[2]  ="receptor"
colnames(selected_receptor_expr)[3]  ="expr"


data$ligand_receptor = paste(data$Ligand, data$Receptor, sep = "_")
head(cellchat@LR$LRsig)
test = subset(cellchat@LR$LRsig, interaction_name %in% data$ligand_receptor)
test = test %>% select(interaction_name)

p <- netVisual_bubble(cellchat, remove.isolate = FALSE, return.data = T,sources.use = c(6,7,8,9,10),targets.use = 11)
communication = p$communication


communication$source_ligand = paste(communication$source, communication$ligand)
selected_ligand_expr$source_ligand = paste(selected_ligand_expr$source, selected_ligand_expr$ligand)
for (i in communication$source_ligand){
communication[which(communication$source_ligand ==i),"source_ligand_expr"] = selected_ligand_expr[which(selected_ligand_expr$source_ligand ==i),"expr"]
}


communication$target_receptor = paste(communication$target, communication$receptor)
selected_receptor_expr$target_receptor= paste(selected_receptor_expr$target,selected_receptor_expr$receptor)
for (i in unique(communication$target_receptor)){
communication[which(communication$target_receptor ==i),"target_receptor_expr"] =unique(selected_receptor_expr[which(selected_receptor_expr$target_receptor ==i),"expr"])
}


communication$exprs = communication$source_ligand_expr + communication$target_receptor_expr
communication$LOG2exprs = log2(communication$exprs/2+0.0001)




df = communication


source("/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230906_Plot/1.Cellchat_plot/visualization.R")


p =  netVisual_bubble_plot(cellchat,communication_df = communication, remove.isolate = FALSE, return.data = T,sources.use = c(6,7,8,9,10),targets.use = 11, pairLR.use = test, font.size = 12)
output_dir = "/public/scRNA_works/works/liuxuan/Project/scRNA/DZOE2023020735_Human/houxu20230906_Plot/1.Cellchat_plot"
ggsave(file.path(output_dir, "PMC.pdf"),plot = p$gg.obj,height = 7, width = 7, limitsize = FALSE)


communication$pval <- gsub("1", "p > 0.05", communication$pval)
communication$pval <- gsub("2", "0.01 < p < 0.05", communication$pval)
communication$pval <- gsub("3", "p < 0.01", communication$pval)


write.table(communication, file.path(output_dir, "PMC_cell_comm_annotation.xls"), sep = "\t", col.names = T, row.names = F, quote = F)

##################挑几个基因测试与表达量趋势是否一致###############################
ID = rownames(data.use)
avgexpr = as.data.frame(data.use)
avgexpr$geneID = ID
avgexpr =avgexpr %>% gather( key = "Barcodes", value = "values", -"geneID")

meta = cellchat@meta %>% select(ident) %>%rownames_to_column("Barcodes")

avgexpr = avgexpr %>%left_join(meta, by = "Barcodes")
avgexpr = avgexpr %>% dplyr::group_by( ident, geneID) %>% dplyr::mutate(exprs = log2(mean(values)+0.0001) )

TEST = avgexpr %>% select(geneID, ident, exprs) %>%  dplyr::distinct()
TEST%>%filter(geneID == "COL1A1"&ident == "iCAF")
#-0.557
TEST%>%filter(geneID == "SDC4"&ident == "Cancer_pre")
#-4.22
TEST%>%filter(geneID == "COL4A5"&ident == "iCAF")
#-5.07


TEST%>%filter(geneID == "LAMB3"&ident == "C5_MSLN")
#-2.63
TEST%>%filter(geneID == "LAMB2"&ident == "C5_MSLN")
#-3.74
TEST%>%filter(geneID == "MDK"&ident == "C5_MSLN")
#-3.33
