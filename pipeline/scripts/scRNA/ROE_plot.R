#计算ROE值
#函数来源
#https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Tissue_distribution_analysis.R
#文献来源
#https://www.sciencedirect.com/science/article/pii/S0092867421000106#figs3


library(dplyr)
options(stringsAsFactors=FALSE)
library(reticulate)


ROIE <- function(crosstab){
## Calculate the Ro/e value from the given crosstab
##
## Args:
#' @crosstab: the contingency table of given distribution
##
## Return:
## The Ro/e matrix
rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
#按照细胞类型计算细胞总数
rowsum.matrix[,1] <- rowSums(crosstab)
colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
colsum.matrix[1,] <- colSums(crosstab)
allsum <- sum(crosstab)
roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
row.names(roie) <- row.names(crosstab)
colnames(roie) <- colnames(crosstab)
return(roie)
}


divMatrix <- function(m1, m2){
## Divide each element in turn in two same dimension matrixes
##
## Args:
#' @m1: the first matrix
#' @m2: the second matrix
##
## Returns:
## a matrix with the same dimension, row names and column names as m1.
## result[i,j] = m1[i,j] / m2[i,j]
dim_m1 <- dim(m1)
dim_m2 <- dim(m2)
if( sum(dim_m1 == dim_m2) == 2 ){
div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
row.names(div.result) <- row.names(m1)
colnames(div.result) <- colnames(m1)


for(i in 1:dim_m1[1]){
for(j in 1:dim_m1[2]){
div.result[i,j] <- m1[i,j] / m2[i,j]
}
}
return(div.result)
}
else{
warning("The dimensions of m1 and m2 are different")
}
}




summary <- table(seurat_ob@meta.data[,c('celltype','group')])
roe <- as.data.frame(ROIE(summary))




#roe就是绘图使用的矩阵
#先做宽变长数据转换
library(reshape2)
df = roe
df$celltype = rownames(df)
df = melt(df,variable.name='group',value.name='ROE',id.vars = "celltype")
setwd("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230629/1.New_celltype/ROE")
df$celltype = as.factor(df$celltype)
df$group = as.factor(df$group)
nlevel=length(unique(df$group))


CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}


p = ggplot( df,  aes(x = celltype, y = ROE,group = group,colour= group)) +geom_line(linetype = "dashed", size = 1.5)+geom_point()+theme_classic()+theme(axis.text.x= element_text(
                size=9,angle=90,hjust=1,vjust = 0.5))+labs(y = "Ro/e")+scale_color_manual(values=CustomCol2(1:nlevel))


ggsave("ROE.png",p)
ggsave("ROE.pdf",p)


write.table(df,file = file.path("./ROE.xls"),
            col.names =T,row.names = F,sep = "\t",quote=F)
