library(dplyr)
library(SingleCellSignalR)
library(circlize)
library(grDevices)
source("visualize_interactions.R")

#导入cellchat的表格
df<- read.delim("communication.xls",sep='\t')
data = df %>% select(source,target,ligand,receptor,annotation,prob,dataset)
source_cell=c("Macrophages","NK")
target_cell=c("Endothelia_cells","CD.PC_DCT","PT")
output_dir <- normalizePath(".")

#设置环的颜色
my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")

new_celltype=c("B_cell","Endothelia_cells","Macrophages","FIB","CD.PC_DCT","PT","Neutrophils","T_cell","NK")
new_celltype_col=my_palette[1:length(new_celltype)]
names(new_celltype_col)=new_celltype

#转换成SingleCellSignalR绘图的格式
for (group in unique(as.vector(data$dataset))){
    if (! file.exists(group)) {
        dir.create(group, recursive= T )
    }
    print(paste0("绘制group:",group))
    data_group <- subset(data,dataset==group)
    signal=list()
    for (i in source_cell) {
        for (j in target_cell){
            if(i!=j){
                name1=i
                name2=j
                signal[[paste0(i,"-",j)]]= subset(data_group,source==i & target==j) %>% select(ligand,receptor,annotation,prob)
                colnames(signal[[paste0(i,"-",j)]])=c(i,j,"annotation","prob")
            }else{
                next
            }
         }
     }
     for (name in names(signal)) {
         if(dim(signal[[name]])[1]>1){
             pdf(file.path(group,paste0(name,'_LR_prob.pdf')))
             visualize_interactions(signal, show.in = name, grid.col=new_celltype_col)
             dev.off()
             write.table(signal[[name]],file.path(group,paste0(name,'_LR_prob.xls')),sep='\t',quote=F,row.names=F)
         }else{
             write.table(signal[[name]],file.path(group,paste0(name,'_LR_prob.xls')),sep='\t',quote=F,row.names=F)
         }
     }
     setwd(paste0(output_dir,"/",group))
     system("for i in  `ls *.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")
     setwd(output_dir)
}
