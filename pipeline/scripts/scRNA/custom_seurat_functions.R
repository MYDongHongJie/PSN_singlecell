plot.clusters.group = function (data = seurat_data,clusters = seurat_clusters,
                                    group = orig.ident,group2 = orig.ident,widths = c(3,1),log =TRUE,
                                    legend.title = "Group",color = 1,xlab = ""){
 ## take an integrated Seurat object, plot distributions over orig.ident
 library(Seurat)
 library(patchwork)
 library(ggplot2)
 library(reshape2)
 library(RColorBrewer)
 library(paletteer)
 mytheme = theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 10,color ="black"),
                 axis.text = element_text(size=10,color = "black"),
                 #axis.line = element_line(color = "black"),
                 #axis.ticks = element_line(color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 panel.grid=element_blank(), # 去网格线
                 # legend.position = "none",
                 legend.text = element_text(size=8),
                 legend.title= element_text(size= 8),
                 # axis.text.x = element_text(angle = 45, hjust=1, vjust=1)
)

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",                                                 
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",                                                 
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",                                                 
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",                                                 
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")                                                 
  return(my_palette[n])
}
colors2use = CustomCol2(1:length(unique(seurat_ob@meta.data[,"sampleid"])))
print(colors2use)
 
 count_table <- table(seurat_ob@meta.data[,"new_celltype"], seurat_ob@meta.data[,"sampleid"])
 count_mtx <- as.data.frame.matrix(count_table)
 count_mtx$cluster <- rownames(count_mtx)
 melt_mtx <- melt(count_mtx)
 melt_mtx$cluster <- as.factor(melt_mtx$cluster)
 
 cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
 
 if("0" %in% cluster_size$cluster){
   sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
}else{
   sorted_labels <- paste(cluster_size$cluster[order(cluster_size$value)])
}
 
 cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
 melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
 
 colnames(melt_mtx)[2] <- "dataset"
   
################### p1  
 if(log){
   p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") +
     theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("") + mytheme
}else{
   p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") +
     theme_bw() + xlab("") + ylab("") + mytheme
}
################### p2    
 ########################### color 1  
 if(color==1){
 if(length(unique(melt_mtx$dataset)) < 21){
   p2 <- ggplot(data,aes(x=cluster,y=value,fill=dataset)) +
     geom_bar(position="fill", stat="identity",) + theme_bw() + coord_flip() +
     scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))+
     ylab("") + xlab("") +
     theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
     scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
}else{
   warning("The color limit is <21")
   p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
   geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
     ylab("") + xlab("") +
   theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
   scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
}
}
 ########################### color 2
 if(color==2){
 if(length(unique(melt_mtx$dataset)) < 9){
   p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
     geom_bar(position="fill", stat="identity",) + theme_bw() + coord_flip() +
     scale_fill_manual(values = colors2use)+
     ylab("") + xlab("") +
     theme(legend.position="top") + guides(fill = guide_legend(title = "new_celltype")) +
     scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
}else{
   warning("The color limit is <9")
   p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
     geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
     scale_fill_manual(values = colors2use)+
     ylab("") + xlab("") +
     theme(legend.position="top") + guides(fill = guide_legend(title = "Group")) +
     scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
}
}
##############p3
    count_table <- table(data@meta.data[,clusters], data@meta.data[,group2])
    count_mtx <- as.data.frame.matrix(count_table)
    count_mtx$cluster <- rownames(count_mtx)
    melt_mtx <- melt(count_mtx)
    melt_mtx$cluster <- as.factor(melt_mtx$cluster)
    
    cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
    
    if("0" %in% cluster_size$cluster){
    sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
    }else{
    sorted_labels <- paste(cluster_size$cluster[order(cluster_size$value)])
    }
    
    cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
    melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"

   p3 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
     geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
     ylab("") + xlab("") +
     theme(legend.position="top") + guides(fill = guide_legend(title = "")) +
     scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme


 wrap_plots(ncol = 3,p2,p3,p1,widths = widths)
}