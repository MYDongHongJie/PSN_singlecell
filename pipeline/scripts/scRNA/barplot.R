library('ggplot2')
setwd("D:/นคื๗/HT2020-19541/")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

gene_data<-read.csv("gene_data.csv",sep = ",")
meta_data<-read.csv("meta_data.csv",sep = ",")
head(gene_data)
colnames(gene_data)
for (i in colnames(gene_data)) {
  meta_data[,i]=""
  meta_data[rownames(gene_data),i]=gene_data[rownames(gene_data),i]
  meta_data[,i]=as.numeric(meta_data[,i])
}
head(meta_data)

w<-which(meta_data$celltype_Ep=="")
meta_data<-meta_data[-w,]
meta_data$celltype_Ep<-factor(meta_data$celltype_Ep,levels = sort(unique(meta_data$celltype_Ep)))
levels(meta_data$group2)<-c("GS","IM","GC")
levels(meta_data$group3)<-c("GS_HP_N","GS_HP_P", "IM_HP_N", "IM_HP_P","GC_HP_N", "GC_HP_P")

for (i in colnames(gene_data[2:8])) {
  for (g in c("group","group2","group3")) {
  group <- ggplot(data = meta_data, mapping=aes(x = interaction(orig.ident,celltype_Ep), 
                                                y = meta_data[,i], fill=celltype_Ep))+
    geom_col(position = 'dodge')+geom_bar(stat = "identity", position = "dodge")+
    facet_wrap(~meta_data[,g])+labs(x= "celltype",y = i)+
    scale_fill_manual(values=CustomCol2(1:length(unique(meta_data$celltype_Ep))))
  
  ggsave(paste0(i,"_",g,"_barplot.png"),group,width = 15,height = 7)
  ggsave(paste0(i,"_",g,"_barplot.pdf"),group)
  }
}

















