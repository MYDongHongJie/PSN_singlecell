suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("magrittr"))
library(reshape2)
### 调色板
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

print("read h5 and csv")
data_ob = OESingleCell::ReadX(input = '/public/scRNA_works/works/chenhaoruo/project/HT2021-18054-6-Human/20221116/new_celltype_T.h5seurat', informat = 'h5seurat',verbose = F)


plotdata=data_ob@meta.data

plotdata$group2 = ""

plotdata[which(plotdata$sampleid %in% c("P435_T1","P435_T2","P435_N")),"group2"]="P435"
plotdata[which(plotdata$sampleid %in% c("P436_T1","P436_T2","P436_N")),"group2"]="P436"
plotdata[which(plotdata$sampleid %in% c("P437_T1","P437_T2","P437_N")),"group2"]="P437"
plotdata[which(plotdata$sampleid %in% c("P438_T1","P438_T2","P438_N","P438_T3","P438_T4")),"group2"]="P438"
plotdata[which(plotdata$sampleid %in% c("p439_T1","p439_T2","P439_N")),"group2"]="P439"
plotdata$sampleid=as.character(plotdata$sampleid)
plotdata[which(plotdata$sampleid %in% c("p439_T1")),"sampleid"]="P439_T1"
plotdata[which(plotdata$sampleid %in% c("p439_T2")),"sampleid"]="P439_T2"

plotdata$sampleid=factor(plotdata$sampleid,levels=c("P435_T1","P436_T1","P437_T1","P438_T1","P439_T1","P435_T2","P436_T2","P437_T2","P438_T2","P439_T2","P435_N","P436_N","P437_N","P438_N","P439_N","P438_T3","P438_T4"))
plotdata$group2=factor(plotdata$group2,levels=c("P435","P436","P437","P438","P439"))
table(plotdata$sampleid)
table(plotdata$group2)

    cluster_ids <- plotdata[,"new_celltype"]
    sample_ids <- plotdata[, "sampleid"]
    counts <- table(cluster_ids, sample_ids)
    df <- melt(t(round(t(counts)/colSums(counts) * 100, 2)), 
        varnames = c("new_celltype", "sampleid"), value.name = "freq")

df[which(df$sampleid %in% c("P435_T1","P435_T2","P435_N")),"group2"]="P435"
df[which(df$sampleid %in% c("P436_T1","P436_T2","P436_N")),"group2"]="P436"
df[which(df$sampleid %in% c("P437_T1","P437_T2","P437_N")),"group2"]="P437"
df[which(df$sampleid %in% c("P438_T1","P438_T2","P438_N","P438_T3","P438_T4")),"group2"]="P438"
df[which(df$sampleid %in% c("P439_T1","P439_T2","P439_N")),"group2"]="P439"
num=1
for(i in unique(df$group2)){
    df[which(df$group2==i),"x"]=num
    num=num+1
}

dfxx = dplyr::arrange(df,x)
num=0
for (sample in unique(dfxx[,"sampleid"])){
    dfxx[which(dfxx[,"sampleid"]==sample),"xx"]=num+dfxx[which(dfxx[,"sampleid"]==sample),"x"]
    num=num+1
}


plot = ggplot(dfxx,aes(x=xx,y=freq,fill=new_celltype))+
  geom_bar(stat="identity",position = "stack")+
#   facet_wrap(~group2)+
  scale_x_continuous(breaks = unique(dfxx$xx),
                     labels = unique(dfxx$sampleid))+
  scale_fill_manual(values = CustomCol2(1:length(unique(dfxx$new_celltype)))) +
 scale_y_continuous(expand = c(0, 
        0), labels = seq(0, 100, 25)) + labs(x = NULL, y = "Proportion [%]") + 
        theme_bw() + theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
            color = NA), panel.border = element_blank(), axis.ticks.x = element_blank(), 
        axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 90, 
            hjust = 1, vjust = 0.5))

ggsave(plot,file="plot1.png",width = 10,height=5)
ggsave(plot,file="plot1.pdf",width = 10,height=5)


    cluster_ids <- plotdata[,"new_celltype"]
    sample_ids <- plotdata[, "sampleid"]
    counts <- table(cluster_ids, sample_ids)
    df <- melt( counts   , 
        varnames = c("new_celltype", "sampleid"), value.name = "counts")
df[which(df$sampleid %in% c("P435_T1","P435_T2","P435_N")),"group2"]="P435"
df[which(df$sampleid %in% c("P436_T1","P436_T2","P436_N")),"group2"]="P436"
df[which(df$sampleid %in% c("P437_T1","P437_T2","P437_N")),"group2"]="P437"
df[which(df$sampleid %in% c("P438_T1","P438_T2","P438_N","P438_T3","P438_T4")),"group2"]="P438"
df[which(df$sampleid %in% c("P439_T1","P439_T2","P439_N")),"group2"]="P439"
num=1
for(i in unique(df$group2)){
    df[which(df$group2==i),"x"]=num
    num=num+1
}

dfxx = dplyr::arrange(df,x)
num=0
for (sample in unique(dfxx[,"sampleid"])){
    dfxx[which(dfxx[,"sampleid"]==sample),"xx"]=num+dfxx[which(dfxx[,"sampleid"]==sample),"x"]
    num=num+1
}

plot = ggplot(dfxx,aes(x=xx,y=counts,fill=new_celltype))+
  geom_bar(stat="identity",position = "stack")+
#   facet_wrap(~group2)+
  scale_x_continuous(breaks = unique(dfxx$xx),
                     labels = unique(dfxx$sampleid))+
  scale_fill_manual(values = CustomCol2(1:length(unique(dfxx$new_celltype)))) +
 labs(x = NULL, y = "Counts") + 
        theme_bw() + theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
            color = NA), panel.border = element_blank(), axis.ticks.x = element_blank(), 
        axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 90, 
            hjust = 1, vjust = 0.5))

ggsave(plot,file="plot2.png",width = 10,height=5)
ggsave(plot,file="plot2.pdf",width = 10,height=5)
