# install.packages("sscVis")
# library("sscVis") #报错 采用phetamap出图
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
#设置图片输出目录
out.prefix <- "./Fig1"
# 读入数据：
data_ob = OESingleCell::ReadX(input = '/public/scRNA_works/works/chenhaoruo/project/HT2021-18054-6-Human/20220822/result/new_celltype_total/new_celltype_total.h5seurat', informat = 'h5seurat',verbose = F)

meta <- data_ob@meta.data

  library(data.table)

##
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

meta=data_ob@meta.data
tab.S0 <- as.data.frame(table(meta[,c("sampleid", "group","new_celltype")]))
tab.S0 <- tidyr::spread(data = tab.S0,key = sampleid,value = Freq)
head(tab.S0)
for (i in 3:19) {
 tab.S0[,i] = tab.S0[,i] / colSums(tab.S0[,c(3:19)])[i-2] *100
}
head(tab.S0)


tab.S0_New =  melt(tab.S0)
head(tab.S0_New)
colnames(tab.S0_New)=c("group","celltype","sample","prop")  #设置行名
tab.S0_New$celltype = as.character(tab.S0_New$celltype)
head(tab.S0_New)
#去0
tab.S0_New = tab.S0_New[which(tab.S0_New$prop!=0),]


# for (i in unique(tab.S0_New$celltype)) {
#   plotdata=tab.S0_New[which(tab.S0_New$celltype==i),]

i="new_celltype"
plotdata=tab.S0_New
    ymax=max(plotdata$prop)
  p3 = ggplot(plotdata,aes(x = group, y = prop)) +
  geom_bar(stat="summary",fun=mean,width = NULL,colour="black",aes(fill=group))+
  stat_summary(geom = "errorbar", colour = "black",
                width = 0.5,position = position_dodge( .9)) +
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  #  geom_point(pch = 20,size = 4,aes(x = group, y = prop , colour=group)) + # ,position = position_dodge(0.9)
  geom_point(pch = 21,size = 4,aes(fill=group)) + # ,position = position_dodge(0.9)
  #  scale_colour_manual(values = c("#1CB4B8", "#EB7369"))+
  labs(y="RO/E",x= NULL,title=NULL) +
   facet_grid(cols =  vars(celltype),scales="free") +
  #scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))+
  theme_bw()+theme(plot.title = element_text(size = 15,color="black",hjust = 0.5),
                      axis.title = element_text(size = 12,color ="black"),
                      axis.text = element_text(size=12,color = "black"),
                      #axis.line = element_line(color = "black"),
                      #axis.ticks = element_line(color = "black"),
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      panel.grid=element_blank(), # 去网格线
                      legend.position = "right",
                      legend.text = element_text(size=12),
                      legend.title= element_text(size= 12),
                    axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) #+ 
                    # stat_compare_means(aes(group = group),
                    #   label = "p.format",
                    #   method = "wilcox.test",
                    #   label.y=ymax+3,
                    #   label.x = 1.5,
                    #   hide.ns = F) +
                    # stat_compare_means(aes(group = group),
                    #   label = "p.signif",
                    #   method = "wilcox.test",
                    #   label.y=ymax+5,
                    #   label.x = 1.5,
                    #   hide.ns = F)
  ggsave(plot=p3,paste0(i,"_ROE_plot_wilcox_test.png"),dpi=1000,width=14,height=10)
  ggsave(plot=p3,paste0(i,"_ROE_plot_wilcox_test.pdf"),width=14,height=10)
# }



                  #  geom_signif(comparisons = list(c("GS","IM"),
                  #                c("GS","GC"),
                  #                c("IM","GC")),step_increase = 0.1,map_signif_level = F,test = t.test)

library(ggraph)
library(ggforce)
library(tidyr)
for (i in unique(tab.S0_New$celltype)) {
    # i="Mast_cells"
    plotdata=tab.S0_New[which(tab.S0_New$celltype==i),]
    plotdata$sample=as.character(plotdata$sample)
    temp=unlist(strsplit(plotdata$sample,split = "_"))
    plotdata$sample2=temp[seq(1,length(temp),by=2)]
    plotdata[which(plotdata$sample2=="p439"),"sample2"]="P439"
    plotdata=plotdata[,c("group","prop","sample2")]
    tempdata=plotdata[which(plotdata$group=="T"),] 
    for(r in rownames(tempdata)){
      if(tempdata[r,"sample2"] %in% plotdata[which(plotdata$group=="N"),"sample2"]){
        tempdata[r,"n_prop"]=plotdata[which(plotdata$group=="N" & plotdata$sample2==tempdata[r,"sample2"]),"prop"]
      }
    }
    line= data.frame(
      x = rep("N",dim(tempdata)[1]),
      y= tempdata[,"n_prop"],
      xend = rep("T",dim(tempdata)[1]),
      yend = tempdata[,"prop"]
    )
    # plotdata=plotdata[order(plotdata$group),]
  # i="new_celltype"
  # plotdata=tab.S0_New
      ymax=max(plotdata$prop)
    p3 = ggplot(plotdata,aes(x = group, y = prop)) +
    geom_boxplot(colour="black",aes(fill=group))+
    # stat_summary(geom = "errorbar", colour = "black",
    #               width = 0.5,position = position_dodge( .9)) +
    scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
    #  geom_point(pch = 20,size = 4,aes(x = group, y = prop , colour=group)) + # ,position = position_dodge(0.9)
    geom_point(pch = 21,size = 4,aes(fill=group)) + # ,position = position_dodge(0.9)
    geom_diagonal(data=line,aes(x = x,y = y, xend = xend, yend = yend),strength = 0.75) +
    #  scale_colour_manual(values = c("#1CB4B8", "#EB7369"))+
    labs(y="RO/E",x= NULL,title=NULL) +
    #  facet_grid(cols =  vars(celltype),scales="free") +
    #scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))+
    theme_bw()+theme(plot.title = element_text(size = 15,color="black",hjust = 0.5),
                        axis.title = element_text(size = 12,color ="black"),
                        axis.text = element_text(size=12,color = "black"),
                        #axis.line = element_line(color = "black"),
                        #axis.ticks = element_line(color = "black"),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid=element_blank(), # 去网格线
                        legend.position = "right",
                        legend.text = element_text(size=12),
                        legend.title= element_text(size= 12),
                      axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) #+ 
                      # stat_compare_means(aes(group = group),
                      #   label = "p.format",
                      #   method = "wilcox.test",
                      #   label.y=ymax+3,
                      #   label.x = 1.5,
                      #   hide.ns = F) +
                      # stat_compare_means(aes(group = group),
                      #   label = "p.signif",
                      #   method = "wilcox.test",
                      #   label.y=ymax+5,
                      #   label.x = 1.5,
                      #   hide.ns = F)
    ggsave(plot=p3,paste0(i,"_ROE_plot_wilcox_test.png"),dpi=1000,width=4,height=5)
    ggsave(plot=p3,paste0(i,"_ROE_plot_wilcox_test.pdf"),width=4,height=5)
}



