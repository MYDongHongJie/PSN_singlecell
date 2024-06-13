
library(ggplot2)
score = read.delim("/public/scRNA_works/works/guokaiqi/project/scRNA/Wangzheng/DOE20224843/houxu-20230214/geneset/geneset_visualization/geneset.addmodeulescore.xls")
plot_data = reshape2::melt(score, id.vars=c('Barcode',"clusters", "sampleid","group"), variable.name = 'geneset', value.name = 'data')
plot_data$group = gsub("1","",plot_data$group)
plot_data$group = gsub("2","",plot_data$group)
          # # vlnplot
          # vln[[geneset]] = ggplot(score,aes_string(x=groupby,y=geneset,fill=groupby))+
            # geom_violin()+
            # labs(x="Clusters",y=paste0(geneset,"_score")) +
            # theme(panel.grid.major =element_blank(),
                  # panel.grid.minor = element_blank(),
                  # panel.background = element_blank(),
                  # axis.line = element_line(color = "black"))+
            # scale_fill_manual(values=CustomCol2(1:length(unique(metadata[,groupby]))))
          # ggsave(file.path(output_dir,paste0(geneset,"_score_vlnplot.pdf")),plot=vln[[geneset]],width=8)
          # ggsave(file.path(output_dir,paste0(geneset,"_score_vlnplot.png")),plot=vln[[geneset]],width=8)

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

output_dir = "/public/scRNA_works/works/guokaiqi/project/scRNA/Wangzheng/DOE20224843/houxu-20230214/geneset/geneset_visualization"
p = ggplot(plot_data,aes_string(x="sampleid",y='data',fill="sampleid")) + geom_violin()+
            labs(x="Clusters",y=paste0("geneset_score")) + facet_grid(geneset~clusters)+
            theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black"))+
            scale_color_manual(values=CustomCol2(1:6))+scale_fill_manual(values=CustomCol2(1:6))
          ggsave(file.path(output_dir,paste0("geneset_score_vlnplot_by_sampleid.pdf")),plot=p,width=10,height=10)
          ggsave(file.path(output_dir,paste0("geneset_score_vlnplot_by_sampleid.png")),plot=p,width=10,height=10)

p = ggplot(plot_data,aes_string(x="sampleid",y='data',fill="sampleid")) + geom_violin()+
            labs(x="Clusters",y=paste0("geneset_score")) + facet_grid(geneset~clusters, scales = "free_y")+
            theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black"))+
            scale_color_manual(values=CustomCol2(1:6))+scale_fill_manual(values=CustomCol2(1:6))
          ggsave(file.path(output_dir,paste0("geneset_score_vlnplot_free_y_by_sampleid.pdf")),plot=p,width=10,height=10)
          ggsave(file.path(output_dir,paste0("geneset_score_vlnplot_free_y_by_sampleid.png")),plot=p,width=10,height=10)

p = ggplot(plot_data,aes_string(x="group",y='data',fill="group")) + geom_violin()+
            labs(x="Clusters",y=paste0("geneset_score")) + facet_grid(geneset~clusters)+
            theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black"))+
            scale_color_manual(values=CustomCol2(1:6))+scale_fill_manual(values=CustomCol2(1:6))
          ggsave(file.path(output_dir,paste0("geneset_score_vlnplot_by_group.pdf")),plot=p,width=10,height=10)
          ggsave(file.path(output_dir,paste0("geneset_score_vlnplot_by_group.png")),plot=p,width=10,height=10)

p = ggplot(plot_data,aes_string(x="group",y='data',fill="group")) + geom_violin()+
            labs(x="Clusters",y=paste0("geneset_score")) + facet_grid(geneset~clusters, scales = "free_y")+
            theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black"))+
            scale_color_manual(values=CustomCol2(1:6))+scale_fill_manual(values=CustomCol2(1:6))
          ggsave(file.path(output_dir,paste0("geneset_score_vlnplot_free_y_by_group.pdf")),plot=p,width=10,height=10)
          ggsave(file.path(output_dir,paste0("geneset_score_vlnplot_free_y_by_group.png")),plot=p,width=10,height=10)

