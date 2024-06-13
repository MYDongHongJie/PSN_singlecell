suppressPackageStartupMessages( library("dplyr") )
suppressPackageStartupMessages( library("ggplot2") )
suppressPackageStartupMessages( library("tidyr") )

Control="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2021-11-10/9.Commu/Cellphone/Control/cell_comm_annotation.xls"
Herts33="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2021-11-10/9.Commu/Cellphone/Herts33/cell_comm_annotation.xls"
LaSota="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2021-11-10/9.Commu/Cellphone/LaSota/cell_comm_annotation.xls"

data1 = read.table(Control, header = T, sep = "\t")
data2 = read.table(Herts33, header = T, sep = "\t")
data3 = read.table(LaSota, header = T, sep = "\t")
##### 筛选 top5 的配体对 #####
data1$id="Control"
data2$id="Herts33"
data3$id="LaSota"

screenvar="expr"
topn=5
opt=list()
opt$pvalue=0.05

all_data=rbind(data1,data2,data3)
## all_data <- all %>% group_by(receptor_cell,ligand_cell,id) %>% filter( pval < opt$pvalue ) %>% top_n(topn, !!ensym(screenvar))

tamp = all_data %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|")

# tamp2 = all %>% unite( "pair", receptor, ligand, sep = "|") 
# top5  <- tamp2 %>% filter( pair %in% tamp$pair)

# top5_2 <- separate(top5,"pair", into= c("receptor","ligand"),sep= "\\|")

#     这个受体配体对做下dotplot图，几个组合并出图，同一个细胞类型，不同组的
# cellphoneDB的结果是2021-11-10

# B作为受体，CD4+，CD8+，Fibroblast作为配体
output_dir="dotplot"
dir.create(output_dir)

# data <- all_data

Control="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2022-01-24/3.cellphome_dotplot/control_cellphonedb.xls"
Herts33="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2022-01-24/3.cellphome_dotplot/Herts33_cellphonedb.xls"
LaSota="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2022-01-24/3.cellphome_dotplot/control_cellphonedb.xls"

data1 = read.table(Control, header = T, sep = "\t")
data2 = read.table(Herts33, header = T, sep = "\t")
data3 = read.table(LaSota, header = T, sep = "\t")

data1 = read.table(Control, header = T, sep = "\t")
data2 = read.table(Herts33, header = T, sep = "\t")
data3 = read.table(LaSota, header = T, sep = "\t")
##### 筛选 top5 的配体对 #####
data1$id="Control"
data2$id="Herts33"
data3$id="LaSota"


sig_data=rbind(data1,data2,data3)
#data_output <- data %>% filter(ligand_cell %in% ligand_cell_list ) %>% filter(receptor_cell %in% receptor_cell_list )
tamp2 = sig_data %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|")

# tamp2 = sig_data %>% unite( "pair", receptor, ligand, sep = "|") 

sig  <- tamp %>% filter( pair %in% tamp2$pair)  %>% filter( clusters %in% tamp2$clusters) 

# data_output <- separate(top5,"pair", into= c("receptor","ligand"),sep= "\\|")

#data = data_output[,c("receptor","ligand","receptor_cell","ligand_cell","pval","expr","id")]
## %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|")
# filter_data = sig  %>% unite("clusters", clusters, id, sep = " " )

filter_data = sig  

palette = c("black", "blue", "yellow", "red")
xsize=10
xangle=45

suppressPackageStartupMessages( library("ggplot2") )
filter_data$pval[filter_data$pval==0] = 0.0009
  filter_data$expr[filter_data$expr==0] = 1
  filter_data$mean = as.numeric(log2(filter_data$expr))

  plot.data = filter_data[,c("pair","clusters","pval","mean","id")]
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean','id')

  my_palette <- colorRampPalette(palette, alpha=TRUE)(n=399)
  dotplot <- ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    facet_wrap(~id) + 
    guides(color = guide_colorbar(order = 1, title = 'Log2 mean (Molecule 1, Molecule 2)')) + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=xsize, colour = "black"),
          axis.text.x = element_text(size = 10, angle = xangle, hjust = 1),
          axis.text.y = element_text(size=xsize, colour = "black"),
          axis.title=element_blank(),
          text = element_text(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

ggdot = dotplot + theme(axis.text.y = element_text(face = "bold")) 


png_width = length(unique(plot.data$clusters))
png_height = length(unique(plot.data$pair))


write.table(sig_data, file.path(output_dir, "cell_comm_annotation.xls"),sep="\t",row.names=F,quote=F)
ggsave(file.path(output_dir, "cell_comm_dotplot.pdf"), limitsize = F,
        plot = ggdot, width = png_width*0.8+2, height = png_height*0.3-1)
ggsave(file.path(output_dir, "cell_comm_dotplot.png"), limitsize = F,
        plot = ggdot, width = png_width*0.8+2, height = png_height*0.3-1)



######################## 不展示配体对在其他细胞类型中的 展示 ######################################
output_dir="dotplot2"

dir.create(output_dir)
Control="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2022-01-24/3.cellphome_dotplot/control_cellphonedb.xls"
Herts33="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2022-01-24/3.cellphome_dotplot/Herts33_cellphonedb.xls"
LaSota="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-16870-chicken/Further_analysis_2022-01-24/3.cellphome_dotplot/control_cellphonedb.xls"

data1 = read.table(Control, header = T, sep = "\t")
data2 = read.table(Herts33, header = T, sep = "\t")
data3 = read.table(LaSota, header = T, sep = "\t")

data1 = read.table(Control, header = T, sep = "\t")
data2 = read.table(Herts33, header = T, sep = "\t")
data3 = read.table(LaSota, header = T, sep = "\t")
##### 筛选 top5 的配体对 #####
data1$id="Control"
data2$id="Herts33"
data3$id="LaSota"


sig_data=rbind(data1,data2,data3)
#data_output <- data %>% filter(ligand_cell %in% ligand_cell_list ) %>% filter(receptor_cell %in% receptor_cell_list )
#tamp2 = sig_data %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|")

# tamp2 = sig_data %>% unite( "pair", receptor, ligand, sep = "|") 

#sig  <- tamp %>% filter( pair %in% tamp2$pair)  %>% filter( clusters %in% tamp2$clusters) 

# data_output <- separate(top5,"pair", into= c("receptor","ligand"),sep= "\\|")

#data = data_output[,c("receptor","ligand","receptor_cell","ligand_cell","pval","expr","id")]
## %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|")
# filter_data = sig_data  %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|") %>% unite("clusters", clusters, id, sep = " " )

filter_data = sig_data  %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|") 

##%>% unite("clusters", clusters, id, sep = " " )

palette = c("black", "blue", "yellow", "red")
xsize=10
xangle=45

suppressPackageStartupMessages( library("ggplot2") )
filter_data$pval[filter_data$pval==0] = 0.0009
  filter_data$expr[filter_data$expr==0] = 1
  filter_data$mean = as.numeric(log2(filter_data$expr))

  plot.data = filter_data[,c("pair","clusters","pval","mean","id")]
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean','id')

  my_palette <- colorRampPalette(palette, alpha=TRUE)(n=399)
  dotplot <- ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    facet_wrap(~id) + 
    guides(color = guide_colorbar(order = 1, title = 'Log2 mean (Molecule 1, Molecule 2)')) + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=xsize, colour = "black"),
          axis.text.x = element_text(size = 10, angle = xangle, hjust = 1),
          axis.text.y = element_text(size=xsize, colour = "black"),
          axis.title=element_blank(),
          text = element_text(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

ggdot = dotplot + theme(axis.text.y = element_text(face = "bold")) 


png_width = length(unique(plot.data$clusters))
png_height = length(unique(plot.data$pair))

write.table(sig_data, file.path(output_dir, "cell_comm_annotation.xls"),sep="\t",row.names=F,quote=F)
ggsave(file.path(output_dir, "cell_comm_dotplot.pdf"), limitsize = F,
        plot = ggdot, width = png_width*0.8+2, height = png_height*0.3-1)
ggsave(file.path(output_dir, "cell_comm_dotplot.png"), limitsize = F,
        plot = ggdot, width = png_width*0.8+2, height = png_height*0.3-1)



