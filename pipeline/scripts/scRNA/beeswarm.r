suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggbeeswarm"))

option_list = list(
    make_option(c("-i", "--input"), type="character",
          help="The table of Common Transcriptome Differential Analysis"),
    make_option(c("-g", "--genelist"),type="character",
          help="The table of Marker Analysis"),
    make_option( c("--outdir","-o"),type="character", default = "./",
          help="the output directory of Clustering results." ),
    make_option( c("--cex","-c"),type="character", default = "0.5",
          help="the size of point" )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if ( is.null(opt$outdir) ){
    output_dir = getwd()
} else {
    if ( file.exists(opt$outdir) ){
        output_dir = opt$outdir
    } else {
        output_dir = opt$outdir
        dir.create(output_dir,recursive=T)
    }
}

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","grey","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
} 

marker <- read.delim(opt$genelist,sep='\t')
diff <- read.delim(opt$input,sep='\t')

diff_select <- diff %>% select(gene_id,log2FoldChange,p.value,Regulation,description,GO_id,GO_term,pathway) %>% rename(gene=gene_id) 

data <- left_join(marker , diff_select,by="gene") %>% dplyr::filter(!is.na(log2FoldChange))
num <-data %>% group_by(Regulation,cluster) %>% summarise(num = n())

data$link<- paste(data$Regulation,data$cluster,sep = "-")
num$link <- paste(num$Regulation,num$cluster,sep = "-")
for (i in unique(num$num) ){
    subset_ob= subset(num, num==i)
    data[which(data$link %in% as.vector(subset_ob$link)) ,"num"] = i
}

data_result <- data %>% select(!link)

group_name = unlist(strsplit( opt$input,"-",perl = T))
write.table(data_result,paste0(group_name[1],"-vs-",group_name[3],"_beeswarm.txt"), col.names=T, row.names=F,sep='\t',quote=F)

data_plot <- data %>% select(cluster,log2FoldChange,Regulation,num,link)
cluster <- unique(data_plot$cluster)
cluster_temp <- c(paste("Down", cluster,sep='-') , paste("Up", cluster,sep='-'))
NA_data <- data.frame(link2=cluster_temp,log2FoldChange=0 ,num=0,link=cluster_temp)
NA_data_separate <- separate(NA_data,col='link2',into=c("Regulation","cluster"),sep="-")
data_plot <- rbind(data[,c("Regulation","cluster","log2FoldChange","num","link")],anti_join(NA_data_separate,data[,c("Regulation","cluster","log2FoldChange","num","link")],by="link"))

Up <- subset(data_result,Regulation=='Up')
ylim = Up %>% subset(log2FoldChange != 'Inf') %>% select(log2FoldChange) %>% max() 
Down <- subset(data_result,Regulation=='Down')

data[which(data$log2FoldChange=='Inf' & data$Regulation == 'Up'),'log2FoldChange'] = 5E+300
data[which(data$log2FoldChange=='-Inf' & data$Regulation == 'Down'),'log2FoldChange'] = -5E+300

p <- ggplot(data_plot, aes(cluster, log2FoldChange, color = cluster)) +
     geom_beeswarm(cex = as.numeric(opt$cex), show.legend = FALSE) + #cex 可以通过点的尺寸定义宽度
     scale_color_manual(values = CustomCol2(1:length(unique(data_plot$cluster)))) + #颜色自定义赋值
     theme(panel.grid = element_blank(), panel.background = element_blank(),
           axis.line = element_line(color = 'black'),
               plot.margin = ggplot2::margin(0, 0, 0, 0.1, unit = "cm"),
           axis.text = ggplot2::element_text(margin = ggplot2::unit(0,"null")),
           axis.title.x = ggplot2::element_text(size = 12),
           axis.title.y = ggplot2::element_text(size = 12),
           axis.text.x = ggplot2::element_text(size = 10, angle = 90,vjust=0.5, hjust=1),
           axis.text.y=ggplot2::element_text(size = 8)) +        #去除默认的背景色、边框
     labs(x = '', y = 'log2FoldChange') + #坐标轴标签  
     geom_hline(aes(yintercept = 0),colour='black',linetype='dashed') + #水平线
        geom_text(data=Up,aes( y = ylim+4, label = num),colour = "red") + 
        geom_text(data=Down,aes( y = ylim+3, label = num),colour = "blue")
ggsave(file.path(output_dir,paste0(group_name[1],"-vs-",group_name[3],"_beeswarm.png",collapse="")),width=8,height=5, plot = p)
ggsave(file.path(output_dir,paste0(group_name[1],"-vs-",group_name[3],"_beeswarm.pdf",collapse="")),width=8,height=5, plot = p)


