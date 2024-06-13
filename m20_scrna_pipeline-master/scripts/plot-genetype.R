# Created by: kun.sun
# Created on: 2023/9/8
sub_genetype<- subparsers$add_parser("genetype", help = "绘制gene type分布条形图")
sub_genetype$add_argument("-m", "--starsolo", type = "character", help =  "STAR result 文件夹路径")
sub_genetype$add_argument("-s", "--prefix", type = "character", help = "sample name id ")
sub_genetype$add_argument("-r","--refgenome",type="character",help="reference genome pathway")
args <- commandArgs(TRUE)
if ("genetype" %in% args) {
  opt <- intial_setting()
  if(opt$sub_name == "genetype" ){
    ditto = c(
         "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
         "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
         "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
         "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
         "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C")
    ##data praper
    counts_mtx <- Seurat::Read10X(glue::glue("{opt$starsolo}/filtered/"),gene.column=1)
    data_ob <- Seurat::CreateSeuratObject(counts = counts_mtx, assay = "RNA")
    refdata<-read.table(glue::glue("{opt$refgenome}/STAR_index/geneInfo.tab"),skip=1,sep="\t")
    colnames(refdata)<-c("geneID","gene_name","gene_type")
    biotypes <- unique(refdata[["gene_type"]])
    num <- Matrix::rowSums(Seurat::GetAssayData(data_ob, assay="RNA",slot= "counts"))
    num <- num[which(num > 0)]
    num <- data.frame(nCounts = num) %>% tibble::rownames_to_column("geneID")
    gene_type_summary <- merge(num,refdata,by.x="geneID")
    #export total data
    write.table(gene_type_summary,glue::glue("{output_dir}/{opt$prefix}_gene_type_totaldata.csv"),sep=",",row.names=F)

    plot_rawdata<-gene_type_summary %>%
      dplyr::group_by(gene_type) %>%
      dplyr::summarise(Counts=length(geneID)) %>%
      as.data.frame %>%
      dplyr::arrange(desc(Counts))

    ## merger IG_/TR_ gene type
    IG_gene<-biotypes[grep("IG_",biotypes)]
    TR_gene<-biotypes[grep("TR_",biotypes)]
    others_counts<-sum(plot_rawdata[which(plot_rawdata$gene_type %in% c(IG_gene,TR_gene)),]$Counts)
    plot_rawdata<- plot_rawdata %>%
      dplyr::filter( ! gene_type %in% c(IG_gene,TR_gene))
    plot_rawdata<-rbind(plot_rawdata,c("others",others_counts))

    ##plot
    plot_rawdata$gene_type<-factor(plot_rawdata$gene_type,levels=rev(plot_rawdata$gene_type))
    plot_rawdata$Counts<-as.numeric(plot_rawdata$Counts)
    sum_colors <- ditto[1:length(unique(plot_rawdata$gene_type))]

    ##export plotdata
    write.table(plot_rawdata,glue::glue("{output_dir}/{opt$prefix}_gene_type_plotdata.csv"),sep=",",row.names=F)

    p <- ggplot(plot_rawdata,mapping=aes(x=Counts,y=gene_type,fill=gene_type)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = Counts)) +
      scale_x_log10(labels = scales::number)
    plot <- p+theme_bw()+
      ylab("Gene type")+
      labs(title = opt$prefix)+
      theme(text=element_text(size=14),axis.text=element_text(size=12),plot.title = element_text(hjust = 0.5))+
      scale_fill_manual(values = sum_colors)+
      guides(fill = "none")

    OESingleCell::save_ggplots(plot, filename = glue::glue("{output_dir}/{opt$prefix}_gene_type"), width = 7, height = 6)
  }
}
