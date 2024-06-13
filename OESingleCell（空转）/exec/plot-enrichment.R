sub_enrich_pt = subparsers$add_parser("enrich_pt", help = "plot for different expression result.")
sub_enrich_pt$add_argument("--background_go", type="character",
             help="the go background file dir.")
sub_enrich_pt$add_argument("--category", type="character",default='/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/enrich_background/category.xls',
             help = "the category address.[default: %(default)s]")
sub_enrich_pt$add_argument("--endnode", type="character",default="/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/endNode3.csv",
             help = "the endnode file address.[default: %(default)s]")
sub_enrich_pt$add_argument("--genelist",type="character",
             help = "the genelist file that must have header 'gene'")
sub_enrich_pt$add_argument("--background_kegg", type="character",
             help = "the kegg background file dir")
sub_enrich_pt$add_argument("--circl", type = "character", default = "TRUE",
             help = "if to plot circlplot.[default: %(default)s]")

args <- commandArgs(TRUE)
if ( "enrich_pt" %in% args ){

  opt<-intial_setting()
  if(opt$sub_name == "enrich_pt"){
    suppressPackageStartupMessages(library('clusterProfiler'))
    suppressPackageStartupMessages(library("GO.db"))
    suppressPackageStartupMessages(library("circlize"))
    suppressPackageStartupMessages(library("ggplot2"))
    ref <- read.table(opt$background_go,sep="\t",quote="")
    ref <- split(ref, 1:nrow(ref))
    #/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/enrich_background/category.xls
    category_raw <- read.table(opt$category,sep="\t",header=1,quote="",row.names = 1)

    cpcols <- c('#6495ED','#8FBC8F','#F4A460')
    names(cpcols) <- c('biological_process','cellular_component','molecular_function')
    term2gene <- lapply(ref,function(x){data=OESingleCell::get_backdata_list(x)[[1]]}) %>% dplyr::bind_rows()
    term2name <- lapply(ref,function(x){data=OESingleCell::get_backdata_list(x)[[2]]}) %>% dplyr::bind_rows()
    allgenenumber <- unique(term2gene$geneId)%>% length()

    ###get ancestor file
    endnode<-read.csv(opt$endnode, header=F) %>% as.matrix()
    xx.bp <- as.list(GOBPANCESTOR)
    bp.id<-names(xx.bp)
    xx.cc <- as.list(GOCCANCESTOR)
    cc.id<-names(xx.cc)
    xx.mf <- as.list(GOMFANCESTOR)
    mf.id<-names(xx.mf)

    genedata <- read.table(opt$genelist,sep="\t",quote="",header=1)
    output_go <- file.path(output_dir,'GO_enrichment')
    if (! file.exists(output_go)) { dir.create(output_go, recursive = T) }
    output_go<-normalizePath( output_go)
    outputfile <- file.path(output_go,"enrichement_go.txt")
    thexplot <- OESingleCell::get_go_restule2(diffresult =genedata,
                                       TERM2GENE=term2gene,
                                       TERM2NAME=term2name,
                                       n=10,#cc,bp,mf各10个，总共30个
                                       output=outputfile,
                                       category_raw=category_raw,
                                       colname='gene')
    write.table( thexplot[,-12],file=file.path(output_go,"top_30_GO_trem.txt"),sep='\t',row=F,quote=F)
    #for barplot
    plotdata <- OESingleCell::plot_go_bar_1(thexplot,cpcols,'top 30 GO trem')

    data_dim <-  dim(ggplot2::ggplot_build(plotdata)$data[[1]])[1]
    OESingleCell::save_ggplots(file.path(output_go,'top_30_GO_trem'),
                                                            plot=plotdata,
                                                            width=ifelse(data_dim>50,data_dim/5,10),
                                                            height=ifelse(data_dim>50,data_dim/7,10),
                                                            dpi=200,bg='white')

    #for circle plot

    if(opt$circl=="TRUE"){
    print('plot circl')
    for (i in names(cpcols)){
           data <-  thexplot[thexplot$Category==i,]
           if(dim(data)[1]>0){
           data <- OESingleCell::get_circl_data(data,i,cpcols)
           lgd <- OESingleCell::legend(data,cpcols)
            pdf(file.path(output_go,paste0(i,'.GO.top.Total.circlize.pdf')),w=20,h=8)
            OESingleCell::cir_plot(data,lgd)
            dev.off()

           }
    }}
    #for heatmapplot
    ancestor.bp <- OESingleCell::get_ancestor(thexplot,xx.bp,bp.id)
    ancestor.cc <- OESingleCell::get_ancestor(thexplot,xx.cc,cc.id)
    ancestor.mf <- OESingleCell::get_ancestor(thexplot,xx.mf,mf.id)

    level2heatmap <- OESingleCell::get_heatmap_Data(
    ancestor.bp=ancestor.bp,
    ancestor.cc=ancestor.cc,
    ancestor.mf=ancestor.mf,
    Endnode=endnode,
    go_out=thexplot)

    p <- OESingleCell::plotheatmap(level2heatmap)
    data_dim_gene <- length(unique(level2heatmap$value))
    data_dim_term <- length(unique(level2heatmap$Description))
    OESingleCell::save_ggplots(file.path(output_go,'GO.level2.stat.heatmap'),
                                                plot=p,
                                                width=ifelse(data_dim_term>50,data_dim_term/4,9),
                                                height=ifelse(data_dim_gene>50,data_dim_gene/7,10),
                                                dpi=200,bg='white')

    diffdatabar <- OESingleCell::get_level2_bardata(level2heatmap,
                                                       endnode,
                                                       stat='total')
    p <- ggplot(diffdatabar ,aes(x=Description,y=percent,fill=Category))+
              geom_bar(stat="identity",width=0.8)+
              scale_fill_manual(values=cpcols)+
              theme_test()+
              xlab(NULL)+
              ylab(expression('percent of genes(%)'))+
              labs(title = 'Gene Ontology Classification')+
              theme(plot.title = element_text(hjust = 0.5))+
                # scale_x_discrete(labels=difflevel2$Description)+
              theme(axis.text=element_text(face="plain"))+
              theme(legend.position='left',
                    legend.justification='top',
                    legend.direction='vertical')+
              theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2,size=14,color=diffdatabar$color))+
              theme(axis.ticks.x=element_blank())+
              scale_y_continuous(expression('percent of genes(%)'),sec.axis = sec_axis(~.*allgenenumber/100,
                                                       name = "genes numbers",
                                                       breaks=seq(0,max(diffdatabar$Freq),20),
                                                       labels = seq(0,max(diffdatabar$Freq),20)),
                                                        expand=c(0,0))
    OESingleCell::save_ggplots(file.path(output_go,'GO.level2.stat'),
                                                            plot=p,
                                                            width=15,
                                                            height=9,
                                                            dpi=200,bg='white')

    heatmap_group <- level2heatmap

    heatmap_group_stage <- level2heatmap[,-c(1,3,6)]%>% reshape2::dcast(.,Category+Description~value)
    heatmap_group_stage2 <- sapply(1:dim(heatmap_group_stage)[1],function(x){
    gene <- colnames(heatmap_group_stage)[c(FALSE,FALSE,heatmap_group_stage[x,-c(1,2)]!=0)]
    gene <- paste(gene,collapse=';')
    })
    heatmap_output <- data.frame('GO_classify1'='#Total_gene',
                                 GO_classify2=NA,
                                 termgene_number=dim(heatmap_group_stage)[2]-2,
                                 Gene=NA)
    num <- sapply(heatmap_group_stage2,function(x){
       num <- length(unlist(strsplit(x,';')))

    })
    heatmap_output2 <- data.frame(heatmap_group_stage[,1],
                                 heatmap_group_stage[,2],
                                 unname(num),
                                 heatmap_group_stage2)
    colnames( heatmap_output2) <-  colnames(heatmap_output)
    heatmap_output3 <- rbind( heatmap_output, heatmap_output2)
    write.table( heatmap_output3,
                file=file.path(output_go,'GO.level2.xls'),
                sep="\t",
                quote=F,
                row=F,
                na="")

    #kegg
    re_kegg <- read.table(opt$background_kegg,sep="\t",quote="")
    re_kegg <- split(re_kegg, 1:nrow(re_kegg))
    term2gene_kegg <- lapply(re_kegg,function(x){data=OESingleCell::get_backdata_list(x)[[1]]}) %>% dplyr::bind_rows()
    term2name_kegg <- lapply(re_kegg,function(x){data=OESingleCell::get_backdata_list(x)[[2]]}) %>% dplyr::bind_rows()
    print(opt$background_kegg)


    output_kegg <- file.path(output_dir,'KEGG_enrichment')
    if (! file.exists(  output_kegg)) { dir.create( output_kegg, recursive = T) }
    output_kegg<-normalizePath( output_kegg)
    outputfile <- file.path(output_kegg,'enrichment_kegg.txt')
    thexplot_kegg <- OESingleCell::get_go_restule2(diffresult =genedata,
                                                       TERM2GENE=term2gene_kegg,
                                                       TERM2NAME=term2name_kegg,
                                                       n=20,
                                                       output=outputfile,
                                                       category_raw=category_raw,
                                                       colname='gene')


    data <- thexplot_kegg
    thetitle <- 'kegg enrichment top20 term'
    data$GeneRatio_num <- data$Count/dim(genedata)[1]
    p <- ggplot(data, showCategory = 20,
          aes(GeneRatio_num, forcats::fct_reorder(Description, GeneRatio_num)))+
          geom_segment(aes(xend=0, yend =  Description))+
          geom_point(aes(color=log10p, size = Count))+
          theme_minimal() +
          xlab("GeneRatio") +
          ylab(NULL) +
          ggtitle(thetitle)+labs(color = expression('-log'[10]*' Pvalue'))+
          scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
          scale_size_continuous(range=c(2, 10))



    data_dim <-  dim(ggplot2::ggplot_build(p)$data[[1]])[1]
    OESingleCell::save_ggplots(file.path(output_kegg,'kegg_enrichment'),
                                                plot=p,
                                                width=ifelse(data_dim>50,data_dim/5,7),
                                                height=ifelse(data_dim>50,data_dim/7,7),
                                                dpi=200,bg='white')


}
}
#scVis -o singlego enrich_pt --background_go "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_go.backgroud.xls" --genelist toegene.txt --background_kegg "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_kegg.backgroud.xls"