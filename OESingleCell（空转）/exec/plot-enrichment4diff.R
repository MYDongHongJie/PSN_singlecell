sub_enrichment_pt = subparsers$add_parser("enrichment_pt", help = "plot for different expression result.")
sub_enrichment_pt$add_argument("--background_go", type="character",
             help="the go background file dir.")
sub_enrichment_pt$add_argument("--category", type="character",default='/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/enrich_background/category.xls',
             help = "the category address.[default: %(default)s]")
sub_enrichment_pt$add_argument("--endnode", type="character",default="/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/endNode3.csv",
             help = "the endnode file address.[default: %(default)s]")
sub_enrichment_pt$add_argument("--file",type="character",
             help = "the different expression result address")
sub_enrichment_pt$add_argument("--background_kegg", type="character",
             help = "the kegg background file dir")
sub_enrichment_pt$add_argument("--circl", type = "character", default = "TRUE",
             help = "if to plot circlplot.[default: %(default)s]")


args <- commandArgs(TRUE)
if ( "enrichment_pt" %in% args ){

  opt<-intial_setting()
  if(opt$sub_name == "enrichment_pt"){
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

  ###读入差异结果
  files_dir <- opt$file
  gene_file <- list.files(opt$file,pattern="[0-9]_anno.xls")%>% file.path(files_dir,.)
  gene_file_name <- list.files(files_dir,pattern="[0-9]_anno.xls")%>%
    lapply(.,function(x){strsplit(x,split = "-diff-")%>% unlist()%>% .[1]})%>%unlist()


  output_go <- file.path(output_dir,"GO_enrichment",gene_file_name) %>% lapply(.,function(x){
      print(x)
      if (! file.exists(x)) {  dir.create(x, recursive = T) }
      x
  })%>% unlist
  enrichment_go <- data.frame(groups='groups',
                            'significant term'='significant term',
                            'pvalue<0.05'='pvalue<0.05',
                            "pvalue<0.01"="pvalue<0.01",
                            "padj<0.05"="padj<0.05",
                            "padj<0.01"="padj<0.01")
  write.table(enrichment_go,file=file.path(paste0(output_dir,'/GO_enrichment'),"enrichment_go.xls"),sep="\t",
              quote=F,
              row.names=F,
              col.names=F,
              na="")


  dif_data <- list()#读入的差异结果，分为all,up,down
  x_plot <- list()##根据差异结果做GO，并整理画图数据。也分为all,up,down
  for (i in 1:length(gene_file)){
      data_all <- read.table(gene_file[i],sep="\t",quote="",header = 1)
      data_up <- dplyr::filter(data_all,up_down=="Up")
      data_down <- dplyr::filter(data_all,up_down=="Down")
      dif_data[[gene_file_name[i]]] <- list(Total=data_all,Up=data_up,Down=data_down)
      ##x_plotp[[1]]所有的，x_plotp[[2]]:up,x_plotp[[3]]:down
      gooutput <- c('-Toal.xls',"-Up.xls","-Down.xls")

      x_plot[[gene_file_name[i]]] <- lapply(c(1:3),function(x){
          thedata <- dif_data[[gene_file_name[i]]][[x]]
          outputfile <- file.path(output_go[i],paste0("enrichment-go-",gene_file_name[i],gooutput[x]))
          thexplot <- OESingleCell::get_go_restule(diffresult =thedata,
                                     TERM2GENE=term2gene,
                                     TERM2NAME=term2name,
                                     n=10,#cc,bp,mf各10个，总共30个
                                     output=outputfile,
                                     category_raw=category_raw,
                                     colname='gene')
          return(thexplot)
      })

      name <- c('GO.top.Total','GO.top.Up','GO.top.Down')
      p <- lapply(c(1:3),function(x){
          data=x_plot[[gene_file_name[i]]][[x]]
          write.table(data[,-12],file=file.path(output_go[i],paste0(name[x],'.xls')),sep='\t',row=F,quote=F)
          title=paste0(gene_file_name[i],"(",c('Total','Up','Down')[x],'):',"Top 30 GO Term")
          plotdata <- OESingleCell::plot_go_bar_1(data,cpcols,title)})

      for(j in 1:3){
          if(!is.null(p[[j]])){
          data_dim <-  dim(ggplot2::ggplot_build(p[[j]])$data[[1]])[1]
          OESingleCell::save_ggplots(file.path(output_go[i],name[j]),
                                                          plot=p[[j]],
                                                          width=ifelse(data_dim>50,data_dim/5,10),
                                                          height=ifelse(data_dim>50,data_dim/7,10),
                                                          dpi=200,bg='white')}

      }

  }

  if(opt$circl=="TRUE"){
  for (j in 1:length(x_plot)){
      name <- c('.GO.top.Total.circlize.pdf','.GO.top.Up.circlize.pdf','.GO.top.Down.circlize.pdf')
      print('plot circl')
      for (i in names(cpcols)){
         ##1:3分别对应total,up,down
         lapply(c(1:3),function(x){
             data <- x_plot[[j]][[x]][x_plot[[j]][[x]]$Category==i,]
             if(dim(data)[1]>0){

             data <- OESingleCell::get_circl_data(data,i,cpcols)
             lgd <- OESingleCell::legend(data,cpcols)
              pdf(file.path(output_go[j],paste0(i,name[x])),w=20,h=8)
              OESingleCell::cir_plot(data,lgd)
              dev.off()


             }

             }
          )


         }

  }}

  ancestor.bp <- list()
  ancestor.cc <- list()
  ancestor.mf <- list()
  for (i in 1:length(x_plot)){
      up_down <- c('total',"up","down")
      for (j in 1:3){##total, up ,down

          ancestor.bp[[names(x_plot)[i]]][[up_down[j]]] <- OESingleCell::get_ancestor(x_plot[[names(x_plot)[i]]][[j]],xx.bp,bp.id)
          ancestor.cc[[names(x_plot)[i]]][[up_down[j]]] <- OESingleCell::get_ancestor(x_plot[[names(x_plot)[i]]][[j]],xx.cc,cc.id)
          ancestor.mf[[names(x_plot)[i]]][[up_down[j]]] <- OESingleCell::get_ancestor(x_plot[[names(x_plot)[i]]][[j]],xx.mf,mf.id)

      }

  }


  level2heatmap <- list()
  for (i in 1:length(x_plot)){
      up_down <- c('total',"up","down")
      for (j in 1:3){
          level2heatmap[[names(x_plot)[i]]][[up_down[j]]] <- OESingleCell::get_heatmap_Data(
              ancestor.bp=ancestor.bp[[i]][[j]],
              ancestor.cc=ancestor.cc[[i]][[j]],
              ancestor.mf=ancestor.mf[[i]][[j]],
              Endnode=endnode,
              go_out=x_plot[[i]][[j]] )

      }

  }

  heatmap <- c('GO.level2.stat.Total.heatmap','GO.level2.stat.Up.heatmap','GO.level2.stat.Down.heatmap')
  for(j in 1:length(level2heatmap)){
      for (i in 1:length(level2heatmap[[j]])){
         p <- OESingleCell::plotheatmap(level2heatmap[[j]][[i]] )
          data_dim_gene <- length(unique(level2heatmap[[j]][[i]]$value))
          data_dim_term <- length(unique(level2heatmap[[j]][[i]]$Description))
         OESingleCell::save_ggplots(file.path(output_go[j],heatmap[i]),
                                                          plot=p,
                                                          width=ifelse(data_dim_term>50,data_dim_term/4,9),
                                                          height=ifelse(data_dim_gene>50,data_dim_gene/7,15),
                                                          dpi=200,bg='white')


      }
  }

  ####level画图数据
  diffdatabar <- list()
  for (j in 1:length(level2heatmap)){
      up_down <- c('total',"up","down")
      for (i in 1:length(level2heatmap[[j]])){
          diffdatabar[[names(level2heatmap)[j]]][[up_down[i]]] <- OESingleCell::get_level2_bardata(
                                                     level2heatmap[[names(level2heatmap)[j]]][[up_down[i]]],
                                                     endnode,
                                                     stat=up_down[i])


      }


  }


  level2file <- c('GO.level2.stat.Total','GO.level2.stat.Up','GO.level2.stat.Down')
  title <- c("(Total): Gene Ontology Classification","(Up): Gene Ontology Classification","(Down): Gene Ontology Classification")
  for (j in 1:length(diffdatabar)){
      title <- c("(Total): Gene Ontology Classification","(Up): Gene Ontology Classification","(Down): Gene Ontology Classification")
      for( i in 1:length(diffdatabar[[j]])){
          data <-  diffdatabar[[j]][[i]]
          p <- ggplot(data,aes(x=Description,y=percent,fill=Category))+
            geom_bar(stat="identity",width=0.8)+
            scale_fill_manual(values=cpcols)+
            theme_test()+
            xlab(NULL)+
            ylab(expression('percent of genes(%)'))+
            labs(title = paste0(names(diffdatabar)[j],title[i]))+
            theme(plot.title = element_text(hjust = 0.5))+
              # scale_x_discrete(labels=difflevel2$Description)+
            theme(axis.text=element_text(face="plain"))+
            theme(legend.position='left',
                  legend.justification='top',
                  legend.direction='vertical')+
            theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2,size=14,color=data$color))+
            theme(axis.ticks.x=element_blank())+
            scale_y_continuous(expression('percent of genes(%)'),sec.axis = sec_axis(~.*allgenenumber/100,
                                                     name = "genes numbers",
                                                     breaks=seq(0,max(data$Freq),20),
                                                     labels = seq(0,max(data$Freq),20)),
                                                      expand=c(0,0))

          OESingleCell::save_ggplots(file.path(output_go[j],level2file[i]),
                                                          plot=p,
                                                          width=15,
                                                          height=9,
                                                          dpi=200,bg='white')





      }



  }


    ##输出GO.level2.stat.Down.xls
    stage <- c('Total','Up','Down')
    for(j in 1:length(level2heatmap)){
        heatmap_group <- level2heatmap[[j]]
        for(i in 1:length(heatmap_group)){
            heatmap_group_stage <- heatmap_group[[i]][,-c(1,3,6)]%>% reshape2::dcast(.,Category+Description~value)
            heatmap_group_stage2 <- sapply(1:dim(heatmap_group_stage)[1],function(x){
            gene <- colnames(heatmap_group_stage)[c(FALSE,FALSE,heatmap_group_stage[x,-c(1,2)]!=0)]
            gene <- paste(gene,collapse=';')
            })
            heatmap_output <- data.frame('GO_classify1'='#Total_gene',
                                         GO_classify2=NA,
                                         assign(paste0("diff_",stage[i]),dim(heatmap_group_stage)[2]-2),
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
                        file=file.path(output_go[j],paste0(level2file[i],'.xls')),
                        sep="\t",
                        quote=F,
                        row=F,
                        na="")

        }

    }



  comp_up_down <- list()
  for (j in 1:length(diffdatabar)){
      comp_up_down[[names(diffdatabar)[j]]] <- rbind(diffdatabar[[j]][[2]],diffdatabar[[j]][[3]])%>%
                                              dplyr::arrange(Category,Description,up_down,dplyr::desc(Freq))
      comp_up_down[[names(diffdatabar)[j]]]$Up_Down <- Hmisc::capitalize(comp_up_down[[names(diffdatabar)[j]]]$up_down) %>%
      factor(.,levels=c("Up","Down"), ordered=TRUE)



  }


  for (j in 1:length(comp_up_down)){
      data <- comp_up_down[[j]]
       data$Description <- factor(data$Description,levels =unique(data$Description))
       data_color <- data$color[match(levels(data$Description),data$Description)]
      title <- paste0("Gene Ontology Classification(",names(comp_up_down),")")
          p <- ggplot(data,mapping=aes(x=Description,y=percent,fill=Up_Down))+
          geom_bar(stat="identity",width=0.7,position = 'dodge')+scale_fill_manual(values=c("#F985A3","#77E5D0"))+
          theme_test()+
          xlab(NULL)+ylab(expression('percent of genes(%)'))+labs(title = title)+theme(plot.title = element_text(hjust = 0.5))+
          # scale_x_discrete(labels=difflevel2$Description)+
          theme(axis.text=element_text(face="plain"),panel.border=element_blank())+
          theme(legend.position='left',
                  legend.justification='top',
                  legend.direction='vertical')+
          theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2,size=8,color=data_color))+
          theme(axis.ticks.x=element_blank())+
          scale_y_continuous('percent',sec.axis = sec_axis(~.*allgenenumber/100,
                                                 name = "genes numbers",
                                                 breaks=seq(0,max(data$Freq),20),
                                                 labels = seq(0,max(data$Freq),20)))+
          theme(axis.line=element_line(size=0.1))+ scale_y_continuous(expand = expansion(mult = c(0, 0)))

      OESingleCell::save_ggplots(file.path(output_go[j],'Up_vs_Down.GO.level2.stat'),
                                                          plot=p,
                                                          width=15,
                                                          height=9,
                                                          dpi=200,bg='white')




  }
##kegg

  re_kegg <- read.table(opt$background_kegg,sep="\t",quote="")
  re_kegg <- split(re_kegg, 1:nrow(re_kegg))
  term2gene_kegg <- lapply(re_kegg,function(x){data=OESingleCell::get_backdata_list(x)[[1]]}) %>% dplyr::bind_rows()
  term2name_kegg <- lapply(re_kegg,function(x){data=OESingleCell::get_backdata_list(x)[[2]]}) %>% dplyr::bind_rows()
  print(opt$background_kegg)

  output_kegg <- file.path(output_dir,'KEGG_enrichment')%>% file.path(.,gene_file_name)%>%
   lapply(.,function(x){
      print(x)
      if (! file.exists(x)) {  dir.create(x, recursive = T) }
      x
  })%>% unlist

  enrichment_kegg <- data.frame(groups='groups',
                            'significant term'='significant term',
                            'pvalue<0.05'='pvalue<0.05',
                            "pvalue<0.01"="pvalue<0.01",
                            "padj<0.05"="padj<0.05",
                            "padj<0.01"="padj<0.01")
  print("writetable")
  write.table(enrichment_kegg,file=file.path(paste0(output_dir,'/KEGG_enrichment'),"enrichment_kegg.xls"),sep="\t",quote=F,row.names=F,col.names=F,na="")


  ###依然使用dif_data
  kegg_plot <- list()##根据差异结果做GO，并整理画图数据。也分为all,up,down
    keggoutput <- c('KEGG.top.Total',"KEGG.top.Up","KEGG.top.Down")
    keggoutput2 <- c('-Total',"-Up","-Down")
    plot_titlt <- c('(Total): KEGG Enrichment top 20','(Up): KEGG Enrichment top 20','(Down): KEGG Enrichment top 20')
  for (i in 1:length(gene_file)){


      kegg_plot[[gene_file_name[i]]] <- lapply(c(1:3),function(x){
          thedata <- dif_data[[gene_file_name[i]]][[x]]
          outputfile <- file.path(output_kegg[i],paste0("enrichment-kegg-",gene_file_name[i],keggoutput2[x],".xls"))
          thexplot <- OESingleCell::get_go_restule(diffresult =thedata,
                                                   TERM2GENE=term2gene_kegg,
                                                   TERM2NAME=term2name_kegg,
                                                   n=20,
                                                   output=outputfile,
                                                   category_raw=category_raw,
                                                   colname='gene')
          return(thexplot)
      })

      for(j in 1:length(kegg_plot[[gene_file_name[i]]])){
          if(!is.null(kegg_plot[[gene_file_name[i]]][[j]])){
            data <- kegg_plot[[gene_file_name[i]]][[j]]
            write.table(data[,-c(3,12)],file=file.path(output_kegg[i],paste0(keggoutput[j],'.xls')),sep='\t',row=F,quote=F)
          data$GeneRatio_num <- data$Count/dim(dif_data[[gene_file_name[i]]][[j]])[1]
              thetitle <- paste0(gene_file_name[i],plot_titlt[j])
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
          OESingleCell::save_ggplots(file.path(output_kegg[i],keggoutput[j]),
                                                          plot=p,
                                                          width=ifelse(data_dim>50,data_dim/5,7),
                                                          height=ifelse(data_dim>50,data_dim/7,7),
                                                          dpi=200,bg='white')}

          }


      }




}
}

#2. go 背景文件，示例：/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_go.backgroud.xls
#3. kegg 背景文件，示例：/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_kegg.backgroud.xls
#/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/enrich_background/category.xls
#endnode="/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/endNode3.csv"
