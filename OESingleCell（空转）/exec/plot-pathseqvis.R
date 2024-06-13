## 针对单个样本绘图
plot_featureplot<- function(
    data_ob,
    features="Bacteria_UMIs",
    sample,
    HE=TRUE,
    features_rename=NULL,
    title=NULL,
    output_dir)
{
    ## step1 获取单个样本数据==
    pdf(NULL)
    data_ob_sub<-subset(data_ob,subset=sampleid==sample)
    data_ob_sub@meta.data$sampleid<-as.character(data_ob_sub@meta.data$sampleid)
    if (!is.null(Seurat::Images(data_ob_sub))) {
      unuse_images <- Seurat::Images(data_ob_sub)[!Seurat::Images(data_ob_sub) %in% (data_ob_sub@meta.data$sampleid %>% unique)]
      if (length(unuse_images) > 0) { data_ob_sub@images[unuse_images] <- NULL }
    }
    if(!is.null(features_rename)){
        if( features_rename %in% colnames(data_ob@meta.data)){data_ob_sub@meta.data[,features_rename]<-NULL}
        data_ob_sub@meta.data <- data_ob_sub@meta.data %>% dplyr::rename(!!features_rename:=features)
        name <- glue::glue("{features}_{features_rename}")
        features <- features_rename
    }else{
        name<-features
    }
    title<-ifelse(is.null(title),"",title)
    legend_length<- ifelse(stringr::str_length(features)<10,0.3,stringr::str_length(features)*0.03)
    dev.off()
    ## step2 绘制空间featureplot图
    pdf(NULL)
    plot_lab <- ggplot() +
             annotate(geom = "text", x = 1, y = 1, label = sample, angle = 90,size=6) +
             coord_cartesian(clip = "off")+
             theme_void()
   if(HE){
    plot_HEimages<-OESingleCell::SpatialPlot(data_ob_sub,
                                             features = "nCount_Spatial",
                                             alpha=0,images=sample,
                                             #image.alpha=1,
                                             cols="spectral",
                                             pt.size.factor = 1.2,
                                             min.cutoff = 0,
                                             ##max.cutoff ="q100",
                                             ncol=1,
                                             theme_void=FALSE,
                                             combine = F)[[1]]+Seurat::NoLegend()+ggtitle("")
   }
    plot_featureplot <-OESingleCell::SpatialPlot(data_ob_sub,
                             features = features,
                             alpha=1,
                             images=sample,
                             #image.alpha=1,
                             cols="spectral",
                             pt.size.factor = 1.2,
                             min.cutoff = 0,
                             ##max.cutoff ="q100",
                             ncol=1,
                             theme_void=FALSE,
                             combine =F)[[1]]+theme(legend.position = 'right' )+ggtitle(glue::glue("{title}"))
    plot_legend <- ggpubr::get_legend( plot_featureplot) %>%ggpubr::as_ggplot()
    if(HE){
         width<- 4.5*(0.1+1+1+legend_length)
        plot<-patchwork::wrap_plots( plot_lab,
                                     plot_HEimages,
                                     plot_featureplot+Seurat::NoLegend(),
                                     plot_legend,
                                     nrow = 1,
                                     ncol = 4,
                                     widths=c(0.1,1,1,legend_length))
    }
    else{
         width<-4.5*(0.1+1+legend_length)
         plot<-patchwork::wrap_plots( plot_lab,
                                     plot_featureplot+Seurat::NoLegend(),
                                     plot_legend,
                                     nrow = 1,
                                     ncol = 3,
                                     widths=c(0.1,1,legend_length))
    }

    OESingleCell::save_ggplots(plot=plot,
                               filename =glue::glue("{output_dir}/spatialfeatureplot_{name}_for_{sample}"),
                               width= width,
                               height=4.5)

         }

plot_vlnplot<-function(
    data_ob,
    features,
    group.by,
    output_dir)
{
   ##绘制小提琴图
     pdf(NULL)
     plot_vlnplot <-Seurat::VlnPlot(data_ob ,
                                    features = features,
                                    group.by= group.by,
                                    cols = OESingleCell::SelectColors(levels(data_ob@meta.data[,group.by]),palette = "blindless"))
     OESingleCell::save_ggplots(plot=plot_vlnplot ,
                                filename =glue::glue("{output_dir}/vlnplot_for_{features}_group_by_{group.by}"),
                                width=1*length(levels(data_ob@meta.data[,group.by]))+2,
                                height=7 )
    dev.off()
    }
plot_top10_piechart_featureplot <-function(
    data_ob,
    sample,
    top10_genus,
    top10_reads_stat,
    top10_umi_stat,
    sum_colors,
    output_dir)
{
    futile.logger::flog.info("1:绘制top10饼图")
    colors<- sum_colors[as.character(top10_umi_stat$pathogen)]
    ## png格式
    png(glue::glue("{output_dir}/piechart/Top10_piechart_for_{sample}.png"),height = 500,width=800+max(stringr::str_length(top10_genus))*0.1*100)
    par(mar = c(0,0,2,0))
    pie(top10_umi_stat$`Relative Abundance (%) by UMI`, radius=0.8,main=sample,col=colors,labels=top10_umi_stat$labels,cex=1.2,cex.main=2)
    legend("right", names(colors), cex=1.5, bty = "n",fill = colors)
    dev.off()
    ## pdf格式
    pdf(glue::glue("{output_dir}/piechart/Top10_piechart_for_{sample}.pdf"),height = 5,width=10+max(stringr::str_length(top10_genus))*0.1)
    par(mar = c(0,0,2,0))
    pie(top10_umi_stat$`Relative Abundance (%) by UMI`, radius=0.8,main=sample,col=colors,labels=top10_umi_stat$labels,cex=1.2,cex.main=2)
    legend("right", names(colors), cex = 1.5,bty = "n", fill = colors)
    dev.off()
    futile.logger::flog.info("2:绘制top10主要微生物的UMI空间分布图")
    data_ob_sub<-subset(data_ob,subset=sampleid==sample)
    for(genus in top10_genus){
       plot_featureplot(data_ob_sub,
                        features=genus,
                        sample=sample,
                        HE=FALSE,
                        features_rename="Bacteria_UMIs",
                        title=genus,
                        output_dir=glue::glue("{output_dir}/featureplot/{sample}"))
    }
    futile.logger::flog.info("3:绘制前10微生物的UMI含量barplot")
    plot_barplot<- ggplot(top10_umi_stat%>%dplyr::mutate(pathogen=factor(pathogen,levels=top10_umi_stat$pathogen%>%rev))%>%dplyr::arrange(desc("pathogen")) ,
                   aes(y= pathogen,x=`Number of UMIs`)) +
                   geom_col(fill = "#F6C6A9") +
                   ggtitle(sample)+
                  # geom_text(aes(label = `Number of UMIs`),hjust =0.5)+
                   ylab("Number of UMIs")+
                   theme_light(base_size = 20)+
                   theme(panel.grid=element_blank(),
                         axis.text = element_text(size = 15),
                         axis.title.x  = element_blank(),
                         plot.title = element_text(hjust = 0.5))
   OESingleCell::save_ggplots(plot=  plot_barplot,
                               filename =glue::glue("{output_dir}/barplot/barplot_UMIs_count_for_{sample}"),
                               width=8 +max(stringr::str_length(top10_genus))*0.1,
                               height=4 )
   dev.off()
    #===================================================================================================================
    futile.logger::flog.info("4:绘制前10微生物的Reads含量barplot")
    plot_barplot<-ggplot(top10_reads_stat%>%dplyr::mutate(pathogen=factor(pathogen,levels=top10_umi_stat$pathogen%>%rev))%>%dplyr::arrange(desc("pathogen")) ,
                  aes(y= pathogen,x=`Number of Reads`)) +
                  geom_col(fill = "#a1aaa6") +
                  ggtitle(sample)+
                  #geom_text(aes(label = `Number of Reads`),hjust =0.5)+
                  ylab("Number of Reads")+
                  theme_light(base_size = 20)+
                  theme(panel.grid=element_blank(),
                        axis.text = element_text(size = 15),
                        axis.title.x  = element_blank(),
                        plot.title = element_text(hjust = 0.5))
   OESingleCell::save_ggplots( plot =  plot_barplot,
                               filename =glue::glue("{output_dir}/barplot/barplot_Transcript_Reads_for_{sample}"),
                               width=8+max(stringr::str_length(top10_genus))*0.1,
                               height=4 )
}

statistic_info<-function(data_ob, input_dir, output_dir)
{
 futile.logger::flog.info("1.获取样本信息列表")
sample_list<-list.files(path = input_dir, pattern ="visium.raw.readnamepath") %>%
                 stringr::str_replace_all(pattern = '.visium.raw.readnamepath',replacement = "")
 futile.logger::flog.info("2.筛选有效reads并保存")
for(sample in sample_list){
    if(!file.exists(glue::glue('{input_dir}/{sample}.visium.filtered.readnamepath'))){
        suppressMessages(
            validate.data<- readr::read_csv( glue::glue('{input_dir}/{sample}.visium.validate.csv'),
                                             col_names = c("barcode_fasta","pathogen"))%>%
                   tidyr::separate(barcode_fasta,sep="\\+",into = c("barcode","fasta")) %>%
                   dplyr::select(!"pathogen"))
        suppressMessages(
            raw.data<-readr::read_tsv( glue::glue('{input_dir}/{sample}.visium.raw.readnamepath'),
                              col_names = c("read","barcode","fasta","evalue","score","pathogen")))

        dplyr::inner_join(raw.data,validate.data,by = c("barcode","fasta"))%>%
                readr::write_tsv(glue::glue('{input_dir}/{sample}.visium.filtered.readnamepath'))
    }
}
 futile.logger::flog.info("3.导入所有样本的UMI、Reads信息并汇总") #==========================================================================
 futile.logger::flog.info("3.1 导入UMI信息，统计并存入seurat对象中")
 suppressMessages(result_genus_raw  <-do.call(dplyr::bind_rows, lapply(sample_list, function(sample)
                readr::read_csv(glue::glue('{input_dir}/{sample}.visium.genus.csv'))%>%
                                     tidyr::separate(barcode,sep="-",into=c("a","b")) %>%
                                     dplyr::mutate(barcode=paste0(sample,"-",a),sampleid=sample)%>%
                                     dplyr::select(!c("a","b"))%>%
                                     # replace(is.na(.), 0)  %>%
                                     # dplyr::rowwise(barcode)%>%
                                     # dplyr::mutate(Bacteria_UMIs= sum(dplyr::c_across(where(is.numeric))))%>%
                                     tibble::column_to_rownames("barcode"))))
 futile.logger::flog.info("获取barcode X genus 矩阵，存入misc ,然后统计每个spot中Bacteria_UMIs，并导入meta.data")
    result_genus<- result_genus_raw %>%
                    tibble::rownames_to_column("barcode")%>%
                    replace(is.na(.), 0)
    data_ob@misc$genus_baceria_list<-result_genus
    result_genus <- result_genus %>%
                    dplyr::rowwise(barcode)%>%
                     dplyr::select(!"sampleid")%>%
                    dplyr::mutate(Bacteria_UMIs= sum(dplyr::c_across(where(is.numeric))))%>%
                    tibble::column_to_rownames("barcode")
    result_genus[result_genus==0] <- NA
    colnames(result_genus)<- gsub("\\(","",colnames(result_genus) )
    colnames(result_genus)<- gsub("\\)","",colnames(result_genus) )
    data_ob<-Seurat::AddMetaData(data_ob,metadata = result_genus)
 futile.logger::flog.info("3.2 导入reads信息，统计并存入seurat对象中")
suppressMessages(result_reads<-do.call(rbind, lapply(sample_list, function(sample)
               readr::read_tsv( glue::glue('{input_dir}/{sample}.visium.filtered.readnamepath'))%>%
                        dplyr::group_by(barcode)%>%
                        dplyr::summarise(freq=dplyr::n())%>%
                         `colnames<-`(c("barcode","Bacteria_reads")) %>%
                        tidyr::separate(barcode,sep="-",into=c("a","b")) %>%
                        dplyr::mutate(barcode=paste0(sample,"-",a))%>%
                        dplyr::select(!c("a","b"))%>%
                        tibble::column_to_rownames("barcode"))))

data_ob<-Seurat::AddMetaData(data_ob,metadata = result_reads)
data_ob@meta.data<- data_ob@meta.data%>% dplyr::mutate(Bacteria_reads=ifelse(is.na(Bacteria_reads),NA,Bacteria_reads))
data_ob@meta.data$sampleid<-factor(data_ob@meta.data$sampleid, levels=data_ob@meta.data$sampleid%>%unique%>%sort)
 futile.logger::flog.info("4.统计整个样本水平微生物含量数据")#==========================================================================
sample_stat<- dplyr::bind_cols(
    data_ob@meta.data$sampleid%>%table()%>%as.data.frame()%>%
        `colnames<-`(c("sampleid","Total Spaceranger Detect Spots")) %>%tibble::column_to_rownames("sampleid"),
    data_ob@meta.data%>% dplyr::filter(!is.na(Bacteria_UMIs))%>%.$sampleid %>% table()%>%as.data.frame()%>%
        `colnames<-`(c("sampleid","Number of Positive Genus Capture Spots")) %>%tibble::column_to_rownames("sampleid"),
    data_ob@meta.data%>% as.data.frame()%>%dplyr::select("sampleid","Bacteria_reads")%>%
        dplyr::group_by(sampleid)%>%dplyr::summarize(`Total  Number of Genus Reads` = sum(Bacteria_reads,na.rm=TRUE))%>%
        tibble::column_to_rownames("sampleid"),
    data_ob@meta.data%>% as.data.frame()%>%dplyr::select("sampleid","Bacteria_UMIs")%>%
        dplyr::group_by(sampleid)%>%dplyr::summarize(`Total  Number of Genus UMIs` = sum(Bacteria_UMIs,na.rm=TRUE))%>%
        tibble::column_to_rownames("sampleid"))
    readr::write_tsv(sample_stat%>%tibble::rownames_to_column("sampleid"),
                     glue::glue("{output_dir}/1.Bacteria_reads_Bacteria_UMIs_distribution/summary_statistic.xls"))
 futile.logger::flog.info("5.统计样本水平微生物含量数据")#==========================================================================
    result_genus_raw <- result_genus_raw %>%
                         tibble::rownames_to_column("barcode") %>%
                        dplyr::filter(barcode %in% rownames(data_ob@meta.data))%>%
                        tibble::column_to_rownames("barcode")
    ## 统计每种微生物中UMI数目
    pathogen_UMI<- result_genus_raw %>%
         replace(is.na(.), 0) %>%
         dplyr::group_by(sampleid)%>%
         tidyr::pivot_longer(!sampleid, names_to = "pathogen",values_to = "Number of UMIs")%>%
         dplyr::group_by(sampleid,pathogen)%>%
         dplyr::summarize(`Number of UMIs` = sum(`Number of UMIs`,na.rm=TRUE))%>%
         dplyr::filter(`Number of UMIs`>=1)
    ###  统计每种微生物中spot数目
    pathogen_spot<- result_genus_raw %>%
    replace(is.na(.), 0) %>%
    tibble::rownames_to_column("barcode")%>%
    tidyr::pivot_longer(!c("sampleid","barcode"), names_to = "pathogen",values_to = "Number of UMIs") %>%
    dplyr::filter(`Number of UMIs`>0)%>%
    dplyr::group_by(sampleid)%>%
    dplyr::count(pathogen)%>%
    dplyr::rename("Number of Positive Capture Spots"=n)
    ###  统计每种微生物中reads数目
    results<-list()
    for( sample in sample_list){
      results[[sample]]<- readr::read_tsv( glue::glue('{input_dir}/{sample}.visium.filtered.readnamepath'))%>%
           tidyr::separate_rows(pathogen, sep=",")%>%
           dplyr::mutate(pathogen=factor(pathogen,levels=c(pathogen_UMI%>%dplyr::arrange(sampleid,desc(`Number of UMIs`))%>%dplyr::filter(sampleid==sample)%>%.$pathogen )))%>%
           dplyr::arrange(desc(pathogen))%>%
           dplyr::ungroup() %>%
           dplyr::distinct(read,.keep_all = T)%>%
           dplyr::mutate(sampleid=sample)
    }
    ### 合并数据
    pathogen_reads<-dplyr:: bind_rows(results)%>%
                    dplyr::group_by(sampleid)%>%
                    dplyr::count(pathogen)%>%
                    dplyr::rename("Number of Reads"=n)
    pathogen_stat<- dplyr::left_join( dplyr::left_join(pathogen_spot,pathogen_UMI,by=c("sampleid","pathogen")),
                                      pathogen_reads,
                                      by=c("sampleid","pathogen")) %>%
                    dplyr::left_join(sample_stat%>%tibble::rownames_to_column("sampleid") ,by="sampleid")%>%
         dplyr::mutate(`Relative Abundance (%) by UMI`=`Number of UMIs`/`Total  Number of Genus UMIs`)  %>%
         dplyr::arrange(desc(`Number of UMIs`))%>%
         dplyr::select(!colnames(sample_stat))

 futile.logger::flog.info("6.统计前10微生物的情况")
    top10_genus_sum <-lapply(sample_list, function(sample )
              pathogen_stat %>%
                dplyr::filter(sampleid==sample) %>%
                dplyr::top_n(10,`Number of UMIs`) %>% .$pathogen)%>%unlist%>%unique

    top10_sum_colors<-OESingleCell::SelectColors(c(top10_genus_sum),"blindless")
    top10_sum_colors['others']<-"grey"

    for(sample in sample_list){
        ## top10_umi_stat
        readr::write_tsv(pathogen_stat%>%dplyr::filter(sampleid==sample)%>%dplyr::select(!"sampleid"),
                         glue::glue("{output_dir}/1.Bacteria_reads_Bacteria_UMIs_distribution/{sample}_pathogen_statistic.xls"))
        top10_genus <-pathogen_stat %>%
                dplyr::filter(sampleid==sample)%>%dplyr::arrange(desc(`Number of UMIs`))%>%
                dplyr::top_n(10,`Number of UMIs`) %>%.$pathogen
        top10_umi_stat<- pathogen_stat %>%
                    dplyr::filter(sampleid==sample)%>%
                    dplyr::mutate(pathogen=ifelse(pathogen%in% top10_genus,pathogen,"others")) %>%
                    dplyr::group_by(pathogen) %>%
                    dplyr::summarise(`Number of UMIs` = sum(`Number of UMIs`)) %>%
                    dplyr::mutate(sampleid= glue::glue("{sample}"))%>%
                    dplyr::left_join(sample_stat%>%tibble::rownames_to_column("sampleid")%>%dplyr::filter(sampleid==sample) ,by="sampleid") %>%
                    dplyr::mutate(`Relative Abundance (%) by UMI`=round(`Number of UMIs`*100/`Total  Number of Genus UMIs`,2))%>%
                    dplyr::select(!colnames(sample_stat))  %>%
                    dplyr::mutate(labels = ifelse(`Relative Abundance (%) by UMI`<1 ,NA,paste0(`Relative Abundance (%) by UMI`,"%")))

        if("others" %in% top10_umi_stat$pathogen){cell_levels <- c(top10_genus,"others")}else{cell_levels <-top10_genus }
        top10_umi_stat <- top10_umi_stat %>%
                          dplyr::mutate(pathogen=factor(pathogen,levels=cell_levels))%>%
                          dplyr::arrange(pathogen) %>%dplyr::select(!c("sampleid"))
        readr::write_tsv(top10_umi_stat,glue::glue("{output_dir}/2.Top10_most_dominant_Bacterial_genera/{sample}_UMIs_statistic_for_top10.xls"))
         ## top10_reads_stat
        top10_reads_stat<- pathogen_stat %>%
                    dplyr::filter(sampleid==sample)%>%
                    dplyr::mutate(pathogen=ifelse(pathogen%in% top10_genus,pathogen,"others")) %>%
                    dplyr::group_by(pathogen) %>%
                    dplyr::summarise(`Number of Reads` = sum(`Number of Reads`))
        top10_reads_stat <- top10_reads_stat %>%  dplyr::mutate(pathogen=factor(pathogen,levels=cell_levels))%>%dplyr::arrange(pathogen)
        readr::write_tsv(top10_reads_stat,glue::glue("{output_dir}/2.Top10_most_dominant_Bacterial_genera/{sample}_Reads_statistic_for_top10.xls"))
        pdf(NULL)
        ## 绘图
        plot_top10_piechart_featureplot(data_ob,
                                        sample,
                                        top10_genus,
                                        top10_reads_stat,
                                        top10_umi_stat,
                                        top10_sum_colors,
                                        glue::glue("{output_dir}/2.Top10_most_dominant_Bacterial_genera"))
        dev.off()
        }
    return(data_ob)
}
docstring <- " example1:\\n\\n\\
  scVis  -i  query.rds   -f rds  -o result/pathseq/pathseq_visualization -d rds   --assay SCT  --dataslot data  pathseqvis --pathseq_dir result/pathseq/pathseq_visualization --groupby  sampleid,clusters  --palette  customecol2   "
sub_pathseq <- subparsers$add_parser("pathseqvis",
    description = docstring,
    formatter_class= 'argparse.RawTextHelpFormatter' ,
    #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
    argument_default = "True",
    help = "针对空间转录组pathseq输出结果的可视化")
sub_pathseq$add_argument(
    "--pathseq_dir",
    type = "character",
    default="result/pathseq/pathseq_visualization",
    help = "pathseq解析结果所在目录.[default: %(default)s]")
sub_pathseq$add_argument(
    "-g",
    "--groupby",
    type = "character",
    default = "sampleid,clusters",
    help = "绘制小提琴图采用分组信息.[default: %(default)s]")
sub_pathseq$add_argument(
    "--palette",
    type = "character",
    default = "customecol2",
    help = paste0(
        "绘制前十微生物分布饼图所采用的调色板,颜色数目默认为所有样本前十微生物数目总和.[default: %(default)s]",
        " Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,",
        " tableau20:20,Buen:17,UKBB:18,TF1:17,paired:12"
    ))
args <- commandArgs(TRUE)
if ("pathseqvis" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "pathseqvis") {
        # ===================================================================================================================
        futile.logger::flog.info("step1:按照提供信息，读入rds文件") # ===============
        # for visualization usually only "data" slot is loaded,which can save a lot of time and memory
        suppressMessages(data_ob <- OESingleCell::ReadX(
            input = opt$input,
            informat = opt$informat,
            assays = assays,
            data.use = dataslots, # "data","scale.data"
            reductions = opt$reduct,
            graphs = FALSE,
            images = FALSE, # no graph object needed here
            verbose = FALSE
        ))
        input_dir<-opt$pathseq_dir
        # ===================================================================================================================
        futile.logger::flog.info("step2:汇总各个spot中的微生物reads及UMI总数，并存入seurat对象中")
        if (!file.exists(glue::glue("{output_dir}/1.Bacteria_reads_Bacteria_UMIs_distribution"))){
             dir.create(glue::glue("{output_dir}/1.Bacteria_reads_Bacteria_UMIs_distribution"), recursive = T)
          }
        if (!file.exists(glue::glue("{output_dir}/2.Top10_most_dominant_Bacterial_genera/piechart"))){
            dir.create(glue::glue("{output_dir}/2.Top10_most_dominant_Bacterial_genera/piechart"), recursive = T)
        }
        data_ob<- statistic_info(data_ob, input_dir,output_dir)
        # ===================================================================================================================
        futile.logger::flog.info("step3.绘制各个spot中的Bacteria_reads、Bacteria_UMIs空间分布图及小提琴图")
        futile.logger::flog.info("绘制空间分布图")
        sample_list<-list.files(path = input_dir, pattern ="visium.raw.readnamepath") %>%
                 stringr::str_replace_all(pattern = '.visium.raw.readnamepath',replacement = "")
        lapply(sample_list, function(sample)
            plot_featureplot(data_ob,
            features="Bacteria_reads",
            sample=sample,
            output_dir=glue::glue("{output_dir}/1.Bacteria_reads_Bacteria_UMIs_distribution/featureplot")))
        lapply(sample_list, function(sample)
            plot_featureplot(data_ob,
            features="Bacteria_UMIs",
            sample=sample,
            output_dir=glue::glue("{output_dir}/1.Bacteria_reads_Bacteria_UMIs_distribution/featureplot")))
        futile.logger::flog.info("绘制小提琴图")
        data_ob_sub<-subset(data_ob,subset=Bacteria_reads >0)
        for(group in unlist(stringr::str_split(opt$groupby,pattern =","))){
            plot_vlnplot(data_ob_sub,
                         feature="Bacteria_reads",
                         group.by=group,
                         output_dir=glue::glue("{output_dir}/1.Bacteria_reads_Bacteria_UMIs_distribution/vlnplot"))
            plot_vlnplot(data_ob_sub,
                         feature="Bacteria_UMIs",
                         group.by=group,
                         output_dir=glue::glue("{output_dir}/1.Bacteria_reads_Bacteria_UMIs_distribution/vlnplot"))
       }
            ## 保存seurat对象
            saveRDS(data_ob, file = glue::glue("{output_dir}/seurat_pathseq.rds"))
           # ===================================================================================================================
            futile.logger::flog.info("step5: output session information")
                write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
                quit()
    }
}
