sub_barcode_rank <- subparsers$add_parser("barcode_rank_reads", help = "绘制barcode拐点图")
sub_barcode_rank$add_argument( "--kraken2_output", type = "character", help = "kraken  output file: result/Kraken2/A416_kraken.output")
sub_barcode_rank$add_argument( "--sc_taxonomy", type = "character", help = "sc_taxonomy file: result/clean_data/A416_sc_taxonomy_10000_bc.report ")
sub_barcode_rank$add_argument("-s", "--prefix", type = "character", help = "sample name id,eg:XHR-bm-ctr-74")
args <- commandArgs(TRUE)
if ("barcode_rank_reads" %in% args) {
  opt <- intial_setting()
  if(opt$sub_name == "barcode_rank_reads" ){
    futile.logger::flog.info("step1:导入Kraken2及sc_taxonomy结果并整理")#=================================================
    sc_taxonomy<-data.table::fread(opt$sc_taxonomy,header = T)
    kraken_output<-data.table::fread(opt$kraken2_output,header = F) %>%
                   #dplyr::filter(V1=="C")%>%
                   tidyr::separate(V2,sep="_",into=c("reads","barcode"))%>%
                   dplyr::count(barcode)
    barcode_read<-dplyr::left_join(kraken_output,sc_taxonomy,by="barcode")%>%
                  dplyr::rename(nCount_Reads='n')%>%
                  dplyr::mutate(kept=ifelse(!is.na(species),"True","False"))%>%
                  dplyr::mutate(kept=factor(kept,levels=c("True","False")))
    futile.logger::flog.info("step2:绘制Rank_Barcode-vs-Gene_count")#===================================================
    barcode_read$rank <- rank(-barcode_read$nCount_Reads)
    plot<- ggplot(barcode_read, aes(x = rank, y = nCount_Reads)) +
            geom_point(shape = 1, aes(colour = kept)) +
            scale_x_log10(labels = scales::number) +
            scale_y_log10(labels = scales::number) +
            scale_colour_manual(values = c("#3558AE","grey")) +
            ggtitle(glue::glue("{opt$prefix}_Reads"))+
            ylab("Total Counts") +
            xlab("Rank Barcode")+
            theme_bw(base_size = 20)+
            theme(panel.grid=element_blank(),legend.position = c(.12, .12),plot.title = element_text(hjust = 0.5) )  +
            theme(legend.background = element_rect(color = "black", linetype = "solid", size = 0.1))
    OESingleCell::save_ggplots(plot,filename = glue::glue("{output_dir}/{opt$prefix}_Barcode_vs_Gene_Counts"),width=8,height=8)
    futile.logger::flog.info("step4：saving session information ") #=============================================
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
}}
