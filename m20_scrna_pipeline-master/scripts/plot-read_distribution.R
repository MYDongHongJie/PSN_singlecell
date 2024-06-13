sub_read_distribution<- subparsers$add_parser("read_distribution", help = "获取read分布饼图")
sub_read_distribution$add_argument("-r", "--rna_matrics", type = "character", help = "picard output rna_matrics.txt")
sub_read_distribution$add_argument("-s", "--prefix", type = "character", help = "sample name id ")
args <- commandArgs(TRUE)
if ("read_distribution" %in% args) {
  opt <- intial_setting()
  if(opt$sub_name == "read_distribution" ){
    futile.logger::flog.info("step1:导入rna_matrics结果并整理数据")#========================================================
    rt <- read.table(opt$rna_matrics, nrows=1, header=T, fill=T)
    df <- data.frame( #"PF Bases" = rt$PF_BASES,
                     #"PF Aligned Bases" = rt$PF_ALIGNED_BASES,
                     "Coding Bases" = rt$PCT_CODING_BASES,
                     "UTR Bases" = rt$PCT_UTR_BASES,
                     "Intronic Bases" = rt$PCT_INTRONIC_BASES,
                     "Intergenic Bases" = rt$PCT_INTERGENIC_BASES,
                     "Ribosomal Bases" = rt$PCT_RIBOSOMAL_BASES,
                     #"Median CV Coverage" = rt$MEDIAN_CV_COVERAGE,
                     #"Median 5' Bias" = rt$MEDIAN_5PRIME_BIAS,
                     #"Median 3' Bias" = rt$MEDIAN_3PRIME_BIAS,
                     # "Median 5' to 3' Bias" = rt$MEDIAN_5PRIME_TO_3PRIME_BIAS,
                     #"% Stranded" = rt$PCT_CORRECT_STRAND_READS * 100,
                     #"% rRNA bases" = (rt$RIBOSOMAL_BASES/rt$PF_BASES)*100,
                     check.names= F)
    df <- as.data.frame(t(df)) %>% tibble::rownames_to_column()
    names(df) <- c("Metric", "Value")
    dfx <- data.frame(" "="Bases",
                      "Coding" = df$Value[which(df$Metric=="Coding Bases")],
                      "UTR" = df$Value[which(df$Metric=="UTR Bases")],
                      "Intronic" = df$Value[which(df$Metric=="Intronic Bases")],
                      "Intergenic" = df$Value[which(df$Metric=="Intergenic Bases")],
                      "Ribosomal" = df$Value[which(df$Metric=="Ribosomal Bases")], check.names=F)%>%
            tibble::column_to_rownames(" ")%>%
            t%>%
            as.data.frame()%>%
            tibble::rownames_to_column("Cate")%>%
            dplyr::arrange(desc(Cate)) %>%
            dplyr::mutate(prop = Bases*100 ) %>%
            dplyr::mutate(label =ifelse(prop<5," ",paste0(as.character(round(prop,1)),"%")))%>%
            dplyr::mutate(ypos = cumsum(prop)- 0.5*prop )
    futile.logger::flog.info("step2:绘制reads分布饼图")#=============================================
    sum_colors <- OESingleCell::SelectColors(dfx$Cate,"ditto")
    plot<-  ggplot(dfx, aes(x="", y=prop, fill=Cate)) +
            geom_bar(stat="identity", width=0.5, color="white") +
            coord_polar("y", start=0) +
            theme_void(base_size = 22) +
            ggtitle(opt$prefix)+
            #theme(legend.position="none") +
            geom_text(aes(y = ypos, label =label ), color = "black", size=3) +
            scale_fill_manual(values = sum_colors) +
            theme(#legend.position = "bottom",
                  legend.title = element_blank(),
                  plot.title = element_text(size=20,hjust=0.5),
                  plot.margin = margin(t = 1 ,  # 顶部边缘距离
                                        r = 1 ,  # 右边边缘距离
                                        b = 1 ,  # 底部边缘距离
                                        l = 40)) # 左边边缘距离
    OESingleCell::save_ggplots(plot,filename =glue::glue("{output_dir}/reads_distribution_for_{opt$prefix}"),height = 8,width = 8)
    futile.logger::flog.info("step3：saving session information ") #=============================================
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
}}
