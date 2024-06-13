args <- commandArgs(TRUE)
if ( "summarize"  %in% args ){
  opt<-intial_setting()
  DATA <- OESingleCell::colData(data_ob)[,c(opt$groups, opt$propby)]  %>%
              dplyr::group_by( .dots= c(opt$groups, opt$proby)) %>%
              dplyr::summarize(count = n()) %>%
              dplyr::mutate(freq = (count / sum(count)) * 100)
  write.table(as.data.frame(DATA), file.path(output_dir,file="clust_cond_freq_info.xls"),
              sep="\t",col.names=T, row.names =F)
  clust_sum_all = OESingleCell::PlotAbundances(data_ob,
                                               prop.by = opt$propby,
                                               group.by = opt$groups,
                                               split.by = opt$facetby,
                                               method = opt$plot,
                                               ncol = ifelse(opt$plot=="pie",4,1),
                                               cols= OESingleCell::SelectColors(1:length(unique(data_ob@meta.data[,opt$propby])))
  )
  ggsave(file.path(output_dir,paste0("groupby-",opt$groups,"_resolution-",resolution,"summary_plot.pdf",collapse="")),plot=clust_sum_all)
  ggsave(file.path(output_dir,paste0("groupby-",opt$groups,"_resolution-",resolution,"summary_plot.png",collapse="")),dpi=1000, plot = clust_sum_all)
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}