sub_sankey <- subparsers$add_parser("sankey", help = "visualize the bunch of features in a heatmap using ComplexHeatmap.")
sub_sankey$add_argument("-s",
                        "--sankeyfile",
                        type = "character",
                        help = "the raw file for sankey plot")
args <- commandArgs(TRUE)
if ("sankey" %in% args) {
  opt <- intial_setting()
  if(opt$sub_name == "sankey" ){
    futile.logger::flog.info("step1 read sankeyfiles ") #================================================
    data<-read.table(opt$sankeyfile,sep="\t",header=1)
    groups<-colnames(data)
    futile.logger::flog.info("step2 prepare for links ") #================================================
    data_long_list<-list()
    data_long_list<-lapply( 1:(length(groups)-1),function(i) {
                        d<-data %>%
                           dplyr::select(groups[i],groups[i+1])%>%
                           table()%>%
                           as.data.frame() %>%
                           dplyr::mutate(group=groups[i])%>%
                           dplyr::filter(Freq>0)
                        colnames(d) <- c("source", "target", "value","group")
                        d
            })
    data_long<-data.frame(do.call(rbind,data_long_list))
    futile.logger::flog.info("step3 prepare for nodes ") #================================================
    nodes_list<-list()
    nodes_list<-lapply( 1:(length(groups)),function(i) {
                 d<-data %>%
                    dplyr::select(groups[i]) %>%
                    dplyr::distinct()%>%
                    dplyr::mutate(group=groups[i]) %>%
                    dplyr::rename(name=groups[i])
                 d
            })
    nodes<-data.frame(do.call(rbind,nodes_list))
    rownames(nodes)<-nodes$name
    color<-OESingleCell::SelectColors(as.character(nodes$name)%>%
                                                        sort()%>%
                                                         unique,"blindless")%>%
            as.data.frame()
    colnames(color)<-c("colors")
    nodes<-cbind(nodes,color)
    rownames(nodes) <- 1:dim(nodes)[1]
    row_numbers<-nodes$group%>%table()%>%max

    data_long$IDsource <- match(data_long$source, nodes$name) - 1
    data_long$IDtarget <- match(data_long$target, nodes$name) - 1


    colors <- paste(sapply(nodes$colors, function(x) {
                    paste0("d3.rgb(",
                           paste(c(col2rgb(x), 0.5), collapse = "," ), ")") }), collapse = ", ")
    colorJS <- paste0('d3.scaleOrdinal([', colors, '])')
    futile.logger::flog.info("step4 sankey plotting ") #================================================
    network <- networkD3::sankeyNetwork(Links = as.data.frame(data_long),
                                        Nodes = nodes,
                                        Source = "IDsource",
                                        Target = "IDtarget",
                                        Value = "value",
                                        NodeID = "name",
                                        height=row_numbers*50,
                                        width=length(groups)*300,
                                       # NodeGroup="group",
                                        LinkGroup="source",
                                        colourScale = colorJS,
                                        sinksRight = FALSE,
                                        nodeWidth = 60,
                                        fontSize = 12,
                                        iterations = 0,
                                        nodePadding = 10)
    networkD3::saveNetwork(network,
                           paste0(output_dir, "/", "cluster_celltype_sankeyNetwork.html"),
                           selfcontained = TRUE)
    webshot::webshot(paste0(output_dir, "/", "cluster_celltype_sankeyNetwork.html"),
                     paste0(output_dir, "/", "cluster_celltype_sankeyNetwork.png"))
    webshot::webshot(paste0(output_dir, "/", "cluster_celltype_sankeyNetwork.html"),
                     paste0(output_dir, "/", "cluster_celltype_sankeyNetwork.pdf"))
    futile.logger::flog.info("step5 saving session information ") #=============================================
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}