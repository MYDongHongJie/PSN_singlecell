##wnn analysis
# Title     : TODO
# Objective : TODO
# Created by: sunku
# Created on: 2021/8/11
# Weighted Nearest Neighbor Analysis
#===================================================
docstring <- "example:\\n\\n\\
sctool -i seurat.h5seurat -f h5seurat -o ./ -d h5seurat --assay ATAC --subassay ATAC --dataslot data wnn --reduction RNA_pca,ATAC_lsi --clusteringuse snn --resolution 0.2 --palette blindless"
sub_wnn <- subparsers$add_parser("wnn",
                                 description = docstring,
                                 formatter_class = "argparse.RawTextHelpFormatter",
                                 argument_default = "True",
                                 help = "Weighted Nearest Neighbor Analysis.")
sub_wnn$add_argument("-d",
                     "--reduction",
                     type = "character",
                     default = NULL,
                     help = "scATAC and scRNA reduction methods,such as ATAC_lsi,RNA_pca")
sub_wnn$add_argument("-r",
                     "--resolution",
                     type = "double",
                     default = 0.4,
                     help = "vaule used to set the resolution of cluster distiguish, use a value above(below)1.0 if you want to obtain a larger(smaller) number of communities.[default: %(default)s]")
sub_wnn$add_argument("-c",
                     "--clusteringuse",
                     type = "character",
                     default = "snn",
                     help = "the only supported clustering algorithms currently:snn,louvain,leiden. For other methods, use subcommand 'clusterx'.[default: %(default)s].")
sub_wnn$add_argument("--palette",
                     type = "character",
                     default = "blindless",
                     help = "the discrete color schema for each cell cluster, blindless as default. Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,tableau20:20, paired:12, purplegray12:12.[default: %(default)s]")

if ("wnn" %in% args) {
  opt <- intial_setting()
  ###step1 load data====================================================================================================
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots,
                                 verbose = F)
  message("Loading Data object Finished!")
  ###step2 wnn reduction================================================================================================
  reductions <- unlist(strsplit(opt$reduction, ",", perl = T))
  reductions1_dim <- ncol(OESingleCell::Embeddings(data_ob, reduction = reductions[1]))
  reductions2_dim <- ncol(OESingleCell::Embeddings(data_ob, reduction = reductions[2]))
  reduction_name <- paste("WNN",
                          unlist(strsplit(reductions, split = "_"))[2],
                          unlist(strsplit(reductions, split = "_"))[4],
                          "umap",
                          sep = "_")
  data_ob <- Seurat::FindMultiModalNeighbors(data_ob,
                                               reduction.list = as.list(reductions),
                                               dims.list = list(1:reductions1_dim, 1:reductions2_dim))
  data_ob <- Seurat::RunUMAP(data_ob,
                             nn.name = "weighted.nn",
                             reduction.name = reduction_name,
                             reduction.key = "WNNumap_")
  clustering.alg <- c("snn" = 1, "louvain" = 2, "slm" = 3, "leidn" = 4)
  data_ob <- Seurat::FindClusters(data_ob,
                                  resolution = opt$resolution,
                                  graph.name = "wsnn",
                                  algorithm = clustering.alg[opt$clusteringuse],
                                  verbose = F)
    #===================================================================================================================
    # step4 save the clustering results
    # ==================================================================================================================
    ## 4.1 renane cluster id start from 1
    new_cluster_id <- Seurat::Idents(data_ob) %>%  as.vector () %>% as.numeric() %>% +1 %>% as.factor()
    names(new_cluster_id) <-Seurat::Idents(data_ob) %>% names()
    Seurat::Idents(data_ob) <- new_cluster_id
    cluster_result_colname <-  paste("wnn","res",opt$resolution,sep = ".")
    ## 4.2  as default the seurat will store the clustering results in the data_object@ident to keep the clutering results using
    # the specified resolution to data_object@metadata for reuse
    data_ob[[cluster_result_colname]] <- Seurat::Idents(data_ob)
    data_ob[["clusters"]] <- Seurat::Idents(data_ob)
    reordered_cell_count_by_cluster <- table(Seurat::Idents(data_ob))
    cell_count_labels <- glue::glue("{names(reordered_cell_count_by_cluster)}-{reordered_cell_count_by_cluster} cells")
    ggdimplot <- Seurat::DimPlot(object = data_ob,
                                 reduction = reduction_name,
                                 dims = c(1, 2),
                                 label = T,
                                 label.size =5,
                                 group.by= cluster_result_colname)+
                 ggplot2::labs(title = "") +
                 ggplot2::scale_colour_manual( values = OESingleCell::SelectColors(1:length(unique(Seurat::Idents(data_ob))),
                                                                                   palette = opt$palette),
                                               breaks=levels(Seurat::Idents(data_ob)),
                                               labels = cell_count_labels)+
                theme(legend.text=element_text(size=14))
    ## image width-height ratio
    dim_1_range<- Seurat::Embeddings(data_ob, reduction = glue::glue(reduction_name))[,1] %>% range()
    dim_2_range<- Seurat::Embeddings(data_ob, reduction = glue::glue(reduction_name))[,2] %>% range()
    width_height_ratio<- (dim_1_range[2]-dim_1_range[1])/(dim_2_range[2]-dim_2_range[1])
    width_height_ratio <-ifelse( width_height_ratio >1.3,1.3,width_height_ratio ) ##避免比例差距过大
    OESingleCell::save_ggplots(glue::glue("{output_dir}/wnn_groupby_cluster_resolution{opt$resolution}_plot"),
                               plot = ggdimplot,
                               width= 8*width_height_ratio ,
                               height=8,
                               dpi = 1000)
    ##4.3  update or save seurat object
    if ( as.logical(opt$update) ){
      SeuratDisk::UpdateH5Seurat(file = opt$input,
                                 object = data_ob,
                                 reduction = reduction_name)
    }else{
      OESingleCell::SaveX(data_ob,
                          output = opt$output,
                          update = FALSE,
                          outformat = opt$outformat,
                          prefix = opt$prefix)
    }
    ##4.4 save session informations
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
