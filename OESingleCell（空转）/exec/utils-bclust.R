sub_cluster <- subparsers$add_parser("bclust", help = "one stop run of data dimension reduction and clustering.")
sub_cluster$add_argument(
    "-t",
    "--components",
    type = "integer",
    default = NULL,
    help = "the appropriate number of dimensions for primary reduction. Recommendations:pca:30,mnn:10[default: %(default)s]"
)
sub_cluster$add_argument(
    "-p",
    "--perplexity",
    type = "integer",
    default = 30,
    help = "The value of the perplexity only availiable for tSNE.[default: %(default)s]"
)
sub_cluster$add_argument(
    "-m",
    "--metadata",
    type = "character",
    default = NULL,
    help = "[OPTIONAL]the additional sample metadata which must include sample id in this assay design.[default: %(default)s]"
)
sub_cluster$add_argument(
    "-c",
    "--clusteringuse",
    type = "character",
    default = "snn",
    help = "the only supported clustering algorithms currently:snn,louvain,leiden. For other methods, use subcommand 'clusterx'.[default: %(default)s]."
)

sub_cluster$add_argument(
    "--leiden_method",
    type = "character",
    default = "matrix",
    help = "Method for running leiden (defaults to matrix which is fast for small datasets). Enable method： igraph  to avoid casting large data to a dense matrix.[default: %(default)s]."
)
sub_cluster$add_argument(
    "-r",
    "--resolution",
    type = "double",
    default = 0.4,
    help = "vaule used to set the resolution of cluster distiguish, use a value above(below)1.0 if you want to obtain a larger(smaller) number of communities.[default: %(default)s]"
)
sub_cluster$add_argument(
    "-s",
    "--pointsize",
    type = "double",
    default = 0.5,
    help = "[OPTIONAL]the point size in the plot..[default: %(default)s]"
)
sub_cluster$add_argument(
    "-b",
    "--batchid",
    type = "character",
    default = "batchid",
    help = "[OPTIONAL]the batch information column name of sample in the metadata.[default: %(default)s]"
)
sub_cluster$add_argument(
    "--reduct1",
    type = "character",
    default = "pca",
    choices = c("ica", "cca", "pca", "lsi", "mnn", "tsne", "pca,harmony", "lsi,harmony","Flt-SNE", "UMAP", "swne", "diffusion", "phate", "fa2"),
    #choices = c("ica", "cca", "pca", "lsi", "mnn", "tsne", "harmony", "Flt-SNE", "UMAP", "swne", "diffusion", "phate", "fa2"),
    help = "the primary dimension reduction results used for clustering. Now ica,cca,pca,lsi,mnn,tsne,harmony,Flt-SNE,UMAP,swne,diffusion,phate,fa2 are supported.[default: %(default)s]"
)
sub_cluster$add_argument(
    "--reduct2",
    type = "character",
    choices=c("tsne","Flt-SNE","umap","swne","diffusion","phate","fa2" ),
    help = "the secondary reduction method usually used to embed community. Now tsne,Flt-SNE,UMAP,swne,diffusion,phate,fa2 are supported. If NULL supplied, NO secondary reduction will be carried out.[default: %(default)s]")
sub_cluster$add_argument("--palette",
     type = "character",
     default = "customecol2",
     choices=c("blindless","customecol2","cold","glasbey","ditto","alphabet","alphabet2","colx22","cyclic","tableau20", "paired", "purplegray12"),
    help = "the discrete color schema for each cell cluster, blindless as default. Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,tableau20:20, paired:12, purplegray12:12.[default: %(default)s]"
)
# sub_cluster$add_argument("--no-update", dest = "update", action="store_false",
#              help="produce new data object on disk. Only availiable for h5seurat input.")
#              help="produce new data object on disk. Only availiable for h5seurat input.")
sub_cluster$add_argument(
    "--rerun",
    default = "TRUE",
    type = "character",
    help = "Rerun the primary reduction in case of parameter changes.[default: %(default)s]"
)

sub_cluster$add_argument(
    "--dosubset",
    default = FALSE,
    type = "character",
    help = "did subset to get data .[default: %(default)s]"#是否已做过subset,如果填写了predicate，也需要填TRUE
)
sub_cluster$add_argument(
    "--assay4SCT",
    default = "RNA",
    type = "character",
    help = "RNA or Spatial for SCT if dosubset is TRUE.[default: %(default)s]"#是否已做过subset,如果填写了predicate，也需要填TRUE
)

# sub_cluster$add_argument("--no-rerun", dest = "rerun", action="store_false",
#              help="NO rerun the primary reduction."
# sub_cluster$add_argument(
#     "--which_cells",
#     type = "character",
#     default = NULL,
#     help = "The subset of cluster ids used for subtyping."
# )
args <- commandArgs(TRUE)
if ("bclust" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "bclust") {
        # =================================================================================================================
        futile.logger::flog.info("step1: read the specified assay and data slot in data object into memory")
        # =================================================================================================================
         suppressMessages(data_ob <- OESingleCell::ReadX(
             input = opt$input,
             informat = opt$informat,
             assays = assays,
             data.use = dataslots,
             verbose = F
        ))
        ## For mnn run ,to save a new h5seurat object,must import all previous data
        # assays<-assays[1]
        Seurat::DefaultAssay(data_ob) <- assays
        futile.logger::flog.info("Loading Data object Finished!")

        if (!is.null(opt$metadata)) {
            futile.logger::flog.info("update the metedata in the data_ob@meta.data with new additional sample metadata")
            additional_metadata <- read.csv(opt$metadata, sep = ",", header = T)
            rownames(additional_metadata) <- additional_metadata$sampleid
            data_ob <- Seurat::UpdatCellMeta(data_ob, metadata = additional_metadata, cell.delm = "-")
        }

        if (is.null(OESingleCell::colData(data_ob)[["clusters"]])) {
            # data_ob = StashIdent(data_ob, save.name = "clusters")
            data_ob[["clusters"]] <- Seurat::Idents(data_ob)
        } else {
            # if it is the first time to run this script on the data object, the clusters here is actually the sample index.
            # After running this script, the cluster id will be overwrite with the actual cell cluster id.
            data_ob <- Seurat::SetIdent(data_ob, value = "sampleid")
        }

        # get the subset of cells used for visualization if necessay subset cells on the expression matrix using logical exprssion on features
        if (!is.null(opt$predicate)) {
            #示例： predicate sampleid %in% c(\"A\",\"B\")
            futile.logger::flog.info("get the subset of cells used for visualization if necessay subset cells on the expression matrix using logical exprssion on features")
            df <- OESingleCell::colData(data_ob)
            desired_cells <- subset(df, eval(parse(text = opt$predicate)))
            data_ob <- data_ob[, rownames(desired_cells)]
        }
        # ===================================================================================================================
        futile.logger::flog.info("step2: Determine statistically significant principal components for clustering")
        # ===================================================================================================================
        scale_data <- Seurat::GetAssayData(data_ob,
                                           assay = Seurat::DefaultAssay(data_ob),
                                           slot = "scale.data")
        reduct1s <- as.character(unlist(stringr::str_split(opt$reduct1, ",")))
        ## Tips: 如果对数据进行了抽提，重新进行sctranscform，重新scale, 获取新的高变基因
        if (opt$dosubset) {
            data_ob <- Seurat::SCTransform(data_ob,
                                           vars.to.regress = "percent.mito",
                                           verbose = FALSE,
                                           assay = opt$assay4SCT) %>%
              Seurat::FindVariableFeatures()
        }
        ##Tips: 第一次降维，考虑到harmony的特殊性（需要基于PCA或者LSI的降维结果），因此会根据reduct1提供结果进行循环处理
        for (reduct1 in c(reduct1s)) {
            ## 如果对象中已存在该降维结果，并且--rerun设置为FALSE,执行以下操作，默认情况下，会重新运行
            if (reduct1 %in% sub(glue::glue(assays, "_"), "", OESingleCell::Reductions(data_ob)) & as.logical(opt$rerun) == F) {
                #===============================================================================================================
                ## 如果前面已经运行过降维，且本次参数rerun设置为False,则采用之前结果，并获取最佳PCS
                # if previous reduction has been calculated ,donot rerun
                if (reduct1 %in% c("pca", "cca", "harmony", "ica", "lsi")) {
                    # if the prevous reduction is primary reduction check the components find the optimal components for secondary reduction
                    optimal_pc <- tryCatch(
                        expr = Seurat::Command(
                            data_ob,
                            command = glue::glue("FindNeighbors.{opt$assay}.{reduct1}"),
                            value = "dims"
                        ),
                        error = function(...) {
                            return(NULL)
                        }
                    )
                    if (is.null(optimal_pc)) {
                        # if previous optimal components not available
                        futile.logger::flog.info("NO previous optimal components is AVAILABLE, the optimal components number be detected automatically.")
                        elb <- Seurat::ElbowPlot(data_ob, reduction = glue::glue("{assays}_{reduct1}"))
                        optimal_pc <- OESingleCell::FindElbow(
                            elb$data,
                            reduction = glue::glue("{assays}_{reduct1}"),
                            method = "slope",
                            ndims = 30
                        )
                        # suppressWarnings({ Misc(data_ob, "optimal_pc") = optimal_pc})
                    }
                } else { # reduct1 is not primary reduction and previouly run without primary reduction
                    # the exception is mnn, it's better to be specified from command line with relative low value
                    optimal_pc <- min(opt$components, ncol(OESingleCell::Embeddings(data_ob, reduction = glue::glue("{assays}_{reduct1}"))))
                }
            } else {
                #===============================================================================================================
                ## 如果前面未运行过降维，或者本次参数rerun设置为True,则重新进行降维，并获取最佳PCS
                # this can be primary reduction or secondary reduction
                futile.logger::flog.info("NO specified primary reduction found or forced to rerun! \n Reduction begins!")
                futile.logger::flog.info(paste0("Running ", reduct1, " Dimension Reduction"))
                dim_outdir <- file.path(output_dir, paste0(reduct1, "_Dimension_Reduction"))
                if (!dir.exists(dim_outdir)) { dir.create(dim_outdir) }
                set.seed(2022)
                data_ob <- OESingleCell::RunDimReduc(
                    data_ob,
                    reduct1.use = ifelse(
                        reduct1 == "harmony",
                        yes =  dplyr::case_when(
                            assays == "ATAC" ~ glue::glue("{assays}_lsi"),
                            assays == "SCT" ~ glue::glue("{assays}_pca"),
                            assays == "RNA" ~ glue::glue("{assays}_pca"),
                            assays == "integrated" ~ glue::glue("{assays}_pca")
                        ),
                        no = reduct1
                    ),
                    reduct2.use = reduct1,
                    feature.use = Seurat::VariableFeatures(data_ob),
                    perplexity = opt$perplexity,
                    assay.use = Seurat::DefaultAssay(data_ob),
                    batch.use = opt$batchid,
                    npcs.use = opt$components,
                    reduction.name = glue::glue("{assays}_{reduct1}")
                )
                ## 兼容老版本脚本未存储rawbc情况
                if (!"rawbc" %in% colnames(data_ob@meta.data)) { data_ob[["rawbc"]] <- data_ob[["orig.ident"]] }

                Seurat::FetchData(data_ob, vars = c("rawbc", paste0(Seurat::Key(data_ob)[glue::glue("{assays}_{reduct1}")], 1:2))) %>%
                  dplyr::rename("Barcode" = "rawbc") %>%
                    readr::write_csv(file = glue::glue("{dim_outdir}/{reduct1}_Dimension_Reduction_coordination.csv"))

                optimal_pc <- opt$components

                if (reduct1 %in% c("pca", "cca", "harmony", "ica", "lsi")) { # if the prevous reduction is primary reduction
                    if (is.null(optimal_pc)) {
                        optimal_pc <- OESingleCell::FindElbow(
                            data_ob,
                            reduction = glue::glue("{assays}_{reduct1}"),
                            method = "slope",
                            ndims = 30
                        )
                    }
                }
            }
            # in case of dimension number out of the actual.
            components2use <- min(optimal_pc, ncol(OESingleCell::Embeddings(data_ob, reduction = glue::glue("{assays}_{reduct1}"))))
            suppressWarnings({
                Seurat::Misc(data_ob, glue::glue("{assays}_{reduct1}_components_num")) <- opt$components
            })
            suppressWarnings({
                Seurat::Misc(data_ob, glue::glue("{assays}_{reduct1}_optimal_pc")) <- components2use
            })
            futile.logger::flog.info(glue::glue("The dimension number for use: {components2use }"))
            # opt$components <- components2use
        }
        # after primary reduction, the scale.data will not be used, so as to delete it for saveing memory
        if (tolower(assays[1]) %in% c("rna", "Spatial", "SCT") &
          !as.logical(opt$update) &
          opt$informat != "h5seurat") {
            data_ob[[Seurat::DefaultAssay(data_ob)]] <- Seurat::SetAssayData(
              data_ob[[Seurat::DefaultAssay(data_ob)]],
              slot = "scale.data",
              new.data = scale_data
            )
        }
        invisible(gc(full = T, verbose = F))
        # ==============================================================================================================
        # Clustering with the reduction results using different clustering methods======================================
        futile.logger::flog.info(glue::glue("Beginning to group the cells using algorithm {opt$clusteringuse} with top {components2use} {reduct1} components!"))
        clustering.alg <- c("snn" = 1, "louvain" = 2, "slm" = 3, "leiden" = 4)
        ## ref:https://satijalab.org/seurat/reference/findclusters
        ## 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm
        data_ob <- Seurat::FindNeighbors(
          data_ob,
          reduction = glue::glue("{assays}_{reduct1}"),
            dims = 1:components2use,
            features = Seurat::VariableFeatures(data_ob),
            nn.eps = 0,
            force.recalc = T,
            verbose = F
        )
        data_ob <- Seurat::FindClusters(
            object = data_ob,
            resolution = opt$resolution,
            method =  opt$leiden_method  ,
            algorithm = clustering.alg[opt$clusteringuse],
            verbose = F
        )
        futile.logger::flog.info("Clustering Finished!")
        # ==================================================================================================================
        futile.logger::flog.info("step3 secondary reduction: used to detect the community in graph if graph-based clustering method used")
        # ==================================================================================================================
        if (!is.null(opt$reduct2)) {
            futile.logger::flog.info(glue::glue("Beginning {opt$reduct2} Dimension Reduction"))
            output_dir <- glue::glue("{output_dir}/{opt$reduct2}_Dimension_Reduction")
            if (!dir.exists(output_dir)) {
                dir.create(output_dir)
            }
            data_ob <- OESingleCell::RunDimReduc(
                data_ob,
                reduct1.use = glue::glue("{assays}_{reduct1}"),
                reduct2.use = opt$reduct2,
                feature.use = Seurat::VariableFeatures(data_ob),
                assay.use = Seurat::DefaultAssay(data_ob),
                batch.use = opt$batchid,
                npcs.use = components2use,
                reduction.name = glue::glue("{assays}_{opt$reduct2}")
            )
            Seurat::FetchData(data_ob, vars = c("rawbc", paste0(Seurat::Key(data_ob)[glue::glue("{assays}_{opt$reduct2}")], 1:2))) %>%
              dplyr::rename("Barcode" = "rawbc") %>%
                readr::write_csv(file = glue::glue("{output_dir}/{opt$reduct2}_Dimension_Reduction_coordination.csv"))

        } else {
            opt$reduct2 <- reduct1
            output_dir <- glue::glue("{output_dir}/{opt$reduct2}_Dimension_Reduction")
        }
        # ===================================================================================================================
        futile.logger::flog.info("step4 save the clustering results")
        # ==================================================================================================================
        futile.logger::flog.info("4.1 renane cluster id start from 1")

        new_cluster_id <- Seurat::Idents(data_ob) %>%
            as.vector() %>%
            as.numeric() %>%
            + ifelse(opt$clusteringuse=="snn",1,0) %>%
            as.factor()
        names(new_cluster_id) <- Seurat::Idents(data_ob) %>% names()
        Seurat::Idents(data_ob) <- new_cluster_id

        cluster_result_colname <- glue::glue("{Seurat::DefaultAssay(data_ob)}.{opt$reduct2}.res.{opt$resolution}")
        ## 4.2  as default the seurat will store the clustering results in the data_object@ident to keep the clutering results using
        # the specified resolution to data_object@metadata for reuse
        # data_ob <- Seurat::StashIdent(object = data_ob, save.name = cluster_result_colname)
        data_ob[[cluster_result_colname]] <- Seurat::Idents(data_ob)
        # data_ob <- Seurat::StashIdent(object = data_ob, save.name = "clusters")
        data_ob[["clusters"]] <- Seurat::Idents(data_ob)
            
        reordered_cell_count_by_cluster <- table(Seurat::Idents(data_ob))
        if(is.null(Seurat::Images(data_ob))){
          cell_count_labels <- glue::glue("{names(reordered_cell_count_by_cluster)}-{reordered_cell_count_by_cluster} cells")
        }else{
          cell_count_labels <- glue::glue("{names(reordered_cell_count_by_cluster)}-{reordered_cell_count_by_cluster} spots")
        }
        ggdimplot <- Seurat::DimPlot(
            object = data_ob,
            reduction = glue::glue("{assays}_{opt$reduct2}"),
            dims = c(1, 2),
            label = T,
            label.size = 5,
            group.by = cluster_result_colname
        ) +
            ggplot2::labs(title = "") +
            ggplot2::scale_colour_manual(
                values = OESingleCell::SelectColors(
                    1:length(unique(Seurat::Idents(data_ob))),
                    palette = opt$palette
                ),
                breaks = levels(Seurat::Idents(data_ob)),
                labels = cell_count_labels
            ) +
            ggplot2::theme(legend.text = ggplot2::element_text(size = 14))
        ## image width-height ratio
        dim_1_range <- Seurat::Embeddings(data_ob, reduction = glue::glue("{assays}_{opt$reduct2}"))[, 1] %>% range()
        dim_2_range <- Seurat::Embeddings(data_ob, reduction = glue::glue("{assays}_{opt$reduct2}"))[, 2] %>% range()
        width_height_ratio <- (dim_1_range[2] - dim_1_range[1]) / (dim_2_range[2] - dim_2_range[1])
        width_height_ratio <-ifelse( width_height_ratio >1.3,1.3,width_height_ratio )
        OESingleCell::save_ggplots(
            glue::glue("{output_dir}/{opt$reduct2}_groupby_cluster_resolution{opt$resolution}_plot"),
            plot = ggdimplot,
            width = 8 * width_height_ratio + max(stringr::str_length(cell_count_labels)) / 5,
            height = 8,
            dpi = 1000
        )

        ## 4.3  update or save seurat object
        if (as.logical(opt$update)) {
            # OESingleCell::SaveX(data_ob, output = opt$input,update = TRUE,
            #                     outformat = opt$outformat, graphs = Seurat::Graphs(data_ob),
            #                     reduction = c(reduct1, reduct2), commands = Seurat::Command(data_ob),
            #                     misc = NULL, tools = NULL )
            SeuratDisk::UpdateH5Seurat(
                file = opt$input,
                object = data_ob,
                reduction = unique(c(reduct1s, opt$reduct2))
            )
        } else {
            # if (!is.null(stor_slot)) {
            #     data_ob$SCT@scale.data <- stor_slot
            #     print("add scale.data")
            # }
            OESingleCell::SaveX(
              data_ob,
              output = opt$output,
              update = FALSE,
              outformat = opt$outformat,
              prefix = opt$prefix
            )
        }
        clustersinf=data.frame(Barcode=data_ob$rawbc,sampleid=data_ob$sampleid,clusters=data_ob$clusters,group=data_ob$group)
        write.table(clustersinf, file = file.path(output_dir, 'clusters_infor.csv'), sep = ",", quote = F, row.names = F)
        # reduct=data_ob@reductions$umap@cell.embeddings
        # rownames(reduct)=data_ob$rawbc
        # write.table(reduct,file=file.path(output_dir,paste0(opt$reduct2,
        #                                '_Dimension_Reduction_coordination.csv')),sep=",",quote=F,row.names=T)

        ## 4.4 save session informations
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}

#示例 sctool -i result/count_qc/filtered_seurat.rds -f rds -o result/cluster_seurat -d rds \
# --assay SCT --image TRUE --update FALSE --prefix singlecell_object.clustering_resolution0.4 bclust --reduct1 pca,harmony --reduct2 umap --batchid batchid \
# --component  30 --clusteringuse snn --resolution 0.4 --rerun  T
# --predicate def_celltype==\"Oligodendrocyte\"
