## reference:  https://satijalab.org/seurat/reference/findallmarkers
## usage help:
docstring <- "reference:  https://satijalab.org/seurat/reference/findallmarkers \\n\\n\\
example1: sctool -i st.rds -f rds -o results --assay SCT findallmarkers -e MAST  -n clusters --anno /data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_annotation.xls \\n\\n\\
example2: sctool -i sc.rds -f rds -o results --assay RNA findallmarkers -e MAST  -n clusters --anno /data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_annotation.xls \\n\\n\\
"
sub_markers  <- subparsers$add_parser("findallmarkers",
                                      description = docstring,
                                      formatter_class = 'argparse.RawTextHelpFormatter',
                                      #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                      argument_default = "True",
                                      help = "find markers between all groups.")
sub_markers$add_argument(
    "--group_by",
    "-n",
    type = "character",
    default = NULL,
    help = "the name of groupping column from clustering result metadata used to find markers. The example can be tsne.2D.res.1."
)
sub_markers$add_argument(
    "--test",
    "-e",
    type = "character",
    default = "wilcox",
    help = "test methods used to cell markers. Options are: wilcox,presto,venice, t, bimod,poisson, negbinom, MAST, DESeq2, DESeq2LRT,limma, edgeR.[default: %(default)s] "
)
sub_markers$add_argument(
    "--anno",
    required = TRUE,
    type = "character",
    help = "[OPTIONAL]Annotation file submitted when species is Other.[default: %(default)s]"
)
sub_markers$add_argument(
    "--strict",
    "-s",
    type = "character",
    default = "FALSE",
    help = "whether to use strict mode, using all the filtering options, to find markers. Notice that, this may result in no markers remain for some clusters.[default: %(default)s] "
)
sub_markers$add_argument(
    "--avg_log2FC",
    "-k",
    type = "double",
    default = 1,
    help = "The average log2FC of the gene UMI count in its cluster against the all other clusters.[default: %(default)s]"
)
sub_markers$add_argument(
    "--pvalue",
    "-p",
    type = "double",
    default = 0.05,
    help = "the P-value of the gene differential expression.[default: %(default)s]"
)
sub_markers$add_argument(
    "--FDR",
    "-q",
    type = "double",
    default = 0.1,
    help = "the FDR of the gene differential expression.[default: %(default)s] "
)
sub_markers$add_argument(
    "--min_pct1",
    "-t",
    type = "double",
    default = 0.5,
    help = "the minimium ratio of cells expressing one specific gene in a cluster.[default: %(default)s]"
)
sub_markers$add_argument(
    "--max_pct2",
    "-T",
    type = "double",
    default = 0.5,
    help = "the maximiium ratio of cells expressing one specific gene in all other clusters.[default: %(default)s]"
)
sub_markers$add_argument(
    "--pct_fold",
    "-c",
    type = "double",
    default = 2,
    help = "the minimiu fold of pct1 for gene in a specific cluster against pct2 for all other cluster.[default: %(default)s]"
)
sub_markers$add_argument(
    "--topn_marker",
    "-N",
    type = "integer",
    default = 10,
    help = "the maximium number of ranked marker genes on the top for each cluster.[default: %(default)s]"
)
sub_markers$add_argument(
    "--sub_obj",
    default = NULL,
    help = "if the inputted rds had been subclustered or subsetted, it should be given by this parameter.[default: %(default)s]"
)
sub_markers$add_argument(
    "--recorrect_umi",
    default = TRUE,
    help = "whether to recorrect umi when run FindMarkers. If the object had run PrepSCTFindMarkers and subset, it should be 'FALSE'.[default: %(default)s]"
)

# ============== Subcmd: findallmarkers Find All markers for each cell group ===============
args <- commandArgs(TRUE)
if ("findallmarkers" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "findallmarkers") {
        # ==============================================================================================================
        futile.logger::flog.info("step1:load query object data and get annotation files")
        suppressMessages(data_ob <- OESingleCell::ReadX(
            input = opt$input,
            informat = opt$informat,
            assays = assays,
            data.use = dataslots, # "data"
            reductions = FALSE, # only the used reduction results needed
            graphs = FALSE, # no graph object needed here
            verbose = FALSE,
            image = opt$image
        ))
        anno <- read.delim(opt$anno, stringsAsFactors = F, sep = "\t", quote = "")
        # ==============================================================================================================
        futile.logger::flog.info("step2:define the group for marker finding ")
        # in default, the FindAllMarkers() function will use the default identity class
        # In order to select the different clustering result to find markers, change this
        # class value to one of the column in the cell annotation metadata table
        suppressWarnings({
            if (is.null(opt$group_by)) {
                data_ob[["clusters"]] <- Seurat::Idents(data_ob)
            } else {
                # data_ob = Seurat::SetIdent( data_ob, value = opt$cluster_name)
                Seurat::Idents(data_ob) <- opt$group_by
            }
        })
        # ==============================================================================================================
        futile.logger::flog.info("step3:Running findallmarkers ")
        # find the differential expressed genes for each clusters
        # FindAllMarkers() is primarily used to find markers, but here it was also used
        # to find differentially expressed genes. The default logfc threshold was set to 0.25
        # ,here we set to 0 as no prefiltering and then manually filter the genes to find markers
        # note that if test.use is "negbinom", "poisson", or "DESeq2", slot will be set to "counts" automatically
        if(assays=="SCT"){
            if((length(data_ob@assays$SCT@SCTModel.list) > 1) & is.null(data_ob@misc$PreSCT)){
            ## https://rdrr.io/cran/Seurat/man/PrepSCTFindMarkers.html
            futile.logger::flog.info("Prepare object to run differential expression on SCT assay ")
            data_ob<-Seurat::PrepSCTFindMarkers(data_ob, assay = "SCT", verbose = TRUE)
            data_ob@misc$PreSCT <- TRUE
            }
        }
        saveRDS(data_ob, glue::glue("{output_dir}/seurat_object_markers.rds"))
        if (!is.null(opt$predicate)) {
            futile.logger::flog.info(glue::glue("step1.5: get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
            #factors <- stringr::str_split(opt$predicate, ",")[[1]]
            df <- slot(data_ob, "meta.data")
            desired_cells <- subset(df, eval(parse(text = opt$predicate)))
            data_ob <- subset(data_ob, cells = rownames(desired_cells))
            data_ob@meta.data$sampleid <- factor(data_ob@meta.data$sampleid, levels = unique(data_ob@meta.data$sampleid))
        }
        if (!is.null(Seurat::Images(data_ob))) {
            unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique())]
            if (length(unuse_images) > 0) {
                data_ob@images[unuse_images] <- NULL
            }
        }
        if(!is.null(opt$sub_obj)){
            futile.logger::flog.info("Get subseted object after PreSCT")
            suppressMessages(data_sub <- OESingleCell::ReadX(
                    input = opt$sub_obj,
                    informat = opt$informat,
                    assays = assays,
                    data.use = dataslots, # "data"
                    verbose = FALSE
            ))
            data_ob <- subset(data_ob, cells = rownames(data_sub@meta.data))
            if (!is.null(Seurat::Images(data_ob))) {
                unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique())]
                if (length(unuse_images) > 0) {
                    data_ob@images[unuse_images] <- NULL
                }
            }
            metadata <- data_ob@meta.data %>% dplyr::mutate(new_groupby = data_sub[[opt$group_by]])
            data_ob <- SeuratObject::AddMetaData(data_ob, metadata$new_groupby, "new_groupby")
            SeuratObject::Idents(data_ob) <- data_ob[["new_groupby"]]
        }
        futile.logger::flog.info("Run FindAllMarkers")
        recorrect_umi <- as.logical(opt$recorrect_umi)
        global_DEGs <- Seurat::FindAllMarkers(
            object = data_ob,
            only.pos = T,
            test.use = opt$test,
            logfc.threshold = 0,
            min.pct = 0.25,
            recorrect_umi = recorrect_umi
        )
        global_DEGs <- global_DEGs  %>%
            dplyr::mutate(gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>%
            dplyr::select(gene, dplyr::everything()) %>% dplyr::rename(`p-value` = p_val) %>% dplyr::rename(`q-value` = p_val_adj)
        ## 进行基因注释
        diff_columns <- colnames(global_DEGs)
        global_DEGs <- OESingleCell::ann_marker(global_DEGs, anno)
        readr::write_tsv(
            global_DEGs,
            file = glue::glue("{output_dir}/all_markers_for_each_cluster_anno.xls"),
        )
        readr::write_tsv(
            global_DEGs %>% dplyr::select(all_of(diff_columns)),
            file = glue::glue("{output_dir}/all_markers_for_each_cluster.xls"),
        )
        # ==============================================================================================================
        futile.logger::flog.info("step4:slect top marker use setting paramters  ")
        # find the significant differential expressed genes for each clusters against all other clustershttp://127.0.0.1:14418/help/library/Venice/html/CreateSignacObject.html
        # feature plot of potential marker gene for each cluster
        # for each gene to be a potential marker,in add to be a significant expressed gene, the gene should account for
        # large proportion in the interested cluster but as small as possiable in the other clusters, that means the pct.1 should
        # be bigger than the pct.2.
        if (as.logical(opt$strict) == T) {
            if (!is.null(opt$pvalue)) {
                topn_markers <- global_DEGs %>%
                    dplyr::group_by(cluster) %>%
                    dplyr::filter(avg_log2FC >= opt$avg_log2FC &
                        `p-value` < opt$pvalue &
                        pct.1 > opt$min_pct1 &
                        pct.2 < opt$max_pct2 &
                        gene_diff > opt$pct_fold) %>%
                    dplyr::arrange(`p-value`, dplyr::desc(avg_log2FC), dplyr::desc(gene_diff)) %>%
                    dplyr::filter(gene_diff > opt$pct_fold) %>%
                    dplyr::top_n(opt$topn_marker, gene_diff)
            } else {
                topn_markers <- global_DEGs %>%
                    dplyr::group_by(cluster) %>%
                    dplyr::filter(avg_log2FC >= opt$avg_log2FC &
                        `q-value` < opt$FDR &
                        pct.1 > opt$min_pct1 &
                        pct.2 < opt$max_pct2 &
                        gene_diff > opt$pct_fold) %>%
                    dplyr::arrange(`p-value`, dplyr::desc(avg_log2FC), dplyr::desc(gene_diff)) %>%
                    dplyr::filter(gene_diff > opt$pct_fold) %>%
                    dplyr::top_n(opt$topn_marker, gene_diff)
            }
        } else {
            topn_markers <- global_DEGs %>% dplyr::group_by(cluster) %>%
                # filter(p_val < opt$pvalue ) %>%
                dplyr::arrange(`p-value`, dplyr::desc(avg_log2FC), dplyr::desc(gene_diff)) %>%
                # filter(gene_diff > opt$pct_fold)  %>%
                dplyr::top_n(opt$topn_marker, gene_diff)
        }
        ##  save topn markers
        readr::write_tsv(
            topn_markers,
            file = glue::glue("{output_dir}/top{opt$topn_marker}_markers_for_each_cluster_anno.xls"),
        )
        readr::write_tsv(
            topn_markers %>% dplyr::select(all_of(diff_columns)),
            file = glue::glue("{output_dir}/top{opt$topn_marker}_markers_for_each_cluster.xls"),
        )
        ## =============================================================================================================
        ## save session informations
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}