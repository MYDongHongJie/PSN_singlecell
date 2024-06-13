# create the parser for subcommand celltyping
# the dataslot here is used to specifiy which data slot to use when celltyping
docstring<- "[Example1(with  custom reference database)]:  \\n\\
sctool -i  seurat.rds \\\\\\n\\
       -o results/singleR_results_allen_mouse_cortex_hippocampus_10x \\\\\\n\\
       -f rds \\\\\\n\\
       --assay RNA \\\\\\n\\
       --dataslot counts \\\\\\n\\
       -d rds \\\\\\n\\
       --prefix seurat  \\\\\\n\\
       --update FALSE  \\\\\\n\\
       -j 6 \\\\\\n\\
      celltyping  --refcustom   allen_mouse_cortex_hippocampus_10x.rds  \\\\\\n\\
                  --refinformat rds  \\\\\\n\\
                  --refassay RNA  \\\\\\n\\
                  --refdataslot counts  \\\\\\n\\
                  --annolevel  subclass \\\\\\n\\
                  --refmarker NULL \\\\\\n\\
                  --demethod  wilcox  \\\\\\n\\
                  --usecluster clusters  \\\\\\n\\
                  -v 0.8   \\\\\\n\\
                  -n 25 \\\\\\n\\
                  --palette   blindless \\\\\\n\\
                  --reduct    umap  \\\\\\n\\
                  --species  mouse_brain  \\\\\\n\\
                  --doplot  TRUE \\\\\\n\\
sctool -i multitome.h5seurat \\\\\\n\\
              -f h5seurat \\\\\\n\\
              -o result/singleR_results_allen_mouse_cortex_hippocampus_10x \\\\\\n\\
              -d h5seurat \\\\\\n\\
              --assay SCT \\\\\\n\\
              --dataslot data,counts \\\\\\n\\
              celltyping    --refbuiltin hpca \\\\\\n\\
                            --annolevel main \\\\\\n\\
                            --refmarker NULL \\\\\\n\\
                            --demethod wilcox \\\\\\n\\
                            --clusterby SCT.umap.res.0.4 \\\\\\n\\
                            -v 0.8 \\\\\\n\\
                            -n 25 \\\\\\n\\
                            --palette blindless \\\\\\n\\
                            --reduct SCT_umap \\\\\\n\\
                            --species human \\\\\\n\\
                            --doplot TRUE"

sub_celltype <- subparsers$add_parser("celltyping",
                                       description = docstring,
                                       formatter_class= 'argparse.RawTextHelpFormatter' ,
                                       #formatter_class= "lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                       argument_default = "True" ,
                                       help = "find the celltype for each cell in data object.")
sub_celltype$add_argument("--refbuiltin", type = "character",default=NULL,
             choices=c("hpca","blueprint_encode","schcl","immgen","mouse.rnaseq","scmca"),
             help = paste0("One or more SummarizedExperiment/SingleCellExperiment data object containing log-transformed matrix used as non-built-in reference.",
                           "human:<hpca,blueprint_encode,schcl> mouse:<immgen,mouse.rnaseq,scmca>"))
sub_celltype$add_argument("--refcustom", type = "character",default=NULL,
             help = "one or more reference data object name in seurat object.")
sub_celltype$add_argument("--refinformat", type = "character", default = NULL,
            help = "The format of data object, the possible choices can be:h5seurat,(seurat)rds,(sce)rds, loom. ")
sub_celltype$add_argument("--refassay", type = "character", default = NULL,
              help = "the main assay in data object to use. When it comes to multimodal assay, this is the assay used to initialize the object, all the other assays will merged into it. ")
sub_celltype$add_argument("--refdataslot", type = "character", default = "data",
          help = "the data slot in the specified assay used for this run, one of 'counts','data' must be provided.")

sub_celltype$add_argument("--annolevel", type = "character", default=NULL,
             help = paste0("cell identities annotation column in reference data for each cell in test data. If set, ",
                           "annotation is performed on the cluster level, otherwise it defaults to per-cell annotation.",
                           "[default: %(default)s]"))
sub_celltype$add_argument("--refmarker", type = "character", default = NULL,
             help = paste0("ref markers for refenence database."))

sub_celltype$add_argument("--demethod", type = "character", default = "classic",
             choices=c("classic","wilcox","t"),
             help = paste0("String specifying how DE genes should be detected between pairs of labels in ref data. ",
                           "Choices can be:classic,wilcox and t.[default: %(default)s]."))

sub_celltype$add_argument("--usecluster", type = "character", default = "FALSE",
             help = "whether to predict identitiy for each cluster of test data rather than single cell.[default: %(default)s]")

sub_celltype$add_argument("--clusterby", type = "character", default = "FALSE",
             help = "select cluster to analysis and plot heatmap.[default: %(default)s]")
sub_celltype$add_argument("-v", "--pointsize", type = "double", default = 1,
             help = "the size for each point in scatter plot.[default: %(default)s]")
sub_celltype$add_argument("-n", "--topn", type = "integer", default = 25,
             help = "the top number of reference cell type score for each test cell to visulize on heatmap.[default: %(default)s]")
sub_celltype$add_argument("--colorschema", type = "character", default = "spectral",
             help = "the default color schema mapping for value in each heatmap cell. The format:low, middle, high.[default: %(default)s]")
sub_celltype$add_argument("--palette",type = "character",default= "blindless",
             help=paste0("the discrete color schema for each cell types, blindless as default. Choices:blindless:69,cold:32",
                          "glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,tableau20:20, paired:12, purplegray12:12[default: %(default)s]"))
sub_celltype$add_argument("--reduct", type = "character", default = "umap",
             help = "[OPTIONAL]the reduction results used to embedding the celltypeing results.[default: %(default)s]")
sub_celltype$add_argument("--species", type = "character", default = "human",
             help = "[OPTIONAL]the organism both test and reference data derived.[default: %(default)s]")
sub_celltype$add_argument("--doplot", type = "character", default = "TRUE",
             help = "visualize the prediction results using the heatmap.[default: %(default)s]")
# ================= Subcmd: Celltyping invoked=========
args <- commandArgs(TRUE)
if ( "celltyping"  %in% args ){
     opt<-intial_setting()
     if(opt$sub_name == "celltyping" ){
       ##===============================================================================================================
       species <- opt$species
       #========================================1.test.sce preparing =====================================================
       ## 1.1. read the specified assay and data slot in data object into memory
        data_ob <- OESingleCell::ReadX(input = opt$input,
                                      informat = opt$informat,
                                      assays = assays,
                                      data.use = dataslots, ##"counts"
                                      reductions = opt$reduct, # only the used reduction results needed
                                      graphs = FALSE, # no graph object needed here
                                      images = FALSE,
                                      verbose = FALSE)
        Seurat::Idents(data_ob) <- data_ob[[opt$clusterby]]
        data_ob[["clusters"]] <- data_ob[[opt$clusterby]]
        ##1.2 subset test data object if specify
        if ( !is.null(opt$predicate) ){
            df <- colData(data_ob)
            desired_cells <- subset(df, eval(parse(text=opt$predicate)))
            data_ob <- data_ob[, rownames(desired_cells)]
        }
        ##1.3 convert to sceformart and logNormcounts
        test.sce <- Seurat::as.SingleCellExperiment(data_ob)
        test.sce <- scuttle::logNormCounts(test.sce)
       #========================================2.ref.sce preparing ======================================================
        # 2.Reference datasets are log2 transformed rather than natural log transformed (Seurat).
        if ( is.null(opt$refbuiltin) & is.null(opt$refcustom) ) stop("NO reference AVAILABE!")
        refd.list <- list()
        ##referce database
        if (!is.null(opt$refcustom)) {
          cref <- unlist(strsplit(opt$refcustom, ",") )
          for( refx in cref ){
            refx <- normalizePath(refx)
            refid <- gsub("\\..[^\\.]*$", "", make.names(basename(refx)))
            refd.list[[refid]] <- OESingleCell::ReadX(input = refx,
                                                      informat = opt$refinformat,
                                                      assays = opt$refassay,  ##"RNA"
                                                      data.use = opt$refdataslot, ##"counts"
                                                      verbose = FALSE)
            refd.list[[refid]]<-Seurat::as.SingleCellExperiment(refd.list[[refid]])
            refd.list[[refid]]<-scuttle::logNormCounts(refd.list[[refid]])
          }
        }

        ## refrence  database from celldex : https://github.com/LTLA/celldex/tree/master/R
         if ( !is.null(opt$refbuiltin) ){
          bref <- unlist(strsplit(opt$refbuiltin, ",") )
          for( refx in bref ){
              refrds = paste0(refx,".rds")
              refd.list[[refx]] = readRDS(file.path("/data/database/celltype_refdata/logNorm_rds", refrds))
          }
        }

        # prepare the reference data set
        if ( length(refd.list) > 1 ){
          print("using combined result from multiple datasets")
          ref.labels <- lapply(refd.list, function(x) factor(x[[opt$annolevel]]) )
        }else{
          ref.labels <- refd.list[[1]][[opt$annolevel]]
        }

        #  if ( length(refd.list) == 1 ){
        #    refd.list = refd.list[[1]]
        #  }
        refdata <- paste0(names(refd.list), collapse = "_")

        #======================================3.Run prediction now====================================================
        ## running singleR
        if(length(refd.list)>1){ref.sce<-refd.list}else{ref.sce<-refd.list[[1]]}
        pred <- SingleR::SingleR(test= test.sce,
                                 ref =  ref.sce,
                                 labels = ref.labels,
                                 genes= ifelse(!is.null(opt$refmarker),"de", opt$refmarkers),
                                 de.method=opt$demethod,
                                 assay.type.test = "logcounts",
                                 assay.type.ref = "logcounts",
                                 aggr.ref=ifelse( !is.null(opt$refcustom),"TRUE","FALSE"),
                                 BPPARAM = BiocParallel::MulticoreParam(workers = future::nbrOfWorkers()))
        ## add predict celltype results to seurat
        raw.metaname <- glue::glue("raw_{opt$annolevel}_celltype")
        data_ob <- Seurat::AddMetaData(data_ob, metadata = pred[["labels"]], col.name = raw.metaname)
        if ( length(refd.list) > 1 ){
          celltyping_cols <- data.frame(labels = pred$labels)
          for ( dx in names(refd.list) ){
            celltyping_cols[[glue::glue("{dx}_labels")]] <- pred[["orig.results"]][[dx]][["labels"]]
          }
          data_ob <- Seurat::AddMetaData(data_ob, metadata = celltyping_cols)
        }

        celltyping_stat <- Seurat::FetchData(data_ob, vars = c("clusters", raw.metaname)) %>%
          dplyr::group_by(clusters, !!rlang::sym(raw.metaname)) %>%
          dplyr::summarize(cell_num = dplyr::n()) %>%
          dplyr::ungroup()
        write.table(celltyping_stat,
                    quote = F,
                    glue::glue("{output_dir}/{opt$species}ref_{refdata}_{opt$annolevel}_celltyping_statistics.xls"),
                    sep = "\t",
                    row.names = F)

        #========================================4.DotPlot =============================================================
        pdf(NULL)
        if ( as.logical(opt$doplot) ){
          ## heatmap
          message("4.Plot celltyping heatmap")
          meta.data <- Seurat::FetchData(data_ob, vars = "clusters") %>%
            tibble::rownames_to_column(var = "cells")
          if ( ncol(data_ob) > 40000){
            message("Screening of 50,000 cells in proportion...")
            meta.data <- meta.data %>%
                         group_by(clusters) %>%
                         sample_frac(40000/dim(data_ob)[2])
          }
          ## set color schema
          colstrings <- unlist(strsplit(opt$colorschema, ","))
          if (length(colstrings)==1 ){
           # heatmap_color_schema <- (grDevices::colorRampPalette(c("#D1147E", "white", "#00A44B")))(100)
            heatmap_color_schema <- OESingleCell::continuous_palette[[opt$colorschema]]
          }else{
            heatmap_color_schema <- (grDevices::colorRampPalette(colstrings))(100)
          }
          ## set picture par
          ncells <- ifelse(dim(meta.data)[1] > 70000, 70000, dim(meta.data)[1]  )
          plot_width <- (1e5 - ncells)/1e4 * ceiling(ncells/30000)
          annotation_colors = list(Clusters=OESingleCell::SelectColors(sort(data_ob[["clusters"]][[1]]), palette = opt$palette))
          if ( length(refd.list) == 1 ){
            ggheat <- SingleR::plotScoreHeatmap(pred,
                                                cells.use = meta.data$cells,
                                                clusters = meta.data$clusters,
                                                max.labels = opt$topn,
                                                show.labels = F,
                                                colors = heatmap_color_schema,
                                                order.by = "clusters",
                                                annotation_colors=annotation_colors)
            width4legendlabs <- max(nchar(names(head(sort(table(as.factor(pred$labels)), decreasing = T), opt$topn))))
            OESingleCell::save_ggplots(glue::glue("{output_dir}/{opt$species}ref_{refdata}_{opt$annolevel}_celltyping_heatmap"),
                                       plot = ggheat,
                                       width = plot_width + width4legendlabs/15 )
          }else{
            refid <- c("combine", names(refd.list))
            for ( inx in 0:length(refd.list) ){
              ggheat <- SingleR::plotScoreHeatmap(pred,
                                                  scores.use = inx,
                                                  cells.use = meta.data$cells,
                                                  clusters = meta.data$clusters,
                                                  max.labels = opt$topn,
                                                  show.labels = F,
                                                  colors = heatmap_color_schema,
                                                  order.by = "clusters")
              width4legendlabs <- max(nchar(names(head(sort(table(as.factor(pred$labels)), decreasing = T), opt$topn))))
              OESingleCell::save_ggplots(filename=glue::glue("{output_dir}/{opt$species}ref_{refid[inx+1]}_{opt$annolevel}_celltyping_heatmap"),
                                         plot = ggheat,
                                         width = plot_width + width4legendlabs/15 )
            }
          }
          ##========================================5. plot raw celltype ===============================================
          pdf(NULL)
          message("5.Plot celltyping embeddings")
          ## setting picture par
          custom_cols <- OESingleCell::SelectColors(data_ob[[raw.metaname]][[1]], palette = opt$palette)
          nlevel <- length(unique(data_ob[[raw.metaname]][[1]]))
          ncols <- ifelse(nlevel > 30, as.integer(nlevel/30)+1, 1)
          ## plot
          ggdim <- Seurat::DimPlot(data_ob,
                                   reduction = opt$reduct,
                                   pt.size = opt$pointsize,
                                   group.by = raw.metaname ) +
                   ggplot2::theme( plot.title = element_text(hjust = 0.5),
                                   legend.key.size = grid::unit(0.9,"lines"),
                                   legend.text = element_text(size = 10)) +
                   ggplot2::scale_colour_manual( values = custom_cols ) +
                   ggplot2::guides(colour = guide_legend(ncol = ncols,
                                                         override.aes = list(size=2)))
          OESingleCell::save_ggplots(filename = glue::glue("{output_dir}/{opt$species}ref_{refdata}_{opt$annolevel}_celltyping_plot"),
                                     plot = ggdim,
                                     width = (plot_width+ (width4legendlabs/13)* ncols),
                                     dpi = 1000)
        }

        ##===================================== 6.TOPN result============================================================
        pdf(NULL)
        # data_ob = Seurat::SetIdent( data_ob, value = "clusters")
        top_celltype <- celltyping_stat %>%
                        dplyr::group_by(clusters) %>%
                        dplyr::top_n(1, cell_num) %>%
                        dplyr::ungroup()
        data_ob[[glue::glue("{refdata}.{opt$annolevel}.celltype")]] <- plyr::mapvalues(x = data_ob[["clusters"]][[1]],
                                                                            from = as.vector(top_celltype[["clusters"]]),
                                                                            to = as.vector(top_celltype[[raw.metaname]]))
        data_ob[["celltype"]] <- data_ob[[paste0(refdata,".",opt$annolevel,".celltype")]] # the final celltyping results column

        full_celltyping_df <- OESingleCell::colData(data_ob) %>%
                              tibble::rownames_to_column(var = "barcode_inuse") %>%
                              dplyr::rename( "Barcode" = "orig.ident" ) %>%
                              dplyr::select( Barcode, everything())
        write.table(full_celltyping_df,
                    quote = F,
                    glue::glue("{output_dir}/{opt$species}ref_{refdata}_{opt$annolevel}_celltyping_results.xls"),
                    sep = "\t",
                    row.names = F)

        simplified_celltyping_df <- OESingleCell::colData(data_ob) %>%
                                    dplyr::rename( "Barcode" = "orig.ident") %>%
                                    dplyr::select( Barcode, sampleid, celltype, clusters,group)

        write.table(simplified_celltyping_df,
                    quote = F,
                    sep =",",
                    row.names = F,
                    glue::glue("{output_dir}/{opt$species}ref_{refdata}_{opt$annolevel}_simplified_celltyping_results.csv"))
        ##==============================================================================================================
        ## obtain the top cluster number for each celltype (for the purpose of coloring)
        # top_celltype[["colors"]] = OESingleCell::SelectColors(top_celltype[["clusters"]], palette = opt$palette)
        colorx <- top_celltype %>%
                    dplyr::mutate(colors = OESingleCell::SelectColors(clusters, palette = opt$palette)) %>%
                    dplyr::select(-c(clusters, cell_num)) %>%
                    dplyr::group_by(.data[[raw.metaname]]) %>%
                    dplyr::slice_head(n = 1) %>%
                    dplyr::ungroup() %>% tibble::deframe()
        nlevel <- length(colorx)
        ncols <- ifelse(nlevel > 30, as.integer(nlevel/30)+1, 1)
        ggdim2 <- Seurat::DimPlot(data_ob,
                                  reduction = opt$reduct,
                                  pt.size = opt$pointsize,
                                  group.by = "celltype" ) +
                  theme( plot.title = element_text(hjust = 0.5)) +
                  scale_colour_manual( values = colorx )
        OESingleCell::save_ggplots(filename=glue::glue("{output_dir}/{opt$species}ref_{refdata}_top.{opt$annolevel}_celltyping_plot"),
                                   plot = ggdim2,
                                   width = (plot_width+ (width4legendlabs/13)* ncols),
                                   dpi = 1000,)
        ##=====================================8. output results =======================================================
        if ( as.logical(opt$update) ){
          SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob, verbose = FALSE )
        }else{
          OESingleCell::SaveX(data_ob,
                              output = opt$output,
                              update = FALSE,
                              outformat = opt$outformat,
                              prefix = opt$prefix)
        }
        ##session information output
        write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
        quit()
     }
}