summary_stat_plot <- function(
  object,
  show_features,
  summx_features,
  output_dir,
  prefix,
  pointsize,
  crop
) {
    pdf(NULL)
    nsamples <- unique(Seurat::FetchData(object, vars = "sampleid")[["sampleid"]])
    for(features in show_features ){
        ## vlnplot======================================================================================================
        plot1 <- Seurat::VlnPlot(
            object,
            features = features,
            group.by = "sampleid",
            # alpha = 0.5,
            pt.size = pointsize
        )+Seurat::NoLegend() + ggtitle("") + theme(axis.title.x = element_blank())
        ## spatialfeatureplot ==========================================================================================
        if (!is.null(Seurat::Images(object))) {
            plot2 <- OESingleCell::SpatialPlot(
              object,
              assay = "Spatial",
              features = features,
              ncol = length(nsamples),
              HE = F,
              cols = "spectral",
              alpha = 1,
              pt.size.factor = 1.2,
              crop = crop
            )
            plot<-patchwork::wrap_plots(plot1,
                                         plot2,
                                         ncol = 2,
                                         nrow=1,
                                         widths = c(length(nsamples)*1.2,length(nsamples)*4.5+stringr::str_length(features)*0.1),
                                         heights=c(4.5+max(stringr::str_length(nsamples))*sin(45)*0.1, 4.5))+
                patchwork::plot_annotation(title = glue::glue('{features}'),
                                           theme = theme(legend.direction ="horizontal",
                                                         legend.position = "top",
                                                         legend.key.width=unit(1,'cm'),
                                                         plot.title = element_text(size = 20,hjust = 0.5,face="bold")))
        }

        OESingleCell::save_ggplots(
            filename = glue::glue("{output_dir}/{prefix}_featureplot_for_{features}"),
            plot = plot,
            width = ifelse(is.null(Seurat::Images(object)),
                           length(nsamples)*1,
                           length(nsamples)*1.2+length(nsamples)*4.5+stringr::str_length(features)*0.1) ,
            height =4.5+max(stringr::str_length(nsamples))*sin(45)*0.1,
            limitsize = FALSE
        )
    }

    ## statistic========================================================================================================
    stat <- Seurat::FetchData(object, vars = c(summx_features, "sampleid")) %>%
        dplyr::group_by(sampleid) %>%
        dplyr::summarise(
            dplyr::across(
                where(is.numeric),
                list(
                    ~ round(mean(.x, na.rm = TRUE), 3),
                    ~ median(.x, na.rm = TRUE)
                )
            ),
            total_cells = dplyr::n()
        ) %>%
        tibble::column_to_rownames(var = "sampleid")
    indicator <- expand.grid(c("mean", "median"), summx_features, prefix)
    if (opt$assay != "RNA"){
        stat_colnames <- c(
          glue::glue("{indicator[[1]]}_{indicator[[2]]}_{indicator[[3]]}"),
          glue::glue("Total_Spots_{prefix}")
        )
    }else {
        stat_colnames <- c(
          glue::glue("{indicator[[1]]}_{indicator[[2]]}_{indicator[[3]]}"),
          glue::glue("Total_cells_{prefix}")
        )
    }
    print(stat_colnames)
    colnames(stat) <- stat_colnames
    return(stat)
}

# on stop of running data quanlity control of
docstring <- "[Example1(For 10x multiome data]:  \\n\\
sctool --input result/2.Count_QC/multitome.h5seurat \\\\\\n\\
              --informat h5seurat \\\\\\n\\
              --output result/2.Count_QC \\\\\\n\\
              --outformat h5seurat \\\\\\n\\
              --assay RNA \\\\\\n\\
              --subassay ATAC \\\\\\n\\
              --dataslot counts \\\\\\n\\
              --prefix  multitome \\\\\\n\\
              --ncores 3 \\\\\\n\\
              --update TRUE \\\\\\n\\
              qc   --vars2regress \"nCount_RNA|nCount_ATAC\" \\\\\\n\\
                   --filters nFeature_RNA,nCount_RNA,percent.mito,atac_fragments,TSS.enrichment,nucleosome_signal \\\\\\n\\
                   --lower NULL,NULL,NULL,1000,2,0 \\\\\\n\\
                   --upper NULL,NULL,NULL,100000,100,2 \\\\\\n\\
                   --cut.1 median \\\\\\n\\
                   --cut.2 mad \\\\\\n\\
                   --nfold 2 \\\\\\n\\
                   --features2filter NULL \\\\\\n\\
                   --mincell4gene  1 \\\\\\n\\
                   --rmdoublets TRUE \\\\\\n\\
                   --method doubletfinder \\\\\\n\\
                   --ident NULL \\\\\\n\\
                   --normmeth \"LogNormalize|logtf-idf\" \\\\\\n\\
                   --nvfeatures 2000 \\\\\\n\\
                   --pointsize 0.1"

sub_qc <- subparsers$add_parser(
    "qc",
    description = docstring,
    formatter_class = "argparse.RawTextHelpFormatter",
    # formatter_class= "lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
    argument_default = "True",
    help = "one stop of data control."
)

sub_qc$add_argument(
    "-r",
    "--vars2regress",
    type = "character",
    default = "nCount_RNA,percent.mito",
    help = paste0(
        "comma separated list of sources of variation to regress out from the data. ",
        "The other options can be percent.ribo, CC.Difference etc.[default: %(default)s]"
    )
)
sub_qc$add_argument(
    "-c",
    "--filters",
    type = "character",
    default = "nFeature_RNA,nCount_RNA,percent.mito",
    help = paste0(
        "the variables used to remove the outliers of cells.",
        "The other options can be percent.ribo etc.[default: %(default)s]"
    )
)
sub_qc$add_argument(
    "-l",
    "--lower",
    type = "character",
    default = "200,2000,0",
    help = "the lower bounds list for the parameters specified by --filters.[default: %(default)s]"
)
sub_qc$add_argument(
    "-L",
    "--upper",
    type = "character",
    default = NULL,
    help = paste0(
        "the upper bounds list for the parameters specified by --filters.If not specified,",
        "it will be determined by the parameter --cut.1, --cut.2 and -n automatically.[default: NULL]"
    )
)

sub_qc$add_argument(
    "--cut.1",
    type = "character",
    default = "median",
    help = paste0(
        "the method to determine the center of an array, ",
        "here refering to the median or mean of the QC metric.[default: %(default)s]"
    )
)
sub_qc$add_argument(
    "--cut.2",
    type = "character",
    default = "mad",
    help = "the method to find the difference unit to the data center, here refering to median, sd or mad.[default: %(default)s]"
)
sub_qc$add_argument(
    "-n",
    "--nfold",
    type = "double",
    default = 2.0,
    help = paste0(
        "the fold of intervals(sd for mean, mad for median, median for median) for parameters ",
        "specified by --filters determine the outlier threshold by +/- n*cut.2.[default: %(default)s]"
    )
)


sub_qc$add_argument(
    "-u",
    "--rlm.xy",
    type = "character",
    default = NULL,
    help = paste0(
        "[OPTIONAL]using the robust linear model of variables, for example nCount_RNA:nFeature_RNA,",
        "to remove the outliers of cells."
    )
)
sub_qc$add_argument(
    "-z",
    "--features2filter",
    type = "character",
    default = NULL,
    help = "[OPTIONAL]the file of specified feature list with 'gene' as header to remove from the features-cell matrix.[default: NULL]"
)
sub_qc$add_argument(
    "--nvfeatures",
    type = "double",
    default = 3000,
    help = "[OPTIONAL]the number of top variable features to keep for RNA-seq or the minimum percentile(0~1), for example 0.95, for ATAC-seq ."
)
sub_qc$add_argument(
    "-x",
    "--mincell4gene",
    type = "double",
    default = 0.01,
    help = paste0(
        "the minimium cell number one gene detected.If the value is less than 1, it is a proportion,",
        "otherwise a integerial cell number.[default: %(default)s]"
    )
)


sub_qc$add_argument(
    "--rmdoublets",
    type = "character",
    default = "FALSE",
    help = "remove the doublet cells.[default: %(default)s]"
)
sub_qc$add_argument(
    "--method",
    type = "character",
    default = NULL,
    # choices=c("scrublet", "doubletfinder"),
    help = "the method used to detect doublets, choices: scrublet, doubletfinder.[default: NULL]"
)
sub_qc$add_argument(
    "--ident",
    type = "character",
    default = "clusters",
    help = "[OPTIONAL]ONLY AVAILABLE for doubletfinder,the name of cell group annotation.[default: NULL]"
)

sub_qc$add_argument(
    "-m",
    "--normmeth",
    type = "character",
    default = "LogNormalize",
    # choices=c("sctransform", "LogNormalize", "scran", "CR", "CLR", "logtf-idf" ,"tf-logidf","logtf-logidf","idf"),
    help = paste0(
        "the method to normalize the raw counts. ",
        "Choices can be: sctransform, LogNormalize, scran, CR, CLR,logtf-idf ,tf-logidf,logtf-logidf,idf etc.",
        "For feature-barcoding/cite-seq, CLR is recommended. For scRNA-seq, sctransform is recommended.",
        "Notice that if sctransform used, the default assay should be changed to SCT automatically in the later analysis.",
        "For ATAC assay, logTF-IDF id recommended."
    )
)

sub_qc$add_argument(
    "-t",
    "--sct_split",
    #type = "character",
    default = NULL,
    help = "[OPTIONAL].do sctransform  using Split object or not ,eg: NULL or sampleid or projectid [default: NULL]"
)
sub_qc$add_argument(
    "-v",
    "--pointsize",
    type = "double",
    default = 1,
    help = "the size for each point in QC metrics vlnplot plot. Set 0 for point-freee vlnplot."
)
sub_qc$add_argument(
    "--crop",
    default = 'TRUE',
    help = "whether to crop in spatialplot, data from cytassist project should be 'TRUE'.[default: %(default)s]"
)

# ========= Subcmd: clustering for dimension reduction and clustering ==========
args <- commandArgs(TRUE)
if ("qc" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "qc") {
        # ===================================================================================================================
        futile.logger::flog.info("step1 read the specified assay and data slot in data object into memory")
        # ===================================================================================================================
        suppressMessages(data_ob <- OESingleCell::ReadX(
            input = opt$input,
            informat = opt$informat,
            assays = assays,
            data.use = dataslots
        ))
        data_ob<- subset(data_ob, subset = nCount_Spatial > 0) #在SCT之前过滤掉总UMI数为0的spot，防止报错
        misc<-data_ob@misc
        # =================================================================================================================
        futile.logger::flog.info("step2 QC by gene,mt and so on")
        futile.logger::flog.info("2.1. parse the parameters for metrics' bounds")
        ## 1.parse norm parameters
        normmeth_list <- sapply(unlist(strsplit(opt$normmeth, "\\|", perl = T)), function(i) strsplit(i, ","))
        names(normmeth_list) <- assays

        ## 2.parse QC parameters

        vars2regress_list <- sapply(
          ifelse(is.null(opt$assay) ,opt$vars2regress,unlist(strsplit(opt$vars2regress, "\\|", perl = T))),
          function(i) ifelse( i== "NULL", list(NULL), strsplit(i, ","))
          #function(i)   strsplit(i, ",")
        )
        names(vars2regress_list) <- assays

        ## lower_threhold
        filter_params <- unlist(strsplit(opt$filters, ",", perl = T))
        lower_threshold <- sapply(
            unlist(strsplit(opt$lower, ",")),
            function(x) ifelse(x == "NULL", NA, as.numeric(x))
        )
        if (length(lower_threshold) != length(filter_params)) {
              stop("The lower threshold setting is not consistant with the parameters in --filters!")
          }
        names(lower_threshold) <- filter_params

        ## upper_threhod
        upper_threshold <- sapply(
            unlist(strsplit(opt$upper, ",")),
            function(x) ifelse(x == "NULL", NA, as.numeric(x))
        )
        if (length(upper_threshold) != length(filter_params)) {
              stop("The upper threshold setting is not consistant with the parameters in --filters!")
          }
        names(upper_threshold) <- filter_params

        ##
        bounds_list <- list()
        for (x in filter_params) bounds_list[[x]] <- c(min = unname(lower_threshold[x]), max = unname(upper_threshold[x]))
        # ===================================================================================================================
        futile.logger::flog.info("2.2. auto calculate for percent.mito cutoff")
        if ((!is.null(bounds_list$percent.mito))) {
            if (is.na(bounds_list$percent.mito["max"])) {
                Q95 <- ceiling(quantile(data_ob$percent.mito, 0.95) / 0.05) * 0.05
                if (Q95 > 0.3) {
                    futile.logger::flog.info("WARNING: The 95th percentile of percent.mito excede 30%. ")
                    Q95 <- 0.3
                }
                futile.logger::flog.info(paste0("The percent.mito cut-off is set as ", Q95))
                bounds_list$percent.mito["min"] <- 0
                bounds_list$percent.mito["max"] <- Q95
            }
        }
        # ===================================================================================================================
        futile.logger::flog.info("2.3. calculate the QC metrics for each cell")
        # The number of genes and UMIs are automatically calculated for every object by Seurat.
        # For non-UMI data, nUMI represents the sum of the non-normalized values within a cell.
        # We calculate the percentage of mitochondrial genes here and store it in percent.mito using AddMetaData.
        # AddMetaData adds columns to cell annotation table, and is a great place to stash QC stats
        # ,since this represents non-transformed and non-log-normalized counts.
        # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
        # calculate the mitochondrial gene derived transcript proportion and add to the metadata
        # assays = ifelse( !is.null(opt$subassay), opt$assay, c(opt$assay, unlist(strsplit(opt$subassay,","))))
        for (assayx in assays) {
            if (assayx %in% c("ADT", "Spatial", "CRISPR")) {
                futile.logger::flog.info(glue::glue("skip this step for {assayx}"))
                next
            }
            # Filter unquanlified cells using different methods
            # We filter out cells that have unique gene counts over 2,500 or less than 200. Note that low.thresholds and
            # high.thresholds are used to define a 'gate'. -Inf and Inf should be used if you don't want a lower or upper
            # threshold. If each filter's lower or upper threshold is specified from the command line the user defined value
            # will be used. Otherwise,the default threshold will be calculated using the distribution significant level
            # specified automatically
            # All the thresholds are specific for each sample, so when the threshold and the standard variance fold are both
            # available from the command line, we will compare the value with the one derived from the distribution lower outliner
            outliers <- OESingleCell::FindOutliers(
                data_ob,
                vars = filter_params,
                var.limit = bounds_list,
                batch = "sampleid",
                type = "both",
                cut.1 = opt$cut.1,
                cut.2 = opt$cut.2,
                n = opt$nfold,
                log = FALSE
            )
            outliercells <- do.call(cbind, outliers)
            metric_outlier <- apply(outliercells, 1, function(x) any(x == T))
            data_ob <- OESingleCell::AddMetaData(
                data_ob,
                metadata = metric_outlier,
                col.name = "is_metric_outlier"
            )
            outlier_variables <- "is_metric_outlier"
            # =============================================================================================================
            if (!is.null(opt$rlm.xy)) {
                futile.logger::flog.info("find outlier cells using roboust linear model between the specified corresponding variable and the dependent variable")
                rlm.xy <- unlist(strsplit(opt$rlm.xy, ":"))
                rlmoutlier <- OESingleCell::FindRlmOutlier(data_ob, y = rlm.xy[2], x = rlm.xy[1]) ## 基于基因或者分子，使用回归分析过滤细胞
                data_ob <- OESingleCell::AddMetaData(
                    data_ob,
                    metadata = rlmoutlier,
                    col.name = "is_rlmoutliers"
                )
                outlier_variables <- c(outlier_variables, "is_rlmoutliers")
            }
            # =============================================================================================================
            futile.logger::flog.info("get is_valid_cell information")
            is_valid_cell <- !apply(Seurat::FetchData(data_ob, vars = outlier_variables), 1, function(x) any(x == T))
            data_ob <- OESingleCell::AddMetaData(data_ob, metadata = is_valid_cell, col.name = "is_valid")
            futile.logger::flog.info("statistics of metadat before and after QC")
            thresholds <- do.call(
                cbind,
                lapply(names(outliers), function(x) {
                    cut <- t(attr(outliers[[x]], "threshold"))
                    colnames(cut) <- paste(x, colnames(cut), sep = "_")
                    cut
                })
            )
        }
        # TO DO
        # add support of median calculation
        # 3 valid digits
        # vlnplot before and after QC
        # =================================================================================================================
        futile.logger::flog.info("2.4 vlnplot and statistics of metadata berfore QC")
        crop <- as.logical(opt$crop)
        statistics_beforeQC <- summary_stat_plot(
            object = data_ob,
            show_features = filter_params,
            summx_features =  unique(c(filter_params,
                                       unlist(vars2regress_list[[assayx]] %>% purrr::keep( ~ !is.null(.)) ))),
            output_dir = output_dir,
            prefix = ifelse(
                assayx %in% c("ADT", "Spatial", "CRISPR"),
                "QC",
                "beforeQC"
            ),
            pointsize = opt$pointsize,
            crop = crop
        )
        # =================================================================================================================
        futile.logger::flog.info("2.5 vlnplot and statistics of metadata after QC")
        if (!assayx %in% c("ADT", "Spatial", "CRISPR")) {
            statistics_afterQC <- summary_stat_plot(
                object = subset(data_ob, subset = is_valid == TRUE),
                show_features = filter_params,
                summx_features = unique(c(filter_params,
                                          unlist(vars2regress_list[[assayx]] %>% purrr::keep( ~ !is.null(.)) ))),
                output_dir = output_dir,
                prefix = "afterQC",
                pointsize = opt$pointsize,
                crop = crop
            )
            cell_statitics <- cbind(statistics_beforeQC, statistics_afterQC, thresholds) ## merge the statitics by columns
        } else {
            futile.logger::flog.info(glue::glue("skip this step for {assayx}"))
            cell_statitics <- statistics_beforeQC
        }
        cell_statitics <- cell_statitics %>%
            tibble::rownames_to_column(var = "sample") %>%
            dplyr::select(sample, dplyr::everything())
        readr::write_tsv(
            cell_statitics,
            glue::glue("{output_dir}/statitics_for_QC.xls")
        )
        # }
        # =================================================================================================================
        futile.logger::flog.info("2.6 filter genes by mincellforgene {opt$mincell4gene}")
        # filter genes by using the minimium cell number one gene is detected this step shoud be run after cell filtering,
        # because there may be satisfied cell number for one gene before filtering but fails after cell filtering.
        ## 1.import genes list to filtered out
        if (!assays %in% c("ADT", "Spatial", "CRISPR")) {
            if (!is.null(opt$features2filter) & file.exists(opt$features2filter)) {
                futile.logger::flog.info(glue::glue("using the {opt$features2filter} input gene column  as remove genes"))
                features2filter <- read.csv(opt$features2filter, sep = ",", header = T)
                features2filter <- as.vector(features2filter$gene)
            } else {
                features2filter <- NULL
            }
            ## 2.filter genes
            data_ob <- OESingleCell::FilterGenes(
                data_ob,
                assay = assays[1],
                min.cells = ifelse(
                    opt$mincell4gene < 1,
                    # the parameter is a percentage
                    round(opt$mincell4gene * ncol(data_ob)),
                    # the parameter is a integer
                    opt$mincell4gene
                ),
                filter.genes = features2filter
            )
            # # 3.now actually filter outlier cells
            # data_ob <- subset(data_ob, subset = is_valid == TRUE )
        } else {
            futile.logger::flog.info(glue::glue("skip this step for {assayx}"))
        }

        if (opt$rmdoublets == "TRUE") {
            # =================================================================================================================
            futile.logger::flog.info("step3.Remove doublets cells")
            # =================================================================================================================
            for (assayx in assays) {
                if (assayx == "RNA") {
                    Seurat::DefaultAssay(data_ob) <- assayx
                    # =================================================================================================================
                    futile.logger::flog.info(glue::glue("step3.1: define double cells ratio for {assayx}"))
                    ## https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
                    def_double_rat <- function(x) {
                        ## CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf /CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevE.pdf
                        double_rat_data <- data.frame(
                            range1 = c(1, 750, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500, 9500, 10500),
                            range2 = c(749, 1499, 2499, 3499, 4499, 5499, 6499, 7499, 8499, 9499, 10499, Inf),
                            dbl_rates = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076, 0.1)
                        )
                        double_rat_data <- double_rat_data %>% dplyr::filter(length(Seurat::Cells(x)) > range1 & length(Seurat::Cells(x)) < range2)
                        return(double_rat_data$dbl_rates)
                    }

                    futile.logger::flog.info(glue::glue("step3.2: select method: {opt$method} dealing for doublets remove")) # ======================
                    is_doublets <- switch(
                        opt$method,
                        "scrublet" = {
                            data_list <- OESingleCell::SplitObject(data_ob, split.by = "sampleid")
                            dbl_all <- future.apply::future_lapply(
                                data_list,
                                function(x) {
                                    dbl_rates <- def_double_rat(x)
                                    x <- OESingleCell::RunScrublet(
                                      x,
                                      doublet.rate = dbl_rates
                                    )
                                    dbl_res <- Seurat::FetchData(x, var = "predicted_doublets")
                                },
                                future.seed = 123
                            )
                        },
                        "doubletfinder" = {
                            cellmeta <- OESingleCell::colData(data_ob)
                            if (!opt$ident %in% colnames(cellmeta)) {
                                futile.logger::flog.info("NO specfied cell identity column FOUND in cell annotation!")
                                data_ob <- Seurat::FindVariableFeatures(
                                    data_ob,
                                    loess.span = 0.3,
                                    clip.max = "auto",
                                    mean.function = "FastExpMean",
                                    dispersion.function = "FastLogVMR",
                                    num.bin = 20,
                                    nfeature = "2000",
                                    binning.method = "equal_width"
                                )
                                data_ob <- Seurat::ScaleData(data_ob, features = rownames(data_ob), verbose = T)
                                data_ob <- Seurat::RunPCA(
                                    data_ob,
                                    npcs = 30,
                                    features = OESingleCell::VariableFeatures(data_ob),
                                    do.print = F,
                                    verbose = F
                                )
                            }
                            dbl_all <- list()
                            data_list <- OESingleCell::SplitObject(data_ob, split.by = "sampleid")
                            for (x in 1:length(data_list)) {
                                dbl_rates <- def_double_rat(x)
                                print("OK!!!!!!")
                                data_list[[x]] <- OESingleCell::RunDoubletFinder(
                                  data_list[[x]],
                                  PCs = 1:10,
                                  doublet.rate = dbl_rates,
                                  identity = opt$ident,
                                  sct = ifelse(
                                    normmeth_list[[assayx]] == "sctransform",
                                    TRUE,
                                    FALSE
                                  )
                                )
                                dbl_all[[x]] <- Seurat::FetchData(data_list[[x]], var = "predicted_doublets")
                            }
                            dbl_all
                        }
                    )
                    data_ob <- OESingleCell::AddMetaData(
                        data_ob,
                        metadata = do.call(rbind, is_doublets)[[1]],
                        col.name = "is_doublets"
                    )
                    # append all other data not loaded before subsetting for h5seurat formated file.
                    # if ( opt$informat == "h5seurat" ) data_ob = SeuratDisk::AppendData(opt$input, data_ob)
                    # TO DO
                    # add support of median calculation
                    # 3 valid digits
                    # futile.logger::flog.info("vlnplot ,statistics of metadata before  rmdoublets")
                    # statistics_before_rmdoublets <-summary_stat_plot(object=data_ob,
                    #                                                     show_features=filter_params,
                    #                                                     summx_features=unique(c(filter_params, unlist(vars2regress_list) )),
                    #                                                     output_dir=output_dir,
                    #                                                     image=FALSE,
                    #                                                     prefix="before_rmdoublets",
                    #                                                     pointsize=opt$pointsize )
                    futile.logger::flog.info("statistics of metadata after rmdoublets")
                    statistics_after_rmdoublets <- summary_stat_plot(
                        object = subset(data_ob, subset = is_doublets == FALSE),
                        show_features = filter_params,
                        summx_features =  unique(c(filter_params,
                                                   unlist(vars2regress_list[[assayx]] %>% purrr::keep( ~ !is.null(.)) ))),
                        output_dir = output_dir,
                        image = FALSE,
                        prefix = "after_rmdoublets",
                        pointsize = opt$pointsize,
                        crop = crop
                    )
                    # if ( opt$QC != "TRUE" ) { ## rmdoublet only
                    #     cell_statitics <- cbind(statistics_before_rmdoublets, statistics_after_rmdoublets)
                    #     cell_statitics <- cell_statitics %>%
                    #                       tibble::rownames_to_column(var = "sample") %>%
                    #                       dplyr::select(sample, dplyr::everything())
                    #     write.table(cell_statitics,
                    #                 file.path(output_dir,"statitics_before_after_rmdoublets.xls"),
                    #                 sep="\t",
                    #                 col.names=T,
                    #                 row.names=F)
                    # } else { ## merge final stats
                    cell_statitics <- cbind(
                        statistics_beforeQC,
                        statistics_afterQC,
                        statistics_after_rmdoublets,
                        thresholds
                    )
                    cell_statitics <- cell_statitics %>%
                        tibble::rownames_to_column(var = "sample") %>%
                        dplyr::select(sample, dplyr::everything())
                    write.table(
                        cell_statitics,
                        file.path(output_dir, "statitics_for_QC-rmdoublets.xls"),
                        sep = "\t",
                        col.names = T,
                        row.names = F
                    )
                }
                data_ob <- subset(data_ob, subset = is_doublets == FALSE)
            }
        }
        # ===================================================================================================================
        futile.logger::flog.info("step4.Normalize the raw counts and scale ")
        # ===================================================================================================================
        for (assayx in assays) {
            ##### ============================================================================================================
            if (assayx == "ATAC") {
                Seurat::DefaultAssay(data_ob) <- assayx
                futile.logger::flog.info(glue::glue("Normalize the raw counts for {assayx} using method {normmeth_list[[assayx]]}"))
                tfidf_method <- c("logtf-idf" = 1, "tf-logidf" = 2, "logtf-logidf" = 3, "idf" = 4)
                data_ob <- Signac::RunTFIDF(
                    data_ob,
                    assay = assayx,
                    method = tfidf_method[normmeth_list[[assayx]]]
                )
                data_ob <- Signac::FindTopFeatures(data_ob, min.cutoff = opt$nvfeatures)
                #vars2regress <- unique(vars2regress_list[[assayx]])
                data_ob <- Seurat::ScaleData(
                    data_ob,
                    features = rownames(data_ob),
                    verbose = T,
                    vars.to.regress = vars2regress_list[[assayx]]
                )
            }
            ##### ============================================================================================================
            if (assayx == "RNA" | assayx == "Spatial") {
                ##### ============================================================================================================
                futile.logger::flog.info(glue::glue("Normalize the raw counts for {assayx} using method {normmeth_list[[assayx]]}"))
                Seurat::DefaultAssay(data_ob) <- assayx
                if (normmeth_list[[assayx]] == "sctransform") {
                    # Results are saved in a new assay (named SCT by default) with counts being (corrected) counts, data being log1p(counts),
                    # scale.data being pearson residuals; sctransform::vst intermediate results are saved in misc slot of new assay.
                    # normalize data with SCTransform()
                    ## https://satijalab.org/seurat/articles/cell_cycle_vignette.html
                    if ("CC.Difference" %in% vars2regress_list[[assayx]]) {
                        s.genes <- Seurat::CaseMatch(search = cc.genes$s.genes, match = rownames(data_ob))
                        g2m.genes <- Seurat::CaseMatch(search = cc.genes$g2m.genes, match = rownames(data_ob))
                        data_ob <- Seurat::CellCycleScoring(data_ob, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
                        cycleS <- Seurat::FetchData(data_ob, vars = c("S.Score", "G2M.Score"))
                        data_ob <- OESingleCell::AddMetaData(
                          data_ob,
                          col.name = "CC.Difference",
                          metadata = cycleS[["S.Score"]] - cycleS[["G2M.Score"]]
                        )
                    }
                    # normalise again but this time including also the cell cycle scores
                    #vars.to.regress_params <- ifelse(vars2regress_list[[assayx]]=="NULL", NULL, vars2regress_list[[assayx]])
                    futile.logger::flog.info(glue::glue("using vars.to.regress is {vars2regress_list[[assayx]]}"))
                    if(is.null(opt$sct_split)){
                        futile.logger::flog.info(glue::glue("sctransform on whole object  and return {opt$nvfeatures} variable features"))
                        data_ob <- Seurat::SCTransform(
                          data_ob,
                          assay = assayx,
                          method = "glmGamPoi", ##  glmGamPoi package which substantially improves the speed of the learning procedure.
                          vars.to.regress = vars2regress_list[[assayx]],
                          variable.features.n = opt$nvfeatures,
                          verbose = FALSE,
                          min_cells = 5,  ##Set filtering of genes to FALSE, default is 5 ,ref:https://github.com/satijalab/sctransform/issues/27，作者不建议将此值设置的太低，https://github.com/satijalab/sctransform/issues/22
                          return.only.var.genes = TRUE)
                    }else{
                        futile.logger::flog.info(glue::glue("sctransform on {opt$sct_split} split object and return {opt$nvfeatures} variable features"))
                        data_ob.list <- Seurat::SplitObject(data_ob, split.by=opt$sct_split)
                        ## remove used images for visium data
                        for(i in 1:length(data_ob.list)){
                            if (!is.null(Seurat::Images(data_ob.list[[i]]))) {
                                unuse_images <- Seurat::Images(data_ob.list[[i]])[!Seurat::Images(data_ob.list[[i]]) %in% (data_ob.list[[i]]@meta.data$sampleid %>% unique())]
                                if (length(unuse_images) > 0) {
                                    data_ob.list[[i]]@images[unuse_images] <- NULL
                                } }
                        }
                        ##
                        data_ob.list <- lapply(X = data_ob.list,
                                               FUN = Seurat::SCTransform,
                                                assay = assayx,
                                                method = "glmGamPoi", ##  glmGamPoi package which substantially improves the speed of the learning procedure.
                                                vars.to.regress =  vars2regress_list[[assayx]]  ,
                                                variable.features.n =opt$nvfeatures,
                                                verbose = FALSE,
                                                min_cells=5,  ##Set filtering of genes to FALSE, default is 5 ,ref:https://github.com/satijalab/sctransform/issues/27
                                                return.only.var.genes = TRUE)
                       var.features <- Seurat::SelectIntegrationFeatures(object.list = data_ob.list, nfeatures = opt$nvfeatures)
                       data_ob <- merge(x = data_ob.list[[1]], y = data_ob.list[2:length(data_ob.list)], merge.data=TRUE)
                       data_ob@misc <- misc
                       Seurat::VariableFeatures(data_ob) <- var.features
                    }
                } else {
                    # When you initialise your Seurat object, both counts and data contain your raw transcripts counts (assuming that's
                    # your raw data). While the matrix stored in counts generally remains the raw data, the data in the data slot will
                    # be normalised when you run NormalizeData()
                    data_ob <- Seurat::NormalizeData(
                        data_ob,
                        assay = assayx,
                        normalization.method = normmeth_list[[assayx]],
                        scale.factor = 10000
                    )
                    # wether to regress out the cell cycle effect
                    # genes.inuse <- rownames(Seurat::GetAssayData(data_ob, assay = assayx, slot = "counts"))
                    if ("CC.Difference" %in% vars2regress_list[[assayx]]) {
                        futile.logger::flog.info(glue::glue("cell cycle effect for {assayx}"))
                        s.genes <- Seurat::CaseMatch(search = cc.genes$s.genes, assay = assayx, match = rownames(data_ob))
                        g2m.genes <- Seurat::CaseMatch(search = cc.genes$g2m.genes, assay = assayx, match = rownames(data_ob))
                        data_ob <- Seurat::CellCycleScoring(
                            data_ob,
                            assay = assayx,
                            s.features = s.genes,
                            g2m.features = g2m.genes,
                            set.ident = F
                        )
                        cycleS <- Seurat::FetchData(data_ob, assay = assayx, vars = c("S.Score", "G2M.Score"))
                        data_ob <- OESingleCell::AddMetaData(
                            data_ob,
                            col.name = "CC.Difference",
                            metadata = cycleS[["S.Score"]] - cycleS[["G2M.Score"]]
                        )
                    }
                    # ===================================================================================================================
                    futile.logger::flog.info(glue::glue("FindVariableFeatures and ScaleData for {assayx}"))
                    #vars2regress <- unique(vars2regress_list[[assayx]])
                    data_ob <- Seurat::FindVariableFeatures(
                        data_ob,
                        assay = assayx,
                        loess.span = 0.3,
                        clip.max = "auto",
                        mean.function = "FastExpMean",
                        dispersion.function = "FastLogVMR",
                        num.bin = 20,
                        nfeature = opt$nvfeatures,
                        binning.method = "equal_width"
                    )
                    # regress out all the specified the varibales
                    data_ob <- Seurat::ScaleData(
                        data_ob,
                        assay = assayx,
                        features = rownames(data_ob),
                        vars.to.regress = vars2regress_list[[assayx]],
                        verbose = T
                    ) # takes some time
                }
            }
        }
        ## save seurat object===============================================================================================
        OESingleCell::SaveX(
            data_ob,
            output = output_dir,
            outformat = opt$outformat,
            prefix = opt$prefix,
            update = FALSE
        )
        ## save session informations
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}
