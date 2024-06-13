# sce: singlecellexperiement object supported by R packages:scater,scran etc.
# h5: hdf5 formatted molecular quantification from cellranger etc.
# loom: loom format, a specialized hdf5 format.
# h5seurat: the hdf5 formated seurat object.
# h5ad: the format support by scanpy from anndata in python.
# aggr: the result from cellranger aggr
# tenx:the directory of cellranger count/aggr output with sampleid named subdirectory.
# xsv: the raw gene count matrix file, it can be very large. Please make sure the data path is Your_Project/sampleid/sampleid.(tsv|csv).
# cellDataSet: the format for monocle.
# Seurat: the  binary seurat object .rds from the clustering results.
# visium: the output from spaceranger ananlysis of 10X Visium spatial transcriptomics data.", '"""'))
# se: summarizedExperient object in rds format
# htseq: the results from htseq-count
# vdj: the results from cellranger vdj
# smart-seq2:
# star: the results from STAR alignment result TO DO
# create the parsers for subcommand create
sub_create <- subparsers$add_parser("create",
                                    help = "make the specified format data object out of the supported input data sources.")
sub_create$add_argument("-s", "--source", type = "character", default = "mtx",
                        help = "The source of data, possible choices:h5ad, h5, mtx, xsv(including BD,tsv and csv), Visium, slideseq.[default: %(default)s]")
sub_create$add_argument("-m", "--metadata", type = "character",
                        help = "the sample metadata which must include sample id in this assay design.")
sub_create$add_argument("-x", "--metrics", type = "character", default = "percent.mito",
                        help = "the additional QC metrics list to calculate in advance.[default %(default)s]")
sub_create$add_argument("--gcolumn", type = "integer", default = 2,
                        help = "Specify which column of genes.tsv or features.tsv to use for feature names.[default %(default)s]")
sub_create$add_argument("--aggr", type = "character", default = "FALSE",
                        help = "[OPTIONAL] wether the results derived from cellranger aggr.[default %(default)s]")
sub_create$add_argument("--feature.meta", type = "character", default = NULL,
                        help = "[OPTIONAL]the reference genome annotaion file from 10X.")
sub_create$add_argument("--chrom.size", type = "character", default = NULL,
                        help = "[OPTIONAL]length of each chromosome of the genome file from 10X.")
sub_create$add_argument("--cell.meta", type = "character", default = "FALSE",
                        help = "[OPTIONAL] logical indication wether there is cell annotation file for, 'singlecell.csv' for in ATAC-seq for example.[default: FALSE]")
sub_create$add_argument("--fragment", type = "character", default = "FALSE",
                        help = "[OPTIONAL] use the fragment annotion file 'fragments.tsv.gz' in outs directory from cellranger-atac count results for each sample, Only valid for scATAC-seq assay.[default FALSE]")
sub_create$add_argument("--transpose", type = "character", default = "FALSE",
                        help = "[OPTIONAL]wether to transpose the count matrix when comes to the cell-feature matrix.[default %(default)s]")
sub_create$add_argument("--gset", type = "character", default = NULL,
                        help = "[OPTIONAL]the gmt formated customized genes set file for QC metrics' calculation.")
sub_create$add_argument("--bin.size", default = NULL,
                        help = "[OPTIONAL]bin size setting for stereo-seq data, must be consist with the input h5ad bin size,if is null, means cellbin.")
args <- commandArgs(TRUE)
if ( "create"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name == "create"){
  data_ob <- switch(tolower(opt$source),
                    "mtx" = { # data come from cellranger count results in mtx format
                      seux <- OESingleCell::CreateX(data.dir = opt$input,
                                                    assay.meta = opt$metadata,
                                                    source = opt$source, assay = opt$assay,
                                                    aggr = as.logical(opt$aggr),
                                                    add.meta = as.logical(opt$cell.meta),
                                                    use.fragment = as.logical(opt$fragment),
                                                    gene.column = opt$gcolumn,
                                                    gtf = opt$feature.meta,
                                                    chrom.size = opt$chrom.size)
                    },
                    "h5" = { # read raw data from cellranger count results in hdf5 format
                      seux <- OESingleCell::CreateX(data.dir = opt$input,
                                                    assay.meta = opt$metadata,
                                                    add.meta = as.logical(opt$cell.meta),
                                                    use.fragment = as.logical(opt$fragment),
                                                    source = opt$source,
                                                    assay = opt$assay,
                                                    aggr = as.logical(opt$aggr),
                                                    gtf = opt$feature.meta,
                                                    chrom.size = opt$chrom.size)
                    },
                    "h5ad" = { # read raw data from cellranger count results in stereo h5ad format
                      seux <- OESingleCell::CreateX(data.dir = opt$input,
                                                    assay.meta = opt$metadata,
                                                    add.meta = as.logical(opt$cell.meta),
                                                    use.fragment = as.logical(opt$fragment),
                                                    source = opt$source,
                                                    assay = opt$assay,
                                                    aggr = as.logical(opt$aggr),
                                                    bin_size = opt$bin.size,
                                                    gtf = opt$feature.meta,
                                                    chrom.size = opt$chrom.size)
                    },
                    "xsv" = { # data from text file
                      seux <- OESingleCell::CreateX(data.dir = opt$input,
                                                    assay.meta = opt$metadata,
                                                    add.meta = as.logical(opt$cell.meta),
                                                    source = opt$source,
                                                    assay = opt$assay,
                                                    aggr = as.logical(opt$aggr),
                                                    gene.column = opt$gcolumn,
                                                    transpose = as.logical(opt$transpose),
                                                    random.barcode = FALSE,
                                                    names.field = 1,
                                                    names.delim = "_")
                    }
  )

  if (!is.null(opt$gset)) {
    futile.logger::flog.info("进行线粒体比例统计等")
    # gmt = GSEABase::getGmt(con=opt$gset)
    gset_list <- GSEABase::geneIds(GSEABase::getGmt(con = opt$gset))
    data_ob <- OESingleCell::CalculateMetrics(data_ob, metrics = NULL, gset = gset_list)
    print(head(data_ob@meta.data))
  }

  # metrics = c("percent.mito", "CC.difference"),
  if (!is.null(opt$metrics)) {
    metrics <- unlist(strsplit(opt$metrics, ","))
    data_ob <- OESingleCell::CalculateMetrics(data_ob, metrics = metrics)
  }

  if (!is.null(opt$feature.meta)) {
    data_ob <- OESingleCell::addFeatureAnno(data_ob,
                                            assay = opt$assay,
                                            gtf = opt$feature.meta)
  }

  if(opt$assay=='Spatial' ){
    if(opt$source=="h5"){
        imagedir <- glue::glue("{opt$input}/{unique(data_ob$sampleid)}/outs/spatial/tissue_hires_image.png")
        scalefactors_dir <- glue::glue("{opt$input}/{unique(data_ob$sampleid)}/outs/spatial/scalefactors_json.json")
        #添加image information
        sampleimage <- lapply(imagedir, png::readPNG)
        names(sampleimage) <- unique(data_ob$sampleid)
        data_ob@misc$hires_image <- sampleimage
        #添加jsonfile
        jsonfile <- lapply(scalefactors_dir,jsonlite::fromJSON)
        names(jsonfile) <- unique(data_ob$sampleid)
        data_ob@misc$scalefactors <- jsonfile
    }
  }
  # Now the h5seurat format is used as the standard data object in disk or better interpre
  if (tolower(opt$outformat) == "h5seurat") {
    futile.logger::flog.info("生成h5seurat对象")
    prefix <- file.path(output_dir, opt$prefix)
    SeuratDisk::SaveH5Seurat(data_ob,
                             filename = glue::glue("{prefix}.h5seurat"),
                             overwrite = TRUE,
                             verbose = FALSE)
  }else {
    futile.logger::flog.info("生成seurat对象")
    OESingleCell::SaveX(data_ob,
                        output = output_dir,
                        outformat = opt$outformat,
                        prefix = opt$prefix,
                        update = FALSE)
  }        ## save session informations
  write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
  quit()
  }
}

# =============== Subcmd: convert, change the data object to target format ============
sub_convert <- subparsers$add_parser("convert",
                                     help = "convert the data object from input format to the target format.")
sub_convert$add_argument("--from", type = "character", default = "h5seurat",
                         help = "the format of the input data object, choices can be: monocle2(CellDataSet), monocle3(cell_data_set), Seurat, sce(SingleCellExperiment), anndata.[default: %(default)s] ")
sub_convert$add_argument("--to", type = "character", default = "",
                         help = "the destination output format of the data object, choices can be: monocle2(CellDataSet), monocle3(cell_data_set), Seurat, sce(SingleCellExperiment), anndata.[default: %(default)s] .")
# sub_convert$add_argument("--seperateby", type = "character", default = NULL,
#      help = "[OPTIONAL]save the expression matrix in  mtx format like cellranger output for each sample." )

args <- commandArgs(TRUE)
if ( "convert"  %in% args ){
  opt<-intial_setting()
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots)
  x_mapping <- c(monocle2 = "CellDataSet",
                 monocle3 = "cell_data_set",
                 sce = "SingleCellExperiment",
                 anndata = "AnnDataR6",
                 Seurat = "Seurat")
  data_obx <- OESingleCell::ConvertX(data_ob,
                                     from = x_mapping[opt$from],
                                     to = x_mapping[opt$to])
  OESingleCell::SaveX(data_obx,
                      output = opt$input,
                      assay = assays[1],
                      outformat = opt$outformat,
                      update = F)
}


# =============== Subcmd: downsample, subset cells randomly while keep heterogencity ============
# create the parsers for subcommand downsample
sub_downsample <- subparsers$add_parser("downsample",
                                        help = "downsample the exppression matrix in case of too many cell and scarcely disequal cell numbers in each sample")
sub_downsample$add_argument("--targetN", type = "integer",
                            help = "the number of cell to keep after downsample.")
sub_downsample$add_argument("--usePCA", type = "character", default = "TRUE",
                            help = "logical parameter indicating wether to use fbpca during before downsample.[default: %(default)s]")
sub_downsample$add_argument("--subsampleby", type = "character",
                            help = "radom subsample by the specified group using the specified ratio,example:'clusters:0.5'.")

args <- commandArgs(TRUE)
if ("downsample" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "downsample") {
      # ================================================================================================================
      futile.logger::flog.info("step1 read the specified assay and data slot in data object into memory")
      # ================================================================================================================
      data_ob <- OESingleCell::ReadX(input = opt$input,
                                     informat = opt$informat,
                                     assays = assays,
                                     data.use = dataslots)
      # ================================================================================================================
      futile.logger::flog.info("step1 read the specified assay and data slot in data object into memory")
      # ================================================================================================================
      data_ob <- OESingleCell::DownSample(  data_ob,
                                            opt$targetN,
                                            slot =dataslots ,
                                            use_PCs = opt$usePCA,
                                            PCs = 100,
                                            k = "auto",
                                            seed = 4563,
                                            replace = FALSE,
                                            alpha = 0.1,
                                            max_iter = 200,
                                            verbose = F)
      ## save seurat object=============================================================================================
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

}}


# =============== Subcmd: fetch, extract desired cell annotation ========================================================
# create the parsers for subcommand fetech
sub_fetch <- subparsers$add_parser("fetch",
                                   help = "Retreives data (feature expression, PCA scores, metrics, etc.)")
# for a set of cells.A table with cells as rows and cellular data as columns returns.
sub_fetch$add_argument("--vars",
                       help = "comma seperated List of all variables to fetch, use keyword ident to pull identity classes.")
sub_fetch$add_argument("--header",
                       help = "[OPTIONAL]change the header of the output table if specified. Notice that the first column is always the cell barcode.")

# Retreives data (feature expression, PCA scores, metrics, etc.) for a set of cells in a Seurat object
args <- commandArgs(TRUE)
if ( "fetch"  %in% args ){
  opt<-intial_setting()
  if (!is.null(opt$vars)) {
    var_list <- unlist(strsplit(opt$vars, ",", perl = T))
  }else {
    var_list <- c("group", "ident")
  }
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots)

  desired_df <- Seurat::FetchData(data_ob, slot = dataslots[1], vars = var_list)
  desired_df$cellbarcode <- rownames(desired_df)
  if (!is.null(opt$header)) {
    colnames(desired_df) <- unlist(strsplit(opt$header, ",", perl = T))
  }
  write.table(desired_df,
              file.path(output_dir, glue::glue("{opt$prefix}_cellinfo.xls")),
              row.names = F,
              col.names = T,
              sep = "\t")
  quit()
}

# subcommand merge is invoked
# create the parsers for subcommand merge
sub_merge <- subparsers$add_parser("merge", help = "merge data object")
sub_merge$add_argument("--toadd", type = "character",
                       help = "one or more of object list to merge into the target object.")
sub_merge$add_argument("--adformat", type = "character",
                       help = "The format of input expression matrix to add")
sub_merge$add_argument("--mergedata", type = "character",
                       help = "logical string T/F to indicate wether to merge the data slot")
args <- commandArgs(TRUE)
if ( "merge"  %in% args ){
  opt<-intial_setting()
  seuy <- lapply(strsplit(opt$toadd, ","), function(seu) {
    seux <- readObject(input = seu,
                       assay = opt$assay,
                       informat = opt$adformat)
    if (seux@version < 3) {
      seux <- UpdateSeuratObject(seux)
    }
  })
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots)
  data_ob <- merge(data_ob, y = seuy, merge.data = T)
  saveObject(data_ob,
             outdir = output_dir,
             splitby = opt$seperateby,
             outformat = opt$outformat,
             prefix = opt$prefix)
  quit()
}


#================================================================================================================
### create the parsers for subcommand multimodal
#================================================================================================================
sub_multimodal <- subparsers$add_parser("multimodal",
                                        help = "make the specified format data object out of cellranger multimodal assay, currently only scRNA+scATAC.")
sub_multimodal$add_argument("-s", "--source", type = "character", default = "h5",
                            help = "The source of multimodal assay data, currently only hdf5 formated files are supported.[default %(default)s]")
sub_multimodal$add_argument("-m", "--metadata", type = "character",
                            help = "the sample metadata which must include sample id in this assay design.")
sub_multimodal$add_argument("-x", "--metrics", type = "character", default = "percent.mito",
                            help = "the additional QC metrics list to calculate in advance.[default %(default)s]")
sub_multimodal$add_argument("--gset", type = "character", default = NULL,
                            help = "[OPTIONAL]the gmt formated customized genes set file for QC metrics' calculation.[default: %(default)s] ")
sub_multimodal$add_argument("--gcolumn", type = "integer", default = 2,
                            help = "Specify which column of genes.tsv or features.tsv to use for feature names.[default %(default)s]")
sub_multimodal$add_argument("--feature.meta", type = "character", default = NULL,
                            help = "[OPTIONAL]the reference genome annotaion file from 10X.[default: %(default)s] ")
sub_multimodal$add_argument("--cell.meta", type = "character", default = "FALSE",
                            help = "[OPTIONAL] logical indication wether there is cell annotation file for, 'singlecell.csv' for in ATAC-seq for example.[default: %(default)s] ")
sub_multimodal$add_argument("--fragment", type = "character", default = "FALSE",
                            help = "[OPTIONAL] use the fragment annotion file 'fragments.tsv.gz' in outs directory from cellranger-atac count results for each sample, Only valid for scATAC-seq assay.[default: %(default)s]")

# ================ Subcmd: multimodal, create data object from multimodal assay ========
args <- commandArgs(TRUE)
if ( "multimodal"  %in% args ){
  opt <- intial_setting()
  if(opt$sub_name == "multimodal" ){
    if(tolower(opt$source)=="h5"){
      data_ob <- OESingleCell::CreateMultiModal(data.dir = opt$input,
                                                assay.meta = opt$metadata,
                                                add.meta = as.logical(opt$cell.meta),
                                                use.fragment = as.logical(opt$fragment),
                                                assay = assays[1],
                                                subassay = assays[2],
                                                source = opt$source,
                                                gtf = opt$feature.meta,
                                                chrom.size = NULL)
    if (!is.null(opt$gset)) {
      # gmt = GSEABase::getGmt(con=opt$gset)
      gset_list <- GSEABase::geneIds(GSEABase::getGmt(con = opt$gset))
      data_ob <- OESingleCell::CalculateMetrics(data_ob, metrics = NULL, gset = gset_list)
    }
    if (!is.null(opt$metrics)) {
      metrics <- unlist(strsplit(opt$metrics, ","))
      data_ob <- OESingleCell::CalculateMetrics(data_ob, metrics = metrics)
    }
    if (!is.null(opt$feature.meta)) {
      data_ob <- OESingleCell::addFeatureAnno(data_ob, assay = opt$assay, gtf = opt$feature.meta)
    }
    # Now the h5seurat format is used as the standard data object in disk or better interpre
    if (tolower(opt$outformat) == "h5seurat") {
      prefix <- file.path(output_dir, opt$prefix)
      SeuratDisk::SaveH5Seurat(data_ob, filename = glue::glue("{prefix}.h5seurat"), overwrite = TRUE, verbose = FALSE)
    }else {
      OESingleCell::SaveX(data_ob, output = output_dir, outformat = opt$outformat, prefix = opt$prefix, update = FALSE)
    }
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
    }
  }
}


# ===================== Subcmd: subset the data object using condditional expression on cells =======
# create the parsers for subcommand subset
sub_subset <- subparsers$add_parser("subset",
                                    help = "filter the cells using the specified conditions.")
sub_subset$add_argument("--geneset", type = "character",
                        help = "features/variables to keep in TAB format.")
sub_subset$add_argument("--feature4subset", type = "character",
                        help = "Logical expression on features in matrix indicating cells to keep.Example: 'MSA4>1'.")
sub_subset$add_argument("--levels", type = "character",
                        help = "The subset of variable levels in variable specified by --group2use.")
sub_subset$add_argument("--group2use", type = "character",
                        help = "the cell names or annotation column names set for the cell annotation table.")
sub_subset$add_argument("--low", type = "double", default = NULL,
                        help = "the lower bound of the continious variable specified by --group2use.[default: %(default)s]")
sub_subset$add_argument("--high", type = "double", default = NULL,
                        help = "the higher bound of the continious variable specified by --group2use.[default: %(default)s]")
sub_subset$add_argument("--invert", type = "character", default = "FALSE",
                        help = "reverse the selection conditions specified by all other conditions.[default: %(default)s]")

args <- commandArgs(TRUE)
if ( "subset"  %in% args ){
  opt<-intial_setting()
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots)
  group2use <- opt$group2use
  levels <- opt$levels
  data_obx <- data_ob
  if (is.null(opt$low) | opt$low == "NULL") {
    opt$low <- -Inf
  }
  if (is.null(opt$high) | opt$high == "NULL") {
    opt$high <- Inf
  }
  if (!is.null(levels)) {
    levels_id <- unlist(strsplit(levels, ",", perl = T))
    if (is.null(group2use)) { group2use <- "idents" }
    if (group2use == "idents") {
      if (typeof(Idents(data_obx)) == "integer") { # the idents is clustering id
        data_obx <- SubsetData(data_obx,
                               ident.use = levels_id,
                               low.threshold = opt$low,
                               high.threshold = opt$high)
      }
    } else { #subset cells on metadata using other cell annotations
      data_obx <- SubsetData(data_obx,
                             subset.name = group2use,
                             accept.value = levels_id,
                             low.threshold = opt$low,
                             high.threshold = opt$high)
    }
    desired_cells <- Seurat::FetchData(data_obx, vars = rownames(data_obx)) %>% as.data.frame
  }else {
    desired_cells <- Seurat::FetchData(data_obx, vars = rownames(data_obx)) %>% as.data.frame
  }
  # subset cells on the expression matrix using logical exprssion on features
  if (!is.null(opt$predicate)) {
    # levels = gsub( ".*\\((.*)\\)", "\\1", opt$cellfilter, perl =T)
    # predicate = gsub( "\\(.*\\)", glue::glue("({levels})"), opt$cellfilter, perl = T)

    df <- OESingleCell::colData(data_obx)
    desired_cells <- subset(df, eval(parse(text = opt$predicate)))
    data_obx <- data_obx[, rownames(desired_cells)]
  }

  if (!is.null(opt$feature4subset)) {
    feature_predicate <- opt$feature4subset
    count_df_T <- Matrix::t(GetAssayData(data_obx, assay = opt$assay, slot = opt$slot))
    desired_cells <- subset(as.data.frame(count_df_T), eval(expr = parse(text = feature_predicate)))
    data_obx <- data_obx[, rownames(desired_cells)]
  }

  desired_cell_id <- colnames(data_obx)

  # keep only specified features
  if (!is.null(opt$geneset)) {
    genes <- read.table(opt$geneset, sep = "\t", header = T)
    desired_features <- as.vector(genes$gene)
  }else {
    desired_features <- NULL
  }

  data_ob <- subset(data_ob,
                    features = desired_features,
                    cells = desired_cell_id,
                    invert = as.logical(opt$invert))
  saveObject(data_ob,
             outdir = output_dir,
             splitby = opt$seperateby,
             outformat = opt$outformat,
             prefix = opt$prefix)
  quit()
}

#=======================================================================================================================
# create the parsers for subcommand split

sub_split <- subparsers$add_parser("split", help = "split the object according to the cell metadata.")
sub_split$add_argument("--splitby", type = "character", default = "ident",
                       help = "the cell annotation column name.If NOT supplied, %(default)s will be used as default!")
args <- commandArgs(TRUE)
if ( "split"  %in% args ){
  opt<-intial_setting()
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat) # all assay must be loaded

  if (!is.null(opt$splitby)) {
    seurat_list <- OESingleCell::SplitObject(data_ob,
                                             split.by = opt$splitby)
    lapply(names(seurat_list), function(x) {
      SeuratDisk::SaveH5Seurat(seurat_list[[x]],
                               filename = file.path(output_dir, glue::glue("{x}.rds")),
                               overwrite = T)
    })
  }
  quit()
}


# =============== Subcmd: update the metedata of cells or features in the data_ob ============
# create the parsers for subcommand update
sub_update <- subparsers$add_parser("update",
                                    help = "update the main data object using annotation from other data object or external info.")
sub_update$add_argument("--admeta", type = "character",
          help = "the additional metadata which includes any groupping design in this assay.")
sub_update$add_argument("--annlevel", type = "character", default = "sample",
          help = "the annotation level,sample or cell is supported.[default: %(default)s] ")
sub_update$add_argument("--recode", type = "character", default = NULL,
          help = "rename the level id in specified groups in the file with each column corresponding to header in the object metadata and column element format as from:to.[default: %(default)s] ")
# sub_update$add_argument("--adfeatures", type = "character",
#           help = "the additional metadata of genes in the object.")

args <- commandArgs(TRUE)
if ( "update"  %in% args ){
  opt<-intial_setting()
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots)
  if (!is.null(opt$admeta)) {
    cellnames <- Cells(data_ob)
    additional_metadata <- read.csv(opt$admeta, sep = ",", header = T)
    if (opt$annlevel == "sample") {
      #the index order is the same as the row index of the assay metadata
      sampleidx <- gsub("(_|-)[ATGC]{16,}.*", "", cellnames, perl = T)
      #integrate the additional metadata from the assay design
      additional_cell_meta <- vector()
      for (colidx in colnames(additional_metadata)) {
        additional_cell_meta <- cbind(additional_cell_meta,
                                      as.vector(additional_metadata[sampleidx, colidx]))
      }
      colnames(additional_cell_meta) <- colnames(additional_metadata)
      rownames(additional_cell_meta) <- cellnames
      additional_cell_meta <- as.data.frame(additional_cell_meta)
      data_ob <- AddMetaData(data_ob, additional_cell_meta)
    }else { # the annotation level is cell, make sure the barcodes are in the same pattern as the ones in the object
      additional_cell_meta <- as.data.frame(additional_metadata)
      # TO DO if the cells in the supplied annotation is different from the cells in seurat object
      data_ob <- AddMetaData(data_ob, metadata = additional_metadata)
    }
  }

  cellmeta <- OESingleCell::colData(data_ob)
  if (!is.null(opt$recode)) {
    # the format should be as follows:
    # seurat_clusters
    # 0:Naive CD4 T
    # 1:Memory CD4 T
    # 2:CD14+ Mono
    # 3:B
    # 4:CD8 T

    recode_df <- read.table(opt$recode, header = T, sep = "\t")
    recode_id <- colnames(recode_df)[1]
    from_to <- recode_df %>%
              tidyr::separate(1, c("from", "to"), sep = ":") %>%
              tibble::deframe()
    new_meta <- dplyr::recode(as.vector(cellmeta[[recode_id]]), !!!from_to)
    data_ob <- OESingleCell::AddMetaData(data_ob,
                                         metadata = new_meta,
                                         col.name = recode_id)
  }

  if (as.logical(opt$update)) {
    SeuratDisk::UpdateH5Seurat(file = opt$input,
                               object = data_ob,
                               verbose = FALSE)
  }else {
    OESingleCell::SaveX(data_ob,
                        output = opt$output,
                        update = FALSE,
                        outformat = opt$outformat,
                        prefix = opt$prefix)
  }

  quit()
}


#=======================================================================================================================
sub_write10x <- subparsers$add_parser("write10x", help = "save the data object into the format from 10X Genomics.")
sub_write10x$add_argument("--split.by", "-x", type = "character", default = NULL,
          help = "the groupping variable of cells used to split the data object.[default: %(default)s]")
sub_write10x$add_argument("--overwrite", type = "character", default = "TRUE",
          help = "wether to overwrite the existed data object.[default: %(default)s]")
sub_write10x$add_argument("--version", "-v", type = "character", default = "3",
          help = "the version of cellranger to format the output to.[default: %(default)s]")
sub_write10x$add_argument("--h5", type = "character", default = "FALSE",
          help = "Wether to save the count matrix in hdf5 format defined by 10X Genomics.[default: %(default)s]")

args <- commandArgs(TRUE)
if ( "write10x"  %in% args ){
  opt<-intial_setting()
  # read the specified assay and data slot in data object into memory
  # only counts dataslot needed here
  is.overwrite <- as.logical(opt$overwrite)
  is.h5 <- as.logical(opt$h5)
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots,
                                 verbose = F)
  message("Loading Data object Finished!")
  if (is.h5) {
    # TO DO
    warning("It is not implemented NOW!")
  }else {
    Write10X(data_ob, assay = assays,
             split.by = opt$split.by,
             version = opt$version,
             path = output_dir,
             overwrite = is.overwrite)
  }
  write_session_info(log_dir,sub_name = opt$sub_name )
  quit()
}


