
sub_diffexp = subparsers$add_parser("diffexp", help = "differential test between specified groups.")
sub_diffexp$add_argument("-x", "--design", type = "character", default = NULL,
          help = "The group design for cell clusters or samples to make differential expression analysis.[default: %(default)s]")
sub_diffexp$add_argument("-M", "--addition_metadata", type="character", default = NULL,
          help="[Optional]additional metadata for each sample which includes sample id and additional sample groupping info.[default: %(default)s]")
sub_diffexp$add_argument("-c","--contrast", type = "character",default = NULL,
          help = "[Optional]levels of a factor used to compare with for final differenetial results. The format is Factor:interesting_levle:reference_level.[default: %(default)s]")
sub_diffexp$add_argument("-k", "--FC", type = "double", default = 1,
          help = "The average FC of the gene UMI count in its cluster against the all other clusters.[default: %(default)s]")
sub_diffexp$add_argument("-p", "--pvalue", type = "double", default = 0.05,
          help = "the P-value of the gene differential expression.[default: %(default)s]",metavar = "P-value")
sub_diffexp$add_argument("-q","--fdr", type = "double", default = NULL,
          help = "the fdr of the gene differential expression.[default: %(default)s]",metavar = "fdr")
sub_diffexp$add_argument("--test","-e", type = "character", default = "wilcox",
          help = "the test methods used to find differential expressed genes between specified group levels.Options are:wilcox,presto,venice,limma,roc,t,bimod,poission negbinom,scde,MAST,DESeq2,DESeq2LRT,edgeR.[default: %(default)s] ")
# sub_diffexp$add_argument( "--species",type = "character", help = "e.g. Human, Mouse or Other.")
sub_diffexp$add_argument( "--anno", type = "character", default = NULL,
       help = "[OPTIONAL]Annotation file submitted when species is Other.[default: %(default)s]")
sub_diffexp$add_argument(
    "--sub_obj",
    default = NULL,
    help = "if the inputted rds had been subclustered or subsetted, it should be given by this parameter.[default: %(default)s]"
)
sub_diffexp$add_argument(
    "--recorrect_umi",
    default = TRUE,
    help = "whether to recorrect umi when run FindMarkers. If the object had run PrepSCTFindMarkers and subset, it should be 'FALSE'.[default: %(default)s]"
)

# =============== Subcmd: diffexp, find the differential expressed genes ============
args <- commandArgs(TRUE)
if ( "diffexp"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name == "diffexp" ){
      #parse the design
      if ( is.null( opt$design ) ){
          warning("NO assay design is PROVIDED\n")
      }else{
      design <- opt$design
      }
      if ( is.null(opt$pvalue) && is.null(opt$fdr) ){
      stop("None of P-value or fdr is AVAILABLE! Filtering can be performed using any one of (-p), (-f) at one time.", call.=FALSE)
      }else if ( !is.null(opt$fdr)){
          fdr <- as.numeric(opt$fdr)
          pvalue <- NULL
      }else{
          pvalue <- as.numeric(opt$pvalue)
          fdr <- NULL
      }

 # read the specified assay and data slot in data object into memory
      data_ob <- OESingleCell::ReadX(input = opt$input,
                                      informat = opt$informat,
                                      assays = assays,
                                      data.use = dataslots,
                                      verbose = F)
      if (assays == "SCT") {
          if ((length(data_ob@assays$SCT@SCTModel.list) > 1) & is.null(data_ob@misc$PreSCT)) {
              ## https://rdrr.io/cran/Seurat/man/PrepSCTFindMarkers.html
              futile.logger::flog.info("Prepare object to run differential expression on SCT assay ")
              data_ob <- Seurat::PrepSCTFindMarkers(data_ob, assay = "SCT", verbose = TRUE)
              data_ob@misc$PreSCT <- TRUE
          }
      }
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
      if (!is.null(opt$sub_obj)) {
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
          data_ob@meta.data <- data_sub@meta.data
      }
      message("Loading Data object Finished!")

      #check the metadata of the assay design,which describe the expriement groupping of each
      #sample in the assay
      if (!is.null(opt$addition_metadata)) { #the additional metadata for each sample
          add_assay_metadata <- read.table(opt$addition_metadata, sep = ",", header = T)
          # cellnames = data_ob@cell.names
          cellnames <- colnames(data_ob) #seurat v3 style
          sampleidx <- gsub("_.*", "", cellnames, perl = T) #the index order is the same as the row index of the assay metadata
          #integrate the additional metadata from the assay design
          additional_cell_meta <- vector()
          for (colidx in colnames(add_assay_metadata)) {
              additional_cell_meta <- cbind(additional_cell_meta, as.vector(add_assay_metadata[sampleidx, colidx]))
          }
          colnames(additional_cell_meta) <- colnames(assay_metadata)
          rownames(additional_cell_meta) <- cellnames
          additional_cell_meta <- as.data.frame(additional_cell_meta)
          data_ob <- AddMetaData(data_ob, additional_cell_meta)

          # assay_metadata = data_ob@meta.data
          assay_metadata <- OESingleCell::colData(data_ob)
      }else {
          assay_metadata <- OESingleCell::colData(data_ob)
      }

      # subset cells on the expression matrix using logical exprssion on features
      if (!is.null(opt$predicate)) {
          desired_cells <- subset(assay_metadata, eval(parse(text = opt$predicate)))
          data_ob <- data_ob[, rownames(desired_cells)]
      }

      #parse the contrast for this assay
      #if the constrast is not specified by the user explicitly from the command line,
      #the final differential result will use last level of the last variable in the
      #design formula over the first level of this variable. The levels for each
      #factor is aphabetly ordered by default,which can be reordered by the user.
      #if available,the contrast string from the user must obey the right format:
      #the_interested_factor:the_interested_levels_in_this_factor:the_reference_level_in_this_factor
      if (is.null(opt$contrast)) { #no contrast is provided
          factors_indesign <- strsplit(opt$design, "[~+ ]+", perl = T)
          last_factor_indesign <- factors_indesign[length(factors_indesign)]
          if (is.null(assay_metadata[, last_factor_indesign])) {
              stop("The factor in design formula does not exist in assay metadata.")
          }
          variable_levels <- levels(assay_metadata[, last_factor_indesign])
          contrast <- paste(last_factor_indesign, variable_levels[length(variable_levels)], variable_levels[1], sep = ":")
      }else {
          contrast <- opt$contrast
      }
      contrasts <- unlist(strsplit(contrast, ":", perl = T))
      all_levels <- as.vector(unique(assay_metadata[, contrasts[1]]))
      if (contrasts[2] == "all" & contrasts[3] != "all") {
          all_levels <- all_levels[-which(all_levels == contrasts[3])] #delete the reference level
          all_comparisions <- paste(contrasts[1], all_levels, contrasts[3], sep = ":")
      }else if (contrasts[2] == "all" & contrasts[3] == "all") {

          # combine_of2 = combn(all_levels,2) #random combination of two group
          # all_comparisions = c( paste(contrasts[1],combine_of2[2,],combine_of2[1,],sep = ":"),
          #                       paste(contrasts[1],combine_of2[1,],combine_of2[2,],sep = ":"))
          all_comparisions <- lapply(all_levels,
                                     function(x) paste(contrasts[1], x, paste0(all_levels[-which(all_levels == x)], collapse = ","), sep = ":"))
          all_comparisions <- unlist(all_comparisions)
      }else if (contrasts[2] != "all" & contrasts[3] == "all") {
          #delete the interested level in the  reference level
          ref_levels <- paste0(all_levels[-which(all_levels == contrasts[2])], collapse = ",")
          all_comparisions <- paste(contrasts[1], contrasts[2], ref_levels, sep = ":")
      }else {
          all_comparisions <- contrast
      }

      anno <- read.delim(opt$anno, stringsAsFactors = F, sep = '\t', quote = "")
      future.apply::future_lapply(all_comparisions, function(contrast) {
          OESingleCell::RunDiffexp(object = data_ob,
                                   anno = anno,
                                   test = opt$test, fdr = fdr,
                                   fc.thres = opt$FC,
                                   recorrect_umi = opt$recorrect_umi,
                                   pval.thres = pvalue,
                                   contrast = contrast,
                                   outputdir = output_dir)
      }, future.seed = 2020)


      write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
      quit()
  }
}