docstring<- " example1(without refmarker):\\n\\n\\
  sctool  -i  query.rds   -f rds  -o results  -d rds   --assay Spatial  --dataslot counts  st_deconv  --refexp ref.rds  --refcelltype  celltype --refassay  RNA \\n\\n\\
example2(with refmarker):\\n\\n\\
  sctool  -i  query.rds   -f rds  -o results  -d rds   --assay Spatial  --dataslot counts  st_deconv  --refexp ref.rds  --refcelltype  celltype --refmarker --refmarker.xls  --refassay  RNA"

sub_st_deconv <- subparsers$add_parser("st_deconv",
                                       description = docstring,
                                       formatter_class= 'argparse.RawTextHelpFormatter' ,
                                       #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                       argument_default = "True" ,
                                       help="Using Spotlight to deconvention Spaital transcriptome cell type")
sub_st_deconv$add_argument("--refexp",
                           type = "character",
                           help = "[REQUIRED]the customized reference expression matrix in seurat format, which may come from the quantification results of microarray, Bulk RNA-seq and scRNA-seq sequencing." )
sub_st_deconv$add_argument("--refcelltype",
                           type = "character",
                           help = "[REQUIRED]the cell type annotation column id in the reference expression matrix." )
sub_st_deconv$add_argument("--refmarker",
                           type = "character",
                           default = "auto",
                           help = "the marker gene annotation results for reference data,or auto detected by findallmarker function in seurat.[default: %(default)s] " )
sub_st_deconv$add_argument("--refassay",
                           type = "character",
                           default = "RNA",
                           help = "[REQUIRED]the default assay for reference data to use for the this run.[default: %(default)s] " )
sub_st_deconv$add_argument( "--min_pct1", type = "double", default = 0.1,
                            help ="the minimium ratio of cells expressing one specific gene in a cluster.[default: %(default)s]")
sub_st_deconv$add_argument( "--avg_log2FC", type = "double", default = 0.25,
                            help = "The average log2FC of the gene UMI count in its cluster against the all other clusters.[default: %(default)s]")
sub_st_deconv$add_argument( "--test", type = "character", default = "wilcox",
                            help = "test methods used to cell markers. Options are: wilcox,presto,venice, t, bimod,poisson, negbinom, MAST, DESeq2, DESeq2LRT,limma, edgeR.[default: %(default)s] ")
                            # 'wilcox' : Wilcoxon rank sum test (default).
                            # 'presto' : performace improved Wilcoxon rank sum test for big data.
                            # 'venice' : package from bioTuning with function VeniceAllMarkers imtating the FindAllMarkers in Seurat but much faster.")
                            # 't' : Student\'s t-test.
                            # 'bimod' : Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013).
                            # 'poisson' : Likelihood ratio test assuming an underlying poisson distribution. Use only for UMI-based datasets.
                            # 'negbinom' : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets.
                            # 'MAST' : GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015).
                            # 'DESeq2' : DE based on a model using the negative binomial distribution (Love et al, Genome Biology, 2014).
# ================ Subcmd: st_deconv, deconvolute the cell type composition for spatial transcriptomics or bulk transcriptomics using scRNA-seq data ========
args <- commandArgs(TRUE)
if ( "st_deconv"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name == "st_deconv"){
    futile.logger::flog.info("step1:read the specified assay and data slot in data object into memory") #===============
    suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays[1],
                                 data.use = dataslots,
                                 verbose = F))
    Seurat::DefaultAssay(data_ob) <- assays[1]
    cellmeta <- OESingleCell::colData(data_ob)

    #### set default ident
    if ("clusters" %in% colnames(cellmeta)) {
      data_ob <- SeuratObject::StashIdent(data_ob, save.name = "clusters")
    } else {
      data_ob <- SeuratObject::SetIdent(data_ob, value = "clusters")
    }

    futile.logger::flog.info("step2:load reference object data")#=======================================================
    suppressMessages(ref_ob <- OESingleCell::ReadX(input = opt$refexp,
                                  informat = "rds",
                                  assays = opt$refassay,
                                  data.use = "data",
                                  verbose = F))
    if (!is.null(opt$refexp)) {
      # the input is a seurat object which may contain more than one sample
      SeuratObject::DefaultAssay(ref_ob) <- opt$refassay
      ref_ob <- SeuratObject::SetIdent(ref_ob, value = opt$refcelltype)
    }
    futile.logger::flog.info("step3: FindAllMarkers for reference object data per celltype")#===========================
    if (opt$refmarker !="auto") {
      refmarkers <- data.table::fread(opt$refmarker, header = T)
    } else {
      print("NO marker genes table AVAILABLE for reference data!
          It will use the results from the functon FindAllMarkers using
          default parameters on the reference seurat object!")
      refmarkers <- Seurat::FindAllMarkers(ref_ob,
                                           assay = opt$refassay,
                                           slot = "data",
                                           test.use = opt$test,
                                           only.pos = T,
                                           logfc.threshold = opt$avg_log2FC,
                                           min.pct = opt$min_pct1 )
    }
    ##====================================================================================================================
    #' This function takes in the H coefficient matrix object from and NMF object and returns a matrix object with the topic profile for each cell type
    #'
    #' @param h Object of class matrix, coefficient matrix from NMF model.
    #' @param train_cell_clust Object of class vector with cluster of the cells used to train the model.
    #' @return This function returns a list where the first element is a matrix with the topic profiles of all possible combinations and the 2nd element is the cell composition of each spot.
    #' @export
    #' @examples
    #'
    topic_profile_per_cluster_nmf <- function(h, train_cell_clust) {
      # Check variables
      if (!is(h, "matrix")) stop("ERROR: h must be a matrix object!")
      if (! is(train_cell_clust, "vector")) stop("ERROR: train_cell_clust must be a vector/list object!")

      h_ds <- data.frame(t(h))
      h_ds[, "clust_vr"] <- train_cell_clust

      ct_topic_profiles <- h_ds %>%
                           dplyr::group_by(clust_vr) %>%
                           dplyr::summarise_all(list(median)) %>%
                           tibble::column_to_rownames("clust_vr") %>%
                           as.matrix()
      ct_topic_profiles_t <- t(ct_topic_profiles)
      colnames(ct_topic_profiles_t) <- gsub("[[:punct:]]|[[:blank:]]", "_", colnames(ct_topic_profiles_t))
      return(ct_topic_profiles_t)
    }
      futile.logger::flog.info("step4:deconvolution the spatial assay")#=====================================================================================
      set.seed(2020)
      spotlight_ls <- SPOTlight::spotlight_deconvolution(se_sc = ref_ob,
                                                         counts_spatial = OESingleCell::GetAssayData(data_ob,
                                                                                                     assay = assays[1],
                                                                                                     slot = dataslots),
                                                         clust_vr = opt$refcelltype,
                                                         cluster_markers = refmarkers,
                                                         cl_n = 100,
                                                         hvg = 3000,
                                                         ntop = NULL,
                                                         transf = "uv",
                                                         method = "nsNMF",
                                                         min_cont = 0.09,
                                                         assay = "RNA",
                                                         slot = "counts" )

  ##====================================================================================================================
  #Assess deconvolution Before even looking at the decomposed spots we can gain insight on how well the model performed
  # by looking at the topic profiles for the cell types.
  # The first thing we can do is look at how specific the topic profiles are for each cell typ
    nmf_mod <- spotlight_ls[[1]]
    decon_mtrx <- spotlight_ls[[2]]
    h <- NMF::coef(nmf_mod[[1]])
    rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
    topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(h = h, train_cell_clust = nmf_mod[[2]])
    p<- topic_profile_plts[[2]] +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                     axis.text = ggplot2::element_text(size = 12))
    OESingleCell::save_ggplots(glue::glue("{output_dir}/0.access_deconvolution/NMF_Topic_profies_by_cell_type"),
                               plot = p,
                               dpi = 1000,
                               limitsize = F)
    p<-topic_profile_plts[[1]] +
      theme(axis.text.x = element_text(angle = 90),
            axis.text = element_text(size = 12))
    OESingleCell::save_ggplots(glue::glue("{output_dir}/0.access_deconvolution/NMF_Topic_proportion_within_cell_types"),
                               plot = p,
                               dpi = 1000,
                               limitsize = F)
    basis_spotlight <- data.frame(NMF::basis(nmf_mod[[1]]))
    colnames(basis_spotlight) <- unique(stringr::str_wrap(nmf_mod[[2]], width = 30))
    readr::write_tsv(basis_spotlight%>%as.data.frame()%>% tibble::rownames_to_column("gene") ,
               glue::glue("{output_dir}/0.access_deconvolution/The_most_important_gene_for_each_topic.xls"))
    futile.logger::flog.info("step5:formart output dataframe and save into data_ob@misc$spotlight_results")#============
    decon_mtx <- spotlight_ls[[2]]
    decon_mtx <- decon_mtx[, colnames(decon_mtx) != "res_ss"]
    ##
    rownames(decon_mtx) <- colnames(data_ob)
    decon_mtx<- decon_mtx %>%
                data.frame()
    ##insert unknown celltype's spot
    # decon_mtx ["unknown"] <- 0
    # decon_mtx[which(rowSums(decon_mtx) == 0), "unknown"] <- 1
    ##remove celltypes without any spot infer
    decon_mtx<- decon_mtx %>%  tibble::rownames_to_column("barcodes")
    # summer <- function(x){ if(is.numeric(x)){ sum(x) > 0 } else { TRUE } }
    # decon_mtx<-decon_mtx[, sapply(decon_mtx,  summer)]
    ##save results to seurat object's mis part
    data_ob@misc$spotlight_results <- decon_mtx
    maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
    data_ob<-Seurat::AddMetaData(data_ob,
                                 metadata= data_ob@misc$spotlight_results %>%
                                   tibble::column_to_rownames("barcodes") %>%
                                   dplyr::mutate(
                                        top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
                                        top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
                                    )
                                )
    futile.logger::flog.info("step6: save seurat object with celltype infor in data_ob@misc$spotlight_results")#========
    if ( tolower(opt$outformat) == "h5seurat" & as.logical(opt$update) ){
      SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob )
    }else{
      OESingleCell::SaveX(data_ob,
                          output = opt$output,
                          update = FALSE,
                          outformat = opt$outformat,
                          prefix = opt$prefix)
    }
    ## save session information
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit() }
}