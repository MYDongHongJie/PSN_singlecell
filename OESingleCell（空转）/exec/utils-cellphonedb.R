sub_cellphonedb <- subparsers$add_parser("cellphonedb",
                                         # description = docstring,
                                         formatter_class = 'argparse.RawTextHelpFormatter',
                                         #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                         argument_default = "True",
                                         help = "cellphonedb analysis ")
# sub_cellphonedb$add_argument("-l", "--celltype", type = "character", default = NULL,
#                              help = paste0("the comma seperated list for desired cell types,to be used with the parameter",
#                                            "--column4cell.If NULL, all cell type will be used.[default: %(default)s]"))
sub_cellphonedb$add_argument("-s", "--species", type = "character", default = "human",
                             help = "the species for the count matrix:human, mouse.[default: %(default)s]")
# sub_cellphonedb$add_argument("--deg",type="character", default=NULL,
#                          help = paste0("using Differentially Expressed Genes (cellphonedb method degs_analysis) as an ",
#                           "alternative to the permutation-based approach (cellphonedb method statistical_analysis),input deg files. [default: %(default)s]"))
# sub_cellphonedb$add_argument("--microenvs",type="character", default=NULL,
#                          help = "[options] for spatial transcriptome, special celltype environment[default: %(default)s]")
sub_cellphonedb$add_argument("-c", "--column4cell", type = "character", default = "celltype",
                             help = "the cell type annotation column in the cell annotation meta data.[default: %(default)s]")
sub_cellphonedb$add_argument("-m", "--cellmeta", type = "character", default = NULL,
                             help = "A table with the cell type annotation for each cell barcode.Only two columns: Cell,cell_type needed.[default: %(default)s].")
sub_cellphonedb$add_argument("-p", "--pvalue", type = "double", default = 0.05,
                             help = "the P value threshold for significant interacting pairs.[default: %(default)s]")
sub_cellphonedb$add_argument("-r", "--threshold", type = "double", default = 0.1,
                             help = "minimium percentage of cells expressing the specific ligand/receptor.[default: %(default)s]")
sub_cellphonedb$add_argument("--iterations", "-n", type = "integer", default = 1000,
                             help = "Number of iterations for the statistical analysis, use 1000 as default.[default: %(default)s]")
sub_cellphonedb$add_argument("--customedb", "-u" ,type="character", default="/data/database/cellphonedb/cellphone.db",
                         help = "the user customized interaction annotation database used to run cellphonedb.If NULL(as default), the cellphoneDB will be used.[default: %(default)s]")
sub_cellphonedb$add_argument("--genetype", type="character", default="gene_name",
                         help = "Type of gene identifiers in the counts data, options can be:ensembl,gene_name,hgnc_symbol.[default: %(default)s]")
args <- commandArgs(TRUE)
if ( "cellphonedb"  %in% args){
  opt<-intial_setting()
  if(opt$sub_name == "cellphonedb" ){
    #### load query object data===========================================================================================
    if ( opt$informat == "cellphonedb" ){
      suppressMessages(cellp <- OESingleCell::ParseCpdb(cpdb = opt$input,
                                                        species = opt$species,
                                                        pvalue = opt$pvalue,
                                                        avgexpr = NULL))
      write.table(cellp$ligrec,
                  file.path(output_dir, "cell_comm_annotation.xls"),
                  sep = "\t",
                  col.names = T,
                  row.names = F,
                  quote = F)
      saveRDS(cellp, file.path(output_dir, "cellphonedb_results.rds") )
      quit()
    }

    if ( opt$informat == "rds" ){
      futile.logger::flog.info("step1:load query object data")#===========================================================
      suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                      informat = opt$informat,
                                                      assays = assays[1],
                                                      data.use = dataslots,
                                                      verbose = F))
      Seurat::DefaultAssay(data_ob) <- assays[1]
      cellmeta <- OESingleCell::colData(data_ob)
      #===================================================================================================================
      futile.logger::flog.info(glue::glue("step2:get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
      if (!is.null(opt$predicate)) {
        df <- OESingleCell::colData(data_ob)
        desired_cells <- subset(df, eval(parse(text = opt$predicate)))
        data_ob <- subset(data_ob, cells = rownames(desired_cells))
        if (!is.null(Seurat::Images(data_ob))) {
          unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique)]
          if (length(unuse_images) > 0) { data_ob@images[unuse_images] <- NULL }
        }
      }
      #===================================================================================================================
      futile.logger::flog.info("step3:get the subset of cells used for visualization if necessay")
        #get the count matrix
        counts <- Seurat::GetAssayData(data_ob, slot = opt$dataslot)
        if ( is.null(opt$column4cell) ){
          metadata <- Seurat::FetchData(data_ob, vars = "ident") %>%
                      tibble::rownames_to_column(var = "Cell") %>%
                      dplyr::rename("cell_type" = ident )
        }else{
          metadata <- Seurat::FetchData(data_ob, vars = opt$column4cell ) %>%
                      tibble::rownames_to_column(var = "Cell")
          colnames(metadata) <- c("Cell", "cell_type")
        }
        rm(data_ob)
        gc()
    }else if ( opt$informat == "xsv" ){
        counts <- vroom::vroom(file.path(opt$input), col_names = T, comment = "#")
        counts <- as.data.frame(counts)
        rownames(counts) <- counts[, 1]
        counts <- counts[, -1]
        if ( is.null( opt$cellmeta) | opt$cellmeta == "null" ){
          stop("Cell celltype annotation is Mendatory for user provided count matrix!")
        }else{
          metadata <- read.table(opt$cellmeta,
                                 sep = "\t",
                                 header = T,
                                 stringsAsFactors = F)
          colnames(metadata) <- c("Cell", "cell_type")
        }
        counts <- counts[, metadata$Cell]
        if ( !is.null(opt$celltype)){
            if ( is.null(opt$column4cell )  || opt$column4cell == "null"){
              column4cell <- "cell_type"
            }else{
              column4cell <- opt$column4cell
            }
            desired_cell_types <- unlist(strsplit(opt$celltype, ",", perl = T))
            desired_cells <- metadata %>%
                             dplyr::filter(!!ensym(column4cell) %in% desired_cell_types ) %>%
                             pull( Cell )
            counts <- counts[, desired_cells]
        }
    }

    if ( opt$species == "mouse" ){
      mouse2human <- homologene::human2mouse(annotables::grch38$symbol) %>%
                     dplyr::select( mouseGene, humanGene )
      genes2use <- intersect(names(Seurat::CaseMatch(rownames(counts), mouse2human$mouseGene)),
                             rownames(counts) )
      counts <- counts[genes2use,]
      rownames(counts) <- plyr::mapvalues(rownames(counts),
                                          from = mouse2human$mouseGene,
                                          to = mouse2human$humanGene,
                                          warn_missing = F)
    }

    tempdir <- tempdir()
    # tempdir <- output_dir
    # using vroom to speed up the count matrix writing
    counts_out <- tibble::rownames_to_column( as.data.frame(counts), var = 'Gene')
    vroom::vroom_write(x = counts_out,
                       path = file.path(tempdir,'counts.tsv'),
                       col_names = T,
                       quote = "none",
                       delim = "\t")
    write.table(metadata,
                file = file.path(tempdir,'metadata.tsv'),
                quote = F,
                col.names = T,
                row.names = F,
                sep = '\t')

    outdir <- file.path(output_dir, "out")
    # run the cellphonedb with the statistical method
    #cellphonedb_bin <- system("which cellphonedb", intern = T)
    #cellphonedb_bin <- '/public/dev_scRNA/software/conda_envs/OESingleCell_fangying/bin/cellphonedb'
    cellphonedb_bin <- "/public/dev_scRNA/software/conda_envs/cellphonedb/bin/cellphonedb"
    # run the cellphonedb with the statistical method
    #cellphonedb_bin <- system("which cellphonedb", intern = T)
    if (!file.exists(cellphonedb_bin)) {
      stop("NO cellphonedb command detected! Please make sure cellphonedb path correctly setted!")
    }

    if (!file.exists(file.path(output_dir, "cellphone.db"))) {
      file.copy(opt$customedb, file.path(output_dir, "cellphone.db"))
    }
    futile.logger::flog.info("step4:run cellphonedb:") #================================================
    # if(is.null(opt$deg)){
    cmd <- glue::glue('module purge && ',
                      '{cellphonedb_bin}  method statistical_analysis  {tempdir}/metadata.tsv  {tempdir}/counts.tsv ',
                      ' --counts-data  {opt$genetype}',
                      '  --iterations={opt$iterations}',
                      ' --threshold={opt$threshold}',
                      ' --threads={opt$ncores}  ',
                      ' --pvalue={opt$pvalue} ',
                      ' --output-path={outdir} ',
                      ' --database {opt$customedb} --verbose',
                      ' && module load /home/weihao/modulefiles/OESingleCell/v_3.0.0_visium_produce_advance'
    )
    #   else{
    #   cmd <- glue::glue('{cellphonedb_bin}  method degs_analysis  {tempdir}/metadata.tsv  {tempdir}/counts.tsv  {opt$deg} ',
    #                     ' --counts-data  {opt$genetype}   ',
    #                     ' --threshold={opt$threshold} ',
    #                     ' --threads={opt$ncores}  ',
    #                     ' --pvalue={opt$pvalue} ',
    #                     ' --output-path={outdir} ',
    #                     ' --database {opt$customedb} --verbose')
    # }
    futile.logger::flog.info(cmd)
  system(cmd)
  count_long <- as.data.frame(counts) %>%
    tibble::rownames_to_column(var = "gene") %>%
    tidyr::gather("Cell", "expr", 2:dim(.)[2])
  count_long$cell_type <- plyr::mapvalues(count_long[["Cell"]],
                                          from = metadata$Cell,
                                          to = as.character(metadata$cell_type),
                                          warn_missing = F)
  avgexpr <- count_long %>%
               dplyr::select(gene, cell_type, expr ) %>%
               dplyr::mutate(gene = toupper(gene)) %>%
               dplyr::group_by( cell_type, gene ) %>%
               dplyr::mutate(exprs = log2(mean(expr)+0.0001) ) %>%
               dplyr::select( gene, cell_type, exprs ) %>%
               dplyr::distinct()
    # parse the output of cellphonedb using human genes
    cellphonedb <- OESingleCell::ParseCpdb( outdir ,
                                            species = opt$species,
                                            pvalue = opt$pvalue,
                                            avgexpr = avgexpr  )
    write.table(cellphonedb$ligrec,
                file.path(output_dir, "cell_comm_annotation.xls"),
                sep = "\t",
                col.names = T,
                row.names = F, quote = F)
    saveRDS(cellphonedb, file.path(output_dir, "cellphonedb_results.rds"))
    ## save session information
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
}