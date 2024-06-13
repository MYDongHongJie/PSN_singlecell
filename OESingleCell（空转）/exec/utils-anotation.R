sub_annotation <- subparsers$add_parser("annotation", help = "add annotation for genelist.")
sub_annotation$add_argument("-g","--genelist", type = "character", default = NULL, help = "The gene list file ")
sub_annotation$add_argument("-a","--anno", type = "character", default = NULL, help = "The gene annotation file: gene_annotation.xls")
# =============== Subcmd: annotation, add annotation to a genelist file =========
args <- commandArgs(TRUE)
if ( "annotation"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name == "annotation" ){
    if (is.null(opt$anno)) {
      stop("Please provide the species information with parameter \"-a\" !")
    }else{
      marker = read.delim(opt$genelist, sep = "\t", stringsAsFactors = F)
      anno = read.delim(opt$anno,stringsAsFactors = F, sep = '\t')
      if (colnames(marker[1]) == "gene") {
        res = dplyr::left_join(marker, anno, by = c("gene" = "id"))
      }
      else {
        res = dplyr::left_join(marker, anno, by = c("GeneID" = "id"))
      }
      res[is.na(res)] = "--"
      write.table(res,
                  paste0(strsplit(opt$genelist, "\\.xls")[[1]][1], "_anno.xls", collapse = ""),
                  sep = '\t',
                  quote = F,
                  row.names = F)
    }
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
}