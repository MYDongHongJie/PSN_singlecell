sub_monocl3 = subparsers$add_parser("monocle3", help = "run monocle3 on the input data object")
sub_monocl3$add_argument("--reduct", type = "character", default = "umap",
                         help = "the reduction method used for trajectory.[default: %(default)s] ")

# ================ Subcmd: monocle3, run monocle3 on the input data object. ========
args <- commandArgs(TRUE)
if ( "monocle3"  %in% args ){
  opt<-intial_setting()
   # read the specified assay and data slot in data object into memory
  data_ob = OESingleCell::ReadX(input = opt$input,
                                informat = opt$informat,
                                assays = assays,
                                data.use = dataslots, # data slot is enough
                                verbose = F)

    # object@assays$ATAC@meta.features = object@assays$ATAC@meta.features[,1:2]
    if ( class(data_ob) == "Seurat" ){
      assayx = Seurat::GetAssay(data_ob, assay = assays[1])
      slot(assayx, "meta.features") = slot(assayx, "meta.features")[,1:2]
      data_ob[[assay]] = assayx
      input_cds = SeuratWrappers::as.cell_data_set(data_ob)
    }else{
      input_cds = data_ob
    }

  input_cds = monocle3::cluster_cells( input_cds, reduction_method = "UMAP")
  input_cds = monocle3::learn_graph( input_cds, use_partition = TRUE)
  input_cds = monocle3::order_cells( input_cds, reduction_method = opt$reduct, root_cells = opt$root_cell)
  if ( class(data_ob) == "Seurat" ){
    data_ob = Seurat::AddMetaData(data_ob,
                                  metadata = input_cds@principal_graph_aux@listData$UMAP$pseudotime,
                                  col.name = "pseudotime" )
    OESingleCell::SaveX(data_ob,
                        output = output_dir,
                        outformat = opt$outformat,
                        prefix = opt$prefix,
                        update = FALSE )
  }else{
    OESingleCell::SaveX(input_cds,
                        output = output_dir,
                        outformat = opt$outformat,
                        prefix = opt$prefix,
                        update = FALSE )
  }
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}