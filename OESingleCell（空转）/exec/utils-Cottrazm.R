#=======================================================================================================================
docstring <- " example1:\\n\\n\\
sctool  -i /public/dev_scRNA/weihao/project/demo/demo_human_prostate_2022_08_10/result/spaceranger -o results -j 30 Cottrazm  --sampleid human_prostate_acinar  --resolution 0.4 "
sub_Cottrazm <- subparsers$add_parser(
  "Cottrazm",
  description = docstring,
  formatter_class = "argparse.RawTextHelpFormatter",
  # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
  argument_default = "True",
  help = "Using Cottrazm to define tumor boundary in ST data"
)
sub_Cottrazm$add_argument(
  "--sampleid",
  type = "character",
  default = NULL,
  help = "the sample id to run Cottrazm.[default %(default)s]"
)
sub_Cottrazm$add_argument(
  "--metadata",
  type = "character",
  default = NULL,
  help = "the sample metadata which should include sample id in this assay design (samples.csv).[default %(default)s]"
)
sub_Cottrazm$add_argument(
  "--resolution",
  type = "double",
  default = 0.4,
  help = "vaule used to set the resolution of cluster distiguish, use a value above(below)1.0 if you want to obtain a larger(smaller) number of communities.[default: %(default)s]"
)

#=======================================================================================================================
args <- commandArgs(TRUE)
if ("Cottrazm" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "Cottrazm") {
      if(!is.null(opt$sampleid)) {
          sampleids <- opt$sampleid
      }else if (!is.null(opt$metadata)) {
          sample_file <- read.csv(opt$metadata, header = TRUE, sep = "\t")
          sampleids <- sample_file$sampleid
      }else {
          stop("Please provide sampleid or metadata which contained sampleid!")
      }
      for (sample in sampleids) {
          futile.logger::flog.info("step1:ST PreProcess......") #=======================================================
          input_dir <- paste0(opt$input, "/", sample, "/outs/")
          output_dir <- paste0(opt$output, "/", sample, "/")
          TumorST <- Cottrazm::STPreProcess(InDir = input_dir, OutDir = output_dir, Sample = sample)

          futile.logger::flog.info("step2:Stlearn Clustering......") #==================================================
          if (!file.exists(glue::glue("{output_dir}/TumorST_{sample}.rds"))) {
              TumorST <- Cottrazm::STModiCluster(InDir = input_dir, Sample = sample, OutDir = output_dir, TumorST = TumorST, res = opt$resolution)
              saveRDS(TumorST, glue::glue("{output_dir}/TumorST_{sample}.rds"))
          }else {
              TumorST <- readRDS(glue::glue("{output_dir}/TumorST_{sample}.rds"))
          }

          futile.logger::flog.info("step3:ST inferCNV......") #=========================================================
          if (!file.exists(glue::glue("{output_dir}/InferCNV/output_Spatial/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings"))) {
              STInferCNV <- Cottrazm::STCNV(TumorST = TumorST, OutDir = output_dir, assay = "Spatial", num_threads = opt$ncores)
          }

          futile.logger::flog.info("step4:STCNVScore......") #==========================================================
          TumorST <- Cottrazm::STCNVScore(TumorST = TumorST, assay = "Spatial", Sample = sample, OutDir = output_dir)

          futile.logger::flog.info("step5:Define Boundary......") #=====================================================
          TumorSTn <- Cottrazm::BoundaryDefine(TumorST = TumorST, MalLabel = c(7,8), OutDir = output_dir, Sample = sample)

          futile.logger::flog.info("step6:BoundaryPlot......") #========================================================
          TumorST_LUC_1 <- Cottrazm::BoundaryPlot(TumorSTn = TumorSTn, TumorST = TumorST, OutDir = output_dir, Sample = sample)
      }
      write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
      quit()
  }
}
