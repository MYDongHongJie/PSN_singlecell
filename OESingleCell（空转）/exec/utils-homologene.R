docstring <- "example1:\\n\\n\\
sctool -i in.h5seurat -f h5seurat -o test -d h5seurat homologene -t 7955 -T 9606 \\n\\n\\
example2: \\n\\n\\
"

#=command line parameters setting=============================
# option_list = list(
#     make_option( c("--RDS", "-v"), type = "character", default = "TRUE",
#         help = "the seurat object saved as R object in RDS format."),
#     make_option( c("--output","-o"),type="character", default = "./",
#         help="the output directory of results.", metavar="character"),
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);


#suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("homologene"))

sub_homo <- subparsers$add_parser("homologene", help = "")
sub_homo$add_argument("-b", "--blast", type = "character", default = NULL,
                help = "[OPTIONAL] gne2bd two colomn balst result. query 1st ref 2nd. ")
sub_homo$add_argument("-t", "--inTaxid", type = "character", default = NULL,
                help = "The taxonomy ID consistenting with input RDS corresponding species ")
sub_homo$add_argument("-T", "--outTaxid", type = "character", default = NULL,
                help = "The taxonomy ID of the species that needs to be converted into homologous genes ")
# ========= Subcmd: clustering for dimension reduction and clustering ==========
args <- commandArgs(TRUE)
if ( "homologene"  %in% args ){
    opt<-intial_setting()
  # read the specified assay and data slot in data object into memory
    data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 reductions = "pca",
                                 assays = assays, # only RNA assay is valid
                                 data.use = "counts",
                                 verbose = F)
    if ( is.null(opt$inTaxid)){
        stop("the inTax id is NOT AVAILABLE!")
    }else{
        inTaxid=opt$inTaxid
    }

    if ( is.null(opt$outTaxid)){
        stop("the outTax id is NOT AVAILABLE!")
    }else{
        outTaxid=opt$outTaxid
    }

    #counts=as.data.frame(data_ob@assays$RNA@counts)
    counts=data_ob@assays$RNA@counts
    expr=data_ob@assays$RNA@data
    gene=row.names(counts)
    ## blast 
    if ( is.null(opt$blast)) {
        # inTax is input species tax id , outTax is the transformed species tax id
        in2out = homologene(gene,inTax=inTaxid,outTax=outTaxid)
    } else {
        in2out = read.table(opt$blast,header=F,sep="\t",quote="")
    }
    ## unique
    index <- duplicated(in2out[,1])
    in2out <-in2out[!index,]
    index<-duplicated(in2out[,2])
    in2out <- in2out[!index,]
    homologene <- in2out[,1:2]
    write.table(homologene,file.path(output_dir,paste0("homologene.",inTaxid,"_2_",outTaxid,".xls",collapse=".")),sep="\t",row.names=F,quote=F)
    counts_filted = counts[as.vector(in2out[,1]),]
    rownames(counts_filted) = as.character(in2out[,2])
    data_ob@assays$RNA@counts=counts_filted
    expr_filtered = expr[as.vector(in2out[,1]),]
    rownames(expr_filtered) = as.character(in2out[,2])
    data_ob@assays$RNA@data=expr_filtered

    #saveRDS(data_ob,file.path(output_dir,paste0(inTaxid,"_2_",outTaxid,".",basename(opt$RDS))))
    OESingleCell::SaveX(data_ob,
                      output = output_dir,
                      outformat = opt$outformat,
                      prefix = paste0(inTaxid,"_2_",outTaxid,".",basename(opt$input)),
                      update = FALSE )

    write_session_info(output_dir,sub_name = parser$parse_args()$sub_name )
    quit()
}


