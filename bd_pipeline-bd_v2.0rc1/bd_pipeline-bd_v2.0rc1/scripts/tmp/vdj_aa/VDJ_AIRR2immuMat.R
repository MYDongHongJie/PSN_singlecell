library(R.utils)
library("optparse")
library("dplyr")

#=command line parameters setting=============================
option_list = list(
    make_option( c( "--inputdir","-i"), type = "character",default = "result/BD_Analysis/",
        help = "the input dir."),
    make_option( c("--output","-o"),type="character", default = "result/VDJ_aggr",
        help="the output directory of Clustering results." ),
    make_option( c("--metadata","-m"),type="character", default = "config/samples.csv",
        help="the sample metadata which must include sample id in this assay design.")
    #make_option( c("--distance","-d"),type="character", default = "euclidean",
    #    help="the way to calculate distance ,eg. euclidean,cosine, rankcor, l2"),
    #make_option( c("--root","-r"),type="character", default = NULL,
    #    help="the root clusters,,eg. clusters:3 , new_celltype:NKT "),
    ##make_option( c("--extraGene", "-x"), type = "character", default = NULL,
    ##    help = "The gene list to visualize specified by the user."),
    ##make_option( c("--n_pcs"), type = "double",  default = NULL,
    ##    help = "[Otional] get fast when calculate cells >500  ."),
    #make_option( c("--colorby", "-c" ), type = "character", default = "clusters,sampleid,group",
    #    help = "[Otional]visualize cells' metadata by coloring cells in different color according to cell grouping list."),
    #make_option( c("--groupby", "-g" ), type = "character", default = NULL,
    #    help = "[Otional]the column name in cell grouping metadata used with which_cell parameter for subseting cell groups."),
    #make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
    #    help = "The subset of cluster ids used for subtyping.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

metadata = read.delim(opt$metadata,sep=",",header=T)

#FUNCTION FOR BD-VDJ FILE
format_BDVDJ <- function(AIRR_df,out_type){
  VDJ_Dominant_AIRR <- read.delim(gzfile(AIRR_df),sep="\t",comment.char = '#')
  header_list = c("duplicate_count" ,"cdr3","cdr3_aa","v_call","d_call" ,"j_call" ,
                  "v_sequence_end" ,"d_sequence_start" ,"d_sequence_end" ,"j_sequence_start" ,"cell_id","junction_aa","junction")
  new_header = c("count","CDR3nt","CDR3aa","V","D","J","Vend","Dstart",
                 "Dend","Jstart","barcode","clonotype_aa","clonotype_nt","frequency")
  sort_header = c("count","frequency","CDR3nt","CDR3aa","V","D","J","Vend","Dstart",
                 "Dend","Jstart","barcode","clonotype_aa","clonotype_nt")
  VDJ_Dominant_AIRR$locus = factor(VDJ_Dominant_AIRR$locus,levels = sort(unique(VDJ_Dominant_AIRR$locus)))
  vdj_type = unique(VDJ_Dominant_AIRR$locus)[grep(x = unique(VDJ_Dominant_AIRR$locus),pattern = paste0("^",out_type))]
  vdj_mat = VDJ_Dominant_AIRR[(VDJ_Dominant_AIRR$locus %in% vdj_type),]

  vdj_df = vdj_mat[header_list]
  vdj_df$frequency <- vdj_df$duplicate_count/sum(vdj_df$duplicate_count)
  names(vdj_df) = new_header
  vdj_df = vdj_df[sort_header]
  vdj_df1 = vdj_df[which(vdj_df$CDR3nt !=""),]
  vdj_df1$CDR3nt <- factor(vdj_df1$CDR3nt,levels = sort(unique(vdj_df1$CDR3nt)))
  nucl_N = length(grep("N", vdj_df1$CDR3nt, ignore.case = FALSE))
  if(nucl_N==0){
      vdj_df2 = vdj_df1
  } else if(nucl_N  >0 & nucl_N <20){
      vdj_df2 = vdj_df1[-grep("N", vdj_df1$CDR3nt, ignore.case = FALSE),]
      print("some N nucleobase in your VDJ files,will remove it")
  } else if(nucl_N >20){
      return(message(("too many N nucleobase in your VDJ files,please check it")))
  }
  return(vdj_df2)
}

VDJtype = data.frame("TR","IG")
names(VDJtype)  = c("TCR","BCR")


#1.format files from AIRR file
for (sample in metadata$sampleid){
	VDJ_Dominant_AIRR = file.path(opt$inputdir,sample,paste0(sample,"_VDJ_Dominant_Contigs_AIRR.tsv.gz"))
	for (i in names(VDJtype)){
		outdir = file.path(opt$output,i,"Clonotypes")
		if ( !dir.exists(outdir) ){
			dir.create(outdir,recursive=T)
		}
		out_type = c(matrix(VDJtype[i][,1]))
		out_file = format_BDVDJ(VDJ_Dominant_AIRR,out_type)
		write.table(out_file,file.path(outdir,paste0(sample,"_",i,".xls")),sep="\t",quote=F,row.names=F,col.names=T)
		}
}
#2.get vdj_metadata
for (i in names(VDJtype)){
	outdir = file.path(opt$output,i,"Clonotypes")
	file_name = c(list.files(outdir,pattern = "*CR.xls"))
	new_names = c("#file_name", colnames(metadata)[1:3] )
	assay_metadata = cbind( file_name, metadata[1:3] )
	colnames(assay_metadata) = new_names
	assay_metadata = assay_metadata %>% rename( sample.id = sampleid)
	write.table( assay_metadata, file.path(outdir, "vdj_metadata.xls"),
				sep="\t", col.names = T, row.names=F, quote = F)
}
