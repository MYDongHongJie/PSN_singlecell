library("optparse")
library("dplyr")

#=command line parameters setting=============================
option_list = list(
    make_option( c( "--inputdir","-i"), type = "character",default = "result/BD_Analysis/",
        help = "the input dir."),
    make_option( c("--output","-o"),type="character", default = "result/VDJ_aggr",
        help="the output directory of Clustering results." ),
    make_option( c("--type","-t"),type="character", default = "TCR,BCR",
        help="the VDJ type." ),
    make_option( c("--metadata","-m"),type="character", default = "config/samples.csv",
        help="the sample metadata which must include sample id in this assay design.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

unique_term <- function(df,term){
    df_result = apply(df[term], 1, function(x) {
            unique(strsplit(as.character(x[1]), "; ")[[1]]) %>%paste(collapse = ",")})
    return(df_result)
}

VDJtype = data.frame("TR","IG")
names(VDJtype)  = c("TCR","BCR")
mytype = opt$type
mytype = unlist(strsplit(mytype, ",", perl =T))

VDJtype = VDJtype[mytype]

assay_metadata = read.delim(opt$metadata,sep=",",header=T)
vdj_path = normalizePath(sub("\\/$","",opt$output,perl=T))
if (all(file.exists(Sys.glob(file.path(vdj_path,assay_metadata$sampleid))) )){
    for (i in names(VDJtype)){
      outdir = file.path(opt$output,i,"Clonotypes")
      clonotypes_tbs = file.path(vdj_path,i,"Clonotypes",paste0(assay_metadata$sampleid,"_",i,".xls"))
      names(clonotypes_tbs) = paste0(assay_metadata$sampleid,"_",i)
      vdj_gene_location = paste0(opt$inputdir,"/",assay_metadata$sampleid,"/",assay_metadata$sampleid,"_VDJ_Dominant_Contigs_AIRR.tsv")
      names(vdj_gene_location) = paste0(assay_metadata$sampleid,"_",i)
        merged_profile = future.apply::future_lapply( names(clonotypes_tbs), function(idx){
        # parse the all_contig_annotations.json file
        clonotypes_df = data.table::fread( clonotypes_tbs[idx]) %>% tibble::as_tibble()
        all_info_vdj = data.table::fread( vdj_gene_location[idx]) %>% tibble::as_tibble()
        names(clonotypes_df) = c("duplicate_count" ,"frequency","cdr3","cdr3_aa","v_call","d_call" ,"j_call" ,
                  "v_sequence_end" ,"d_sequence_start" ,"d_sequence_end" ,"j_sequence_start" ,"cell_id","junction_aa","junction")
        filter_df = merge(clonotypes_df,all_info_vdj,all.x=T,by=c("duplicate_count","cdr3","cdr3_aa","v_call","d_call" ,"j_call" ,
		          "v_sequence_end" ,"d_sequence_start" ,"d_sequence_end" ,"j_sequence_start" ,"cell_id","junction_aa","junction"))
        filter_df_select <- filter_df %>% replace(filter_df =="", "NA") %>% mutate(sampleid = idx)%>% mutate(barcode = paste0(idx,"_",cell_id))%>% 
                                    mutate(cdr3s_aa = paste0(locus,":",cdr3_aa)) %>% mutate(cdr3s_nt = paste0(locus,":",cdr3)) %>% 
                                    mutate(cdr3s_region = paste0(locus,":",cdr3_aa)) %>%  tidyr::unite( "v_segment", c("v_call","v_sequence_start", "v_sequence_end"), sep = ":") %>% tidyr::unite( "d_segment", c("d_call","d_sequence_start", "d_sequence_end"), sep = ":") %>% tidyr::unite( "j_segment", c("j_call","j_sequence_start", "j_sequence_end"), sep = ":") %>% tidyr::unite( "segments", c("v_segment","d_segment", "j_segment"), sep = ",") %>% mutate( chainx = locus,read_count = paste0( locus, ":", consensus_count),umi_count = paste0( locus, ":" , duplicate_count) ) 
        productive_summ = filter_df_select %>% dplyr::select( barcode, sampleid, segments, locus,cdr3s_aa, cdr3s_nt,read_count, consensus_count, umi_count,duplicate_count)
        })
      names(merged_profile )=names(clonotypes_tbs)
      merged_profile_all_sample = as.data.frame(do.call(rbind, merged_profile))
      #按照 barcode 和 aa 进行整合
	  merge_result <- merged_profile_all_sample %>%group_by(barcode) %>%
                              summarize_all(paste, collapse = "; ")
      clonotype_count <-  merge_result%>% 
							  dplyr::group_by(cdr3s_aa) %>% dplyr::tally() %>% dplyr::arrange(desc(n))
      counts_sum = merged_profile_all_sample[c("barcode","consensus_count","duplicate_count")] %>% 
	                          group_by(barcode) %>%  summarize_all(sum, na.rm = TRUE) %>% 
							  rename(umis4cell = (duplicate_count),reads4cell = (consensus_count))
      clonotype_count$clonotype_id =  paste("clonotype",1:dim(clonotype_count)[1],sep="")
      integrated_profile1 = dplyr::inner_join(merge_result, clonotype_count,by = "cdr3s_aa") %>%
          dplyr::arrange( desc(n)) %>% dplyr::mutate(is_paired = ifelse(grepl(";", cdr3s_aa), "TRUE", "FALSE" )) %>%
          dplyr::distinct()
	  integrated_profile = dplyr::inner_join(integrated_profile1, counts_sum,by = "barcode") 
	  names(integrated_profile)[4] = "chain"
	  # 对每行字符串进行拆分，并对拆分结果进行去重
      integrated_profile['sampleid'] = unique_term(integrated_profile,"sampleid")
      write.table( as.data.frame(integrated_profile), file.path(outdir, "merged_vdj_contig_annotation.xls"), sep="\t", quote = F, col.names = T, row.names=F)
	  #get counts per clonotype
	   methodx="umi"
       countby = switch( tolower(methodx),
            "umi" = "umis4cell",
            "reads" = "reads4cell",
            "cells" = "id4cell" )
      clonotype_table = integrated_profile %>%
            dplyr::select_at( c("sampleid", "clonotype_id", countby)) %>%
            dplyr::rename("counts" = {{countby}}) %>%
            dplyr::group_by(sampleid,clonotype_id) %>%
            dplyr::summarise(counts = sum(counts)) %>%
            dplyr::ungroup() %>%
            tidyr::spread(sampleid, counts,fill = 0) %>%
            dplyr::mutate(total = rowSums(.[,2:dim(.)[2]])) %>%
            arrange(desc(total)) %>% dplyr::select(-total)
      write.table(clonotype_table, file.path(outdir, glue::glue("clonotype_count_by_{methodx}_for_each_sample.xls")),
                 sep="\t", quote = F, col.names = T, row.names=F)
    }
}else{
    stop("NO BD VDJ output directory available! Please make sure the samples' ID in the metadata are the same as
        the output directories' name of the Cellranger VDJ results.")
}

#merge per sample vdj info
for (i in names(VDJtype)){
    outdir = file.path(opt$output,i,"Clonotypes")
    files <- list.files(path = vdj_path, pattern = "*_VDJ_Summary.csv", full.names = TRUE)
    result <- data.frame()
    for (file_path in files) {
      df <- read.csv(file_path)
      result <- rbind(result, df)
    }
    sub_df = result[result$Chain_Category %in% i,]
    write.table(sub_df, file.path(outdir, glue::glue("Metrics_{i}_for_each_sample.xls")),
                 sep="\t", quote = F, col.names = T, row.names=F)
}
