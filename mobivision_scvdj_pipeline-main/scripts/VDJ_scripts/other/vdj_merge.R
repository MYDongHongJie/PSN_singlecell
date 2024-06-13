# module load  OESingleCell/3.0.a
# the command line api of data analysis of biological experiments
if (!requireNamespace('optparse', quietly = TRUE)) {
  stop("Please install argparse to parse the command line paramters!")
}
suppressPackageStartupMessages( library("dplyr") )
suppressPackageStartupMessages( library("tidyr") )
suppressPackageStartupMessages( library("optparse")  )

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i"), type = "character", default = "result/mobivision",
                 help = "the input directory: mobivision result directory"),
    make_option( c("--metadata","-m"), type ="character",default = "config/samples.csv",
                help="the samples' metadata."),
    make_option( c("--output","-o"),type="character", default = NULL,
                help="the output directory of results.", metavar="character"),
    make_option( c("--type", "-t"), type="character", default = "TCR",
                 help = "vdj type, TCR or BCR")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================

if ( is.null(opt$type) ){
    print("NO vdj type specified,TCR will be used!")
    type = "TCR"
}else{
    type = opt$type
}

if ( type == "TCR" ){
    type_path = "TCR"
}else if( type == "BCR"){
    type_path = "BCR"
}else{
    stop("NO right vdj type found, please check vdj type whether TCR or BCR")
}

if ( !is.null(opt$metadata) ){
  assay_metadata = read.csv(opt$metadata,sep=",",header =T )
  rownames(assay_metadata) = paste0(assay_metadata$sampleid,"_",type)
}else{
  stop("NO sample metadata found!")
}


if ( !is.null(opt$input) ){
  vdj_path = normalizePath(sub("\\/$","",opt$input,perl=T))
}else{
  stop("NO input dir found!")
}
output_dir=opt$output
output_dir=file.path(opt$output,"Clonotypes")
# setting the output directory

  
if ( !file.exists(output_dir) ){
    dir.create(output_dir,recursive = T)
}

output_dir = normalizePath(output_dir )


###################### excute ####################
if (all(file.exists(Sys.glob(file.path(vdj_path,assay_metadata$sampleid))) )){
  # # glob all the clonotype annotation file for all the samples
  clonotypes_tbs = file.path(vdj_path,paste0(assay_metadata$sampleid,"_",type_path),"outs",paste0(assay_metadata$sampleid,"_",type_path,"_clonotypes.csv"))
  names(clonotypes_tbs) = paste0(assay_metadata$sampleid,"_",type)
  filtered_contig_tbs = file.path(vdj_path,paste0(assay_metadata$sampleid,"_",type_path),"outs",paste0(assay_metadata$sampleid,"_",type_path,"_filtered_contig_annotations.csv"))
  names(filtered_contig_tbs) = paste0(assay_metadata$sampleid,"_",type)
  vdj_gene_location = file.path(vdj_path,paste0(assay_metadata$sampleid,"_",type_path),"outs",paste0(assay_metadata$sampleid,"_",type_path,"_airr_rearrangement.tsv"))
  names(vdj_gene_location) = paste0(assay_metadata$sampleid,"_",type)
}else{
  stop("NO Cellranger VDJ output directory available! Please make sure the samples' ID in the metadata are the same as
      the output directories' name of the Cellranger VDJ results.")
}

  ### clonotypes.csv, all_contig_annotations.json

# loess(yCol ~ {2:xCol}, data = data, subset = subset, weights = weights, span = span,
# degree = degree, na.action = {8:na.action})
merged_profile = future.apply::future_lapply( names(clonotypes_tbs), function(idx){
    # parse the all_contig_annotations.json file
    clonotypes_df = data.table::fread( clonotypes_tbs[idx]) %>% tibble::as_tibble()
    filtered_contig_df = data.table::fread( filtered_contig_tbs[idx]) %>% tibble::as_tibble() 
    location_df = data.table::fread( vdj_gene_location[idx]) %>% tibble::as_tibble() %>% dplyr::select( c("sequence_id","v_sequence_start", "v_sequence_end","d_sequence_start","d_sequence_end","j_sequence_start","j_sequence_end") ) %>% dplyr::rename( contig_id = sequence_id)

    filtered_contig_df[which(filtered_contig_df$d_gene %in% "") ,"d_gene" ] = NA 
    productive_summ = filtered_contig_df %>% dplyr::rename( "clonotype_id" = "raw_clonotype_id"  )  %>% 
    dplyr::inner_join( clonotypes_df, by = "clonotype_id")  %>%
    dplyr::select( -c("proportion", "frequency") ) %>% dplyr::inner_join( location_df, by = "contig_id")  %>%  tidyr::unite( "v_segment", c("v_gene","v_sequence_start", "v_sequence_end"), sep = ":") %>% tidyr::unite( "d_segment", c("d_gene","d_sequence_start", "d_sequence_end"), sep = ":") %>% tidyr::unite( "j_segment", c("j_gene","j_sequence_start", "j_sequence_end"), sep = ":") %>% tidyr::unite( "segments", c("v_segment","d_segment", "j_segment"), sep = ",") %>% 

    dplyr::mutate( sampleid = idx ) %>% 
    dplyr::mutate( barcode = paste(sampleid, gsub("-\\d+","",barcode), sep = "-"))
    productive_summ
  })


names(merged_profile )=names(clonotypes_tbs)

  # merge all the annotation from each to one data.frame
  merged_profile = as.data.frame(do.call(rbind, merged_profile))

# unification of two immunological sequences' annotation into one line in the same cell if exits.
  integrated_profile =  merged_profile %>%
    dplyr::mutate( chainx = chain,
            read_count = paste0( chain, ":", reads),
            umi_count = paste0( chain, ":" , umis) ) %>%
    dplyr::select( barcode, sampleid,
    contig_id, segments, chain,cdr3s_aa, cdr3s_nt,read_count, reads, umi_count,umis, clonotype_id ) %>%
    dplyr::group_by( barcode ) %>%
    dplyr::summarize( # merge the paired VDJ annotation for each chain into one line
      sampleid = unique(sampleid),
      cdr3s_aa = unique(cdr3s_aa),
      cdr3s_nt = unique(cdr3s_nt),
      # clonotype_name = unique(clonotype_id),
      #cdr3s_region = paste(unique(cdr3s_region), collapse=";"),
      segments= paste(segments, collapse=";"),
      #CDS_region = paste(unique(CDS_region), collapse=";"),
      vdj_read_count = paste(unique(read_count), collapse=";"),
      vdj_umi_count = paste(unique(umi_count), collapse=";"),
      umis4cell = sum(unique(umis)), # have bugs when double chain have same umi counts
      reads4cell = sum(unique(reads)) ) %>% dplyr::ungroup()


clonotype_count = integrated_profile %>% dplyr::group_by(cdr3s_aa) %>% dplyr::tally() %>% dplyr::arrange(desc(n))
  clonotype_count$clonotype_id =  paste("clonotype",1:dim(clonotype_count)[1],sep="")
  integrated_profile = dplyr::inner_join(integrated_profile, clonotype_count,by = "cdr3s_aa") %>%
    dplyr::arrange( desc(n)) %>% dplyr::mutate(is_paired = ifelse(grepl(";", cdr3s_aa), "TRUE", "FALSE" )) %>%
    dplyr::distinct()


  write.table( as.data.frame(integrated_profile), file.path(output_dir, "merged_vdj_contig_annotation.xls"), sep="\t", quote = F, col.names = T, row.names=F)


`%!in%` <- Negate(`%in%`) #定义%!in%
additional_cell_meta = vector()
  sampleidx =  gsub("(_|-)[ATCG]{16,}.*", "", integrated_profile$barcode,perl=T) #the index order is the same as the row index of the assay metadata
  for ( colidx in colnames(assay_metadata) ){
    if ( colidx %!in% c("sampleid","group") ) {
        additional_cell_meta = cbind(additional_cell_meta,
            as.vector( assay_metadata[sampleidx, colidx]) )
    }else{
        additional_cell_meta = cbind(additional_cell_meta,
            paste0(as.vector(assay_metadata[sampleidx, colidx]),"_",type) )
    }
  }
  colnames(additional_cell_meta) = colnames(assay_metadata)
  rownames(additional_cell_meta) = integrated_profile$barcode
  additional_cell_meta = additional_cell_meta %>% as.data.frame() %>%
                        tibble::rownames_to_column(var = "barcode")
  integrated_profile = dplyr::inner_join(additional_cell_meta, integrated_profile,
                                  by = c("barcode" = "barcode", "sampleid" = "sampleid"), keep = F)



  # count the umi, reads or cells for each group
  opt$quant.use="umi"
  #for ( methodx in unlist( strsplit(opt$quant.use,",") ) ){
    methodx="umi"
    countby = switch( tolower(methodx),
      "umi" = "umis4cell",
      "reads" = "reads4cell",
      "cells" = "id4cell"
    )


    clonotype_table = integrated_profile %>%
      dplyr::select_at( c("sampleid", "clonotype_id", countby)) %>%
      dplyr::rename("counts" = {{countby}}) %>%
      dplyr::group_by(sampleid,clonotype_id) %>%
      dplyr::summarise(counts = sum(counts)) %>%
      dplyr::ungroup() %>%
      tidyr::spread(sampleid, counts,fill = 0) %>%
      dplyr::mutate(total = rowSums(.[,2:dim(.)[2]])) %>%
      arrange(desc(total)) %>% dplyr::select(-total)

    write.table(clonotype_table, file.path(output_dir, glue::glue("clonotype_count_by_{methodx}_for_each_sample.xls")),
                 sep="\t", quote = F, col.names = T, row.names=F)



  opt$vdjformat="vdjtools"
  # change the format clonotype table to the specified
  future.apply::future_lapply(names(clonotypes_tbs), function(idx){
    if ( opt$vdjformat == "vdjtools" ){
      parsed_vdj = merged_profile %>%
        dplyr::filter( sampleid == idx ) %>%
        dplyr::mutate( frequency = umis/sum(umis)) %>%
        tidyr::separate(segments,
                 c("V", "V.start", "V.end", "D", "D.start", "D.end", "J", "J.start", "J.end"), sep= ":|,") %>%
        dplyr::select( umis,frequency,cdr3_nt,cdr3,V,D,J, V.end, D.start, D.end, J.start, barcode,cdr3s_aa, cdr3s_nt) %>%
        dplyr::rename( count = umis, CDR3nt = cdr3_nt, CDR3aa = cdr3,
                Vend = V.end, Dstart = D.start, Dend = D.end, Jstart = J.start,
                clonotype_nt = cdr3s_nt, clonotype_aa = cdr3s_aa ) %>%
        dplyr::arrange(dplyr::desc(count))
      parsed_vdj[is.na( parsed_vdj) ]= "."
    }else if ( opt$vdjformat == "immunarch" ){
      parsed_vdj = merged_profile %>%
        dplyr::filter( sampleid == idx ) %>%
        dplyr::mutate( frequency = umi_count/sum(umi_count)) %>%
        tidyr::separate(segments,
                 c("V", "V.start", "V.end", "D", "D.start", "D.end", "J", "J.start", "J.end"), sep= ":|,") %>%
        dplyr::select( umi_count,frequency,cdr3_nt,cdr3_aa,V,D,J, V.end, D.start, D.end, J.start, contig_id, barcode, cdr3s_nt,cdr3s_aa) %>%
        dplyr::rename( count = umi_count, CDR3nt = cdr3_nt, CDR3aa = cdr3_aa, Vend = V.end, Dstart = D.start, Dend = D.end, Jstart = J.start,clonotype_nt = cdr3s_nt, clonotype_aa = cdr3s_aa  ) %>%
        dplyr::arrange(desc(count))
      parsed_vdj[is.na( parsed_vdj) ]= "."
    }
    write.table( parsed_vdj, file.path(output_dir, paste0( idx, ".xls", collapse = "")),
                 sep="\t", col.names = T, row.names=F, quote = F)
  })

  file_name = paste( rownames(assay_metadata), ".xls", sep = "")
  new_names = c("#file_name", colnames(assay_metadata) )
  assay_metadata = cbind( file_name, assay_metadata )
  colnames(assay_metadata) = new_names
  assay_metadata = assay_metadata %>% rename( sample.id = sampleid)
  for ( i in assay_metadata$sample.id ){
    assay_metadata[which(assay_metadata$sample.id %in% i ), "sample.id"] = paste0(i, "_", type)
  }
  for ( i in assay_metadata$group ){
    assay_metadata[which(assay_metadata$group %in% i ), "group"] = paste0(i, "_", type)
  }

  write.table( assay_metadata, file.path(output_dir, "vdj_metadata.xls"),
               sep="\t", col.names = T, row.names=F, quote = F)


########## 