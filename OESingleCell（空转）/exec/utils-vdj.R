docstring <- " example:\\n\\n\\
  sctool  -i  query.rds   -f rds  -o results  -d rds   --assay Spatial  --dataslot counts  st_deconv  --refexp ref.rds  --refcelltype  celltype --refassay  RNA "
sub_vdj <- subparsers$add_parser("vdj",
                                 description = docstring,
                                 # formatter_class= 'argparse.RawTextHelpFormatter' ,
                                 formatter_class = "lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                 argument_default = "True",
                                 help = "parse the output of cellranger vdj to integrate.")
sub_vdj$add_argument("--vdj", type = "character", default = "./",
                     help = "the directory of 10X cellranger vdj results directory with each sample as subdirectory or the directory with file filtered_contig_annotations.csv,clonotypes.csv, all_contig_annotations.json,all_contig.fasta.[default: %(default)s] ")
sub_vdj$add_argument("--metadata", "-m", type = "character", default = NULL,
                     help = "the sample metadata which must include sample id in this assay design.[default: %(default)s] ")
sub_vdj$add_argument("--type", "-t", type = "character", default = NULL,
                     help = "the type of immunological cell receptor: TCR or BCR.[default: %(default)s] ")
sub_vdj$add_argument("--seurat", "-s", type = "character", default = NULL,
                     help = "add the clonotype annotation to the seurat as an assay.[default: %(default)s] ")
sub_vdj$add_argument("--quant.use", "-q", type = "character", default = "umi",
                     help = "count number for each clonotype in each sample, choices can be: cells: the number of cell for each clonotype in each sample; umi: the UMI number for each clonotype in each sample; reads：the reads for each clonotype in each sample. [default: %(default)s] ")
sub_vdj$add_argument("--maketse", "-e", type = "character", default = "FALSE",
                     help = "whether to construct the TreeSummarizedExperiment object to store data for later diversity analysis.[default: %(default)s] ")
sub_vdj$add_argument("--bcprefix", "-x", type = "character", default = "sampleid",
                     help = "prefix the cell barcode with sample id or row index in metadata.tsv. Options can be sampleid or index.[default: %(default)s]")
sub_vdj$add_argument("--vdjformat", "-v", type = "character", default = "vdjtools",
                     help = "the output format of the V(D)J clonotypes table. Options can be vdjtools or immunarch.[default: %(default)s] ")

# =============== Subcmd: vdj, parsing the vdj output from cellranger =========
args <- commandArgs(TRUE)
if ("vdj" %in% args) {
  opt <- intial_setting()
  suppressPackageStartupMessages( library("dplyr") )
  suppressPackageStartupMessages( library("tidyr") )
  # suppressPackageStartupMessages( library("dplyr") )
  # suppressPackageStartupMessages( library("tidyr") )
  if (!is.null(opt$metadata)) {
    assay_metadata <- read.csv(opt$metadata, sep = ",", header = T)
    rownames(assay_metadata) <- assay_metadata$sampleid
  }else {
    stop("NO sample metadata found!")
  }

  vdj_path <- normalizePath(sub("\\/$", "", opt$vdj, perl = T))
  if (all(file.exists(Sys.glob(file.path(vdj_path, assay_metadata$sampleid))))) {
    # # glob all the clonotype annotation file for all the samples
    if ("subtype" %in% colnames(assay_metadata)) {
      clonotypes_tbs <- file.path(vdj_path,
                                  assay_metadata$sampleid,
                                  "outs",
                                  "per_sample_outs",
                                  assay_metadata$sampleid,
                                  assay_metadata$subtype,
                                  "clonotypes.csv")
      json4samples <- file.path(vdj_path,
                                assay_metadata$sampleid,
                                "outs",
                                "multi",
                                assay_metadata$subtype,
                                "all_contig_annotations.json")
    }else {
      clonotypes_tbs <- file.path(vdj_path, assay_metadata$sampleid, "outs", "clonotypes.csv")
      json4samples <- file.path(vdj_path, assay_metadata$sampleid, "outs", "all_contig_annotations.json")

    }
    names(clonotypes_tbs) <- assay_metadata$sampleid
    names(json4samples) <- assay_metadata$sampleid
    print(clonotypes_tbs)
    print(json4samples)
  }else {
    stop("NO Cellranger VDJ output directory available! Please make sure the samples' ID in the metadata are the same as
        the output directories' name of the Cellranger VDJ results.")
  }
  # join the clonotype.csv and filtered_contig_annotations.csv table using the clonotype id
  # concatenate the table to relabel the clonotype id
  merged_profile <- future.apply::future_lapply(names(json4samples), function(idx) {
    # parse the all_contig_annotations.json file
    clonotypes_df <- data.table::fread(clonotypes_tbs[idx]) %>% tibble::as_tibble() # read in the clonotypes.csv
    productive_summ <- jsonlite::fromJSON(json4samples[idx]) %>% # all contig annotation nested data.frame
                        dplyr::filter(is_cell != "FALSE" &
                                        high_confidence != "FALSE" &
                                        !is.null(cdr3) &
                                        productive == TRUE) %>%
                        dplyr::mutate(raw_clonotype_id = info$raw_clonotype_id) %>%
                        # dplyr::select(-c("info", "filtered", "primer_annotations", "frame", "clonotype", "quals", "is_cell", "productive", "high_confidence")) %>%
                        dplyr::select(-c("info", "filtered", "frame", "clonotype", "quals", "is_cell", "productive", "high_confidence")) %>%
                        tidyr::drop_na(raw_clonotype_id) %>%
                        tidyr::drop_na(cdr3)
    anno_list <- future.apply::future_lapply(productive_summ$annotations, function(anno) {
                 anno <-  anno %>%
                          dplyr::select(-c("annotation_length", "annotation_match_end",
                                           # "annotation_match_start","cigar","mismatches", "score" )) %>%
                                           "annotation_match_start", "cigar", "score")) %>%
                          dplyr::mutate(region_type = feature$region_type,  # extract the immune allel gene segment annotation
                                        region = feature$gene_name, chain = feature$chain) %>%
                          dplyr::filter(!region_type %in% c("5'UTR", "C-REGION")) %>%
                          dplyr::select(-feature)

                if (!"D-REGION" %in% unique(anno[["region_type"]])) {
                  D_region <- c(
                    contig_match_end = "NA",
                    contig_match_start = "NA",
                    region_type = "D-REGION",
                    region = "NA",
                    chain = unique(anno[["chain"]])[1]
                  )
                  anno <- rbind(anno, D_region[colnames(anno)])
                }
                anno_row <- anno %>%
                  dplyr::mutate(region_type = gsub("L?-REGION\\+?", "", region_type, perl = T)) %>%
                  dplyr::arrange(factor(region_type, levels = c("V", "D", "J"))) %>%
                  tidyr::unite("range", c("contig_match_start", "contig_match_end"), sep = ":") %>%
                  tidyr::unite("segments", c("region", "range"), sep = ":") %>%
                  dplyr::summarize(segments = paste(segments, collapse = ","), chain = unique(chain)[1])
    })
    anno_df <- do.call(rbind, anno_list)
    productive_summ <- cbind(productive_summ, anno_df) %>%
                       dplyr::select(-annotations) %>%
                       dplyr::rename("clonotype_id" = "raw_clonotype_id",
                                    "contig_id" = "contig_name", "cdr3_aa" = "cdr3", "cdr3_nt" = "cdr3_seq") %>%
                       dplyr::inner_join(clonotypes_df, by = "clonotype_id") %>%
                       dplyr::select(-c("aa_sequence", "proportion", "frequency")) %>%
                       dplyr::mutate(sampleid = idx) %>%
                       dplyr::mutate(barcode = paste(sampleid, gsub("-\\d+", "", barcode), sep = "-"))
                       productive_summ
  })

  # merge all the annotation from each to one data.frame
  merged_profile <- as.data.frame(do.call(rbind, merged_profile))

  # unification of two immunological sequences' annotation into one line in the same cell if exits.
  integrated_profile <- merged_profile %>%
                        dplyr::mutate(umis = umi_count,
                                      reads = read_count,
                                      chainx = chain,
                                      read_count = paste0(chain, ":", read_count),
                                      umi_count = paste0(chain, ":", umi_count),
                                      cdr3s_region = paste0(chain, ":", cdr3_start, "-", cdr3_stop),
                                      CDS_region = paste0(chain, ":", start_codon_pos, "-", stop_codon_pos)) %>%
                        dplyr::select(barcode, sampleid, contig_id, chain, cdr3s_aa, cdr3s_nt, cdr3s_region, chainx,
                                      segments, CDS_region, read_count, reads, umi_count, umis) %>%
                        dplyr::group_by(barcode) %>%
                        dplyr::summarize( # merge the paired VDJ annotation for each chain into one line
                                        sampleid = unique(sampleid),
                                        cdr3s_aa = unique(cdr3s_aa),
                                        cdr3s_nt = unique(cdr3s_nt),
                                        cdr3s_region = paste(unique(cdr3s_region), collapse = ";"),
                                        segments = paste(segments, collapse = ";"),
                                        CDS_region = paste(unique(CDS_region), collapse = ";"),
                                        vdj_read_count = paste(unique(read_count), collapse = ";"),
                                        vdj_umi_count = paste(unique(umi_count), collapse = ";"),
                                        umis4cell = sum(unique(umis)), # have bugs when double chain have same umi counts
                                        reads4cell = sum(unique(reads))) %>%
                        dplyr::ungroup()
  #
  clonotype_count <-  integrated_profile %>%
                      dplyr::group_by(cdr3s_aa) %>%
                      dplyr::tally() %>%
                      dplyr::arrange(desc(n))
  clonotype_count$clonotype_id <- paste("clonotype", 1:dim(clonotype_count)[1], sep = "")
  integrated_profile <- dplyr::inner_join(integrated_profile, clonotype_count, by = "cdr3s_aa") %>%
                        dplyr::arrange(desc(n)) %>%
                        dplyr::mutate(is_paired = ifelse(grepl(";", cdr3s_aa), "TRUE", "FALSE")) %>%
                        dplyr::distinct() %>%
                        dplyr::mutate(id4cell = 1)
  write.table(integrated_profile, file.path(output_dir, "merged_vdj_contig_annotation.csv"),
              sep = "\t", quote = F, col.names = T, row.names = F)


  additional_cell_meta <- vector()
  sampleidx <- gsub("(_|-)[ATCG]{16,}.*", "", integrated_profile$barcode, perl = T) #the index order is the same as the row index of the assay metadata
  for (colidx in colnames(assay_metadata)) {
    additional_cell_meta <- cbind(additional_cell_meta,
                                  as.vector(assay_metadata[sampleidx, colidx]))
  }
  colnames(additional_cell_meta) <- colnames(assay_metadata)
  rownames(additional_cell_meta) <- integrated_profile$barcode
  additional_cell_meta <- additional_cell_meta %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "barcode")
  integrated_profile <- dplyr::inner_join(additional_cell_meta, integrated_profile,
                                          by = c("barcode" = "barcode", "sampleid" = "sampleid"), keep = F)

  if (!is.null(opt$seurat)) {
    data_ob <- OESingleCell::ReadX(input = opt$input,
                                   informat = opt$informat,
                                   assays = assays, # only RNA assay is valid
                                  data.use = "counts", # counts slot is enough
                                  verbose = F)
    seurat_ob <- OESingleCell::AddMetaData(data_ob, metadata = integrated_profile)
    SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob)
    # saveRDSMC(seurat_ob, opt$seurat)
  }
  # count the umi, reads or cells for each group
  for (methodx in unlist(strsplit(opt$quant.use, ","))) {
    countby <- switch(tolower(methodx),
                      "umi" = "umis4cell",
                      "reads" = "reads4cell",
                      "cells" = "id4cell"
    )

    clonotype_table <- integrated_profile %>%
      dplyr::select_at(c("sampleid", "clonotype_id", countby)) %>%
      dplyr::rename("counts" = { { countby } }) %>%
      dplyr::group_by(sampleid, clonotype_id) %>%
      dplyr::summarise(counts = sum(counts)) %>%
      dplyr::ungroup() %>%
      tidyr::spread(sampleid, counts, fill = 0) %>%
      dplyr::mutate(total = rowSums(.[, 2:dim(.)[2]])) %>%
      dplyr::arrange(desc(total)) %>%
      dplyr::select(-total)

    write.table(clonotype_table, file.path(output_dir, glue::glue("clonotype_count_by_{methodx}_for_each_sample.xls")),
                sep = "\t", quote = F, col.names = T, row.names = F)

    if (as.logical(opt$maketse)) {
      # make the TreeSummarizedExperiment as a data container of the sample level
      # clonotype annotation for diversity ananlysis
      tse_counts <- clonotype_table %>%
        tibble::column_to_rownames(var = "clonotype_id") %>%
        as.matrix()
      clonotype_anno <- integrated_profile %>%
        dplyr::select(clonotype_id, cdr3s_aa) %>%
        unique() %>%
        dplyr::as_tibble() %>%
        tibble::column_to_rownames(var = "clonotype_id")
      tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(counts = tse_counts),
        colData = assay_metadata[colnames(tse_counts),],
        rowData = clonotype_anno[rownames(tse_counts),])
      saveRDS(tse, file.path(output_dir, glue::glue("vdj.tse.by.{methodx}.rds")))
    }
  }

  # change the format clonotype table to the specified
  future.apply::future_lapply(names(json4samples), function(idx) {
    if (opt$vdjformat == "vdjtools") {
      parsed_vdj <- merged_profile %>%
        dplyr::filter(sampleid == idx) %>%
        dplyr::mutate(frequency = umi_count / sum(umi_count)) %>%
        tidyr::separate(segments,
                        c("V", "V.start", "V.end", "D", "D.start", "D.end", "J", "J.start", "J.end"), sep = ":|,") %>%
        dplyr::select(umi_count, frequency, cdr3_nt, cdr3_aa, V, D, J, V.end, D.start, D.end, J.start, barcode, cdr3s_aa, cdr3s_nt) %>%
        dplyr::rename(count = umi_count, CDR3nt = cdr3_nt, CDR3aa = cdr3_aa,
                      Vend = V.end, Dstart = D.start, Dend = D.end, Jstart = J.start,
                      clonotype_nt = cdr3s_nt, clonotype_aa = cdr3s_aa) %>%
        dplyr::arrange(dplyr::desc(count))
      parsed_vdj[is.na(parsed_vdj)] <- "."
    }else if (opt$vdjformat == "immunarch") {
      parsed_vdj <- merged_profile %>%
        dplyr::filter(sampleid == idx) %>%
        dplyr::mutate(frequency = umi_count / sum(umi_count)) %>%
        tidyr::separate(segments,
                        c("V", "V.start", "V.end", "D", "D.start", "D.end", "J", "J.start", "J.end"), sep = ":|,") %>%
        dplyr::select(umi_count, frequency, cdr3_nt, cdr3_aa, V, D, J, V.end, D.start, D.end, J.start, contig_id, barcode, cdr3s_nt, cdr3s_aa) %>%
        dplyr::rename(count = umi_count, CDR3nt = cdr3_nt, CDR3aa = cdr3_aa, Vend = V.end, Dstart = D.start, Dend = D.end, Jstart = J.start, clonotype_nt = cdr3s_nt, clonotype_aa = cdr3s_aa) %>%
        dplyr::arrange(desc(count))
      parsed_vdj[is.na(parsed_vdj)] <- "."
    }
    write.table(parsed_vdj, file.path(output_dir, paste0(idx, ".xls", collapse = "")),
                sep = "\t", col.names = T, row.names = F, quote = F)
  })

  file_name <- paste(rownames(assay_metadata), ".xls", sep = "")
  new_names <- c("#file_name", colnames(assay_metadata))
  assay_metadata <- cbind(file_name, assay_metadata)
  colnames(assay_metadata) <- new_names
  assay_metadata <- assay_metadata %>% rename(sample.id = sampleid)
  write.table(assay_metadata, file.path(output_dir, "vdj_metadata.xls"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
  quit()
}