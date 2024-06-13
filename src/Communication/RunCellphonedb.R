#!/usr/bin/env Rscript
# this is the script used to infer the cell-cell communication
# using the high expressed genes or differential expressed genes
# in each cell. The key point is to find the Ligand-Receptor pairs
# in the selected genes.

##' parse the cellphonedb results
##'
##' @param cpdb the path of the cellphonedb results , containing pvalues.txt, means.txt etc.
##' @param species the species name, Only "human" and "mouse" are supported currenrly.
##' @param pvalue the P value threshlod for filtering
##' @param avgexpr the expressoin of each gene in each group
##'
##' @export
ParseCpdb <- function(
  cpdb,
  species = "human",
  pvalue = 0.05,
  avgexpr = NULL
){
  if ( species == "mouse" ){
    mouse2human <- homologene::human2mouse(annotables::grch38$symbol) %>% dplyr::select( mouseGene, humanGene )
  }
  raw_pval <- read.table(glue::glue("{cpdb}/pvalues.txt"), header=T,
                        stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  raw_pval$interacting_pair = apply(raw_pval,1, function(x){
                      if(grepl("^complex:", x[3], perl = T)[1]){
                        pt_a = gsub("^complex:", "", x[3])
                        x[2]= gsub(glue::glue("{pt_a}_"), glue::glue("{pt_a}|"), x[2], ignore.case = T)
                      }else if ( grepl("^complex:", x[4], perl = T)[1] ){
                        pt_b = gsub("^complex:", "", x[4])
                        x[2] = gsub(glue::glue("_{pt_b}"), glue::glue("|{pt_b}"), x[2])
                      } else if ( grepl("^simple:",x[3], perl = T) && grepl("^simple:", x[4],perl =T)){
                        x[2] = gsub("^([a-zA-Z0-9-]+)_([a-zA-Z0-9-]+)$", "\\1|\\2", x[2])
                      }
    })

  raw_means <- read.table(glue::glue("{cpdb}/means.txt"), header=T,
                         stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  raw_means$interacting_pair = apply(raw_means,1, function(x){
                      if(grepl("^complex:", x[3], perl = T)[1]){
                        pt_a = gsub("^complex:", "", x[3])
                        x[2]= gsub(glue::glue("{pt_a}_"), glue::glue("{pt_a}|"), x[2])
                      }else if ( grepl("^complex:", x[4], perl = T)[1] ){
                        pt_b = gsub("^complex:", "", x[4])
                        x[2] = gsub(glue::glue("_{pt_b}"), glue::glue("|{pt_b}"), x[2])
                      } else if ( grepl("^simple:",x[3], perl = T) && grepl("^simple:", x[4],perl =T)){
                        x[2] = gsub("^([a-zA-Z0-9-]+)_([a-zA-Z0-9-]+)$", "\\1|\\2", x[2])
                      }
    })
  raw_sig_means <- read.table(glue::glue("{cpdb}/significant_means.txt"),
                         header=T, stringsAsFactors = F, sep="\t",
                         comment.char = '', check.names=F)
  raw_sig_means$interacting_pair = apply(raw_sig_means,1, function(x){
                      if(grepl("^complex:", x[3], perl = T)[1]){
                        pt_a = gsub("^complex:", "", x[3])
                        x[2]= gsub(glue::glue("{pt_a}_"), glue::glue("{pt_a}|"), x[2])
                      }else if ( grepl("^complex:", x[4], perl = T)[1] ){
                        pt_b = gsub("^complex:", "", x[4])
                        x[2] = gsub(glue::glue("_{pt_b}"), glue::glue("|{pt_b}"), x[2])
                      } else if ( grepl("^simple:",x[3], perl = T) && grepl("^simple:", x[4],perl =T)){
                        x[2] = gsub("^([a-zA-Z0-9-]+)_([a-zA-Z0-9-]+)$", "\\1|\\2", x[2])
                      }
    })
  deconvoluted <- read.table(glue::glue("{cpdb}/deconvoluted.txt"),
                            header=T, stringsAsFactors = F,
                            sep="\t", comment.char = '', check.names=F)
  if ( species == "mouse" ){
  # all the calculation is based on the human gene, so when it come to mouse,
  # the genes in results should be changed to mouse's before plotting.
      all_interact_pair = unlist(lapply(strsplit(raw_pval$interacting_pair, "\\|", perl = T),
                                    function(m){
                                      for( i in seq_along(m) ){
                                        if ( m[i] %in% mouse2human$humanGene ){
                                          m[i] = mouse2human[which( mouse2human$humanGene == m[i] )[1],"mouseGene"]
                                        }else{
                                          mx = unlist(strsplit(m[i], "_| ", perl = T))
                                          for ( j in seq_along(mx) ){
                                            if ( mx[j] %in% mouse2human$humanGene ){
                                              mx[j] = mouse2human[which( mouse2human$humanGene == mx[j] )[1],"mouseGene"]
                                            }else{
                                              mx[j] = paste0("_", mx[j], collapse = "")
                                            }
                                          }
                                          m[i] = paste0("_", mx, collapse = "")
                                          m[i] = gsub("_+(complex|receptor|\\d)", " \\1", gsub("_+([0-9a-zA-Z]+)",  "_\\1", gsub("^_+", "", m[i])))
                                        }
                                      }
                                      m = paste(m, collapse = "|")
                                }
                            ))
      # pvalues.txt and means.txt have the same interacting gene pairs
      raw_pval$interacting_pair <- all_interact_pair
      raw_means$interacting_pair <- all_interact_pair
      raw_sig_means$interacting_pair <- unlist(lapply(strsplit(raw_sig_means$interacting_pair, "\\|", perl = T),
                                    function(m){
                                      for( i in seq_along(m) ){
                                        if ( m[i] %in% mouse2human$humanGene ){
                                          m[i] = mouse2human[which( mouse2human$humanGene == m[i] )[1],"mouseGene"]
                                        }else{
                                          mx = unlist(strsplit(m[i], "_| ", perl = T))
                                          for ( j in seq_along(mx) ){
                                            if ( mx[j] %in% mouse2human$humanGene ){
                                              mx[j] = mouse2human[which( mouse2human$humanGene == mx[j] )[1],"mouseGene"]
                                            }else{
                                              mx[j] = paste0("_", mx[j], collapse = "")
                                            }
                                          }
                                          m[i] = paste0("_", mx, collapse = "")
                                          m[i] = gsub("_+(complex|receptor|\\d)", " \\1", gsub("_+([0-9a-zA-Z]+)",  "_\\1", gsub("^_+", "", m[i])))
                                        }
                                      }
                                      m = paste(m, collapse = "|")
                                }
                            ))
  }

  desired_cols  <- colnames(raw_pval)
  sig_gene_pair <- raw_pval %>% dplyr::select(desired_cols) %>%
                  dplyr::select(2,12:dim(.)[2]) %>%
                  tidyr::gather( "cell_pair","pval", 2:dim(.)[2]) %>%
                  distinct() %>% filter( pval < pvalue ) %>% pull(interacting_pair)

  desired_means <- raw_means %>%
                dplyr::filter( interacting_pair %in% sig_gene_pair ) %>%
                dplyr::select( desired_cols)
  desired_pval <- raw_pval %>%
                dplyr::filter( interacting_pair %in% sig_gene_pair ) %>%
                dplyr::select( desired_cols )

  desired_sig_pvalx <- desired_pval %>%
                        dplyr::select( interacting_pair, receptor_a, receptor_b, secreted, is_integrin, 12:dim(.)[2] ) %>%
                        tidyr::gather( "cell_pair","pval", 6:dim(.)[2]) %>% distinct()
  desired_sig_meansx <- desired_means %>%
                        dplyr::select( interacting_pair, 12:dim(.)[2] ) %>%
                        tidyr::gather( "cell_pair","expr", 2:dim(.)[2]) %>% distinct()
  merged_df <- dplyr::full_join(desired_sig_meansx, desired_sig_pvalx,
            by =c("interacting_pair" = "interacting_pair", "cell_pair" = "cell_pair") ) %>%
            tidyr::separate( cell_pair, into = c("part_a_cell", "part_b_cell"), sep = "\\|") %>%
            tidyr::separate(interacting_pair, into = c("part_a_gene", "part_b_gene"), sep ="\\|")
   lr_df <- future.apply::future_lapply(1:nrow(merged_df), function(i){
        if ( merged_df[i,"receptor_a"]== "False" & merged_df[i,"receptor_b"] == "True" ){
             a <- merged_df[i,c(2,1,4,3,5,8,9,10)]
            names(a) <- c("receptor", "ligand", "receptor_cell", "ligand_cell", "expr","secreted","is_integrin","pval")
        }else{
            a <- merged_df[i,c(1,2,3,4,5,8,9,10)]
            names(a) <- c("receptor", "ligand", "receptor_cell", "ligand_cell", "expr","secreted","is_integrin","pval")
        }
       a
    })
    lr_df <- do.call(rbind, lr_df)
  if ( !is.null(avgexpr) ){
    avg <- avgexpr$exprs
    names(avg) <- paste(avgexpr$gene, avgexpr$cell_type, sep = "")
    lr_df$receptor_expr <- unlist(
        future.apply::future_lapply( seq_along(lr_df$receptor), function(x){
          receptors = toupper(unlist( strsplit(lr_df$receptor[x], "(_| )+", perl = T)))
          rgene = intersect( receptors , avgexpr$gene )
          if ( length(rgene) == 0 ){
            a = 0
          }else{
            a = sum(avg[paste(rgene, lr_df$receptor_cell[x], sep = "")])/length(rgene)
          }
        }
      )
    )
    lr_df$ligand_expr <- unlist(
        future.apply::future_lapply( seq_along(lr_df$ligand), function(x){
          ligands = toupper(unlist( strsplit(lr_df$ligand[x], "(_| )+", perl = T)))
          lgene = intersect( ligands, avgexpr$gene )
          if ( length(lgene) == 0 ){
            a = 0
          }else{
            a = sum(avg[paste(lgene, lr_df$ligand_cell[x], sep = "")])/length(lgene)
          }
        }
      )
    )
  }

  cellphonedb <- list(raw.pvalues=raw_pval, raw.means = raw_means,
                      ligrec = lr_df,desired_means = desired_means, desired_pval = desired_pval,
                     raw.sigmean = raw_sig_means, deconvoluted = deconvoluted)
  return(cellphonedb)
}


suppressWarnings({
    suppressPackageStartupMessages( library("annotables") )
    suppressPackageStartupMessages( library("future") )
    suppressPackageStartupMessages( library("future.apply") )
    suppressPackageStartupMessages( library("homologene") )
    suppressPackageStartupMessages( library("dplyr") )
    suppressPackageStartupMessages( library("optparse") )
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("OESingleCell") )
})

#=command line parameters setting=============================
option_list = list(
        make_option( c("--input", "-i" ), type = "character",
             help = "The input exprssion matrix in several possible format."),
        make_option( c("--informat", "-f" ), type = "character", default = "seurat",
             help = "The indication of type of input expression matrix, the possible type can be:
                    cellphonedb: the cellphonedb results directory.
                    seurat: the seurat object with celltype annotation.
                    xsv: the raw count matrix table in txt."),
        make_option( c("--celltype", "-l" ), type = "character", default = NULL,
             help = "the comma seperated list for desired cell types,to be used with the parameter --column4cell.
                    If NULL, all cell type will be used."),
        make_option( c("--species", "-s"), type="character", default="human",
            help="the species for the count matrix:human, mouse."),
        make_option( c("--column4cell", "-c"), type="character", default="celltype",
            help="the cell type annotation column in the cell annotation meta data."),
        make_option( c("--cellmeta", "-m" ), type = "character", default = NULL,
             help = "A table with the cell type annotation for each cell barcode.Only two columns: Cell,cell_type needed."),
        make_option( c("--slot", "-d"), type="character", default="data",
            help="the expression data type used as input, options can be raw or (normalized) data."),
        make_option( c("--assay", "-a"), type="character", default="RNA",
            help="the assay used to calulation in case of multimodal data."),
        make_option(c("--downsample", "-e"),
            type = "character", default = "30000",
            help = "the downsample number of cells "),
        make_option( c("--pvalue", "-p"), type="double", default=0.05,
            help="the P value threshold for significant interacting pairs."),
        make_option( c("--threshold", "-r"), type="double", default= 0.1,
            help="minimium percentage of cells expressing the specific ligand/receptor."),
        make_option( c("--threads", "-t"), type="integer", default=10,
            help="the threads use to run cellphonedb, 10 as default."),
        make_option( c("--iteratoins", "-n"), type="integer", default=1000,
            help="Number of iterations for the statistical analysis, use 1000 as default"),
        # make_option( c("--customedb", "-u"), type="character", default=NULL,
        #    help="the user customized interaction annotation database used to run cellphonedb. 
        #          If NULL(as default), the cellphoneDB will be used."),
        make_option( c("--genetype"), type="character", default="gene_name",
            help="Type of gene identifiers in the counts data, options can be:ensembl,gene_name,hgnc_symbol"),
        make_option( c("--output","-o"),type="character", default = "./CellCommunication",
            help="the output directory of Clustering results." ),
        make_option( c("--var2use", "-q" ), type = "character", default = NULL,
            help = "[OPTIONAL]The column name in cell metadata used as identity
                  of each cell combined with levels4var."),
        make_option( c("--levels4var", "-v" ), type = "character", default = NULL,
            help = "[OPTIONAL] subset of factor levels for the specified factor by --var2use.")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir,recursive=T)
    }
}
output_dir = normalizePath(output_dir )

# if ( (!is.null(opt$customedb)) && opt$customedb != "null" ){
#   customedb = NULL
# }else{
#  customedb = opt$customedb
# }

if ( is.null(opt$threads) ){
  threads = 10
}else{
  threads = opt$threads
}

if ( is.null(opt$iterations) ){
  iterations = 1000
}else{
  iterations = opt$iterations
}

if ( is.null(opt$threshold) ){
  threshold = 0.1
}else{
  threshold = opt$threshold
}

if ( is.null(opt$pvalue) ){
  pvalue = 0.05
}else{
  pvalue = opt$pvalue
}

if ( is.null( opt$genetype) ){
  genetype = "gene_name"
}else{
  genetype = opt$genetype
}
if ( is.null( opt$downsample) ){
  downsample = 30000
}else{
  downsample = opt$downsample
}
# setting the cores for parallization
options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
plan("multicore", workers = min(availableCores(), opt$threads)) # parallization using specified CPUs start from here
#=================================================================================
# read in the expression matrix 
#=================================================================================
if ( opt$informat == "cellphonedb" ){
    cellp = ParseCpdb(cpdb = opt$input, species = opt$species, pvalue = opt$pvalue, avgexpr = NULL)
    write.table(cellp$ligrec, file.path(output_dir, "cell_comm_annotation.xls"), sep = "\t", col.names = T, row.names = F, quote = F)
    saveRDS(cellp, file.path(output_dir, "cellphonedb_results.rds") )
    quit()
}

if ( opt$informat == "seurat" ){
    # the input is a seurat object produced by previous analysis
    seurat_ob = OESingleCell::readRDSMC( opt$input, cores = opt$threads )
    # if the input seurat object version is less than 3, upgrade it to version 3
    # if the input seurat object version is less than 3, upgrade it to version 3
    if ( seurat_ob@version < 3 ){
        seurat_ob = Seurat::UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
    }
    if ( !is.null(opt$assay) ){ DefaultAssay(seurat_ob) = opt$assay }

    #get the subset of cells used for running RunCellphonedb
    if ( !is.null(opt$levels4var)){
        if ( is.null(opt$var2use ) ){
            print("NO cell identity column name AVAILABLE! The group will be used as default.")
            ident2use = "group"
        }else{
            ident2use = opt$var2use
        }
        cluster_list = unlist(strsplit( opt$levels4var,",",perl = T))
        seurat_ob = SubsetData( seurat_ob, subset.name = ident2use, accept.value = cluster_list)
        seurat_ob@meta.data[[ident2use]]=factor(seurat_ob@meta.data[[ident2use]],levels = sort(unique(seurat_ob@meta.data[[ident2use]])))
    }


    if (ncol(seurat_ob) > 40000) {
    # library(sampling) need install,instead with manual function
        ratio <- as.numeric(downsample) / ncol(seurat_ob)
        metadata_temp <- as.data.frame(seurat_ob@meta.data)
        # strata(metadata_temp,stratanames="clusters",ratio,description=FALSE)
        cells_sample <- c()
        for (i in unique(seurat_ob$clusters)) {
          cells_temp <- rownames(metadata_temp)[which(metadata_temp$clusters == i)]
          cells_temp_sample <- sample(cells_temp, ceiling(length(cells_temp) * ratio), replace = FALSE, prob = NULL)
          cells_sample <- append(cells_sample, cells_temp_sample)
        }
        seurat_ob <- subset(seurat_ob, cells = cells_sample)
    }
    print(dim(seurat_ob))




    #get the subset of cells used for visualization if necessay
    if ( !is.null(opt$celltype)){
        if ( is.null(opt$column4cell ) ){
            print("NO cell type column name AVAILABLE! The clusters annotation will be used as default.")
            ident2use = "clusters"
        }else{
            ident2use = opt$column4cell
        }
        cluster_list = unlist(strsplit( opt$celltype,",",perl = T))
        seurat_ob = SubsetData(seurat_ob,cells =
                           OldWhichCells( seurat_ob, subset.name= ident2use, accept.value = cluster_list))
    }
    
    #get the count matrix
    counts <- GetAssayData(seurat_ob, slot = opt$slot)
    
    if ( is.null(opt$column4cell) ){
      metadata <- Seurat::FetchData(seurat_ob, vars = "ident") %>%
                    tibble::rownames_to_column(var = "Cell") %>%
                    dplyr::rename("cell_type" = ident )
    }else{
      metadata <- Seurat::FetchData(seurat_ob, vars = opt$column4cell ) %>%
                    tibble::rownames_to_column(var = "Cell")
      colnames(metadata) = c("Cell","cell_type")
    }
    rm(seurat_ob)
    gc()
}else if ( opt$informat == "xsv" ){
    counts = vroom::vroom(file.path(opt$input), col_names = T, comment = "#")
    counts = as.data.frame(counts)
    rownames(counts) = counts[,1]
    counts = counts[,-1]
    if ( is.null( opt$cellmeta) | opt$cellmeta == "null" ){
      stop("Cell celltype annotation is Mendatory for user provided count matrix!")
    }else{
      metadata = read.table(opt$cellmeta, sep = "\t", header = T, stringsAsFactors = F)
      colnames(metadata) = c("Cell", "cell_type")
    }
    counts = counts[,metadata$Cell]
    
    if ( !is.null(opt$celltype)){
        if ( is.null(opt$column4cell )  || opt$column4cell == "null"){
          column4cell = "cell_type"
        }else{
          column4cell = opt$column4cell
        }
        desired_cell_types = unlist(strsplit( opt$celltype,",",perl = T))
        desired_cells = metadata %>% dplyr::filter( !!ensym(column4cell) %in% desired_cell_types ) %>% 
                                pull( Cell )
        counts = counts[, desired_cells]
    }
    
    
}

# read the data, it should be a cell(row)-gene(column) expression matrix in data.frame format
# the expression matrix should integrate the cell type annotation stored in the matrix as 
# cell_type.
if ( opt$species == "mouse" ){
  mouse2human <- homologene::human2mouse(annotables::grch38$symbol) %>% dplyr::select( mouseGene, humanGene )
  genes2use <- intersect(names(Seurat::CaseMatch(rownames(counts), mouse2human$mouseGene)), rownames(counts) )
  counts <- counts[genes2use,]

  rownames(counts) = plyr::mapvalues(rownames(counts),
                                     from = mouse2human$mouseGene,
                                     to = mouse2human$humanGene, warn_missing = F)
}

# tempdir <- tempdir()
tempdir <- output_dir
# using vroom to speed up the count matrix writing
counts_out <- tibble::rownames_to_column( as.data.frame(counts), var = 'Gene')
vroom::vroom_write(x = counts_out, path = file.path(tempdir,'counts.tsv'), col_names = T, quote = "none", delim = "\t")
write.table(metadata, file = file.path(tempdir,'metadata.tsv'), 
            quote = F, col.names = T, row.names = F, sep = '\t')

outdir <- file.path(output_dir, "out")
# run the cellphonedb with the statistical method
cellphonedb_bin = system("which cellphonedb", intern = T)
if ( !file.exists(cellphonedb_bin) ){
    stop("NO cellphonedb command detected! Please make sure cellphonedb path correctly setted!")
}

if(!file.exists(file.path(output_dir, "cellphone.db"))){
  file.copy("/data/database/cellphonedb/cellphone.db", file.path(output_dir, "cellphone.db"))
}

setwd(output_dir)
system(glue::glue("{cellphonedb_bin} method statistical-analysis {tempdir}/metadata.tsv {tempdir}/counts.tsv --counts-data {genetype} --iterations={iterations} --threshold={threshold} --threads={threads} --pvalue={pvalue} --output-path={outdir}  --database cellphone.db --verbose"))

count_long = as.data.frame(counts) %>% tibble::rownames_to_column(var = "gene" ) %>%
    tidyr::gather("Cell","expr", 2:dim(.)[2])
count_long$cell_type = plyr::mapvalues(count_long[["Cell"]],
                            from = metadata$Cell,
                            to = as.character(metadata$cell_type), warn_missing = F)
avgexpr = count_long %>% dplyr::select( gene, cell_type, expr ) %>%
    dplyr::mutate(gene = toupper(gene)) %>%
    dplyr::group_by( cell_type, gene ) %>%
    dplyr::mutate(exprs = log2(mean(expr)+0.0001) ) %>%
    dplyr::select( gene, cell_type, exprs ) %>%
    dplyr::distinct()

# parse the output of cellphonedb using human genes
cellphonedb <- ParseCpdb( outdir , species = opt$species, pvalue = pvalue, avgexpr = avgexpr  )
write.table(cellphonedb$ligrec, file.path(output_dir, "cell_comm_annotation.xls"), sep = "\t", col.names = T, row.names = F, quote = F)
saveRDS(cellphonedb, file.path(output_dir, "cellphonedb_results.rds"))

if(file.exists(file.path(output_dir, "counts.tsv"))){
  file.remove(file.path(output_dir, "counts.tsv"))
}
if(file.exists(file.path(output_dir, "metadata.tsv"))){
  file.remove(file.path(output_dir, "metadata.tsv"))
}
