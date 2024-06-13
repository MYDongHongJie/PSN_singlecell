#!/usr/bin/env Rscript
#================================================================================
# function definition
#================================================================================
saveRDSMC <- function(object, file, threads= 10) {
    threads <- min(threads, future::availableCores() )
    message("using ", threads, " threads for compression.")
    #con <- pipe(paste0("xz -T", threads, " -9 -f > ", file), "wb")
    con <- pipe(paste0("pigz -p ", threads, " -9 -f > ", file), "wb")
    saveRDS(object, file = con)
    on.exit(if(exists("con")) close(con))
}
  
addFeatureAnno <- function(
object,
assay = "RNA",
gtf
){
    # read in the reference genome file if provided
    suppressPackageStartupMessages( library("rtracklayer") )
    gtf =  readGFF(gtf) # gtf is a data.frame of gene annotation
    gtf_df = gtf %>% rename(chromosome = seqid, gene = gene_name)
    gtf_df$chromosome = paste0("chr", gtf_df$chromosome)
    if ( assay == "RNA" ){
    	features_loc = file.path(matrix_path,"features.tsv.gz")
    	gene_loc = file.path(matrix_path[1],"genes.tsv")
    	prev_ver_3 = file.exists(gene_loc)
    	feature.names = read.delim( file = ifelse( prev_ver_3, gene_loc, features_loc), header = F, stringsAsFactors = F)
    	gtf_df$gene = make.unique( gtf_df$gene )
    	rownames(gtf_df) = gtf_df$gene
    	gtf_df = gtf_df[feature.names[,1],] # use the row order in the gene.tsv or feature.tsv.gz
    	object[[assay]]= AddMetaData(seurat_ob[[assay]] , metadata = gtf_df )
    	Misc(object, slot = "FeatureAnnotation") = gtf_df
    }else{
    	gtf_df$gene = make.unique( gtf_df$gene )
    	rownames(gtf_df) = gtf_df$gene
    	Misc(object, slot = "FeatureAnnotation") = gtf_df
    }
    
    return(object)
}

# read cellranger-atac aggr results or recounts
Read10X_aggr <- function(
  data.dir,
  sample.names, # vector of sample names
  split = FALSE,
  gene.column = 2,
  names.delm = "-",
  unique.features = TRUE
) {
  raw.counts <- Read10X(data.dir = data.dir, gene.column = gene.column,
                        unique.features = unique.features)
  rownames(x = raw.counts) <- gsub(pattern = '_', replacement = '-', x = rownames(x = raw.counts))
  raw.bc <- colnames(x = raw.counts)
  lib.id <- stringr::str_extract(string = raw.bc, pattern = '[0-9]*$')
  if (length(x = sample.names) != dim(x = table(lib.id))) {
    stop("The number of sample names dosen't match the number of libraries")
  }
  trimed.bc <- gsub(pattern = '[0-9]*$', replacement = '', x = raw.bc)
  trimed.bc <- gsub(pattern = '-', replacement = '', x = trimed.bc)
  bc <- data.frame(raw.bc = raw.bc,
                   lib.id = as.numeric(lib.id),
                   trimed.bc = trimed.bc)
  bc[['sample']] <- sample.names[bc[['lib.id']]]
  bc[['new.bc']] <- paste0(bc[['sample']], names.delm, bc[['trimed.bc']])
  colnames(x = raw.counts) <- bc[['new.bc']]
  if (split) {
    bc[['sample']] <- factor(x = bc[['sample']], levels = sample.names)
    bc <- split(x = bc, f = bc[, 'sample'], drop = FALSE)
    new.counts <- lapply(X = bc, function(x) {
      raw.counts[, x[['new.bc']], drop = FALSE]
    })
    return(new.counts)
  } else {
    return(raw.counts)
  }
}



    # save the object in specified format
saveObject <- function(
    object,
    outdir = NULL,
    outformat = "seurat",
    assay = "RNA",
    prefix = "seurat",
    splitby = NULL
){
    outformat = tolower(outformat)
    # subcommand create is invoked
    if ( is.null(outformat)  ){
        outformat = "seurat"
    }
    if ( !is.null( outdir ) ){
        output_dir = outdir
    }else{
        warnings("NO output directory AVAIABLE!")
        output_dir = getwd()
    }

    if ( outformat == "seurat" ){
        saveRDSMC(object, file.path(output_dir, paste0(prefix,"_seurat.rds") ))
    }
    if( outformat == "loom" ){
        seurat_loom = as.loom(object,
        filename= file.path(output_dir, paste0(prefix, ".loom")), overwrite=T)
        seurat_loom$close_all()
    }
    if ( outformat == "sce" ){ # output in SingleCellExperiment format
        seurat_sce = as.SingleCellExperiment( object )
        saveRDSMC(seurat_sce, file.path(output_dir, paste0(prefix, "_sce.rds")))
    }
    if ( outformat == "celldataset" ){
        cellDataSet = as.CellDataSet( object )
        saveRDSMC( cellDataSet, file.path(output_dir, paste0(prefix, "_celldataset.rds")))
    }

    #only support the same specie for each sample here
    if ( outformat == "tenx" ){
        # save the expression matrix seperately by the specifed conditions
        # usually by sampleid
        if ( !is.null(splitby) ){
            seurat_ob_list = SplitObject(object, split.by = splitby )
        }else{
            seurat_ob_list = list(Cellranger = object)
        }
        lapply( names(seurat_ob_list), function(namex){
            objectx = seurat_ob_list[[namex]]
            sampleid = as.vector(unique(objectx@meta.data$sampleid))
            suppressPackageStartupMessages( library("DropletUtils"))
            if ( is.null(splitby) & length(sampleid) > 1){ # there are more than one sample in the object
                write10xCounts(
                    GetAssayData(objectx, assay = assay, slot = "counts"),
                    path = file.path( output_dir, namex),
                    barcodes = colnames(objectx), gene.id = rownames(objectx),
                    gene.symbol=rownames(objectx), gene.type="Gene Expression",
                    overwrite=FALSE, type="sparse", version="3")
            }else{ # only one sample in the object, so the sample id prefix for each cell should be removed
                barcodes = data.frame(barcode=colnames(objectx))
                barcodes = apply(barcodes,1,function(x)gsub( "^.*_", "",x, perl = T ))
                write10xCounts(
                    GetAssayData(objectx, assay = assay, slot = "counts"),
                    path=   file.path( output_dir, namex),
                    barcodes = as.vector(barcodes), gene.id = rownames(objectx),
                    gene.symbol=rownames(objectx), gene.type="Gene Expression",
                    overwrite=FALSE, type="sparse", version="3")
            }
        })
    }
}

# read in the expression matrix in possible format
# input the input expression matrix in several format
# informat the format of the input expression matrix
# metadata the sample metadata used to produce the metadata for each cell
# assay the assay for this run
# peak.meta the peak count for each type of region for each cell, used in atac-seq
# fragment the unique fragment in each region
readObject <- function(
    input = NULL,
    informat = NULL,
    metadata = NULL,
    assay = "RNA",
    peak.meta = NULL,
    fragment = NULL,
    aggr = FALSE,
    ...
){
    informat = tolower(informat)
    if ( informat == "seurat" ){# the input is a seurat object which may contain more than one sample
        input = normalizePath(sub("\\/$","",input,perl=T))
        seurat_ob = readRDSMC( input, cores = 10 ) # parallel reading
        if ( seurat_ob@version < 3){
            seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
        }

        # if the metadata is avaiable, update the seurat@meta.data using the input metadata
        if ( !is.null(metadata) ){
            additional_metadata = read.csv(metadata,sep=",",header =T, colClasses = c("character") ,row.names="sampleid")
            cellnames = Cells(seurat_ob)
            sampleidx =  gsub("(_|-)[ATCG]{16,}.*","",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
            #integrate the additional metadata from the assay design
            additional_cell_meta = vector()
            for ( colidx in colnames(additional_metadata) ){
                additional_cell_meta = cbind(additional_cell_meta, as.vector(additional_metadata[sampleidx, colidx]))
            }
            colnames(additional_cell_meta) = colnames(additional_metadata)
            rownames(additional_cell_meta) = cellnames
            additional_cell_meta = as.data.frame(additional_cell_meta)
            seurat_ob = AddMetaData( seurat_ob, additional_cell_meta)
        }
        return
    }

    if ( informat == "tenx" ){ #10X cellranger output regardless of V2 or V3
        if ( !is.null(metadata) ){
            assay_metadata = read.csv(metadata,sep=",",header =T, colClasses = c("character") )
            rownames(assay_metadata) = assay_metadata$sampleid
        }
        # read in the 10X data produced by cellrnager V3.X or Cellranger V2.X
        # Initialize the Seurat object with the raw (non-normalized data).
        # Directory containing the matrix.mtx.gz, feature.tsv.gz, and barcodes.tsv.gz files provided by 10X Cellranger V3
        # or matrix.mtx, gene.tsv, barcodes.tsv produced by Cellranger V2.
        # A vector or named vector can be given in order to load several data directories of samples simutailiously.
        # If a named vector is given, the cell barcode names will be prefixed with the name.
        tenx_path = normalizePath(sub("\\/$","",input,perl=T))
        mtx_exists = length(Sys.glob(file.path(tenx_path,"*","outs","filtered_*","barcodes*"))) > 0
        if ( mtx_exists ){ #cellrnager V3.X
            matrix_path = unlist(lapply(assay_metadata$sampleid,function(x)
                                        unique(dirname(Sys.glob(file.path(tenx_path,x,"outs","filtered_*","barcodes*")) ))))
        }else{ #cellranger V2.X
            matrix_path = unlist(lapply(assay_metadata$sampleid,function(x)
                                        unique(dirname(Sys.glob(file.path(tenx_path,x,"outs","filtered_*","*", "barcodes*")) ))))
        }
        # gzfiles = Sys.glob(file.path(tenx_path,"*",outs_filtered_path,"*.gz"))
        # lapply(gzfiles,function(x)gunzip(x, destname = gsub("[.]gz$", "", x), remove = F, overwrite = T))
        # lapply(Sys.glob(file.path(tenx_path,"*",outs_filtered_path,"features.tsv")),
        #                 function(x)file.rename(x,gsub("features\\.tsv","genes\\.tsv",x)))
        matrix_path = unique(matrix_path)
        if ( aggr == T ){
            if ( assay == "ATAC" ){
                countMatrixSparse = Read10X_aggr(matrix_path, sample.names = assay_metadata$sampleid, gene.column = 1)
            }else{
                countMatrixSparse = Read10X_aggr(matrix_path, sample.names = assay_metadata$sampleid, gene.column = 2)
            }
        }else{
            names(matrix_path) = assay_metadata$sampleid
            if ( assay == "ATAC" ){
                countMatrixSparse <- Read10X(matrix_path, gene.column = 1)
            }else{
                countMatrixSparse <- Read10X(matrix_path, gene.column = 2)
            }
        }
        if ( class(countMatrixSparse) == "list" ){ # there is different library_type in this assay
		   #construct the seurat object using the meta data above
           if ( !is.null(countMatrixSparse[["Gene Expression"]]) ){
		       seurat_ob <- CreateSeuratObject( countMatrixSparse[["Gene Expression"]], names.field = 2, assay = "RNA", names.delim = "-" )
           }
           if ( !is.null(countMatrixSparse[["Antibody Capture"]]) ){
			   rownames(countMatrixSparse[["Antibody Capture"]]) <- 
							   gsub("_[control_]*TotalSeqB", "", rownames(countMatrixSparse[["Antibody Capture"]]))
               seurat_ob[["ADT"]] <- CreateAssayObject(counts = countMatrixSparse[["Antibody Capture"]][, colnames(seurat_ob)])
           }
           if ( !is.null(countMatrixSparse[["CRISPR Guide Capture"]]) ){
			   rownames(countMatrixSparse[["CRISPR Guide Capture"]]) <- 
							   gsub("_[control_]*TotalSeqB", "", rownames(countMatrixSparse[["CRISPR Guide Capture"]]))
               seurat_ob[["CRISPR"]] <- CreateAssayObject(counts = countMatrixSparse[["CRISPR Guide Capture"]][, colnames(seurat_ob)])
           }
           
        }else{
		   #construct the seurat object using the meta data above
		   seurat_ob <- CreateSeuratObject( countMatrixSparse, names.field = 2, assay = assay, names.delim = "-" )
        }

        rm(countMatrixSparse)

        # cell barcode pattern: sample_index_xxxxxxxxxxxxxxxx V2
        # cell barcode pattern: sample_index_xxxxxxxxxxxxxxxxxx V3
        cellnames = Cells(seurat_ob)
        sampleidx =  gsub("(_|-)[ATGC]{16,18}.*","",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
        index = which( assay_metadata$sampleid == rownames(assay_metadata))
        names(index) = rownames(assay_metadata)
        cellindex = index[sampleidx]
        bare_cb = gsub( "^.*(_|-)", "", cellnames, perl = T)
        indexed_cb = paste( bare_cb, cellindex, sep= "-")
        #integrate the metadata from the assay design
        cell_meta = vector( )
        for ( colidx in colnames(assay_metadata) ){
            cell_meta= cbind(cell_meta, as.vector(assay_metadata[sampleidx, colidx]))
        }

        colnames(cell_meta) = colnames(assay_metadata)
        rownames(cell_meta) = cellnames
        cell_meta = as.data.frame(cell_meta)
        cell_meta$orig.ident = indexed_cb
        seurat_ob = AddMetaData(seurat_ob, metadata = cell_meta)


        if ( assay == "ATAC" ){
            if ( !is.null(peak.meta) ){
                peakmeta <- read.table( file = peak.meta, header = TRUE, row.names = 1 )
                peakmeta <- peakmeta[seurat_ob$orig.ident,]
                rownames(peakmeta) <- seurat_ob@meta.data %>%
                                        tibble::rownames_to_column(var = "barcode") %>%
                                        filter( orig.ident %in% rownames(peakmeta) ) %>%
                                        pull( barcode )
                seurat_ob = AddMetaData(seurat_ob, peakmeta )
            }

            if ( !is.null(fragment) ){
                # keep only the barcodes identified as cells in Seurat object
                cellbc = Cells(seurat_ob)
                seurat_ob = RenameCells(seurat_ob, new.names = seurat_ob@meta.data$orig.ident)
                FilterFragments( fragment.path = fragment, cells = colnames(seurat_ob),
                output.path =file.path(output_dir, "fragments.tsv" ) )
                seurat_ob <- SetFragments( seurat_ob, file = file.path(output_dir,"fragments.tsv.bgz" ) )
                seurat_ob <- NucleosomeSignal(seurat_ob)
                seurat_ob = RenameCells(seurat_ob, new.names = cellbc)
            }
        }
    }

    #read data from the count matrix file
    #the input matrix can be stored as sparse format or text format
    # the metadata should include the count matrix file name in on column named as "sampletsv"
    # each count matrix should only be put into a directory named as their samples' name
    if ( informat == "xsv" ){
        if ( !is.null(metadata) ){
            assay_metadata = read.csv(metadata,sep=",",header =T, colClasses = c("character") )
            rownames(assay_metadata) = assay_metadata$sampleid
        }
        tenx_path = normalizePath(sub("\\/$","",input,perl=T))
        matrix_path = file.path( tenx_path, assay_metadata$sampleid,assay_metadata$sampletsv)
        # matrix_path = apply( assay_metadata,1, function(samplex) file.path(tenx_path,samplex["sampletsv"]))
        names(matrix_path) = assay_metadata$sampleid
        # TO DO mkdir
        # tranform the tsv to mtx and then read in the expression matrix using Read10X
        options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
        plan("multicore", workers = length(matrix_path))
        seu_list = future.apply::future_lapply(seq_along(matrix_path), function(idx){
          # use fread to speed up reading count matrix in case of its size
          mycountmatrix = vroom::vroom(file.path(matrix_path[idx]), col_names = T, comment = "#")
          mtx_dir = dirname(matrix_path[idx])
          names(mtx_dir) = names(matrix_path)[idx]
          if (  tolower(opt$transpose) == "true" ){ # count matrix from BD single cell workflow is cell-gene
            mycountmatrix = as.matrix(mycountmatrix[,-1])
            # set the cell barcode as the random DNA string to substitude the digit index
            # in case of count matrix from BD platform
            rownames(mycountmatrix) = unlist(lapply(1:dim(mycountmatrix)[1],
                                          function(x) paste(sample(LETTERS[c(1,3,7,20)], size=27, replace=TRUE), collapse = "")))
            countMatrixSparse = Matrix::Matrix(t(mycountmatrix), sparse=T)
          }else{
            countMatrixSparse = Matrix(as.matrix(mycountmatrix), sparse=T)
          }
          Matrix::writeMM(countMatrixSparse,file.path(mtx_dir,"matrix.mtx"))
          gene_id = data.frame(gene_id=row.names(countMatrixSparse),gene_name=row.names(countMatrixSparse))
          write.table(gene_id, file.path(mtx_dir,"genes.tsv"),sep="\t",col.names=F,row.names=F,quote=F)
          barcodes = data.frame(barcode=colnames(countMatrixSparse))
          barcodes = apply(barcodes,1,function(x)gsub( "^.*_", "",x, perl = T ))
          write.table( barcodes, file.path(mtx_dir,"barcodes.tsv"), col.names=F, row.names=F,quote=F)
          rm(mycountmatrix) #remove the original matrix to reduce memory usage
          countMatrixSparse <- Read10X(mtx_dir)
          #construct the seurat object using the meta data above
          seu <- CreateSeuratObject( countMatrixSparse, names.field = 2, assay = assay,names.delim = "-" )
        })
        
        seurat_ob = merge(seu_list[[1]], seu_list[2:length(seu_list)])
        assay_metadata = assay_metadata %>% select(-sampletsv)
        #add in the metadata of all cells
        #preserve the sample group infomation into the seurat object,the order of the sample
        #is the same as the order in the metadata list
        cellnames = Cells(seurat_ob)
        #the barcodes of cells in first sample don't prefix with sample index number
        #so we add it for convenince
        sampleidx =  gsub("_[ATGC]{16,}.*","",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
        #integrate the metadata from the assay design
        cell_meta = vector( )
        for ( colidx in colnames(assay_metadata) ){
          cell_meta= cbind(cell_meta, as.vector(assay_metadata[sampleidx, colidx]))
        }
        colnames(cell_meta) = colnames(assay_metadata)
        rownames(cell_meta) = cellnames
        cell_meta = as.data.frame(cell_meta)
        cell_meta$orig.ident = cellnames
        #construct the seurat object using the meta data above
        seurat_ob <- AddMetaData( seurat_ob, metadata = cell_meta )
    }



    if ( informat == "visium" ){
        if ( !is.null(metadata) ){
            assay_metadata = read.csv(metadata,sep=",",header =T, colClasses = c("character") )
            rownames(assay_metadata) = assay_metadata$sampleid
        }
        tenx_path = normalizePath(sub("\\/$","",input,perl=T))
        matrix_path = unlist(lapply(assay_metadata$sampleid,function(x)
                            unique(dirname(Sys.glob(file.path(tenx_path,x,"outs","*.h5")) ))))

        names(matrix_path) = assay_metadata$sampleid
        # tranform the tsv to mtx and then read in the expression matrix using Read10X
        options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
        plan("multicore", workers = length(matrix_path))
        seu_list = future.apply::future_lapply(seq_along(matrix_path), function(idx){
            # use fread to speed up reading count matrix in case of its size
            seu_ob = Load10X_Spatial(data.dir=matrix_path[idx],assay = "Spatial", 
                                     filter.matrix=T, slice= paste0("slice",idx))
        })

        seurat_ob = merge( seu_list[[1]], seu_list[2:length(seu_list)] )
        #add in the metadata of all cells
        #preserve the sample group infomation into the seurat object,the order of the sample
        #is the same as the order in the metadata list
        cellnames = Cells(seurat_ob)
        #the barcodes of cells in first sample don't prefix with sample index number
        #so we add it for convenince
        raw_CB =  gsub("([ATGC]{16,}).*","\\1",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
        sampleidx =  as.numeric(gsub("[ATGC]{16,}(-1)?_","",cellnames,perl=T))
        #integrate the metadata from the assay design
        cell_meta = vector( )
        for ( colidx in colnames(assay_metadata) ){
            cell_meta= cbind(cell_meta, as.vector(assay_metadata[sampleidx, colidx]))
        }
        colnames(cell_meta) = colnames(assay_metadata)
        rownames(cell_meta) = cellnames
        cell_meta = as.data.frame(cell_meta)
        cell_meta$orig.ident = cell_meta$sampleid
        #construct the seurat object using the meta data above
        seurat_ob = AddMetaData(seurat_ob, metadata=cell_meta)
        new_CB = paste(cell_meta$sampleid, raw_CB, sep = "-")
        seurat_ob = RenameCells(seurat_ob, new.names = new_CB)
    }


    if ( informat == "loom" ){
        suppressPackageStartupMessages( library("loomR") )
        seurat_loom = connect( filename = input, mode = "r" )
        seurat_ob = as.Seurat(seurat_loom)
    }

    if ( informat == "sce" ){
        suppressPackageStartupMessages( library("SingleCellExperiment") )
        suppressPackageStartupMessages( library("scater") )
        seurat_sce = readRDSMC( input, cores = 10 )
        seurat_ob = as.Seurat(seurat_sce)
    }

    if ( informat == "h5ad"){
        seurat_ob = ReadH5AD(input, assay = assay, overwrite = T)
    }

    # Read count matrix from 10X CellRanger hdf5 file.
    # This can be used to read both scATAC-seq and scRNA-seq matrices.
    if ( informat == "h5" ){
        if ( !is.null(metadata) ){
            assay_metadata = read.csv(metadata,sep=",",header =T, colClasses = c("character") )
        }
        countMatrixSparse = Read10X_h5( input, use.names = T)
        #construct the seurat object using the meta data above
        if ( class(countMatrixSparse) == "list" ){
		   #construct the seurat object using the meta data above
           if ( !is.null(countMatrixSparse[["Gene Expression"]]) ){
		       seurat_ob <- CreateSeuratObject( countMatrixSparse[["Gene Expression"]], names.field = 2, assay = "RNA", names.delim = "-" )
           }
           if ( !is.null(countMatrixSparse[["Antibody Capture"]]) ){
			   rownames(countMatrixSparse[["Antibody Capture"]]) <- 
							   gsub("_[control_]*TotalSeqB", "", rownames(countMatrixSparse[["Antibody Capture"]]))
               seurat_ob[["ADT"]] <- CreateAssayObject(counts = countMatrixSparse[["Antibody Capture"]][, colnames(seurat_ob)])
           }
           if ( !is.null(countMatrixSparse[["CRISPR Guide Capture"]]) ){
			   rownames(countMatrixSparse[["CRISPR Guide Capture"]]) <- 
							   gsub("_[control_]*TotalSeqB", "", rownames(countMatrixSparse[["CRISPR Guide Capture"]]))
               seurat_ob[["CRISPR"]] <- CreateAssayObject(counts = countMatrixSparse[["CRISPR Guide Capture"]][, colnames(seurat_ob)])
           }
           
        }else{
		   #construct the seurat object using the meta data above
		   seurat_ob <- CreateSeuratObject( countMatrixSparse, names.field = 2, assay = assay, names.delim = "-" )
		   if ( assay == "ATAC" ){
		        if ( !is.null(peak.meta) ){
					peakmeta <- read.csv( file = peak.meta, header = TRUE, row.names = 1 )
					seurat_ob = AddMetaData(seurat_ob, peakmeta )
				}

				if ( !is.null(fragment) ){
					# keep only the barcodes identified as cells in Seurat object
                    suppressPackageStartupMessages( library("Signac") )
					FilterFragments( fragment.path = fragment,
									 cells = colnames(seurat_ob),
									 output.path =file.path(output_dir, "fragments.tsv" ) )
					seurat_ob <- SetFragments( seurat_ob, file = file.path(output_dir,"fragments.tsv.bgz" ) )
					seurat_ob <- NucleosomeSignal(seurat_ob)
				}
			}
        }
        cellnames = Cells(seurat_ob)
        # the index order is the same as the row index of the cellranger aggr library metadata, so make sure
        # it is the same order of the sample metadata before running
        sampleidx =  gsub("(_|-)[ATGC]{16,}-?","",cellnames,perl=T)
        #integrate the metadata from the assay design
        cell_meta = vector( )
        for ( colidx in colnames(assay_metadata) ){
            cell_meta= cbind(cell_meta, as.vector(assay_metadata[sampleidx, colidx]))
        }
        colnames(cell_meta) = colnames(assay_metadata)
        barcodes =  gsub("-.*$","",cellnames,perl=T)
        barcodes = paste( cell_meta[,"sampleid"],barcodes, sep = "-")
        cell_meta = as.data.frame(cell_meta)
        rownames(cell_meta) = cellnames
        cell_meta$new_barcode = barcodes
        cell_meta$orig.ident = cell_meta$sampleid
        seurat_ob = AddMetaData( seurat_ob, cell_meta)
        rm(countMatrixSparse)
    }

    DefaultAssay(seurat_ob) = assay
    return( seurat_ob )
}

#================================================================================
# package loading
#================================================================================
suppressWarnings({
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("Matrix") )
    suppressPackageStartupMessages( library("OESingleCell") )
    suppressPackageStartupMessages( library("docopt") )
    suppressPackageStartupMessages( library("dplyr") )
    suppressPackageStartupMessages( library("future") )
    suppressPackageStartupMessages( library("future.apply") )
    suppressPackageStartupMessages( library("tibble") )
})

# =command line parameters setting=============================
'usage:sctool.R
     sctool.R create [options]
     sctool.R subset [options]
     sctool.R removeEmptyDrops [options]
     sctool.R transform [options]
     sctool.R merge  [options]
     sctool.R fetch  [options]
     sctool.R integrate [options]
     sctool.R update [options]
     sctool.R avgexp [options]
     sctool.R summary [options]
     sctool.R split [options]
     sctool.R downsample [options]
     sctool.R export2cbbrowser [options]
     sctool.R visfeature [options]
     sctool.R impute [option]

options:
    global_options:
        -i <input>, --input <input>                   # The input exprssion matrix in several possible format.
        -f <informat>, --informat <informat>          # The indication of type of input expression matrix, the possible type can be:
                                                        sce: singlecellexperiement object supported by R packages:scater,scran etc.
                                                        h5: hdf5 formatted molecular quantification from cellranger etc.
                                                        loom: loom format, a specialized hdf5 format.
                                                        h5ad: the format support by scanpy from anndata in python.
                                                        tenx:the directory of cellranger count/aggr output with sampleid named subdirectory.
                                                        xsv: the raw gene count matrix file, it can be very large. Please make sure the
                                                              data path is Your_Project/sampleid/sampleid.tsv/csv.
                                                        cellDataSet: the format for monocle.
                                                        Seurat: the seurat object from the clustering results.
                                                        visium: the output from spaceranger ananlysis of 10X Visium spatial transcriptomics data.
        -o <outdir>, --outdir <outdir>                # the output directory of results.
        -d <outformat>, --outformat <outformat>       # the output format of the expression matrix, the choice can be:
                                                        sce:singlecellexperiement object supported by R packages:scater,scran etc.
                                                        loom: loom format.
                                                        h5ad: the format support by scanpy from anndata in python.
                                                        tenx:the result of cellranger V3 count/aggr with sampleid as its subdirectory.
                                                        Seurat: the seurat object from the clustering results.
                                                        cellDataSet: the data object for monocle.
        -p <prefix>, --prefix <prefix>                # the prefix of output file without file extension.[default: seurat]
        -k, <seperateby>, --seperateby <seperateby>   # [OPTIONAL]save the expression matrix in  mtx format like cellranger output
                                                        for each group,usually sampleid. If not availiable, all samples will be saved
                                                        in one global object.[default: NULL]
        --slot <slot>                                 # [OPTIONAL]the expression matrix type used for calculation,
                                                        options can be "counts", "data", "scale.data".[default: counts]
        --assay <assay>                               # [OPTIONAL]the assay object to use, which is avaiable for multimod data.[default: RNA]
    subcommand:
        create        # crate the desired expression matrix in the specified format from specific format of input expression.
            options:
                -m <metadata>, --metadata <metadata>   # the sample metadata which must include sample id in this assay design.
                --gtf <gtf>                            # [OPTIONAL]the reference genome annotaion file from 10X.
                --atac_cellmeta <atac_cellmeta>        # [OPTIONAL] the peak annotation file "singlecell.csv" for each cell in ATAC-seq.
                --fragment <fragment>                  # [OPTIONAL] the fragment annotion file "fragments.tsv.gz" in ATAC-seq.
                --transpose <transpose>                # [OPTIONAL] transform the count table before processing the cell-feature matrix,
                                                          FALSE as default.[default: FALSE]
                --aggr <aggr>                          # [OPTIONAL] the results derived from cellranger aggr.[default: FALSE]
        transform     # tranform the specified cell annotation from one to another.
            options:
                --variable <variable>                  # the variable to be transformed.
                --from <from>                          # the original level id list in specified variable
                --to <to>                              # the adjusted level id list with order specified by --from.
        subset        # filter the cells using the specified additions.
            options:
                --geneset <geneset>                    # features/variables to keep in TAB format.
                --feature4subset <feature4subset>      # Logical expression on features in matrix indicating cells to keep.Example: "MSA4>1".
                --cellfilter <cellfilter>              # Logical expression indicating cells to keep. Example: "clusters %in% c(1,2,3)".
                --levels <levels>                      # The subset of variable levels in variable specified by --group2use.
                --group2use <group2use>                # the cell names or annotation column names set for the cell annotation table.
                --low <low>                            # the lower bound of the continious variable specified by --group2use.[default: NULL]
                --high <high>                          # the higher bound of the continious variable specified by --group2use.[default: NULL]
                --invert <invert>                      # reverse the selection conditions specified by all other conditions.[default: F]
        fetch         # Retreives data (feature expression, PCA scores, metrics, etc.) for a set of cells.A table with cells as rows and cellular data as columns returns.
            options:
                --vars <vars>                          # comma seperated List of all variables to fetch, use keyword ident to pull identity classes.
                --header <header>                      # [OPTIONAL]change the header of the output table if specified. Notice that the first column is always the cell barcode.
        avgexp        # get the pseudobulk expression of each gene by the groupping variable specified.
            options:
                --avgby <avgby>                        # cell groupping factor for calculate the average expression of cells in each group
                --features <features>                  # features/variables to keep in TAB format.
        merge         # add another sample to the existing object.
            options:
                --toadd <toadd>                        # one or more of object list to merge into the target object.
                --adformat <adformat>                  # The format of input expression matrix to add
                --mergedata <mergedata>                # logical string T/F to indicate wether to merge the data slot
        split         # split the object according to the cell metadata
            options:
                --splitby <splitby>                    # the cell annotation column name.
        downsample    # downsample the exppression matrix in case of too many cell and scarcely disequal cell numbers in each sample
            options:
                --targetN <targetN>                    # the number of cell to keep after downsample.
                --filters <filters>                    # the conditional expression for subsetting the cells using seurat@meta.data slot for downsample. Example: "clusters %in% c(1,2,3)"
                --doinvert <doinvert>                  # whether to get complement set of specified cells with --filters[default: FALSE].
        removeEmptyDrops  # remove the empty droplet using DropletUitls
            options:
                --fdr <fdr>                            # the significance level for a cell to be identified as an empty droplet.[default: 0.01]
                --iters <iter>                         # the iteration to be used in calculation.[defalut: 1000]
        update        # using the additional cell annotation to updata the cell annotation.
            options:
                --admeta <admeta>                      # the additional metadata which includes any groupping design in this assay.
                --adfeatures <adfeature>               # the additional metadata of genes in the object.
                --annlevel <annlevel>                  # the annotation level,sample or cell is supported[default: sample].
        export2cbbrowser  #export the object to the specified port of a local preinstalled cellbrowser for visualizing
            options:
                --project <project>                    # name of the dataset. Default to the Seurat object project name.
                --reducts <reducts>                    # the comma seperated redcution names list to export to cellbrowser.
                --markers <markers>                    # the marker genes file
                --clusterid <clusterid>                # name of the metadata field of cell clusters result to use.
                --cb.dir <cb.dir>                      # path to directory to ceate UCSC cellbrowser static website content root.
                                                       # e.g. an index.html, .json file, etc. These file can be copied to any webserver.
                                                       # if this is specified, the cellbrowser package has to be accessible from R
                                                       # via reticulate.
                --port <port>                          # on which port to run UCSC cellbrowser webserver after export.
' -> doc
opt = docopt(doc)

if ( is.null(opt$outdir) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( dir.exists(opt$outdir) ){
        output_dir = opt$outdir
    }else{
        output_dir = opt$outdir
        dir.create(output_dir, recursive = T)
    }
}

# #####################################################################
if ( opt$create ){
    if ( tolower(opt$assay) == "atac" ){
        suppressPackageStartupMessages(library("Signac"))
    }
    if ( opt$aggr == "FALSE" ){
        is.aggr = FALSE
    }else{
        is.aggr = TRUE
    }
    if ( is.null(opt$seperateby) | opt$seperateby == "NULL"){
        seperateby = NULL
    }else{
        seperateby = opt$seperateby
    }
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, informat = opt$informat,
                                assay = opt$assay, metadata = opt$metadata,
                                peak.meta = opt$atac_cellmeta, fragment = opt$fragment,
                                aggr = is.aggr
        )
        if ( !is.null(opt$gtf) ){
            seurat_ob = addFeatureAnno( seurat_ob, assay = opt$assay,gtf = opt$gtf  )
        }
        saveObject(seurat_ob, outdir = output_dir,
                splitby = seperateby, outformat = opt$outformat, prefix = opt$prefix)
    }
}

if ( opt$subset ){ # subset the object using different conditions
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                    informat = opt$informat, metadata = opt$metadata )
    }else{
        stop("NO input data FOUND!")
    }
    group2use = opt$group2use
    levels = opt$levels
    seurat_obx = seurat_ob
    if ( is.null(opt$low) | opt$low == "NULL" ){
        opt$low = -Inf
    }
    if ( is.null(opt$high) | opt$high == "NULL" ){
        opt$high = Inf
    }
    if ( !is.null(levels)){
        levels_id = unlist(strsplit( levels,",",perl = T))
        if ( is.null(group2use) ){ group2use = "idents" }
        if ( group2use == "idents" ){
            if ( typeof(Idents(seurat_obx)) == "integer" ){ # the idents is clustering id
                seurat_obx = SubsetData(seurat_obx, ident.use = levels_id,
                                low.threshold = opt$low, high.threshold = opt$high)
            }
        } else{ #subset cells on metadata using other cell annotations
            seurat_obx = SubsetData(seurat_obx, subset.name = group2use,
                                   accept.value = levels_id,
                                   low.threshold = opt$low, high.threshold = opt$high
                      )
        }
        desired_cells = FetchData(seurat_obx, vars = rownames(seurat_obx)) %>% as.data.frame
    }else{
        desired_cells = FetchData(seurat_obx, vars = rownames(seurat_obx)) %>% as.data.frame
    }
    # subset cells on the expression matrix using logical exprssion on features
    if ( !is.null(opt$cellfilter) ){
        # levels = gsub( ".*\\((.*)\\)", "\\1", opt$cellfilter, perl =T)
        # predicate = gsub( "\\(.*\\)", glue::glue("({levels})"), opt$cellfilter, perl = T)
        predicate = opt$cellfilter

        df = seurat_obx@meta.data
        desired_cells= subset(df, eval( parse(text=predicate)))
        seurat_obx = SubsetData( seurat_obx, cells = rownames(desired_cells) )
    }

    if ( !is.null(opt$feature4subset) ){
        feature_predicate = opt$feature4subset
        count_df_T = Matrix::t(GetAssayData(seurat_obx, assay = opt$assay, slot = opt$slot))
        desired_cells = subset( as.data.frame(count_df_T), eval(expr = parse(text=feature_predicate)))
        seurat_obx = SubsetData( seurat_obx, cells = rownames(desired_cells) )
    }

    desired_cell_id = colnames(seurat_obx)


    # keep only specified features
    if ( !is.null( opt$geneset) ){
        genes = read.table(opt$geneset, sep = "\t", header = T)
        desired_features = as.vector( genes$gene)
    }else{
        desired_features = NULL
    }


    to_invert = opt$invert
    seurat_ob = subset(seurat_ob, features = desired_features,
                        cells = desired_cell_id, invert = to_invert)
    saveObject(seurat_ob, outdir = output_dir, splitby = opt$seperateby,
                outformat = opt$outformat, prefix = opt$prefix)
    quit()
}

# subcommand merge is invoked
if ( opt$merge ){
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                                informat = opt$informat)
        if ( seurat_ob@version < 3 ){
            seurat_ob = UpdateSeuratObject(seurat_ob)
        }
    }else{
        stop("NO input data FOUND!")
    }

    seuy = lapply(strsplit(opt$toadd, ","), function(seu){
        seux = readObject( input = seu, assay = opt$assay, informat = opt$adformat)
        if ( seux@version < 3 ){
            seux = UpdateSeuratObject(seux)
        }
    })
    seurat_ob = merge( seurat_ob, y = seuy, merge.data = T)
    if ( !is.null(levels)){
      id=group2use
      seurat_ob@meta.data[,id]=factor(seurat_ob@meta.data[,id],levels=unique(sort(seurat_ob@meta.data[,id])))
    }
    saveObject(seurat_ob, outdir = output_dir, splitby = opt$seperateby,
            outformat = opt$outformat, prefix = opt$prefix)
    quit()
}


#updata the metedata of cells or features in the seurat_ob
if ( opt$update ){
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                                informat = opt$informat)
    }else{
        stop("NO input data FOUND!")
    }
    if ( !is.null(opt$admeta) ){
        cellnames = Cells(seurat_ob)
        additional_metadata = read.csv(opt$admeta,sep=",",header =T ,row.names=1)
        if ( opt$annlevel == "sample" ){
            #the index order is the same as the row index of the assay metadata
            sampleidx =  gsub("(_|-)[ATGC]{16,}.*","",cellnames,perl=T)
            #integrate the additional metadata from the assay design
            additional_cell_meta = vector()
            for ( colidx in colnames(additional_metadata) ){
                additional_cell_meta = cbind(additional_cell_meta,
                as.vector(additional_metadata[sampleidx, colidx]))
            }
            colnames(additional_cell_meta) = colnames(additional_metadata)
            rownames(additional_cell_meta) = cellnames
            additional_cell_meta = as.data.frame(additional_cell_meta)
            seurat_ob = AddMetaData( seurat_ob, additional_cell_meta)
        }else{# the annotation level is cell, make sure the barcodes are in the same pattern as the ones in the object
            # barcodes =  gsub("-|_.*$","",additional_metadata$barcode,perl=T) # remove suffix of barcodes
            # sampleidx =  gsub("(_|-)[ATGC]{16,18}.*","",additional_metadata$barcode,perl=T) #the index order is the same as the row index of the assay metadata
            # samples = unique(as.vector(seurat_ob$sampleid))
            # sampleid = samples[as.numeric(sampleidx)]
            # additional_metadata$barcode = paste( samples,barcodes, sep = "-")
            # rownames(additional_metadata) = additional_metadata$barcode
            additional_cell_meta = as.data.frame(additional_metadata)
            # TO DO if the cells in the supplied annotation is different from the cells in seurat object
            seurat_ob = AddMetaData( seurat_ob, metadata = additional_metadata)
        }
    }else{
            stop("NO metadata AVAIABLE!")
    }
    saveObject(seurat_ob, outdir = output_dir, splitby = opt$seperateby,
                outformat = opt$outformat, prefix = opt$prefix)
    quit()
}

if ( opt$split ){
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                            informat = opt$informat, metadata = opt$metadata )
    } else {
        stop("NO input data FOUND!")
    }
    if ( !is.null( opt$splitby) ){
        if (opt$outformat == "tenx" ) {
            saveObject( seurat_ob, outdir = output_dir,splitby = opt$splitby,outformat = "tenx")
        } else {
            seurat_list = SplitObject( seurat_ob, split.by = opt$splitby )
            lapply( names(seurat_list), function(x){
                saveObject( seurat_list[[x]], outdir = output_dir,outformat = opt$outformat, prefix = x)})
        }
    } else {
        stop("please specify --splitby when using split subcommand.")
    }
    quit()
}

if ( opt$avgexp ){
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                            informat = opt$informat, metadata = opt$metadata )
    }else{
        stop("NO input data FOUND!")
    }
    if ( !is.null( opt$features) ){
        genes = read.table(opt$features, sep = "\t", header = T)
        desired_features = as.vector( genes$gene)
    }else{
        desired_features = NULL
    }
    if ( !is.null(opt$avgby) ){
        avgby = opt$avgby
    }else{
        avgby = "clusters"
    }
    seurat_ob = SetIdent( seurat_ob, value = avgby)
    avg_mtx = AverageExpression(seurat_ob, assays = opt$assay, features = desired_features, slot = opt$slot)
    write.table(avg_mtx, file.path(output_dir,paste0("average_expression_matrix_by_",avgby,".xls")),
                sep = "\t", col.names = T,row.names = T)
    # output metadata either
}

if ( opt$removeEmptyDrops ){
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                            informat = opt$informat, metadata = opt$metadata )
    }else{
        stop("NO input data FOUND!")
    }

    # remove the empty droplets using dropletuitls if necessary
    emptyRemoved = RemoveEmptyDrops(
                  fdrThreshold=as.numeric(opt$fdr),
                  emptyDropNIters= as.numeric(opt$iters),
                  raw_count = GetAssayData(seurat_ob, slot = "counts") )
    subset_se = subset( seurat_ob, cells = colnames(emptyRemoved) )
    # TO DO 
    # make a UMI rank plot
    saveObject(seurat_ob, outdir = output_dir,
            splitby = seperateby, outformat = opt$outformat, prefix = opt$prefix)
    quit()
}

# Retreives data (feature expression, PCA scores, metrics, etc.) for a set of cells in a Seurat object
if ( opt$fetch ){ # subcommand fetch is invoked
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                                informat = opt$informat )
    }else{
        stop("NO input data FOUND!")
    }

    if ( is.null(opt$slot) ){
        slot = "counts"
    }else{
        slot = opt$slot
    }
    if ( !is.null(opt$vars) ){
        var_list = unlist( strsplit(opt$vars, ",", perl = T) )
    }else{
        var_list = c("group", "ident")
    }
    desired_df = FetchData(seurat_ob, slot = opt$slot, vars = var_list )
    desired_df = desired_df %>% rownames_to_column(var = "CellBarcode")
    if ( !is.null(opt$header) ){
        colnames(desired_df) = unlist( strsplit(opt$header, ",", perl = T) )
    }
    write.table(desired_df, file.path(output_dir,
                glue::glue("{opt$prefix}_cellinfo.xls")),
                row.names = F,col.names = T, sep = "\t")
    quit()

}

if ( opt$transform ){
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
        informat = opt$informat, metadata = opt$metadata )
    }else{
        stop("NO input data FOUND!")
    }
    # seurat_ob = RenameCells(seurat_ob, new.names = barcodes )
    seurat_ob = RenameCells()

}

if ( opt$downsample ){  #### TO DO
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                            informat = opt$informat, metadata = opt$metadata )
    }else{
        stop("NO input data FOUND!")
    }
    if (  tolower(opt$filters) == "null" ) filters = NULL
    if (  tolower(opt$doinvert) == "true" ) {
        invert = TRUE
    }else{
        invert = FALSE
    }
    if (  tolower(opt$targetN) == "null" ) {
        targetN = 10000
    }else{
        targetN = as.numeric(opt$targetN)
    }

    if ( !is.null( filters ) ){
        target_cells = rownames(subset(seurat_ob@meta.data, eval( parse(text= filters))))
        if ( invert ) {
            target_cells = setdiff( Cells(seurat_ob), target_cells )
        }
        targetN = min(length(target_cells), targetN)
        down_refcells = Cells(DownSample(SubsetData( seurat_ob, cells = target_cells ),
                             N = targetN, slot = "data", verbose = F))
        malignant_cells = setdiff( Cells(seurat_ob), target_cells )
        down_seurat = subset(seurat_ob, cells = c( down_refcells, malignant_cells))
    }else{
        targetN = min(dim(seurat_ob)[2], targetN)
        down_seurat = DownSample(seurat_ob, N = targetN, slot = "data", verbose = F)
    }

    saveObject(down_seurat, outdir = output_dir,
            splitby = seperateby, outformat = opt$outformat, prefix = opt$prefix)
    quit()
}

if ( opt$export2cbbrowser ){ ### TO DO
    if ( !is.null(opt$informat) & !is.null(opt$input) ){
        seurat_ob = readObject( input = opt$input, assay = opt$assay,
                                informat = opt$informat, metadata = opt$metadata )
    }else{
        stop("NO input data FOUND!")
    }

    ExportToCellbrowser( object, dir = output_dir,
                dataset.name = ifelse(is.null(opt$project),names((object)), opt$project),
                reductions = ifelse(is.null(opt$reducts),names(object@reductions), opt$reducts ),#it's not a transcriptable writing
                markers.file = NULL,
                cluster.field = "clusters",
                port = ifelse( is.null(opt$port), "127.0.0.1", opt$port ),
                cb.dir = ifelse( is.null(opt$cb.dir), NULL, opt$cb.dir) )
}
