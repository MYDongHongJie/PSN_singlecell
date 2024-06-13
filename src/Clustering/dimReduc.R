#!/usr/bin/env Rscript

# this is the script used only for data matrix dimension reduction
# it was intened for single cell data analysis mainly, but in the 
# future it will be applied to any data matrix.
# suppressPackageStartupMessages( library("scRobustPCA"))
suppressWarnings({
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("Matrix") )
    suppressPackageStartupMessages( library("optparse") )
    suppressPackageStartupMessages( library("OESingleCell") )
    suppressPackageStartupMessages( library("BiocParallel") )
    suppressPackageStartupMessages( library("Signac") )
})

#=================================================================================
# function defintion
#=================================================================================
# do dimension reduction using the specified method
RunDimReduc <- function(
  object,
  feature.use = NULL,
  reduct2.use = NULL,
  reduct1.use = NULL,
  perplexity = 30,
  batch.use = NULL,
  npcs.use = NULL,
  assay.use = NULL ,
  ...
){

    if ( !is.null(assay.use) ){
        DefaultAssay(object) = assay.use
    }
  if ( is.null( feature.use ) ){
    feature.use = VariableFeatures(object)
  }
  
  npcs.use = ifelse(is.null(npcs.use),yes = 30,no = npcs.use)
  
  if ( tolower(reduct2.use) %in% c("cca", "mnn", "harmony") ){
      if ( is.null( batch.use) ){
        warning("NO batch information specified!The batchid will be used as default!")
        batch.use = "batchid"
      }
      # check the batchid in the metadata column names
      if ( ! batch.use %in% colnames(object@meta.data) ){
        stop("NO specified batch column found in the meta.data slot of seurat object!")
      }
  }

    if ( is.null(reduct1.use) ){
        print("NO primary reduction method specified, PCA will be used as default!")
        reduct1.use = "pca"
    }

    if ( is.null(reduct2.use) ){
        reduct2.use = reduct1.use
    }
  
  object = switch (tolower(reduct2.use),
                      'ica' = RunICA(object, nics = npcs.use, 
                                     features = feature.use, 
                                     verbose = F, ...),
                      'cca' = do.call("RunCCA",c(SplitObject(object, split.by = batch.use ),
                                        list(features = feature.use, 
                                             renormalize =F, rescale = F ))),
                      'pca' = RunPCA(object, npcs = npcs.use,
                                     features = feature.use, 
                                     do.print = F, verbose = F, ...),
                      'lsi' = RunSVD( object, n = npcs.use, features = feature.use,
                                     reduction.key = 'LSI_', reduction.name = 'lsi',
                                     verbose = F, ...),
                      # 'lsi' = RunLSI(object, n = npcs.use,
                      #                 features = feature.use,
                      #                 do.print = F, verbose = F, ...),
                      'lda' = RunLDA(object,  slot = "counts",
                                      topics = 30, method = "Z-score"),
                      'lsa' = RunLSA(object, npcs = npcs.use ),
                      'tsne' = RunTSNE(object, reduction = reduct1.use,
                                       dim.embed = 3, tsne.method = "Rtsne",
                                      features = feature.use, perplexity = perplexity,
                                      dims = 1:npcs.use, num_threads = 10,
                                      max_iter = 2000, check_duplicates = F, ... ),
                      'flt-sne' = RunTSNE(object, reduction = reduct1.use,
                                       dim.embed = 3, tsne.method = "FIt-SNE",
                                      features = feature.use, dims = 1:npcs.use,
                                       fasttsne_path = system("which fast_tsne"),
                                       nthreads = 10, max_iter = 2000, 
                                      check_duplicates = F, ... ),
                      'rpca' = RunRobPCA(object, pc.genes = feature.use,
                                         use.modified.pcscores = T,
                                          do.print = F, npcs = npcs.use, ...),
                      'svd'  = RunSVD( object, n = npcs.use, features = feature.use,
                                        verbose = F, ...),
                      'swne' = RunSWNE(object, reduction.use = reduct1.use,
                                       dims.use = 1:npcs.use,
                                       proj.method = "sammon", 
                                       batch = as.vector(object[[batch.use]]),
                                       var.genes = feature.use,
                                       ncores = 20, return.format = "seurat"),
                      'harmony' = RunHarmony(object, batch.use, 
                                            theta = 2, plot_convergence = F, nclust = 50,
                                            max.iter.cluster = 20, max.iter.harmony = 5,
                                            ...),
                      'mnn' = RunMnn.Seurat(object, features = feature.use, assay = assay.use,
                                        batch = batch.use, ...),
                      'umap' = RunUMAP(object,dims = 1:npcs.use,verbose = F,
                                      reduction = reduct1.use, n.components = 3,...),
                                      #features = feature.use, ...),
                      # Kevin R. Moon et.al(2019).Visualizing Structure and Transitions for Biological Data Exploration
                      'phate' = RunPHATE.Seurat(object, assay = assay.use,
                                        features = feature.use, npcs = npcs.use, ...),
                      # Mathieu Jacomy et.al(2014). ForceAtlas2, a Continuous Graph Layout Algorithm for Handy Network
                      # Visualization Designed for the Gephi Software
                      'fa2' = RunForceAtlas2.Seurat(object, reduction = reduct1.use,
                                         features = feature.use, npcs = 2:npcs.use, ...),
                      'diffusion' = RunDiffusion.Seurat(object, features = feature.use,
                                                     dims = 1:npcs.use,
                                                     reduction = reduct1.use, ...) )

  message(paste0(reduct2.use, " Dimension Reduction finished"))
  return( object )
}

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
                 help = "The input exprssion matrix in several possible format."),
    make_option( c("--informat", "-f" ), type = "character", default = "seurat",
                 help = "The indication of type of input expression matrix, the possible type can be:
                        seurat: the seurat object from the clustering results."),
    make_option( c("--components", "-t"), type="integer", default=20,
        help="the appropriate number of statistically significant components to use for clustering,
                 which derived from the JackStraw result."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of Clustering results." ),
    make_option( c("--metadata", "-m" ), type="character", default = NULL,
        help="[OPTIONAL]the additional sample metadata which must include sample id in this assay design."),
    make_option( c("--batchid", "-b"), type = "character", default = NULL,
         help = "[OPTIONAL]the batch information column name of sample in the metadata."),
    make_option( c("--perplexity", "-p"), type="integer", default=30,
        help="The value of the perplexity used for tSNE"),
    make_option( c("--assay" ), type="character", default = "RNA",
        help="[OPTIONAL]the assay used to calulation in case of multimodal data."),
    make_option( c("--forcerun" ), type="logical", default = FALSE,
    help="[OPTIONAL]Rerun the dimension reduction regardless of previous run."),
    make_option( c("--reduct1.use","-d"),type = "character",default=NULL,
         help="[OPTIONAL]the primary reduction methods whose results can be used as input for secondary reduction.
                The supported methods can be one of the followings:
                  ica, Independent Compononet Analysis(ICA).
                  cca, Canonical Correlation Analysis(CCA), can be used to remove batch effect.
                  pca, principal component analysis(PCA).
                  rpca, robust PCA, an improved PCA.
                  lsi, Latent Semantic Indexing(LSI) on binary count matrix, mianly used for scATAC-seq.
                  lsa, Latent Semantic Analysis,mianly used for scATAC-seq.
                  svd, Singular Value Decomposition
                  mnn, Mutual Nearest Neighbours Analysis(MNN) after PCA, as a dimension reduction method with batch removing.
                  harmony,Scalable integration of single cell RNAseq data for batch correction and meta analysis,
                  swne,
                  fa2,Forceatlas2, an force layout method used in the method used in the SPRING.
                  tsne,t-Stochastic Neighborhood Embedding
                  flt-tsne, FFT-accelerated Interpolation-based t-SNE
                  umap, Uniform Manifold Approximation and Projection.
                  phate,
                  diffusion, Diffusion map."),
    make_option( c("--reduct2.use","-D"),type = "character",default="tsne",
         help="one or more of the scondary dimension reduction methods used for this run.
                the options can be tsne, flt-sne, umap, duffusion, swne. That is to say not all the primary redcution
                method can used as scondary one."),
    make_option( c("--ident2use", "-q" ), type = "character", default = NULL,
         help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
            help = "[OPTIONAL]The subset of cluster ids used for subtyping."),
    make_option( c("--overwrite" ), type = "logical", default = TRUE,
         help = "[OPTIONAL]overwrite the input seurat object in situe or output in the output directory.[default: TRUE]")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# =================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$reduct1.use )){
  print("NO first level dimension reduction methods specified,the default PCA will be used!")
  reduct1.use = "pca"
}else{
  reduct1.use = tolower(opt$reduct1.use)
}

if ( !is.null(opt$reduct2.use) ){
    reduct2.use = tolower(unlist(strsplit(opt$reduct2.use,",",  perl = T)))
}else{
    # if no secondary reduction method specified, use primary reduction instead
    reduct2.use = reduct1.use
}

if ( !is.null(opt$batchid) ){
  batchid = opt$batchid
}

if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir)
    }
}
output_dir = normalizePath(output_dir )

if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( tolower(opt$informat) == "seurat" ){
        # the input is a seurat object produced by previous analysis
        seurat_ob = readRDSMC( opt$input )
        # if the input seurat object version is less than 3, upgrade it to version 3
        if ( seurat_ob@version < 3 ){
            seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
        }
    }else if ( tolower(opt$informat)  == "sce" ){
      # the input is a SingleCellExperiment object
        sc_object = readRDSMC( opt$input )
        seurat_ob = as.Seurat(sc_object)
    }

    #check the value of perplexity
    suppressWarnings({
        if ( is.null(opt$perplexity) ){
            if ( is.null(Misc(seurat_ob, "perplexity")) ){
                print( "NO previous perplexity is AVAILABLE, 30 will be used as default.")
                Misc(seurat_ob, "perplexity") = 30
            }else{
                perplexity = Misc(seurat_ob, "perplexity")
            }
        }else{
            perplexity = opt$perplexity
            Misc(seurat_ob, "perplexity") = opt$perplexity
        }
    })

    suppressWarnings({
        #check the components to use
        if ( is.null(opt$components) ){
            if ( is.null(Misc(seurat_ob, "components_num")) ){
                print( "NO previous components calculation is AVAILABLE, 20 will be used as default.")
                Misc(seurat_ob, "components_num") = 20
            }else{
                component2use = Misc(seurat_ob, "components_num")
            }
        }else{
            Misc(seurat_ob, "components_num") = opt$components
            component2use = opt$components
        }
    })
}


# change the default assay for reduction if necessary
if ( !is.null( opt$assay) ){
    DefaultAssay(seurat_ob) = opt$assay
}else{
    DefaultAssay(seurat_ob) = "RNA"
}

#updata the metedata in the seurat_ob@meta.data with new additional sample metadata
if ( !is.null(opt$metadata) ){
    additional_metadata = read.csv(opt$metadata,sep=",",header =T )
    rownames(additional_metadata) = additional_metadata$sampleid

    cellnames = Cells(seurat_ob)
    sampleidx =  gsub("_[ATGC]{16,}.*","",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
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
}

#get the subset of cells used for visualization if necessay
if ( !is.null(opt$which_cells)){
    if ( is.null(opt$ident2use ) ){
        print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
        ident2use = "clusters"
    }else{
        ident2use = opt$ident2use
    }
    cluster_list = unlist(strsplit(opt$which_cells,",", perl = T))
    #cells used as interested cells
    seurat_ob = SubsetData(seurat_ob,cells = 
                       OldWhichCells( seurat_ob, subset.name= ident2use, accept.value = cluster_list))
}

# =================================================================================
# dimension reduction 
# =================================================================================
# check out wether the specified primary reduction has beed calculated
if ( !is.null(reduct1.use )){ # primary reduction is specified
  # if no previous reduction found, catch the error and rerun the
  # reduction without stop with error
    if ( !reduct1.use %in% names(Key(seurat_ob)) | opt$forcerun ){
        print( "NO specified primary reduction found in the object! Rerun begins!")
        message(paste0("Beginning ", reduct1.use, " Dimension Reduction"))
        dim_outdir = file.path(output_dir,paste0(reduct1.use, "_Dimension_Reduction"))
        if ( !dir.exists(dim_outdir) ){
            dir.create(dim_outdir)
        }

        seurat_ob = RunDimReduc(seurat_ob, reduct1.use = reduct1.use ,
                                reduct2.use = reduct1.use,
                                feature.use = VariableFeatures(seurat_ob),
                                perplexity = perplexity,
                                assay.use = DefaultAssay(seurat_ob),
                                batch.use = batchid, npcs.use = component2use )
        # pdf(file.path(dim_outdir, "genes_highly_associated_with_each_dimension_heatmap.pdf"), width = 10, height = 8)
        # for( idx in 1:min(15,ncol(Embeddings(seurat_ob, reduction = reduct1.use))) ){
        #     DimHeatmap(object = seurat_ob, reduction = reduct1.use, dims = idx, cells = 500 )
        # }
        # dev.off()
        reduct1_coord = FetchData(seurat_ob,
                            var = c("orig.ident", paste0( Key(seurat_ob)[reduct1], 1:2))) %>%
                            dplyr::rename( "Barcode" = "orig.ident")
        write.table( reduct1_coord, file.path(dim_outdir, paste0(reduct1.use, "_Dimension_Reduction_coordination.csv")),
                    sep = ",", col.names = T, row.names = F, quote = F)
    }else{
        print("The specified reduction results have beed calculated!Skip it!")
    }
}

# different reduction method from primary reduction specified or
# forced to rerun
if ( !reduct2.use %in% names(Key(seurat_ob)) | (reduct2.use != reduct1.use & opt$forcerun) ){
    for ( reduction2 in reduct2.use ){
        message(paste0("Beginning ", reduction2, " Dimension Reduction"))
        output_dir = file.path(output_dir,paste0(reduction2, "_Dimension_Reduction"))
        if ( !dir.exists(output_dir) ){
            dir.create(output_dir)
        }
        seurat_ob = RunDimReduc(seurat_ob, reduct1.use = reduct1.use , reduct2.use = reduction2,
                                feature.use = VariableFeatures(seurat_ob),
                                assay.use = DefaultAssay(seurat_ob),
                                batch.use = batchid, npcs.use = component2use)
        reduct1_coord = FetchData(seurat_ob,
                            var = c("orig.ident", paste0( Key(seurat_ob)[reduction2], 1:2))) %>%
                            dplyr::rename( "Barcode" = "orig.ident")
        write.table( reduct1_coord, file.path(dim_outdir, paste0(reduction2, "_Dimension_Reduction_coordination.csv")),
                    sep = ",", col.names = T, row.names = F, quote = F)
    }
}

if ( opt$overwrite ){
    saveRDSMC(seurat_ob,opt$input )
}else{
    saveRDSMC( seurat_ob, file.path( output_dir, basename(opt$input)) )
}
