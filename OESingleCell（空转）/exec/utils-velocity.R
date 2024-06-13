# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2021/8/11
docstring<- "example:\\n\\n\\
    sctool  -i singlecell_object.clustering_resolution0.4.rds  -f rds    -o velocity_results   -d rds  velocity  -l SG_1L1D.loom,SG_E8D.loom,SG_1LM.loom -s 0.5   -r umap -- mnn "
sub_velocity<- subparsers$add_parser(
    "velocity",
    description = docstring,
    formatter_class= 'argparse.RawTextHelpFormatter' ,
    #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
    argument_default = "True",
    help = "plot for visium celltype parse results")
sub_velocity$add_argument(
    "-l",
    "--loom",
    type = "character",
    default = NULL,
    help = "[REQUIRED]The comma seperated loom file list from run velocyto.[default: %(default)s] " )
sub_velocity$add_argument(
    "-s",
    "--pointsize",
    type = "double",
    default = 0.5,
    help = "[REQUIRED][optional]the point size in the plot.[default: %(default)s] " )
sub_velocity$add_argument(
    "--metadata",
    "-m" ,
    type = "character",
    default = NULL,
    help = "[optional]the sample metadata which must include sample id in this assay design.")
sub_velocity$add_argument(
    "--reduct1",
    type = "character",
    default = "pca",
    help = "[optional] the (primary) reduction methods used, here is to calculate cell distance.[default: %(default)s]")
sub_velocity$add_argument(
    "--reduct.vis",
    "-r",
    type = "character",
    default = "tsne",
    help = "[optional]the results of computed reduction methods used as input for secondary reduction. choice can be ica,cca,pca,rpca,mnn,tsne,umap etc.[default: %(default)s] ")
sub_velocity$add_argument(
    "--ident2use",
    "-q" ,
    type = "character",
    default = NULL,
    help = "[optional]the column name in cell metadata used as identity of each cell combined with which_cell.[default: %(default)s] ")
sub_velocity$add_argument(
    "--which_cells",
    "-u",
    type = "character",
    default = NULL,
    help = "[optional]the subset of cluster ids used for subtyping.[default: %(default)s] ")
sub_velocity$add_argument(
    "--replot" ,
    type = "logical",
    default = FALSE,
    help = "[optional] replot or not .[default: %(default)s] ")
sub_velocity$add_argument(
    "--color" ,
    type = "character",
    default = "clusters",
    help = "[optional] the metadata colname to color .[default: %(default)s] ")
#=======================================================================================================================
args <- commandArgs(TRUE)
if ( "velocity"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name == "velocity"){
    # suppressPackageStartupMessages(library("SeuratWrappers"))
    # suppressPackageStartupMessages(library("velocyto.R"))
    #===================================================================================================================
    futile.logger::flog.info("step1:read the specified assay and data slot in data object into memory")
    seurat_ob <- OESingleCell::ReadX(input = opt$input,
                                      informat = opt$informat,
                                      assays = assays[1],
                                      data.use = dataslots,
                                      verbose = F)
    Seurat::DefaultAssay(seurat_ob) <- assays[1]
   # futile.logger::flog.logger(glue::glue("Showning reductions information in input object:{reduct}",
    #                                      reduct=sapply(names(seurat_ob@reductions, paste0, collapse=""))))
    #===================================================================================================================
    futile.logger::flog.info("step2:check reduction information")
    if ( !opt$reduct1 %in% names(seurat_ob@reductions) ) stop(paste0("Primary reduction ",opt$reduct1," is not present in the imput file." ))
    if ( !opt$reduct.vis %in% names(seurat_ob@reductions) ) stop(paste0("Secondary reduction ",opt$reduct.use," is not present in the imput file." ))

    #===================================================================================================================
    futile.logger::flog.info("step3:update the metedata in the seurat_ob@meta.data with new additional sample metadata")
    if ( !is.null(opt$metadata) ){
        additional_metadata <- read.csv(opt$metadata, sep=",", header =T )
        rownames(additional_metadata) <- additional_metadata$sampleid
        cellnames <- Cells(seurat_ob)
        #make the index order is the same as the row index of the assay metadata
        eampleidx <- gsub("_|-[ATGC]{16,}.*", "", cellnames, perl=T)
        #integrate the additional metadata from the assay design
        additional_cell_meta <- vector( )
        for ( colidx in colnames(additional_metadata) ){
            additional_cell_meta <- cbind(additional_cell_meta, as.vector(additional_metadata[sampleidx, colidx]))
        }
        colnames(additional_cell_meta) <- colnames(additional_metadata)
        rownames(additional_cell_meta) <- cellnames
        additional_cell_meta <- as.data.frame(additional_cell_meta)
        seurat_ob <- Seurat::AddMetaData(seurat_ob, metadata = additional_cell_meta)
    }
    #===================================================================================================================
    futile.logger::flog.info("step3:get the subset of cells used for visualization if necessay")
    if ( !is.null(opt$which_cells)){
        if ( is.null(opt$ident2use ) ){
            print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
            ident2use <- "clusters"
        }else{
            ident2use <- opt$ident2use
        }
        cluster_list <- unlist(strsplit(opt$which_cells, ",", perl = T))
        #cells used as interested cells
        expr <- Seurat::FetchData(object = seurat_ob, vars = ident2use)%>%
                tibble::rownames_to_column("barcodes") %>%
                dplyr::filter( .data[[ident2use]]  %in% cluster_list )
        seurat_ob <-subset(seurat_ob,cells = expr$barcodes)
        # seurat_ob = subset(seurat_ob,
        # seurat_ob = subset(seurat_ob,
        #                    cells == Seurat::WhichCells( object, subset.name= ident2use, accept.value = cluster_list))
    }
    #===================================================================================================================
    futile.logger::flog.info("step4:plotting")
    if (opt$replot == F) {
        # This is generated from the Velocyto python command line tool.
        # You need a loom file before you can proceed
        # ldat <- read.loom.matrices( opt$loom )
        ldat_seurat_list <- lapply(unlist(strsplit(opt$loom, ",")), function(loomx){
            ldat <- SeuratWrappers::ReadVelocity(loomx)
            ldat_seurat <- Seurat::as.Seurat(ldat)
            # This is a little bit of foo-magic that needs to be adjusted on a per-sample
            # basis depending on the cell names and how you ran the pipeline. Each cell
            # stored in the loom object and seurat have an ID, make sure these are the same.
            # What this step does is essentially this:
            # > head(Cells(ldat_seurat))
            # [1] "possorted_genome_bam_XL2S3:AAAGATGCATACTACGx"
            # [2] "possorted_genome_bam_XL2S3:ACCTTTATCTTTAGTCx"
            # [3] "possorted_genome_bam_XL2S3:AAGGAGCCACGCATCGx"

            cellids <- gsub(":", "-", Seurat::Cells(ldat_seurat))##或者是_
            cellids <- gsub("x$", "", cellids) # trim the suffix "x" at the end of the cell names
            # Now the names in emat and nmat will match up to the cell names used in my seurat object
            # > head(Cells(ldat_seurat))
            # [1] "possorted_genome_bam_XL2S3_AAAGATGCATACTACG"
            # [2] "possorted_genome_bam_XL2S3_ACCTTTATCTTTAGTC"
            # [3] "possorted_genome_bam_XL2S3_AAGGAGCCACGCATCG"
            ldat_seurat <- Seurat::RenameCells(ldat_seurat, new.names = cellids)
        })
        if (length(ldat_seurat_list) >1 ) {
            ldat_seurat <- merge(ldat_seurat_list[[1]],
                                                 y = ldat_seurat_list[2:length(ldat_seurat_list)])
                                                 #merge.data = T)
        }else {
            ldat_seurat <- ldat_seurat_list[[1]]
        }

        ldat_seurat <- subset(ldat_seurat, cells = Seurat::Cells(seurat_ob))
        seurat_ob[["spliced"]] <- ldat_seurat[["spliced"]]
        seurat_ob[["unspliced"]] <- ldat_seurat[["unspliced"]]

        # Gather the spliced and unspliced estimates
        # emat <- ldat$spliced
        # nmat <- ldat$unspliced
        seurat_ob <- SeuratWrappers::RunVelocity(seurat_ob,
                                                 deltaT = 1,
                                                 kCells = 25,
                                                 fit.quantile = 0.02,
                                                 reduction = opt$reduct1)
        #saveRDSMC(seurat_ob, file.path(output_dir, "seurat_Velocity.rds"))
        ##save results to seurat object's mis part
        if ( tolower(opt$outformat) == "h5seurat" & as.logical(opt$update) ){
            SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob )
         }else{
            OESingleCell::SaveX(seurat_ob,
                                output = opt$output,
                                update = FALSE,
                                outformat = opt$outformat,
                                prefix = opt$prefix)
        }
        # Estimate the cell-cell distances
        # cell.dist <- as.dist(1-armaCor(t( emb )))

        # Main velocity estimation
        # rvel.cd <- gene.relative.velocity.estimates(emat,nmat, deltaT=2,
        #                                             kCells=10, cell.dist=cell.dist,
        #                                             fit.quantile=fit.quantile, n.cores=24)

        # This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
        # gg <- TSNEPlot(object)
        # ggplot_build(gg)$data
        # colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
        # names(colors) <- rownames(emb)
        # ident.colors <- (scales::hue_pal())(n = length(x = levels(x = seurat_ob)))
         Seurat::Idents(seurat_ob) <- seurat_ob[[opt$color]]
         ident.colors <- OESingleCell::SelectColors(1:length(unique(Seurat::Idents(seurat_ob))),palette ="customecol2" )
         names(ident.colors) <- levels(seurat_ob)
         cell.colors <- ident.colors[Seurat::Idents(seurat_ob)]
         names(cell.colors) <- colnames(seurat_ob)

        # take embedding from the Seurat data object
        # NOTE: This assumes you have a seurat data object loaded
        # into R memory prior to using this script. STOP and rerun seurat
        # pipeline if you do not have this loaded. In my case, my seurat object is simply myData
        emb <- Seurat::Embeddings(seurat_ob, reduction = opt$reduct.vis)
        arrow.scale <- (range(emb[, 1])[2]-range(emb[, 1])[1])/40*4 # grid.n=40. Here the max arrow.scale is twice the length of one grid cell
        pdf(file.path(output_dir,"velocity_plot.pdf"), width = 10, height = 10)
        velocyto.R::show.velocity.on.embedding.cor(emb,
                                                   vel = Seurat::Tool(object = seurat_ob, slot = "SeuratWrappers::RunVelocity"),
                                                   n = 200,
                                                   scale = "sqrt",
                                                   min.arrow.size = 0,
                                                   cell.colors = velocyto.R::ac(x = cell.colors, alpha = 0.8),
                                                   cex = opt$pointsize ,
                                                   arrow.scale = arrow.scale,
                                                   show.grid.flow = TRUE,
                                                   min.grid.cell.mass = 0.5,
                                                   grid.n = 40,
                                                   arrow.lwd = 1,
                                                   do.par = F,
                                                   n.cores=1,
                                                   cell.border.alpha = 0)
        dev.off() }
      ## save session information
      write_session_info(output_dir,sub_name = parser$parse_args()$sub_name )
      quit()
  }
}