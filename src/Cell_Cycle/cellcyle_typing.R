#!/usr/bin/env Rscript
rm(list=ls())
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("pheatmap") )
suppressPackageStartupMessages( library("dplyr") )
suppressPackageStartupMessages( library("OESingleCell") )
suppressPackageStartupMessages( library("ggplot2") )
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

#=command line parameters setting=============================
option_list = list(
    make_option( c("--INEXPRESS", "-i" ), type = "character",
                 help = "The input exprssion matrix in several possible format."),
    make_option( c("--INFORMAT", "-f" ), type = "character", default = "tenx",
                 help = "The indication of type of input expression matrix, the possible type can be:
                        seurat: the seurat object from the clustering results."),
    make_option( c("--species","-s"),type="character",
        help="the species of the sample" ),
    make_option( c("--splitby","-b"),type="character",default = "sampleid",
        help="[OPTIONAL]visualize cells in seperate plot split by this groupping variable" ),
    make_option( c("--bartype","-t"),type="character",default = "stack",
        help="[OPTIONAL]barplot type is stack or dodge" ),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of Clustering results." ),
    make_option( c("--method","-m"),type="character", default = NULL,
        help="the method use to calculate cellcyle." ),
    make_option( c("--ident2use", "-q" ), type = "character", default = NULL,
        help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
        help = "[OPTIONAL]The subset of levels used for subtyping."),
    make_option( c("--proportion", "-v"), type = "double", default = NULL,
                 help = "the proportion of cells used to subsample from each cluster."),
    make_option( c("--pt.sze"), type = "character", default = 1,
                help = "the point size setting in the dimension reduction plot."),
    make_option( c("--cycle", "-k" ), type = "character",
         help = "the cell cycle gene markers for human and mouse"),
    make_option( c("--reduct"), type = "character", default = "tsne",
        help = "the dimension reduction result used to visualize the cells"),
    make_option( c("--groupby", "-g"), type = "character", default = "clusters",
                 help = "groupping cells on heatmap.")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================

if ( is.null(opt$method) ){
    warning("NO methods provides, dropbead will be used." )
    method = "dropbead"
}else{
    if (!opt$method %in% c("dropbead","scran")) stop("Please select a method from one of the followings: dropbead, scran.")
    method = opt$method
}

cellcycle_name=paste0(method,"_CellCycle")
groupbydesign = opt$groupbydesign
if ( is.null(opt$species) ){
  stop("NO SPECIE SPECIFIED!" )
}else{
  spe = tolower(opt$species)
}

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

if ( !is.null( opt$splitby) & opt$splitby != "NULL" ){
    facetby = unlist( strsplit( opt$splitby, ",", perl =T ) )
}else{
    facetby = NULL
}

reduct = opt$reduct
if ( is.null( opt$groupby ) ){
    print("The groupping information is not specified. The sampleid will be used for plot")
    cellfeature2vis = "clusters"
}else{
    cellfeature2vis = unlist(strsplit(opt$groupby, ",", perl = T))
}



output_dir = normalizePath(output_dir )

# ####################################################################
#read in the 10X data in different format
# ####################################################################
#whatever format the input file has, all the input expression matrix
#should be transformed to the seurat object finally.
if ( !is.null(opt$INEXPRESS) ){
    # the input is a seurat object which may contain more than one sample
    if ( opt$INFORMAT == "seurat" ){
        seurat_ob = readRDSMC( opt$INEXPRESS )

        findcluster_record = Command(seurat_ob, command = "FindClusters")
        resolution = findcluster_record$resolution

        if ( is.null(seurat_ob@meta.data$clusters) ){
            seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
        }else{
            # if it is the first time to run this script on a seurat object, the
            # clusters here is actually the sample index but cell cluster ids.
            # After running this script, the cluster id will be overwrite with
            # the actual cell cluster id.
            seurat_ob = SetIdent( seurat_ob, value = "clusters")
        }

        #get the subset of cells used for visualization
        if ( !is.null(opt$which_cells)){
            cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
            seurat_ob = SubsetData(seurat_ob,
                                    cells = OldWhichCells(seurat_ob,
                                                subset.name = opt$ident2use,
                                                accept.value = cluster_list ))
        }
    }
}

if (dim(seurat_ob)[2] < 500){
    pointsize=1.5
}else{
    pointsize=as.numeric(opt$pt.sze)
}

# ####################################################################
## Cell cycle phase with Dropbead
## Point towards cell cycle gene list supplied with Dropbead
# ####################################################################
if (method=="dropbead") {
    suppressPackageStartupMessages( library("dropbead") )
    if ( is.null(opt$cycle) ){
      stop("NO cell cycle marker gene file is AVAILABLE!")
    }else{
      cell_cycle_marker = opt$cycle
    }
    # random subsample the cells in each cluster for heatmap visualization in case of big data
    if ( !is.null(opt$proportion) ){
        sampled_cells = seurat_ob@meta.data %>% mutate( barcode = rownames(seurat_ob@meta.data)) %>%
            group_by(clusters) %>% sample_frac( size = opt$proportion, replace=F)
        seurat_subset = SubsetData( seurat_ob, cells = sampled_cells$barcode )
    
        ## Create a new single species object
        SingleSS_ob <- new("SingleSpeciesSample",
                            species1= spe, ## Change to mouse if working with mouse data!
                            cells= rownames( seurat_subset@meta.data ),
                            genes=rownames(GetAssayData(seurat_subset,slot="data")),
                            dge=as.data.frame(as.matrix(GetAssayData(seurat_subset,slot="data"))))
    }else{
        ## Create a new single species object
        SingleSS_ob <- new("SingleSpeciesSample",
                        species1= spe, ## Change to mouse if working with mouse data!
                        cells= rownames( seurat_ob@meta.data ),
                        genes=rownames(GetAssayData(seurat_ob,slot="data")),
                        dge=as.data.frame(as.matrix(GetAssayData(seurat_ob,slot="data"))))
    }

    phases <- assignCellCyclePhases(SingleSS_ob,gene.file= cell_cycle_marker, do.plot = F) 
    if ( opt$groupby %in% names(seurat_ob@meta.data)) {
        # order1 get original annotation 
        cells=SingleSS_ob@cells[as.integer(rownames(phases))]
        anno=seurat_ob@meta.data[cells,opt$groupby]
        # reorder 
        phases_reorder=phases[order(anno),]
        cells=SingleSS_ob@cells[as.integer(rownames(phases_reorder))]
        anno=as.character(seurat_ob@meta.data[cells,opt$groupby])
        
        annotation_col=data.frame(clusters=anno,row.names=rownames(phases_reorder))
        
        col=CustomCol2(1:length(unique(seurat_ob@meta.data[,opt$groupby])))

        #names(col)=as.character(1:length(unique(seurat_ob@meta.data[,opt$groupby])))
        names(col) <- unique(annotation_col$clusters)
        annotation_colors=list(clusters=col)
        names(annotation_colors) = opt$groupby
        names(annotation_col) = opt$groupby
        # visualize on heatmap 
        pdf(file.path(output_dir,paste0("cell_cycle_annotation_heatmap_orderedby_",opt$groupby,".pdf")),width=20,height = 20)
        pheatmap::pheatmap(as.matrix(t(phases_reorder)),cluster_row=F,cluster_col=F,fontsize = 20,
                           color = colorRampPalette(c("#3794bf", "#FFFFFF", "#cc4140"))(100),
                           show_colnames = F,annotation_col=annotation_col,annotation_colors=annotation_colors)
        dev.off()
    } else {
       warning("groupby not found in metadata")
    }
    pdf(file.path(output_dir,paste0("cell_cycle_annotation_heatmap.pdf")),width=20,height = 20)
    pheatmap::pheatmap(as.matrix(t(phases)),treeheight_row = 0,treeheight_col = 0,fontsize = 20,
                       color = colorRampPalette(c("#3794bf", "#FFFFFF", "#cc4140"))(100),
                       show_colnames = F)
    dev.off()
    sorted_phases <- phases %>% arrange(rownames(phases))
    ## For each cell, check which phase has highest score and assign this to the cell
    ## Code by akrun (https://stackoverflow.com/questions/36274867/getting-column-names-for-max-value-in-each-row)
    j1 <- max.col(sorted_phases, "first")
    value <- sorted_phases[cbind(1:nrow(sorted_phases), j1)]
    cluster <- names(sorted_phases)[j1]
    res <- data.frame(CellCycle = cluster)
    res <- cbind(sorted_phases,res)
    colnames(res) <- paste0(method,"_", colnames(res) )
    rownames(res) <- rownames(seurat_ob@meta.data)
    ## Add the cell cycle phase to the Seurat object as Metadata
    seurat_ob <- AddMetaData(seurat_ob, metadata=res, col.name =colnames(res) ) 
    seurat_ob@meta.data[,cellcycle_name]= factor(seurat_ob@meta.data[,cellcycle_name],levels=colnames(phases))

} else if (method=="scran"){
    if (spe %in% c("human","mouse")) {
        ref.pairs <- readRDS(system.file("exdata", paste0(spe,"_cycle_markers.rds"), package="scran"))
        a <- read.table(paste0("/home/dongjiaoyang/database/genename2id/",spe,".txt"),header=F)
        a=tibble::column_to_rownames(as.data.frame(a),"V1")
        a[,1]=as.character(a[,1])
    } else stop("species not supported.")
    ref=list(
        G1=data.frame(first=a[ref.pairs$G1[,1],],second=a[ref.pairs$G1[,2],],stringsAsFactors=F),
        S=data.frame(first=a[ref.pairs$S[,1],],second=a[ref.pairs$S[,2],],stringsAsFactors=F),
        G2M=data.frame(first=a[ref.pairs$G2M[,1],],second=a[ref.pairs$G2M[,2],],stringsAsFactors=F)
    )

    #=================================================================================
    # Classifying a new dataset:
    #=================================================================================
    if( length(intersect(rownames(seurat_ob),ref[[1]]$first))==0 ) stop("Gene names not matched. please check the species.")
    sce_ob<- as.SingleCellExperiment(seurat_ob) #to sce
    assignments <- scran::cyclone(sce_ob, ref)
    scores= assignments$normalized.scores %>% mutate(barcodes=rownames(seurat_ob@meta.data), phase=assignments$phases) %>% select(barcodes,everything()) %>% tibble::column_to_rownames("barcodes")
    colnames(scores) <- paste0(method,"_",colnames(scores))
    colnames(scores)[4] <- cellcycle_name
    seurat_ob <- AddMetaData(seurat_ob, metadata=scores, col.name =colnames(scores) )
    # names(assignments$phases) = rownames(seurat_ob@meta.data)
    # rownames(assignments$normalized.scores) = rownames(seurat_ob@meta.data)
    # rownames(assignments$scores) = rownames(seurat_ob@meta.data)
    # saveRDS(assignments,"scran_result.rds")

}

#saveRDSMC(seurat_ob,opt$INEXPRESS)

#=================================================================================
# visualizataion
#=================================================================================
nlevel = length(unique(seurat_ob@meta.data[,cellcycle_name]))
ggdim = DimPlot(object = seurat_ob,
               dims = c(1,2),
               reduction = reduct,
               pt.size = pointsize,
               group.by = cellcycle_name ) +
               ggplot2::scale_colour_manual( values = CustomCol2(1:nlevel))

ggplot2::ggsave(file.path(output_dir,paste0(reduct, "_resolution",
    resolution,"_visualize_CellCycle_plot.pdf",collapse="")), plot = ggdim)
ggplot2::ggsave(file.path(output_dir,paste0(reduct, "_resolution",
    resolution,"_visualize_CellCycle_plot.png",collapse="")),dpi=1000,plot = ggdim)

meta.data = seurat_ob@meta.data %>%
                mutate( cell_barcode = rownames(seurat_ob@meta.data)) %>%
                dplyr::select(cell_barcode,everything())
write.table(meta.data, file.path(output_dir, "cell_cycle_annotation_result.xls"),
            col.names =T, row.names =F,sep="\t",quote=F )

## split visualizataion ###
# output_dir = file.path(output_dir, paste0( "visualize_cluster_by_", cellcycle_name, collapse = ""))
# if ( !file.exists(output_dir) ){
#     dir.create(output_dir,recursive=T)
# }
if ( !is.null(facetby) ){
    for ( facetbyx in facetby ){
        DATA <- as.data.frame( meta.data[,c(facetbyx, cellcycle_name)]) %>%
                    table() %>% as.data.frame()  %>% dplyr::rename(cell_number=Freq) %>% 
                    dplyr::arrange(get(facetbyx)) %>% 
                    dplyr::group_by(.dots= facetbyx) %>% 
                    dplyr::mutate(freq =(cell_number / sum(cell_number)) * 100) %>% 
                    as.data.frame()

        write.table(DATA,
            file.path(output_dir,file=paste0(facetbyx,"_clust_cond_freq_info.xls")),
        sep="\t",col.names=T, row.names =F,quote=F)
        # clust_sum_all = PlotAbundances(seurat_ob, prop.by = cellcycle_name , group.by = facetbyx, method = "barplot",
        #     cols= CustomCol2(1:nlevel)
        # )
        ### barplot 
        sum_all <- ggplot(DATA,aes_string(x=facetbyx,y="freq",fill=cellcycle_name)) + 
                    geom_bar(stat="identity", position = opt$bartype ) +
                    #geom_bar(stat ="identity",position ="stack") +
                    labs(x=" ",y=paste0("Proportion [%]")) +
                    theme(panel.grid.major =element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.line = element_line(color = "black"))+
                    scale_y_continuous(expand = c(0,0)) + ## //这个可以去掉与X轴间隙
                    #scale_x_discrete(expand = c(0.5,0)) +
                    scale_fill_manual(values=CustomCol2(1:nlevel)) +
                    guides(fill=guide_legend(title="Cellcycle"))

        ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_resolution-",resolution,"_summary_barplot.pdf",collapse="")),plot=sum_all)
        ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_resolution-",resolution,"_summary_barplot.png",collapse="")),dpi=1000, plot = sum_all)

        ### pie plot ###
        sum_all = ggplot(DATA, aes(x = 0, y = freq)) +
                    geom_bar(aes_string(fill = cellcycle_name), width = 1, stat = "identity") +
                    theme_void() + 
                    scale_fill_manual(values = CustomCol2(1:nlevel)) + 
                    theme(axis.title.y = element_blank(),axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                    coord_polar(theta = "y", start = 0) + 
                    facet_wrap(eval(expr(~!!ensym(facetbyx))),ncol = 4) + 
                    theme(strip.background = element_rect(fill = NA,color = NA), strip.text = element_text(face = "bold")) +
                    guides(fill=guide_legend(title="Cellcycle"))

        ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_resolution-",resolution,"summary_pieplot.pdf",collapse="")),plot=sum_all)
        ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_resolution-",resolution,"summary_pieplot.png",collapse="")),dpi=1000, plot = sum_all)

        ## 降维聚类图 
        seurat_ob = SetIdent( seurat_ob, value = cellcycle_name)
        nrow = ceiling(length(unique(seurat_ob@meta.data[,facetbyx]))/2)
        groupvis_split = DimPlot(object = seurat_ob,
                        dims = c(1,2),reduction = opt$reduct,
                        pt.size = pointsize, ncol = 2,
                        group.by = cellcycle_name, split.by = facetbyx)+
                        theme( plot.title = element_text(hjust = 0.5))
        groupvis_split = groupvis_split + scale_colour_manual( values = CustomCol2(1:nlevel))
        ggsave(file.path(output_dir,paste0("splitby-",facetbyx,"_resolution",resolution,"_split_plot.pdf",collapse="")),
        limitsize = FALSE, plot = groupvis_split, height = 6*nrow, width = 14)
        ggsave(file.path(output_dir,paste0("splitby-",facetbyx,"_resolution",resolution,"_split_plot.png",collapse="")),
        limitsize = FALSE, plot = groupvis_split, width = 10, height = 4*nrow )
    }
}
if(!file.exists(file.path(output_dir, "Scran细胞周期分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/Scran细胞周期分析说明.docx",
  file.path(output_dir, "Scran细胞周期分析说明.docx"))
}



