#!/usr/bin/env Rscript

#' A function to draw a heatmap of top cluster marker genes with
#' subgroup labels. The function uses the "scale.data" slot to make the heatmap.
#' @param object A seurat objected with scaled data and cluster information
#' @param slot the transformation of the count matrix, options can be:"counts", "data", "scale.data".
#' @param assay the assay to use for pulling expression matrix.
#' @param cells The names of the cells to use. If NULL all cells will be used
#' @param features the features to use in the desired matrix.
#' @param show.legend whether to show the expression heatmap legend
#' @param cell.annotation If given, the name vector of variables in the metadata of the seurat object
#' @param gene.annotation A data.frame with chromosome annotation for each gene in the rownames.
#' @param group.colors a list of color schemas specifications with elements named with the strings in the cell.annotation ,
#'        specified according to the HeatmapMetadata format. Use "clusters" slot to specify cluster colors.
#' @param col.split.by A character string to specify the groupping variable to split the cells in rows.
#' @param row.split.by A character string to specify the groupping variable to split genes on the column.
#' @param group.order
#' @param duplicate A logical vector to specify wether to remove duplicated genes
#' @param x.size
#' @param colors the color schema mapped to the value of data matrix.
#' @param raster A logical vector
#' @import ComplexHeatmap Seurat
#' @export
CnvHeatmap <- function(
    object,
    slot = "data",
    assay = "CNV",
    group.by = "cnv",
    features = NULL,
    cells = NULL,
    cell.annotation = NULL,
    gene.annotation = NULL ,
    group.colors = list(ident = "blindless"), # named list of color schema for each annotation variable
    disp.range = "auto",
    col.split.by = NULL,
    row.split.by = NULL,
    group.order = NULL,
    duplicate = FALSE,
    show.legend = TRUE,
    x.size = 5.5,
    colors = NULL, #circlize::colorRamp2(c(-1.5, 0, 1.5), c("#FF08BD","#000000","#FFF42F"))
    raster = TRUE,
    old_color= FALSE,
    ...
){
    if ( is.null(cells) ){
        cells_use <- colnames(object)
    }else{
        cells_use <- cells
    }
    if ( is.null(assay) ){
        assay <- Seurat::DefaultAssay(object)
    }else{
        assay <- assay
    }
    Seurat::DefaultAssay(object) <- assay
    if ( is.null(features) ){
        features <- rownames(object)
    }else{
        features <- features
    }
    # make sure features are present
    possible.features <- unique(rownames(object))
    if (any(!features %in% possible.features)) {
        bad.features <- features[!features %in% possible.features]
        genes_use <- features[features %in% possible.features]
        if(length(genes_use) == 0) {
            stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
        }
        warning("The following features were omitted as they were not found in the ", slot,
        " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
    }else{
        genes_use = features
    }


    # fetch the cell annotation
    if(!is.null(cell.annotation)) {
        if( all(cell.annotation %in% c(group.by,colnames(object[[]]))) ) {
            cell.anno.df <- FetchData(object, vars = unique(c(group.by,cell.annotation)), cells = cells_use)
        }else{
            subgroup.miss <- setdiff(cell.annotation, colnames(object[[]]))
            stop(glue::glue("specified groups: {subgroup.miss} not found in the metadata"))
        }
    }else{
        cell.anno.df <- FetchData(object, vars = group.by, cells = cells_use)
    }

    # Reorder cells by (given) group
    if (!is.null(group.order)){
        for ( x in names(group.order)){
            cell.anno.df[,x] <- factor( cell.anno.df[,x], levels = group.order[[x]] )
            }
    }else{
        #for ( x in unique(c(group.by,cell.annotation))){
        #   cell.anno.df[,x] <- factor( cell.anno.df[,x], levels = unique(sort(object@meta.data[,x])) )
        #}
    }
    cell.anno.df <- cell.anno.df %>% tibble::rownames_to_column(var="barcodes") %>% arrange(!!! rlang::syms(colnames(cell.anno.df))) %>%  tibble::column_to_rownames("barcodes")
    cell.barcode.order <- rownames(cell.anno.df)
    data2use = GetAssayData(object, slot=slot)[genes_use, cells_use]

    if ( (length(disp.range) == 1) & (disp.range[1] == "auto") ) {
        # examine distribution of data that's off-center, since much of the center could
        # correspond to a mass of data that has been wiped out during noise reduction
        x.center <- mean(data2use)
        quantiles = quantile(data2use[ data2use!= x.center], c(0.01, 0.99))

        # determine max distance from the center.
        delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
        low_threshold = x.center - delta
        high_threshold = x.center + delta
        disp.range = c(low_threshold, high_threshold)
    } else {
        # use defined values
        low_threshold = disp.range[1]
        high_threshold = disp.range[2]

        if (low_threshold > x.center | high_threshold < x.center | low_threshold >= high_threshold) {
            stop(paste("Error, problem with relative values of x.range: ", x.range, ", and x.center: ", x.center))
        }
    }

    data2use[data2use < low_threshold] <- low_threshold
    data2use[data2use > high_threshold] <- high_threshold
    data2use<-data2use[, cell.barcode.order]

    # set up the color palette for each annotation variable
    color_schema = list()
    for ( x in colnames(cell.anno.df) ){
        nlevels <- length(unique(cell.anno.df[,x]))
        color_schema[[x]] <- SelectColors(cell.anno.df, value = x, n = nlevels, palette = group.colors[[x]])
    }
    # color_schema <- lapply(colnames(cell.anno.df), function(x){
    #     nlevels <- length(unique(cell.anno.df[,x]))
    #     cluster_cols <- SelectColors(object, value = x, n = nlevels, palette = group.colors[[x]])
    # })
    # names(color_schema) = colnames(cell.anno.df)

    # set the order of legend
    annotation_legend_params = list()
    for ( x in colnames(cell.anno.df) ){
        if ( !is.null( group.order[[x]]) ){
            annotation_legend_params[[x]]  = list( at = group.order[[x]],
                                                   labels = group.order[[x]])
        }else{
            annotation_legend_params[[x]]  = list( at =sort(unique(cell.anno.df[,x])),
                                                   labels = sort(unique(cell.anno.df[,x])) )
        }
    }

    cellAnnotation = ComplexHeatmap::HeatmapAnnotation(
                                            df = cell.anno.df,
                                            col = color_schema,
                                            annotation_legend_param = annotation_legend_params,
                                            which = "row",
                                            show_annotation_name = F)
    # rowdata_split
    if(!is.null(row.split.by)){
        if(length(row.split.by) > 1) stop("rowdata_split can only be one column")
        if(!row.split.by %in% colnames(object[[]])) {
            stop("Requested split not in cell annotations.")
        }
        rowdata_split <- FetchData(object, vars = row.split.by, cells = cells_use )
    }else{
        rowdata_split <- row.split.by
    }

    # get the vector of per-cell cluster names
    if ( !is.null(gene.annotation) ){
        gene.annotation$chr <- gsub("^chr","",gene.annotation$chr)
        gene.anno <- gene.annotation[intersect(rownames(gene.annotation),genes_use),"chr", drop =F]
        gene.anno$chr <- sort(as.numeric(gene.anno$chr))
        if(duplicate == FALSE){
            gene.anno$order <- 1:dim(gene.anno)[1]
            gene.anno$gene <- rownames(gene.anno)
            marker.dup<-names(which(table(gene.anno$gene) >1))
            for(i in marker.dup){
                dup<-gene.anno[which(gene.anno$gene == i),]
                del.col<-dup[-which.max(dup$avg_logFC),"order"]
                gene.anno <- gene.anno[-which(gene.anno$order %in% unlist(del.col)),]
            }
        }
        data2use<-data2use[rownames(gene.anno), ]

        # coldata_split
        if(!is.null(col.split.by)){
            if(length(col.split.by) > 1) stop("rowdata_split can only be one row")

            if(!col.split.by %in% colnames(gene.annotation) ){
                stop("Requested split not in gene annotations.")
            }
            coldata_split <- gene.anno[,col.split.by, drop =F]
        }else{
            coldata_split <- NULL
        }

    }else{
        geneAnnotation <- NULL
    }
    # draw the heatmap
    if (old_color){
        ha <- ComplexHeatmap::Heatmap(t(data2use),name="CNV level",
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_order = cell.barcode.order,
        row_split = rowdata_split,
        row_gap = unit(0, "mm"),
        # column setting
        column_order = rownames(gene.anno),
        column_split = coldata_split,
        column_gap = unit(0, "mm"),
        column_title =unique(gene.anno$chr),
        column_title_side="bottom",
        border = T,
        # top_annotation = geneAnnotation ,
        left_annotation = cellAnnotation,
        use_raster=raster,
        raster_quality=4,
        show_heatmap_legend = TRUE,
        show_row_names =F,
        show_column_names = F
        )
    }else{
      ha <- ComplexHeatmap::Heatmap(t(data2use),name="CNV level",
        cluster_rows = FALSE, cluster_columns = FALSE,
        col = colorRampPalette(c("darkblue", "white", "darkred"))(56),
        row_order = cell.barcode.order,
        row_split = rowdata_split,
        row_gap = unit(0, "mm"),
        # column setting
        column_order = rownames(gene.anno),
        column_split = coldata_split,
        column_gap = unit(0, "mm"),
        column_title =unique(gene.anno$chr),
        column_title_side="bottom",
        border = T,
        # top_annotation = geneAnnotation ,
        left_annotation = cellAnnotation,
        use_raster=raster,
        raster_quality=4,
        show_heatmap_legend = TRUE,
        show_row_names =F,
        show_column_names = F)
      }   
    return(ha)
}

#' Creates a list of unique color values used for plotting
#'
#' @return A named vector of unique hexedecimal color values, either generated from a preselected
#'         vector of 20 unique colors, or from a sequence of colors in hsv colorspace.
#'
#' @param seurat.obj A singular preprocessed Seurat object
#' @param gradient Setting to TRUE will use a sequence of hsv colors instead of 20 unique colors,
#'                 useful for comparisons of more than 20 cell types.
#' @param value The Seurat metadata slot to generate colors for. Defaults to "celltype".
#'
#' @import SingleCellExperiment
#' @import Seurat
#'
#' @seealso \code{\link{as.SingleCellExperimentList}}
#' @seealso \code{\link{ExtractGenes}}
#' @seealso \code{\link{DecoderVariance}}
#' @seealso \code{\link{MeanDecoderVariance}}
#' @seealso \code{\link{GetCharMetadata}}
#'
#' @examples
#' DimPlot(object = seurat.obj,
#'         reduction = "tsne",
#'         cols = SelectColors(seurat.obj),
#'         group.by = "celltype",
#'         label = TRUE,
#'         repel = TRUE)
#'
#' @export
SelectColors <- function(
    object = NULL,
    palette = "blindless",
    value = "celltype",
    n = NULL
){
    if ( !is.null(object) ){
        if ( class(object) == "data.frame" ){
            colid <- ifelse( is.null(value), colnames(object)[1], value )
            if (is.factor(object[[colid]])) {
                names= levels(object[[colid]])
            } else {
                names <- unique(object[[colid]])
            }
        }
        if ( is.factor(object) ){
            names <- levels(object)
        }
        if ( is.vector(object) ){
            names <- unique(object)
        }
        n = length(names)
    }else if ( !is.null(n) ) {
        names = NULL
    }

    colors2pick = switch(palette,
    ##ref: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
    ditto = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
    "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
    "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
    "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"),
    CustomCol2 = c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17"),
    seurat = hcl( h = seq(15, 375, length = n+1), l = 65, c = 100),
    col50 =  c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
    "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
    "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
    "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
    "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
    "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
    "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
    "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
    "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
    "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c"),
    paired = brewer.pal(n = n, 'Paired'),
    colx22 = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
    '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
    '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
    '#000075', '#808080', '#4f34ff', '#f340F0'),
    jet = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
    "#FF7F00", "red", "#7F0000" ),
    tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
    "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
    "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
    "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
    tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
    "#CDCC5D", "#6DCCDA"),
    colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
    "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
    "#FFBC79", "#CFCFCF"),
    trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
    "#FFC156", "#69B764", "#F26C64", "#FFDD71",
    "#9FCD99"),
    purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
    "#5F5A41", "#B4B19B", "#995688", "#D898BA",
    "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
    bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
    "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
    "#BD0A36", "#F4737A"),
    greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
    "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
    "#39737C", "#86B4A9", "#82853B", "#CCC94D"),
    cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
    "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
    "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
    "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"),
    CustomCol = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
    blindless = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
    "#BF5B17", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4",
    "#B3DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3", "#BC80BD",
    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
    "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
    "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8",
    "#984EA3",  "#FFFF33", "#A65628", "#F781BF", "#999999","#FFED6F",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9","#666666"),
    )

    if ( is.null(n) ){
        colors_use <- colors2pick
    }else{
        colors_use <- colors2pick[1:n]
    }
    if ( !is.null(names) ){ names(colors_use) <- names }
    return(colors_use)
}


#' @description Add meta.data about CNAs to a Seurat object from an infercnv_obj
#'
#' @title add_to_seurat()
#' @param seurat_obj Seurat object to add meta.data to (default: NULL)
#' @param infercnv_output_path Path to the output folder of the infercnv run to use
#' @param nclones How many of the largest CNA (in number of genes) to get.
#' @param bp_tolerance How many bp of tolerance to have around feature start/end positions for top_n largest CNVs.
#' @return seurat_obj
#' @export
#'
add_to_seurat_hmm <- function(
    seurat_obj = NULL,
    infercnv_output_path,
    nclones = 4,
    quantile = 0.75
) {
    if (!file.exists(file.path(infercnv_output_path, "run.final.infercnv_obj"))) {
        stop(glue::glue("Could not find \"run.final.infercnv_obj\" file at {infercnv_output_path}!"))
    }
    infercnv_obj = readRDS(file.path(infercnv_output_path, "run.final.infercnv_obj"))

    if (is.null(seurat_obj)) {
        stop("No Seurat object provided, will only write metadata matrix.")
    }else{
        # the cells in Seurat object must be subset of cells in infercnv_obj
        seurat_obj = seurat_obj[, colnames(infercnv_obj@expr.data)]
        if (!all( Cells(seurat_obj) %in%  colnames(infercnv_obj@expr.data)) ){
            stop( "Found cells in Seurat object is not in the infercnv object.
            Please check it out and rerun!")
        }
        desired_cells <- intersect( Cells(seurat_obj),  colnames(infercnv_obj@expr.data))
        seurat_obj <- Seurat::SubsetData(seurat_obj, cells = desired_cells)
    }

    # get and visualize the dendrogram
    #obs_hcl <- infercnv_obj@tumor_subclusters$hc[names(infercnv_obj@observation_grouped_cell_indices)]
     obs_hcl <- names(infercnv_obj@observation_grouped_cell_indices)
    if ( length(obs_hcl) >1 ){
        tmp_dendrogram <- phylogram::as.dendrogram(obs_hcl[[1]])
        for ( i in 2:length(obs_hcl) ){
            tmp_dendrogram <- merge(tmp_dendrogram, phylogram::as.dendrogram(obs_hcl[[i]]))
        }
        obs_hcl <- as.hclust(tmp_dendrogram)
        tumor_ordered_names <- obs_hcl$labels[obs_hcl$order]
        # obs_dendrogram <- phylogram::as.dendrogram(obs_hclust)
        #g <- dendextend::cutree(obs_hcl, k=nclones) #named group results vector
        g <- dendextend::cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=nclones)
        g <- g[tumor_ordered_names]
    }else{
        #tumor_ordered_names <- obs_hcl[[1]]$labels[obs_hcl[[1]]$order]
        tumor_ordered_names <- infercnv_obj@tumor_subclusters$hc$all_observations$labels[infercnv_obj@tumor_subclusters$hc$all_observations$order]
        # cut the dendrogram
        # obs_dendrogram <- phylogram::as.dendrogram(obs_hcl)
        #g <- dendextend::cutree(obs_hcl[[1]], k=nclones) #named group results vector
        g <- dendextend::cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=nclones)
        g <- g[tumor_ordered_names]
    }
    # tp <- color_branches(dend, k = NCLONES,groupLabels = T)

    normal_cell_bc <- setdiff(colnames(infercnv_obj@expr.data), tumor_ordered_names)
    groups <- data.frame(barcode=c(names(g), normal_cell_bc), cnv_group=c(as.numeric(g),rep("normal", length(normal_cell_bc))),
    stringsAsFactors = F) %>% tibble::column_to_rownames(var = "barcode")

    # use the copy number vairable genes to calculate the cnv level not all genes
    msq <-  rowMeans(infercnv_obj@expr.data^2)
    desired_genes <- names(msq)[msq>=quantile(msq,quantile)]
    infercnv_level <- apply(as.data.frame(t(infercnv_obj@expr.data[desired_genes,])), 1, function(x) {
        x[is.na(x)] <- 0
        return(sum(x))
    })
    infercnv_level <- round(scales::rescale(infercnv_level / nrow(infercnv_obj@expr.data[desired_genes,]), c(1, 100)), 0)
    infercnv_level <- infercnv_level[Cells(seurat_obj)]
    # filter the genes on chromes like MT or other contig
    if( startsWith(as.character(infercnv_obj@gene_order$chr[1]), "chrchr") ){
		infercnv_obj@gene_order$chr = factor( sub("chr","",infercnv_obj@gene_order$chr), levels= sub("chr","",levels(infercnv_obj@gene_order$chr))  )
    }
    genes <- infercnv_obj@gene_order %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter( chr %in% paste0("chr", 1:40)) %>%
    tibble::column_to_rownames(var = "gene")
    

    seurat_obj[["CNV"]] <- Seurat::CreateAssayObject(data = infercnv_obj@expr.data[rownames(genes), Cells(seurat_obj)])
    seurat_obj$cnv_level <- infercnv_level
    seurat_obj$cnv_group <- groups[Cells(seurat_obj), "cnv_group"]
    seurat_obj[["CNV"]] <- AddMetaData(seurat_obj[["CNV"]], metadata = genes)

    return(seurat_obj)
}


add_to_seurat <- function(
    seurat_obj = NULL,
    infercnv_output_path,
    nclones = 4,
    quantile = 0.75
) {
    if (!file.exists(file.path(infercnv_output_path, "run.final.infercnv_obj"))) {
        stop(glue::glue("Could not find \"run.final.infercnv_obj\" file at {infercnv_output_path}!"))
    }
    infercnv_obj = readRDS(file.path(infercnv_output_path, "run.final.infercnv_obj"))

    if (is.null(seurat_obj)) {
        stop("No Seurat object provided, will only write metadata matrix.")
    }else{
        # the cells in Seurat object must be subset of cells in infercnv_obj
        if (!all( Cells(seurat_obj) %in%  colnames(infercnv_obj@expr.data)) ){
            print( "Found cells in Seurat object is not in the infercnv object!")
        }
        desired_cells <- intersect( Cells(seurat_obj),  colnames(infercnv_obj@expr.data))
        seurat_obj <- Seurat::SubsetData(seurat_obj, cells = desired_cells)
    }

    # get and visualize the dendrogram
    obs_hcl <- infercnv_obj@tumor_subclusters$hc[names(infercnv_obj@observation_grouped_cell_indices)]
    if ( length(obs_hcl) >1 ){
        tmp_dendrogram <- phylogram::as.dendrogram(obs_hcl[[1]])
        for ( i in 2:length(obs_hcl) ){
            tmp_dendrogram <- merge(tmp_dendrogram, phylogram::as.dendrogram(obs_hcl[[i]]))
        }
        obs_hcl <- as.hclust(tmp_dendrogram)
        tumor_ordered_names <- obs_hcl$labels[obs_hcl$order]
        # obs_dendrogram <- phylogram::as.dendrogram(obs_hclust)
        g <- dendextend::cutree(obs_hcl, k=nclones) #named group results vector
        g <- g[tumor_ordered_names]
    }else{
        tumor_ordered_names <- obs_hcl[[1]]$labels[obs_hcl[[1]]$order]
        # cut the dendrogram
        # obs_dendrogram <- phylogram::as.dendrogram(obs_hcl)
        g <- dendextend::cutree(obs_hcl[[1]], k=nclones) #named group results vector
        g <- g[tumor_ordered_names]
    }
    # tp <- color_branches(dend, k = NCLONES,groupLabels = T)

    normal_cell_bc <- setdiff(colnames(infercnv_obj@expr.data), tumor_ordered_names)
    groups <- data.frame(barcode=c(names(g), normal_cell_bc), cnv_group=c(as.numeric(g),rep("normal", length(normal_cell_bc))),
    stringsAsFactors = F) %>% tibble::column_to_rownames(var = "barcode")

    # use the copy number vairable genes to calculate the cnv level not all genes
    msq <-  rowMeans(infercnv_obj@expr.data^2)
    desired_genes <- names(msq)[msq>=quantile(msq,quantile)]
    infercnv_level <- apply(as.data.frame(t(infercnv_obj@expr.data[desired_genes,])), 1, function(x) {
        x[is.na(x)] <- 0
        return(sum(x))
    })
    infercnv_level <- round(scales::rescale(infercnv_level / nrow(infercnv_obj@expr.data[desired_genes,]), c(1, 100)), 0)
    infercnv_level <- infercnv_level[Cells(seurat_obj)]
    # filter the genes on chromes like MT or other contig
    if( startsWith(as.character(infercnv_obj@gene_order$chr[1]), "chrchr") ){
		infercnv_obj@gene_order$chr = factor( sub("chr","",infercnv_obj@gene_order$chr), levels= sub("chr","",levels(infercnv_obj@gene_order$chr))  )
    }
    genes <- infercnv_obj@gene_order %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter( chr %in% paste0("chr", 1:40)) %>%
    tibble::column_to_rownames(var = "gene")
    

    seurat_obj[["CNV"]] <- Seurat::CreateAssayObject(data = infercnv_obj@expr.data[rownames(genes), Cells(seurat_obj)])
    seurat_obj$cnv_level <- infercnv_level
    seurat_obj$cnv_group <- groups[Cells(seurat_obj), "cnv_group"]
    seurat_obj[["CNV"]] <- AddMetaData(seurat_obj[["CNV"]], metadata = genes)

    return(seurat_obj)
}


suppressWarnings({
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("Matrix") )
    suppressPackageStartupMessages( library("optparse") )
    suppressPackageStartupMessages( library("ggplot2") )
    suppressPackageStartupMessages(library("gridExtra"))
    suppressPackageStartupMessages( library("ComplexHeatmap") )
    suppressPackageStartupMessages( library("OESingleCell"))
    suppressPackageStartupMessages( library("dplyr"))
    suppressPackageStartupMessages( library("infercnv"))
    suppressPackageStartupMessages( library("tibble"))
    suppressPackageStartupMessages( library("RColorBrewer"))
    suppressPackageStartupMessages( library("dittoSeq"))
})

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i"), type = "character",
        help = "the processed seurat object saved as R object in RDS format."),
    make_option( c("--informat", "-f" ), type = "character", default = "seurat",
        help = "The indication of type of input expression matrix, the possible type can be:
                                seurat: the seurat object from the clustering results.
                                [default \"%default\"]"),
    make_option( c("--infercnv_output_path", "-l"), type = "character",
        help = "[OPTIONAL]the infercnv result directory." ),
    make_option( c("--groupby", "-g"), type = "character",
        help = "the groupping of cells in the metadata of seurat object." ),
    make_option( c("--gorder", "-e"), type = "character", default = NULL,
        help = "[OPTIONAL]comma seperated list of levels for the specified parameter '--groupby' to order the cells.
                The format can be gvariable1:l1,l2,l3;gvariable2:l1,l2,l3." ),
    make_option( c("--subset", "-s"), type = "character",
        help = "Logical expression indicating cells to keep. Example: clusters %in% c(1,2,3)"),
    make_option( c("--invert", "-v"), type = "logical", default = FALSE,
        help = "reverse the selection conditions specified by all other conditions.[default \"%default\"]"),
    make_option( c("--output", "-o"), type = "character", default = "./",
        help = "the output directory of results." ),
    make_option( c("--colormapping", "-m"), type = "character",
        help = "The color mapping for groupping column of cells set by the parameters '--groupby'.
                The exmaple format is variable1:colorschema1,variable2:colorschema2. The supported color schemas can be:
                 blindless, col50, ditto, paired, CustomCol2." ),
    make_option( c("--vlnlegend"), type = "logical", default = F,
        help = "whether to show the legend of the vlnplot." ),
    make_option( c("--splitby","-y"), type = "character", default = NULL,
        help = "the column used to split the vlnplot/featureplot." ),
    make_option( c("--vismethod","-p"), type = "character", default = "all",
        help = "plots to visualize. choose from all, heatmap, vlnplot, featureplot." ),
    make_option( c("--pointsize"), type = "double", default =0.8 ,
         help = "the point size in the plot."),
    make_option( c("--reduct"), type = "character", default = "tsne",
               help = "the reduction coordination to visualize the  results.[default:tsne]"),
    make_option( c("--nclust", "-n" ), type="integer", default = 5,
        help="the number of groups to use seperate the cancer cells according to
              the previewed cnv heatmap from infercnv results.[default \"%default\"]"),
    make_option( c("--hmm" ), type="logical", default = FALSE,
        help="If origin CNV use hmm or not"),
    make_option( c("--colors" ), type="logical", default = FALSE,
        help="If this item is set to TRUE, use the old color scheme")
);
opt_parser = OptionParser(option_list=option_list, usage = "Rscripit infercnv_vis.R -i RDS -f seurat -l ./ -g gvariable1,gvariable2 -o ./ -m variable1:colorschema1,variable2:colorschema2");
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
        dir.create(output_dir)
    }
}
output_dir = normalizePath(output_dir )

# set group order
if ( !is.null(opt$gorder)){
    group.order =  list()
    for ( x in  unlist(strsplit(opt$gorder, ";", perl =T))){
        m = unlist( strsplit( x, ":|,", perl =T))
        group.order[[m[1]]] = m[-1]
    }
}else{
    group.order = NULL
}

# set the groupby level
if ( !is.null( opt$groupby) ){
    cell.annos = unlist( strsplit( opt$groupby, ",", perl =T))
    groupby = unlist( strsplit( opt$groupby, ",", perl =T))[1]
}else{
    stop("please provide the groupping of cells in the metadata of seurat object.")
}

# set the color schema of cell annotation bar
if ( !is.null( opt$colormapping) ){
    group_colors = list()
    for ( x in  unlist(strsplit(opt$colormapping, ",", perl =T))){
        m = unlist( strsplit( x, ":", perl =T))
        group_colors[[m[1]]] = m[-1]
    }
}else{
    group_colors= 1:length(cell.annos)
    names(group_colors) = cell.annos
    group_colors = as.list(group_colors)
}

if( opt$vismethod == "all" ){
    vismethods = c("heatmap","vlnplot","featureplot")
}else{
    vismethods = unlist(strsplit(opt$vismethod,","))
}

if ( is.null(opt$pointsize) ){
    if (dim(seurat_ob)[2] < 500){
        pointsize = 1.5
    } else pointsize = 0.8
} else {
pointsize = opt$pointsize
}



if ( is.null( opt$reduct ) ){
    reduct.use = "tsne"
} else {
    reduct.use = opt$reduct
}

#=================================================================================
# readin the results/rds
#=================================================================================
if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){
        # the input is a seurat object produced by previous analysis
        seurat_ob = OESingleCell::readRDSMC( opt$input, cores = 10)
        # if the input seurat object version is less than 3, upgrade it to version 3
        if ( seurat_ob@version < 3 ){
            seurat_ob = Seurat::UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
        }


        # subset cells on the expression matrix using logical exprssion on cell metadata
        if ( !is.null(opt$subset) ){
            predicate = opt$subset
            df = seurat_ob@meta.data
            desired_cells= subset(df, eval( parse(text=predicate)))
            seurat_ob = subset(seurat_ob, cells = rownames(desired_cells), invert = opt$invert)
        }
    }
}

if ( !is.null(opt$infercnv_output_path) ){
    # integrate the infecnv result to the seurat object as cell annotations in the metadata slot
    # singlecell_ob = add_to_seurat(seurat_ob,
    #                             infercnv_output_path=output_dir,
    #                             top_n=10 # How many of the largest CNA (in number of genes) to get.
    #                          )
    if ( opt$hmm){
    seurat_ob = add_to_seurat_hmm(seurat_ob,
                              infercnv_output_path = opt$infercnv_output_path,
                              nclones = opt$nclust )
    }else{
    seurat_ob = add_to_seurat(seurat_ob,
                              infercnv_output_path = opt$infercnv_output_path,
                              nclones = opt$nclust )
    }
    OESingleCell::saveRDSMC( seurat_ob, file.path(output_dir, "cnv_seurat.rds"))
}else{
    if ( !"CNV" %in% Seurat::Assays(seurat_ob) ){
       stop( "NO previous CNV assay results integration FOUND!")
    }
}

cnv_result = seurat_ob@meta.data %>% tibble::rownames_to_column(var = "barcode_inuse") %>%
            dplyr::rename( "Barcode" =  ifelse("rawbc" %in% colnames(seurat_ob@meta.data), "rawbc","orig.ident")) %>%
             select(barcode_inuse, Barcode, sampleid, clusters, group, cnv_level, cnv_group)
write.table(cnv_result, file.path( output_dir, "cnv_result.xls"),
            quote = F,sep ="\t",row.names =F)

if ( !is.null(opt$splitby)){
    facetbyx = unlist(strsplit(opt$splitby,",",perl=T))
    if (facetbyx %in% colnames(seurat_ob@meta.data)) {
        facetbyx = facetbyx
    } else {
        stop("can not find -splitby in the metadata")
    }
}else{
    splitby = NULL
} 
#=================================================================================
# visulization
#=================================================================================
for ( vismethod in vismethods ){
    if ( vismethod == "heatmap" ){ # heatmap
        gene_order_f =  seurat_ob[["CNV"]]@meta.features
        pdf( file.path( output_dir, "cnv_heatmap.pdf"), width = 18, height = 12)
        print(CnvHeatmap(seurat_ob, group.by = groupby,
                        group.order = group.order,
                        features = rownames(gene_order_f), cell.annotation = cell.annos,
                        gene.annotation = gene_order_f[,1,drop=F],row.split.by = NULL,
                        col.split.by = "chr",group.colors = group_colors,old_color= opt$colors))
        dev.off()
    }

    if ( vismethod == "vlnplot" ){ ## dittoPlot
        for (group_by in cell.annos){
            nlevel=length(unique(seurat_ob@meta.data[,group_by] ))
            vlnplot = dittoPlot(seurat_ob, var= "cnv_level", group.by=group_by, plots =c("vlnplot", "boxplot"), legend.show = F ) +
                      scale_fill_manual(values = SelectColors(object=NULL,palette = group_colors[[group_by]], value = group_by ,n = nlevel))
            ggsave(file.path( output_dir,paste0("cnv_vlnplot_groupby_",group_by,".pdf",collapse = "") ), 
                   plot = vlnplot ,
                   width = nlevel*2^0.5,,limitsize = FALSE)
            ggsave(file.path( output_dir,paste0("cnv_vlnplot_groupby_",group_by,".png",collapse = "") ), 
                   plot = vlnplot ,
                   width = nlevel*2^0.5,,limitsize = FALSE)
            
            if ( !is.null(opt$splitby)){
                for ( facetby in facetbyx ){
                    if ( facetby != group_by) {
                        nlevel=length(unique(seurat_ob@meta.data[,group_by] ))
                        nlevel_splitby = ifelse(is.null(facetby),1,length(unique(seurat_ob@meta.data[,facetby] )))
                        vlnplot = dittoPlot(seurat_ob, var= "cnv_level", group.by=group_by, plots =c("vlnplot", "boxplot"), legend.show = F ,split.by = facetby,split.ncol=2) +
                                scale_fill_manual(values = SelectColors(object=NULL,palette = group_colors[[group_by]], value =group_by ,n = nlevel))+
                                theme( plot.title = element_text(hjust = 10),axis.title.x = element_text(size =20),axis.title.y = element_text(size = 20),axis.text.x=element_text(size = 20),axis.text.y=element_text(size = 10))
                        ggsave(file.path( output_dir,paste0("cnv_vlnplot_groupby_",group_by,ifelse(is.null(facetby),"",paste0("_splitby_",facetby)),".pdf",collapse = "") ), 
                            plot = vlnplot ,
                            width = ifelse(nlevel==2, nlevel*3,nlevel*2) ,height= nlevel*nlevel_splitby^0.5+3,,limitsize = FALSE)
                        ggsave(file.path( output_dir,paste0("cnv_vlnplot_groupby_",group_by,ifelse(is.null(facetby),"",paste0("_splitby_",facetby)),".png",collapse = "") ), 
                            plot = vlnplot ,
                            width = ifelse(nlevel==2, nlevel*3,nlevel*2) ,height= nlevel*nlevel_splitby^0.5+3,,limitsize = FALSE)
                    }
                }
            }
        }
    }

    if ( vismethod == "featureplot" ){ ## featureplot
        ggfeature =  FeaturePlot(seurat_ob,features = "cnv_level", reduction= reduct.use,
                                 cols = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100),
                                 pt.size = pointsize, order= T ) +
                                 theme( plot.title = element_text(hjust = 0.5))
        ggsave(file.path( output_dir,paste0("cnv_level_on_",reduct.use,".pdf")),ggfeature,limitsize = FALSE)
        ggsave(file.path( output_dir,paste0("cnv_level_on_",reduct.use,".png")),ggfeature, dpi=1000,limitsize = FALSE )
    }
}
##splitby---featureplot
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
gene_range <- range(FetchData(object = seurat_ob, vars = "cnv_level"))
min.cutoff <- gene_range[1]
max.cutoff <- gene_range[2]
if ( !is.null(opt$splitby)){
    for ( facetby in facetbyx ){
        nlevel = length(unique(seurat_ob@meta.data[,facetby]))
        split <- SplitObject(seurat_ob,split.by=facetby)
        if( length(levels(seurat_ob@meta.data[,facetby])) == length(unique(seurat_ob@meta.data[,facetby])) ) {
            order=levels(seurat_ob@meta.data[,facetby]) 
        }else{
            order=sort(unique(seurat_ob@meta.data[,facetby])) 
        }
        suppressMessages({
            ggfeature = lapply( order,
                    function(x) FeaturePlot(split[[ x ]],
                        features = "cnv_level",
                        reduction= reduct.use,
                        ncol = 2, pt.size = pointsize, order= T  ) +  
                        theme( plot.title = element_text(hjust = 0.5) ) + 
                        labs(title=x ) )
            ggfeature <- lapply(ggfeature, function(x) x
                        + scale_colour_gradientn(colors= myPalette(100),
                                limits = c(min.cutoff,max.cutoff)) )
            plot <- do.call(ggpubr::ggarrange,
                c(ggfeature, list(ncol = 2,
                            nrow = ceiling(length(ggfeature) / 2),
                            common.legend = TRUE,
                            legend = "right",
                            align = "none")))
            }) 

        nrow=ceiling(length(ggfeature) / 2) 
        ggsave(file.path(output_dir,paste0("cnv_featureplot_on_",reduct.use,ifelse(is.null(facetby),"",paste0("_splitby_",facetby)),".pdf",collapse = "")),plot=plot,width = ifelse(length(ggfeature)==1,yes = 5, no = 10),height = nrow*5,limitsize = FALSE)

        ggsave(file.path(output_dir,paste0("cnv_featureplot_on_",reduct.use,ifelse(is.null(facetby),"",paste0("_splitby_",facetby)),".png",collapse = "")),plot=plot, width = ifelse(length(ggfeature)==1,yes = 5, no = 10),height = nrow*5,limitsize = FALSE )

    }
}


if(!file.exists(file.path(output_dir, "inferCNV分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/inferCNV分析说明.docx",
  file.path(output_dir, "inferCNV分析说明.docx"))
}
