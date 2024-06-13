# Title     : TODO
# Objective : TODO
# Created by: Pioneer
# Created on: 2020/3/27
suppressWarnings({
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(phylogram))
    suppressPackageStartupMessages(library(dendextend))
    suppressPackageStartupMessages(library(gplots))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(infercnv))
    suppressPackageStartupMessages(library(DECIPHER))
    suppressPackageStartupMessages(library("future"))
    suppressPackageStartupMessages(library("Seurat"))
    suppressPackageStartupMessages(library("OESingleCell"))
})


#=command line parameters setting=============================
option_list = list(
make_option(c("--input", "-i"), type = "character",
help = "The input exprssion matrix in several possible format."),
make_option(c("--informat", "-f"), type = "character", default = "seurat",
help = "The indication of type of input expression matrix. The possible type can be:
                        seurat: the seurat object from the clustering results.
                 [default:%default]"),
make_option(c("--dendrogram", "-d"), type = "character", default = "infercnv.21_denoised.observations_dendrogram.txt",
help = "infercnv.21_denoised.observations_dendrogram.txt . [default:%default]"),
make_option(c("--observations", "-s"), type = "character", default = "infercnv.21_denoised.observations.txt",
help = "infercnv.21_denoised.observations.txt . [default:%default]"),
make_option(c("--nclones", "-n"), type = "integer", default = "5",
help = "clones number for cut . [default:%default]"),
make_option(c("--output", "-o"), type = "character", default = "./",
help = "the output directory of QC results. [default:%default]", metavar = "outputdir")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
print(opt)
#=================================================================================
#parse the command line parameters
#=================================================================================
input <- opt$input
output <- opt$output
dendrogram <- opt$dendrogram
nclones <- opt$nclones
observations <- opt$observations
# input<-"seurat_object.clustering_resolution0.4.rds"
# output<-"test"
# dendrogram <-"infercnv.21_denoised.observations_dendrogram.txt"
# nclones<-5
# observations<-"infercnv.21_denoised.observations.txt"
#check output directory ,and if not ,mkdir
if (file.exists(output) == "FALSE") {dir.create(output)}
output_dir = normalizePath(output)
seurat_ob = readRDSMC(input, cores = availableCores())
if (seurat_ob@version < 3) {
    seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
}
#=======================================================================================================================
# get and visualize the dendrogram
dend <- as.hclust(ReadDendrogram(file = dendrogram,))
# cut the dendrogram
# tp <- color_branches(dend, k = nclones,groupLabels = T)
# plot(tp,leaflab = 'none',horiz = F) # set horiz = TRUE to plot dendrogram as in infercnv heatmap
g <- cutree(dend, k = nclones)

#barcodes <- unlist(lapply(names(g),clean_barcode))
groups <- data.frame(row.names = gsub(" ", "_", names(g)), #BARCODE=barcodes,
CLONE = as.numeric(g),
stringsAsFactors = F)
cell_id <- rownames(groups)
seurat_subset = SubsetData(seurat_ob, cells = cell_id)
groupsx = groups[Cells(seurat_subset),]
seurat_subset$cnv_groups = groupsx
#=======================================================================================================================
infercnv_output <- as.data.frame(t(read.table(observations, check.names = FALSE, stringsAsFactor = F)))
#infercnv1 <-as.data.frame(vroom::vroom("infercnv.21_denoised.observations.txt",delim=" "),row.names =1)
create_infercnv_level <- function(df) {
    library(scales)
    # sum total CNV levels for each cell:
    infercnv_level <- apply(df, 1, function(x) {
        x[is.na(x)] <- 0
        return(sum(x))
    })
    # normalise to number of genes in plot:
    infercnv_level <- round(rescale(infercnv_level / ncol(df), c(1, 100)), 0)
    # create infercnv_level data frame:
    infercnv_level <- data.frame(cell_id = gsub("X_", "", rownames(df)), infercnv_level = infercnv_level)
    return(infercnv_level)
}
infercnv_level <- create_infercnv_level(infercnv_output)
# resort infercnv_level by cluster groups cell index
infercnv_x = infercnv_level[Cells(seurat_subset),] %>%
    filter(cell_id %in% cell_id) %>%
    tibble::column_to_rownames(var = "cell_id")
#add infercnv_level  in seurat_ob
seurat_subset = AddMetaData(seurat_subset, metadata = infercnv_x)

#=======================================================================================================================
plot <- VlnPlot(seurat_subset, features = "infercnv_level", group.by = "cnv_groups", aphla = 0.5, pt.size = 0.1) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = NA)) +
    scale_fill_discrete(name = "CNV Cluster") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "CNV Cluster", y = "CNV level") +
    theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), plot.title = element_text(size = 18), axis.text = element_text(size = 12, face = "bold"))
ggsave(file.path(output_dir, "cnv_level_vlnplot.png"), plot = plot, dpi = 1000)
ggsave(file.path(output_dir, "cnv_level_vlnplot.pdf"), plot = plot)

saveRDSMC(seurat_subset, file.path(dirname(output_dir), "infercnvlevel.rds"), threads = availableCores())

