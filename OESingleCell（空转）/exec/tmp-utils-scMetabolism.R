#!/usr/bin/Rscript
##fix table information and change default color of heatmap

library("optparse")
option_list = list(
  make_option(c("--RDS", "-v"), type = "character", default = NULL,
              help = "the seurat object saved as R object in RDS format."),
  make_option(c("--output", "-o"), type = "character", default = "./",
              help = "the output directory of results.", metavar = "character"),
  make_option(c("--assays", "-a"), type = "character", default = "SCT",
              help = "Default assays for seurat object,default:SCT", metavar = "character"),
  make_option(c("--collapseby", "-q"), type = "character",
              help = "The variable level of heatmap. e.g. clusters or sampleid or celltype"),
  make_option(c("--colors", "-c"), type = "character", default = NULL,
              help = "colors choise for Heatmap picture : redwhiteblue, redblackgreen ,yellowblackblue "),
  make_option(c("--topn", "-n"), type = "integer", default = 10,
              help = "the number of top Term for heatmap"),
  make_option(c("--pvalue", "-p"), type = "double", default = 0.05,
              help = "the pvalue "),
  make_option(c("--type", "-t"), type = "character", default = "REACTOME ",
              help = "KEGG or REACTOME for your plot."),
  make_option(c("--gmt", "-g"), type = "character", default = NULL,
              help = "[OPTION:]input your KEGG or REACTOME gmtfile"),
  #    make_option( c("--contrast", "-d"),type = "character",default = NULL,
  #            help = "[Required]levels of a factor used to compare with for final differenetial results.
  #                    The format is Factor:interesting_level:reference_level."),
  make_option(c("--cpu"), type = "integer", default = 2,
              help = "the number of cpu."),
  #    make_option( c("--reduct","-r"), type = "character", default = "tsne",
  #        help = "the previous calculated reduction result used in the featureplot,."),
  make_option(c("--method", "-m"), type = "character", default = "VISION",
              help = "method supports VISION, AUCell, ssGSEA"),
  make_option(c("--predicate"), type = "character", default = NULL,
              help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell.")
);
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
#
# suppressWarnings({
#     suppressPackageStartupMessages(library("VISION"))
#     suppressPackageStartupMessages( library("Seurat"))
#     suppressPackageStartupMessages(library("dplyr"))
#     suppressPackageStartupMessages(library("gridExtra"))
#     suppressPackageStartupMessages(library("ggplot2"))
#     suppressPackageStartupMessages(library("tibble"))
#     suppressPackageStartupMessages(library("RColorBrewer"))
#     suppressPackageStartupMessages(library("wesanderson"))
#     suppressPackageStartupMessages(library(pheatmap))
#     suppressPackageStartupMessages(library(scMetabolism))
#     suppressPackageStartupMessages(library(rsvd))
# })

#1.check file and param
# if ( is.null(opt$RDS) ){
#     stop("the seurat object is NOT AVAILABLE!")
# }else{
#     seurat_ob = readRDS(opt$RDS)
#     dim(seurat_ob)
# }
# if ( file.exists(opt$output)){
#     output_dir = opt$output
# }else{
#     output_dir = opt$output
#     dir.create(output_dir, recursive = TRUE)
# }
#if ( seurat_ob@version < 4){
#    my_test = Seurat::UpdateSeuratObject(object=seurat_ob) #make sure the seurat object match with the latest seurat package
#} else {
# my_test = seurat_ob
#}
# my_test@assays
# if ( !is.null(opt$colors) ){
#     if ( opt$colors == "redwhiteblue" ){
#         palette <- colorRampPalette(c("RoyalBlue2", "White", "Red2"))(n=256)
#     }else if ( opt$colors == "redblackgreen" ){
#         palette <- colorRampPalette(c("Green", "Black", "Red"))(n=256)
#     }else if ( opt$colors == "yellowblackblue" ){
#         palette <- colorRampPalette(c("Blue", "Black", "Yellow"))(n=256)
#     }
# }else {
#     palette <- colorRampPalette(c("steelblue", "White", "firebrick2"))(n=256)
# }

# if (!is.null(opt$predicate)) {
#     df <- my_test@meta.data
#     desired_cells <- subset(df, eval(parse(text = opt$predicate)))
#     my_test <- my_test[, rownames(desired_cells)]
# }
# print(dim(my_test))
# print(table(my_test$group))

#2.some default module
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))

CustomCol2 <- function(n) {
  my_palette = c(
    "#7fc97f", "#beaed4", "#fdc086", "#386cb0", "#f0027f", "#a34e3b", "#666666", "#1b9e77", "#d95f02", "#7570b3",
    "#d01b2a", "#43acde", "#efbd25", "#492d73", "#cd4275", "#2f8400", "#d9d73b", "#aed4ff", "#ecb9e5", "#813139",
    "#743fd2", "#434b7e", "#e6908e", "#214a00", "#ef6100", "#7d9974", "#e63c66", "#cf48c7", "#ffe40a", "#a76e93",
    "#d9874a", "#adc64a", "#5466df", "#d544a1", "#54d665", "#5e99c7", "#006874", "#d2ad2c", "#b5d7a5", "#9e8442",
    "#4e1737", "#e482a7", "#6f451d", "#2ccfe4", "#ae6174", "#a666be", "#a32b2b", "#ffff99", "#3fdacb", "#bf5b17")
  return(my_palette[n])
}

center_scale <- function(data) {
  data = data
  center_data = scale(data, center = T, scale = F)
  return(center_data)
}

SeuratObject::DefaultAssay(my_test) = opt$assays
print("OK1")
#3.check gtmFile
if (is.null(opt$gmt)) {
  print("Loading default Human gmtFile")
  countexp.Seurat <- scMetabolism::sc.metabolism.Seurat(obj = my_test, assays = opt$assays, method = opt$method, imputation = F, ncores = opt$cpu, metabolism.type = opt$type)
} else {
  print("Loading your gmtFile")
  countexp.Seurat <- scMetabolism::sc.metabolism.Seurat(obj = my_test, assays = opt$assays, method = opt$method, imputation = F, ncores = opt$cpu, metabolism.type = opt$type, gmtFile = opt$gmt)
}

print("OK2")
metabolism.matrix <- as.data.frame(countexp.Seurat@assays$METABOLISM$score)
metabolism_score = tibble::rownames_to_column(metabolism.matrix, var = opt$type)
write.table(metabolism_score, quote = F, sep = "\t", row.names = F, col.names = T,
            file.path(output_dir, paste0("metabolism_", opt$type, "_score.xls", collapse = "")))

#4.save score

if (opt$type == "KEGG") {
  my_test[['KEGG']] = CreateAssayObject(counts = metabolism.matrix)
  count = as.data.frame(my_test@assays$KEGG@counts)
}else if (opt$type == "REACTOME") {
  my_test[['REACTOME']] = CreateAssayObject(counts = metabolism.matrix)
  count = as.data.frame(my_test@assays$REACTOME@counts)
}
filename = as.character(lapply(strsplit(basename(opt$RDS), split = ".rds", fixed = TRUE), head, 1))
#saveRDS(my_test,paste0(filename,".metabolism_",opt$type,"_score.rds",collapse = ""))

heatmax = max(metabolism.matrix)
heatmin = min(metabolism.matrix)

meta.data = my_test@meta.data
collapseby = opt$collapseby
meta.data$id = rownames(meta.data)
collapsed_count = vector()
if (!collapseby %in% colnames(meta.data)) {
  stop("NO specified column found!")
}

collapsed_group = meta.data %>%
  group_by(.dots = collapseby) %>%
  do(data.frame(cellid = paste0(.$id, collapse = ",")))
if (collapseby == "clusters")  collapsed_group$clusters = paste("cluster", collapsed_group$clusters, sep = "_")

for (cells in collapsed_group$cellid) {
  samplex = unlist(strsplit(cells, ",", perl = T))
  collapsed_count = cbind(collapsed_count, rowMeans(count[, samplex, drop = F]))
}
collapsed_count = as.matrix(collapsed_count)
collapsed_group = as.data.frame(collapsed_group)
colnames(collapsed_count) = as.matrix(collapsed_group[, 1])

ave_score = tibble::rownames_to_column(as.data.frame(collapsed_count), var = opt$type)
write.table(ave_score, quote = F, sep = "\t", row.names = F, col.names = T,
            file.path(output_dir, paste0("average_", opt$type, "_score.xls", collapse = "")))

#5.plot
#prepare data
if (dim(collapsed_count)[2] == 2) {
  scale_type = "none"
  heatmap_data = as.data.frame(t(apply(collapsed_count, 1, center_scale)))
  names(heatmap_data) = colnames(collapsed_count)
}else {
  scale_type = "row"
  heatmap_data = collapsed_count
}

title = paste0("metabolism_", opt$type)
heatmap_plot = pheatmap(heatmap_data, main = title, treeheight_row = 7,
                        lwd = 1, scale = scale_type, border = F, border_color = "white",
                        color = palette, show_rownames = T, show_colnames = T, , fontsize_col = 7, fontsize_row = 7,
                        cluster_rows = T, cluster_cols = F, angle_col = 45,
                        cellheight = 8, cellwidth = 9
)
ggsave(paste0(output_dir, "/average_", opt$type, "_heatmap.pdf"), plot = heatmap_plot, height = nrow(heatmap_data) / 7 + 1, width = 9)
ggsave(paste0(output_dir, "/average_", opt$type, "_heatmap.png"), plot = heatmap_plot, height = nrow(heatmap_data) / 7 + 1, width = 9, dpi = 600, limitsize = F)
#ggsave(paste0(output_dir,"/average_",opt$type,"_heatmap.pdf"),plot=heatmap_plot,width=(dim(collapsed_count)[2]))
#ggsave(paste0(output_dir,"/average_",opt$type,"_heatmap.png"),plot=heatmap_plot,width=(dim(collapsed_count)[2]),dpi = 1000 ,limitsize = F)

##6.Find_diff
assay_metadata = my_test@meta.data
#if ( is.null(opt$contrast ) ){ #no contrast is provided
#    stop("NO contrast string is AVAILABLE!")
#}else{
#    contrast = opt$contrast
#}
contrast = paste0(collapseby, ":all:all")
contrasts = unlist(strsplit(contrast, ":", perl = T))
all_levels = as.vector(unique(assay_metadata[, contrasts[1]]))
if (contrasts[2] == "all" & contrasts[3] != "all") {
  case_levels = all_levels[-which(all_levels == contrasts[3])] #delete the reference level
  all_comparisions = paste(contrasts[1], case_levels, contrasts[3], sep = ":")
}else if (contrasts[2] != "all" & contrasts[3] == "all") {
  ref_levels = all_levels[-which(all_levels == contrasts[2])] #delete the interested level
  all_comparisions = paste(contrasts[1], contrasts[2], ref_levels, sep = ":")
}else if (contrasts[2] == "all" & contrasts[3] == "all") {
  all_comparisions = lapply(all_levels,
                            function(x) paste(contrasts[1], x, paste0(all_levels[-which(all_levels == x)], collapse = ","), sep = ":"))
  all_comparisions = unlist(all_comparisions)
}else {
  all_comparisions = contrast
}
diff_outdir = paste(output_dir, "/1.Diff_Term", sep = "")
dir.create(diff_outdir, recursive = TRUE)

if (contrasts[2] == "all" & contrasts[3] == "all") {
  # find the significant differential pathway for each Interested_group against all other Interested_group
  Diff_results = c()
  for (contrastx in all_comparisions) {
    contrastsx = unlist(strsplit(contrastx, ':', perl = T))
    DEG_META_tmp = FindMarkers(my_test,
                               ident.1 = as.character(unlist(strsplit(contrastsx[2], ",", perl = T))),
                               ident.2 = as.character(unlist(strsplit(contrastsx[3], ",", perl = T))),
                               test.use = "bimod",
                               min.pct = 0,
                               group.by = contrastsx[1],
                               logfc.threshold = -Inf,
                               assay = opt$type,
                               slot = "data")
    DEG_META_tmp$cluster = as.character(contrastsx[2])
    DEG_META = DEG_META_tmp %>% rownames_to_column(var = "geneset")
    Diff_results = rbind(Diff_results, DEG_META)
  }
  #	Diff_results = Diff_results %>% select(geneset, avg_log2FC, pct.1 ,pct.2, p_val_adj,cluster)
  colnames(Diff_results) = c("geneset", "pval", "avg_log2FC", "pct_1", "pct_2", "padj", "Interested_group")
  #Diff_results = Diff_results = Diff_results[, c("geneset", "pct_1", "pct_2","pval",  "padj","avg_log2FC", "Interested_group")]
  Diff_results = Diff_results = Diff_results[, c("geneset", "pval", "padj", "avg_log2FC", "Interested_group")]
  Diff_results = Diff_results %>% dplyr::select(geneset, everything())
  topn_markers = Diff_results %>%
    group_by(Interested_group) %>%
    filter(pval < opt$pvalue) %>%
    arrange(pval, desc(avg_log2FC)) %>%
    # filter(gene_diff > pct_fold_cutoff)  %>%
    top_n(opt$topn, avg_log2FC)

  dim(Diff_results)
  dim(topn_markers)
  write.table(Diff_results, file = file.path(diff_outdir, paste0(opt$type, "_all_results.xls", collapse = "")), quote = F, sep = "\t", col.names = T, row.names = F)
  write.table(topn_markers, file = file.path(diff_outdir, paste0("top", opt$topn, "_", opt$type, "_Significant_results.xls", collapse = "")), quote = F, sep = "\t", col.names = T, row.names = F)

  #####plot
  sort_topn = topn_markers[order(topn_markers$Interested_group),]
  sort_topn$geneset = factor(sort_topn$geneset, levels = unique(sort_topn$geneset))
  plot_list = as.matrix(unique(sort_topn$geneset))
  plot_data = heatmap_data[as.vector(plot_list),]

  topn_score = tibble::rownames_to_column(as.data.frame(plot_data), var = opt$type)
  write.table(topn_score, quote = F, sep = "\t", row.names = F, col.names = T,
              file.path(diff_outdir, paste0("Significant_", opt$type, "_score.xls", collapse = "")))

  heatmap_plot2 = pheatmap(plot_data, main = "Significant_Term", treeheight_row = 7,
                           lwd = 1, scale = scale_type,
                           border_color = "white", color = palette, show_rownames = T, show_colnames = T,
                           cluster_rows = F, cluster_cols = F, angle_col = 45, fontsize_col = 9,
                           fontsize_row = 8, cellheight = 12, cellwidth = 12
  )
  ggsave(paste0(diff_outdir, "/Significant_", opt$type, "_heatmap.pdf"), plot = heatmap_plot2, width = (dim(plot_data)[2] + 3.5), height = (dim(plot_data)[2] + 3),)
  ggsave(paste0(diff_outdir, "/Significant_", opt$type, "_heatmap.png"), plot = heatmap_plot2, width = (dim(plot_data)[2] + 3.5), height = (dim(plot_data)[2] + 3), dpi = 600, limitsize = F)
}

###vlnplot
if (is.null(opt$gmt)) {
  for (pathway in plot_list) {
    pathway_name = gsub('_-_|_/_|__|,_', '_', gsub("[[:space:]]", "_", gsub("[(.*)]", "", pathway)))
    p1 = VlnPlot(my_test, assay = opt$type, features = pathway, group.by = collapseby, pt.size = 0,) +
      NoLegend() +
      geom_boxplot(width = .2)
    ggsave(paste0(diff_outdir, "/", opt$type, "_", pathway_name, "_violin_plot.png"), plot = p1, dpi = 300, limitsize = F)
    ggsave(paste0(diff_outdir, "/", opt$type, "_", pathway_name, "_violin_plot.pdf"), plot = p1, limitsize = F)
  }
} else {
  for (pathway in plot_list) {
    pathway_name = gsub('_-_|_/_|__|,_', '_', gsub("[[:space:]]", "_", pathway))
    p1 = VlnPlot(my_test, assay = opt$type, features = pathway, group.by = collapseby, pt.size = 0,) +
      NoLegend() +
      geom_boxplot(width = .2)
    ggsave(paste0(diff_outdir, "/", opt$type, "_", pathway_name, "_violin_plot.png"), plot = p1, dpi = 300, limitsize = F)
    ggsave(paste0(diff_outdir, "/", opt$type, "_", pathway_name, "_violin_plot.pdf"), plot = p1, limitsize = F)
  }
}

if (!file.exists(file.path(output_dir, "scMetabolism代谢分析说明.docx"))) {
  file.copy("/public/scRNA/works/Documents/scMetabolism代谢分析说明.docx",
            file.path(output_dir, "scMetabolism代谢分析说明.docx")) }
