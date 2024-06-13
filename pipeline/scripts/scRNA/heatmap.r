 diffGene_path = normalizePath("/public/scRNA/works/tangxuan/HT2019-6928-human/houxu20210725/2.sub_Fibroblasts_supply/Diffexp/clusters_5-vs-6-diff-pval-0.05-FC-1.5_anno.xls")
diffGene_name = gsub(".*/","",diffGene_path)
    if(grepl("-diff-",diffGene_name)){
        name = gsub("-diff-.*","",diffGene_name)
    }else if(grepl("-all_",diffGene_name)){
        name = gsub("-all_.*","",diffGene_name)
    }
markers2vis = read.delim( "/public/scRNA/works/tangxuan/HT2019-6928-human/houxu20210725/2.sub_Fibroblasts_supply/Diffexp/clusters_5-vs-6-diff-pval-0.05-FC-1.5_anno.xls", sep="\t", header = T,quote="")
    markers2vis = subset(markers2vis,!grepl("^(mt-|Rps|Rpl|MT-|RPS|RPL)",markers2vis$gene))
    up = filter(markers2vis,FoldChange > 1) %>% arrange(desc(log2FoldChange ))  %>% top_n(200,log2FoldChange )
    down = filter(markers2vis,FoldChange < 1) %>% arrange(log2FoldChange )  %>% top_n(as.numeric(paste0("-",200)),log2FoldChange )
    topn_markers = rbind(up, down)
    write.table(topn_markers,file.path("./",paste0("top", "200", "_", name, "_genes.xls", collapse = "")),quote = F,row.names = FALSE,sep="\t")
seurat_ob = readRDS("/public/scRNA/works/tangxuan/HT2019-6928-human/sub_celltype/sub_Fibro_singlecell_object.clustering_resolution0.4.rds")
seurat_ob = SubsetData( seurat_ob, subset.name = "clusters", accept.value = c("5","6"))
seurat_ob@meta.data[["clusters"]]=factor(seurat_ob@meta.data[["clusters"]],levels = sort(unique(seurat_ob@meta.data[["clusters"]])))

col_anno = seurat_ob@meta.data[,c("clusters","sampleid","group")] %>% rownames_to_column(var="barcodes") %>% arrange(group) %>% arrange(sampleid) %>% arrange(clusters) %>% column_to_rownames(var="barcodes")
plot_data = seurat_ob@assays$RNA@scale.data[unique(as.vector(topn_markers$gene)),rownames(col_anno)]
library(ComplexHeatmap)
library(circlize)

#top5 gene
markers2vis = read.delim( "/public/scRNA/works/tangxuan/HT2019-6928-human/houxu20210725/2.sub_Fibroblasts_supply/Diffexp/clusters_5-vs-6-diff-pval-0.05-FC-1.5_anno.xls", sep="\t", header = T,quote="")
    markers2vis = subset(markers2vis,!grepl("^(mt-|Rps|Rpl|MT-|RPS|RPL)",markers2vis$gene))
    up_5 = filter(markers2vis,FoldChange > 1) %>% arrange(desc(log2FoldChange ))  %>% top_n(5,log2FoldChange )
    down_5 = filter(markers2vis,FoldChange < 1) %>% arrange(log2FoldChange )  %>% top_n(as.numeric(paste0("-",5)),log2FoldChange )
    topn_markers_5 = rbind(up_5, down_5)

index = as.vector(which(rownames(plot_data) %in% as.character(topn_markers_5$gene)))
index = c(index,c(45,53,85,132,166))
genelist = rownames(plot_data[index,])

lab = rowAnnotation(ano = anno_mark(at = index,labels = genelist,labels_gp = gpar(fontsize = 8)))
col_fun = colorRamp2(c(-1,0,5), c("#406AA8","white","#D91216"))

sampleid=as.vector(CustomCol2(3:12))
names(sampleid) = as.vector(sort(unique(col_anno$sampleid)))
group = as.vector(CustomCol2(13:14))
names(group) = as.vector(sort(unique(col_anno$group)))
ha = HeatmapAnnotation(clusters = as.character(col_anno$clusters),
                        sampleid = as.character(col_anno$sampleid),
                        group = as.character(col_anno$group),
    col = list(clusters = c("5" = CustomCol2(1), "6" = CustomCol2(2)),
                sampleid = sampleid,
                group = group))


#palette <- colorRampPalette(c("#406AA8", "white", "#D91216"))(n=299)
pdf("heatmap.pdf",width=6,height=6)
Heatmap(plot_data,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    right_annotation = lab,
    col = col_fun,
    name=" ",
    top_annotation = ha)
dev.off()

##arrange B2N
test = plot_data[c(1:100,166,101:165,167:366),]
index = as.vector(which(rownames(test) %in% as.character(topn_markers_5$gene)))
index = c(index,c(45,53,85,101,133))
genelist = rownames(test[index,])

lab = rowAnnotation(ano = anno_mark(at = index,labels = genelist,labels_gp = gpar(fontsize = 8)))
col_fun = colorRamp2(c(-1,0,5), c("#406AA8","white","#D91216"))

sampleid=as.vector(CustomCol2(3:12))
names(sampleid) = as.vector(sort(unique(col_anno$sampleid)))
group = as.vector(CustomCol2(13:14))
names(group) = as.vector(sort(unique(col_anno$group)))
ha = HeatmapAnnotation(clusters = as.character(col_anno$clusters),
                        sampleid = as.character(col_anno$sampleid),
                        group = as.character(col_anno$group),
    col = list(clusters = c("5" = CustomCol2(1), "6" = CustomCol2(2)),
                sampleid = sampleid,
                group = group))


#palette <- colorRampPalette(c("#406AA8", "white", "#D91216"))(n=299)
pdf("heatmap2.pdf",width=6,height=6)
Heatmap(test,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    right_annotation = lab,
    col = col_fun,
    name=" ",
    top_annotation = ha)
dev.off()

