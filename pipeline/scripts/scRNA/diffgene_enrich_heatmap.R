library(Seurat)
library(ggplot2)

S20207=readRDS("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-20207-qingxie/houxu-20230831/rds/sub_celltype2023-08-31.rds") 
sp=SplitObject(S20207,split.by="celltype")
embeddings = S20207@reductions$umap@cell.embeddings

color1=c('#F768A1', '#AE017E', '#7A0177', '#49006A')
p = DimPlot(sp[["Gill 1"]], reduction = "umap",pt.size = 0.6,group.by="sub_celltype")+
        ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_colour_manual( values =color1) + coord_fixed(ratio = 1)+
        ylim(c(min(embeddings[,2]),max(embeddings[,2]))) + xlim(c(min(embeddings[,1]),max(embeddings[,1])))+
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),axis.line=element_blank())
ggsave("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-20207-qingxie/houxu-20230831/Gill_1.png",width=7.2)
ggsave("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-20207-qingxie/houxu-20230831/Gill_1.pdf",width=7.2)


color1=c('#4EB3D3', '#084081')
p = DimPlot(sp[["Gill 3"]], reduction = "umap",pt.size = 0.6,group.by="sub_celltype")+
        ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_colour_manual( values =color1) + coord_fixed(ratio = 1) + 
        ylim(c(min(embeddings[,2]),max(embeddings[,2]))) + xlim(c(min(embeddings[,1]),max(embeddings[,1])))+
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),axis.line=element_blank())
ggsave("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-20207-qingxie/houxu-20230831/Gill_3.png",width=7.2)
ggsave("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-20207-qingxie/houxu-20230831/Gill_3.pdf",width=7.2)

color1=c('#FE9929', '#EC7014', '#CC4C02', '#993404', '#662506')
p = DimPlot(sp[["Gill 5"]], reduction = "umap",pt.size = 0.6,group.by="sub_celltype")+
        ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_colour_manual( values =color1) + coord_fixed(ratio = 1)+
        ylim(c(min(embeddings[,2]),max(embeddings[,2]))) + xlim(c(min(embeddings[,1]),max(embeddings[,1])))+
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),axis.line=element_blank())
ggsave("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-20207-qingxie/houxu-20230831/Gill_5.png",width=7.2)
ggsave("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-20207-qingxie/houxu-20230831/Gill_5.pdf",width=7.2)

############################
celltype=c("bulk", "Gill_1-1", "Gill_1-2", "Gill_1-3", "Gill_1-4", "Gill_2", "Gill_3-1", "Gill_3-2", "Gill_4", "Gill_5-1", "Gill_5-2", "Gill_5-3", "Gill_5-4", "Gill_5-5")
diff_sum = NULL
for (i in celltype){
    diff_stat = read.delim(paste0("Diffexp/",i,"/diffexp_results_stat.xls"),sep="\t")
    diff_stat$celltype = i
    diff_sum = rbind(diff_sum, diff_stat)
}
diff_sum$diff = paste0(diff_sum$Case,"-vs-",diff_sum$Control)

for (i in unique(diff_sum$diff)){
    plot_data = subset(diff_sum, diff == i)
    c = setdiff(celltype, plot_data$celltype)
    plot_data = reshape2::melt(plot_data, id.vars=c("diff","celltype"), measure.vars = c("Up_diff","Down_diff"),variable.name = "UP_DOWN",value.name = "NUM")
    complementary = data.frame(diff = rep(unique(plot_data$diff),length(c)*2), celltype = rep(c,2) ,UP_DOWN = c(rep("Up_diff",length(c)),rep("Down_diff",length(c))),NUM = 0 )
    plot_data = rbind(plot_data, complementary)

    p = ggplot(plot_data,aes(x=celltype , y=NUM ,fill = UP_DOWN)) +
        geom_bar(position = position_dodge(), stat = "identity", width = 0.6, size=1.2) +
        geom_text(aes(label=NUM), size=2, vjust=-0.5, position = position_dodge(0.6)) +
        theme_classic() + scale_y_continuous(expand=c(0,0),limits = c(0,ceiling(max(plot_data$NUM) / 100) * 100)) + 
        scale_fill_manual( values = c("#FFA500","#1873CC")) +
        theme(legend.position=c(0.85,0.85),
              legend.background = element_rect(fill = 'transparent'),
              legend.title = element_blank(),
              axis.text = element_text(color = "black"))
    ggsave(paste0("Diffexp/",i,"_barplot.png"),width=8)
    ggsave(paste0("Diffexp/",i,"_barplot.pdf"),width=8)

}
i="CS(sampleid)-vs-AS(sampleid)"
i="AS(sampleid)-vs-IS(sampleid)"
    plot_data = subset(diff_sum, diff == i)
    c = setdiff(celltype, plot_data$celltype)
    plot_data = reshape2::melt(plot_data, id.vars=c("diff","celltype"), measure.vars = c("Up_diff","Down_diff"),variable.name = "UP_DOWN",value.name = "NUM")
    complementary = data.frame(diff = rep(unique(plot_data$diff),length(c)*2), celltype = rep(c,2) ,UP_DOWN = c(rep("Up_diff",length(c)),rep("Down_diff",length(c))),NUM = 0 )
    plot_data = rbind(plot_data, complementary)

    p = ggplot(plot_data,aes(x=celltype , y=NUM ,fill = UP_DOWN)) +
        geom_bar(position = position_dodge(), stat = "identity", width = 0.6, size=1.2) +
        geom_text(aes(label=NUM), size=2, vjust=-0.5, position = position_dodge(0.6)) +
        theme_classic() + scale_y_continuous(expand=c(0,0),limits = c(0,ceiling(max(plot_data$NUM) / 100) * 100),
                                             labels = seq(0,ceiling(max(plot_data$NUM) / 100) * 100, 100),
                                             breaks = seq(0,ceiling(max(plot_data$NUM) / 100) * 100, 100)) + 
        scale_fill_manual( values = c("#FFA500","#1873CC")) +
        theme(legend.position=c(0.85,0.85),
              legend.background = element_rect(fill = 'transparent'),
              legend.title = element_blank(),
              axis.text = element_text(color = "black"))
    ggsave(paste0("Diffexp/",i,"_barplot.png"),width=8)
    ggsave(paste0("Diffexp/",i,"_barplot.pdf"),width=8)
###########################################################
library(dplyr)
library(ggplot2)
library(pheatmap)
#########
polychrome = c(
  "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE",
  "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD",
  "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
  "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0",
  "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5",
  "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB"
)
diff = unique(gsub("\\(sampleid\\)","",diff_sum$diff ))
for (i in diff){

enrich = NULL
terms = NULL
top10 = NULL
for (c in celltype){
file=paste0("Diffexp/",c,"/enrichment/KEGG_enrichment/sampleid_",i,"/enrichment-kegg-sampleid_",i,"-Total.xls")
if (file.exists(file)){
enrich1 = read.delim(file, head=TRUE, sep="\t", quote="")
enrich1 <- enrich1[which(enrich1["ListHits"]>2), ]
if(nrow(enrich1)==0){
    print(paste0(c," enrich1 KEGG terms is equal to zero, program exit!"))
} else { 
enrich1$celltype = c
enrich1[, "Term"] <- sub("^path:", "", paste(enrich1[,1], enrich1[,2], sep=": "))
enrich = rbind(enrich, enrich1)
top1 <- head(enrich1[order(enrich1[,"q.value"]),], 10)
top10 = rbind(top10, top1)
terms = c(terms, as.vector(top1$Term))
}
}}
top10 = subset(enrich, Term %in% terms)

output_dir = "Diffexp/"
write.table(top10, file.path(output_dir, paste0(i,"_KEGG.top.10.xls")), sep="\t",
	quote=FALSE, col.names=TRUE, row.names=FALSE)
#top10[, "term"] <- sub("^path:", "", paste(top10[,1], top10[,2], sep=": "))
top10$celltype = factor(top10$celltype, levels = celltype)
sp = split(top10, top10$celltype)
pval_list = lapply(sp, function(x){
                      x = subset(x, select = c("Term","q.value"))
                      return(x)
                      })
plot_data = Reduce(function(x,y) merge(x,y,by = "Term", all = T), pval_list)
plot_data = tibble::column_to_rownames(plot_data, var = "Term")
plot_data = plot_data[unique(terms),]
colnames(plot_data) = names(sp) 
# plot_data = -log10(plot_data)
annotation_row = data.frame(
  Term = top10$Term,
  Class = top10$Classification_level2
)
annotation_row = unique(annotation_row)
rownames(annotation_row) = NULL
annotation_row = tibble::column_to_rownames(annotation_row, var = "Term")
anno_colors = polychrome[1:length(unique(annotation_row$Class))]
names(anno_colors) = sort(unique(annotation_row$Class))
annotation_colors = list(anno_colors)
names(annotation_colors) = "Class"
######################
col = c(colorRampPalette(c('red', 'white'))(5), colorRampPalette(c('white', 'green'))(ceiling(quantile(plot_data,na.rm=T)[4])/0.01))
bk = c(seq(0,0.05,5), seq(0.05,ceiling(quantile(plot_data,na.rm=T)[4]),by=0.01))
p = pheatmap( plot_data,
                    scale = "none",
                    color = col,
                    show_colnames = T,
                    cluster_cols = F,
                    cluster_rows = F,
                    border_color = "black",
                    fontsize_row = 10 ,
		            cellheight = 14, cellwidth = 30, 
                    annotation_row = annotation_row, annotation_colors = annotation_colors,
		            legend_breaks=c(0,0.05, seq(0.05,ceiling(quantile(plot_data,na.rm=T)[4]),ceiling(quantile(plot_data,na.rm=T)[4])/5)),
                    breaks=bk,display_numbers = TRUE , number_format = "%.1e")
ggsave(file.path(output_dir, paste0(i,"_kegg.top", 10, "_enrich_term.pdf")), plot=p, width=15, height=16)
ggsave(file.path(output_dir, paste0(i,"_kegg.top", 10, "_enrich_term.png")), plot=p, width=15, height=16, dpi=700)

}
###############

dotplot=read.delim("dotplot.txt")
S20207@meta.data$new_celltype = factor(S20207@meta.data$new_celltype ,levels=sort(unique(S20207@meta.data$new_celltype )))
S20207 = SetIdent( S20207, value = "new_celltype" )
ggdots = Seurat::DotPlot(object = S20207, features = rev(dotplot$gene))  +
                              Seurat::RotatedAxis() + 
                              ggplot2::scale_colour_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100) ) +
                              theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),axis.line=element_blank(),
                                    plot.margin = unit(c(0,1,0,5),'mm'), axis.text =element_text(size=10))+
                              coord_fixed(ratio = 1) 

ggsave("featureplot_vlnplot/marker_gene_dotplot.png",width=10,bg="white") 
ggsave("featureplot_vlnplot/marker_gene_dotplot.pdf",width=10,bg="white") 
####################
dotplot=read.delim("dotplot2.txt")
S20207_sub = SubsetData(S20207, cells = S20207@meta.data$sub_celltype %in% c("Gill 1-1", "Gill 1-2", "Gill 1-3", "Gill 1-4", "Gill 3-1", "Gill 3-2", "Gill 5-1", "Gill 5-2", "Gill 5-3", "Gill 5-4", "Gill 5-5")) 
S20207_sub@meta.data$sub_celltype = factor(S20207_sub@meta.data$sub_celltype ,levels=sort(unique(S20207_sub@meta.data$sub_celltype )))
S20207_sub = SetIdent( S20207_sub, value = "sub_celltype" )

ggdots = Seurat::DotPlot(object = S20207_sub, features = unique(rev(as.vector(dotplot$genes))))  +
                              Seurat::RotatedAxis() + 
                              ggplot2::scale_colour_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100) ) +
                              theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),axis.line=element_blank(),
                                    plot.margin = unit(c(0,1,0,5),'mm'), axis.text =element_text(size=10))+
                              coord_fixed(ratio = 1) 

ggsave("featureplot_vlnplot/marker_gene_dotplot2.png",width=10,bg="white") 
ggsave("featureplot_vlnplot/marker_gene_dotplot2.pdf",width=10,bg="white") 
############
groupby = "new_celltype"
root_dir = "Marker"
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

markers2vis = read.delim("Marker/all_markers_for_each_cluster.xls",sep="\t")
markers2vis[["cluster"]] = factor(markers2vis[["cluster"]] , levels = sort(unique(S20207@meta.data[["new_celltype"]])))
topn_markers  = markers2vis %>% group_by(cluster) %>%
        arrange(p.value,desc(avg_log2FC),desc(gene_diff)) %>%
        arrange(cluster) %>% mutate(folder_suffix = paste0("cluster",cluster)) %>% select(cluster,gene,folder_suffix)

markers2vis4heatmap = unique(as.vector(topn_markers$gene))
subseted_seurat = S20207

subseted_seurat@meta.data[,groupby] = as.factor(subseted_seurat@meta.data[,groupby])
colors2use = CustomCol2(1:length(unique(subseted_seurat@meta.data[,groupby])))
        if (length(markers2vis4heatmap) > 135){
            sz = 4-log(length(markers2vis4heatmap)/100)
            heig = 5+log2(length(markers2vis4heatmap)/10)
            wid = 7.5
        }else if (length(markers2vis4heatmap) < 75){
            sz = 6-log2(length(markers2vis4heatmap)/80);heig = 7;wid = 7
        }else{
            sz = 4-log2(length(markers2vis4heatmap)/120);heig = 7;wid = 7
        }

        ggheat = DoHeatmap( object = subseted_seurat,
                            features = markers2vis4heatmap,
                            group.colors = colors2use,
                            group.by = groupby, group.bar = T, label = F) +
                            theme(axis.text.y = element_text(size = sz, face = "bold"))
                            # group.cex = 10, cex.row = 4,
                            # slim.col.label = T, group.label.rot = F)
        ggheat + guides(fill = guide_colorbar( title.position = "top", order = 1), color = guide_legend(order = 2, override.aes = list(alpha = 1)))
        ggsave(file.path(root_dir,paste0("top","marker_gene_heatmap.pdf", collapse = "_")),height = heig, width = wid)
        ggsave(file.path(root_dir, paste0("top", "marker_gene_heatmap.png", collapse = "_")),height = heig, width = wid, dpi = 1000 ,limitsize = F)
###############

sampled_cellmeta = S20207@meta.data %>% rownames_to_column() %>%
                group_by( .dots= groupby ) %>%
                sample_n( 200,replace = F) %>% column_to_rownames()
subseted_seurat = SubsetData(S20207, cells = rownames(sampled_cellmeta))

ggheat = DoHeatmap( object = subseted_seurat,
                            features = markers2vis4heatmap,
                            group.colors = colors2use,
                            group.by = groupby, group.bar = T, label = F) +
                            theme(axis.text.y = element_text(size = sz, face = "bold"))
                            # group.cex = 10, cex.row = 4,
                            # slim.col.label = T, group.label.rot = F)
ggheat + guides(fill = guide_colorbar( title.position = "top", order = 1), color = guide_legend(order = 2, override.aes = list(alpha = 1)))
ggsave(file.path(root_dir,paste0("top","marker_gene_heatmap_sub.pdf", collapse = "_")),height = heig, width = wid+2)
ggsave(file.path(root_dir, paste0("top", "marker_gene_heatmap_sub.png", collapse = "_")),height = heig, width = wid+2, dpi = 1000 ,limitsize = F)
#####################################################
diff = unique(gsub("\\(sampleid\\)","",diff_sum$diff ))
"NS-vs-CS" "NS-vs-AS" "NS-vs-IS" "CS-vs-AS" "CS-vs-IS" "AS-vs-IS"
for (i in celltype){
    for ( base_prefix in diff){
    outputdir = paste0("Diffexp/",i)
    if (file.exists(file.path(outputdir,paste0("sampleid_",base_prefix,"-all_diffexp_genes",".xls")))){
    res = read.delim(file.path(outputdir,paste0("sampleid_",base_prefix,"-all_diffexp_genes",".xls")), sep="\t")
    #filter the genes by the threshold specified by the user
    res_Significant = dplyr::filter(res, p.value < 0.05 , abs(log2FoldChange) > log2(2))
    if ( nrow(res_Significant) == 0 ){
        warning(paste0("NO Significant Differential Genes Identified for ", contrast,sep =""))
        return(0)
    }
    res_Significant[which(res_Significant$log2FoldChange  > 0), "up_down"] <- "Up"
    res_Significant[which(res_Significant$log2FoldChange  < 0), "up_down"] <- "Down"
    #write out the desired differential gene expression analysis results
    write.table(res_Significant,
                    file.path(outputdir,
                    paste0("sampleid_",base_prefix,"-diff-","pval-0.05-FC-2.xls")),
                    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="")
}}}



