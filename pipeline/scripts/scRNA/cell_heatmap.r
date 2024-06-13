library(Seurat)
library(dplyr)

seurat_ob = readRDS("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-22420-mihoutao/houxu-20220605/rds/data_ob_v3.rds")
# diff = read.delim("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-22420-mihoutao/HT2022-10204-1_add_MT/result/Diffexp/group_BAM1-vs-H1-all_diffexp_genes_anno.xls",sep="\t")
diff = read.delim("/public/scRNA_works/works/guokaiqi/project/scRNA/HT2021-22420-mihoutao/HT2022-10204-1_add_MT/result/Diffexp/group_BAM1-vs-H1-diff-pval-0.05-FC-1.5_anno.xls",sep="\t")



count = as.matrix(GetAssayData(seurat_ob, slot = "data")[diff$gene,])
meta.data = seurat_ob@meta.data
collapseby = "new_celltype_clusters"

meta.data$id = rownames(meta.data)
collapsed_count = vector()

collapsed_group = meta.data %>% group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))

for ( cells in collapsed_group$cellid ){
    samplex = unlist(strsplit(cells, ",", perl =T))
    collapsed_count= cbind(collapsed_count,rowMeans( count[,samplex,drop=F] ))
}
collapsed_count = as.matrix( collapsed_count )
collapsed_group = as.data.frame(collapsed_group)
colnames(collapsed_count) = as.matrix(collapsed_group[,1])

max_counts<-apply(collapsed_count, 1, function(x){colnames(collapsed_count)[which.max(x)]})

data = tibble::rownames_to_column(as.data.frame(collapsed_count),"gene")
data$max_counts = max_counts

data = left_join(data,diff,by ="gene")

# write.table(data, file.path("group_BAM1-vs-H1-all_diffexp_genes_anno.xls"),quote = F, row.names = F, sep = "\t")
write.table(data, file.path("group_BAM1-vs-H1-diff-pval-0.05-FC-1.5_anno.xls"),quote = F, row.names = F, sep = "\t")
##########################
genes=read.delim("genelist.txt")
diff=read.delim("group_BAM1-vs-H1-diff-pval-0.05-FC-1.5_anno.xls",sep="\t")
markers2vis=subset(diff,gene %in% genes$genes)
up=subset(markers2vis,up_down=="Up") %>% dplyr::arrange(desc(log2FoldChange ))
down=subset(markers2vis,up_down=="Down") %>% dplyr::arrange(log2FoldChange) 
topn_markers = rbind(up, down) 
markers2vis4heatmap = unique(as.vector(topn_markers$gene))
ggheat = Seurat::DoHeatmap( object = data_ob,  features = markers2vis4heatmap,group.by = "group",, group.bar = T, label = F, draw.lines = F)   
heatmap_data = ggheat$data
library(reshape2)
plot_data=dcast(heatmap_data,Feature~heatmap_data$Cell,value.var = "Expression")

metadata=data_ob@meta.data[,c("clusters","group","new_celltype_clusters")]
# > head(metadata)
                    # clusters group new_celltype_clusters
# H1-AAACCCAAGCTGACAG        5    H1      Mesophyll_cell_5
# H1-AAACCCACACTACGGC        4    H1            Guard_cell
# H1-AAACCCAGTCAACCTA        3    H1             Epidermis
# H1-AAACCCATCAATCAGC        5    H1      Mesophyll_cell_5
# H1-AAACGAAAGACAGCGT        4    H1            Guard_cell
# H1-AAACGAACAAAGGAGA        3    H1             Epidermis
meta=metadata %>% group_by(new_celltype_clusters,group)

library(pheatmap)
library(ggplot2)
meta$new_celltype_clusters=factor(meta$new_celltype_clusters,levels=sort(unique(meta$new_celltype_clusters)))
meta$new_celltype_clusters=factor(meta$new_celltype_clusters,levels=c("Mesophyll_cell_1","Mesophyll_cell_2", "Mesophyll_cell_5","Phloem_companion_cell","Phloem","Guard_cell","Epidermis","xylem"))
meta= meta[order(meta$new_celltype_clusters,meta$group),]
plot_data = plot_data[markers2vis4heatmap,rownames(meta)]
aa = as.vector(table(meta$new_celltype_clusters))
bb = ""
for (i in 1:length(aa)){bb = c(bb,sum(aa[1:i])) }
bb=bb[-c(1,9)]
bb=as.numeric(bb)

ann_colors=list( group=c(BAM1="#3fdacb", H1="firebrick"), 
                 new_celltype_clusters=c(Mesophyll_cell_1="#7fc97f",Mesophyll_cell_2="#beaed4",Mesophyll_cell_5="#fdc086",Epidermis="#386cb0",
                                         Guard_cell="#f0027f",xylem ="#a34e3b",Phloem="#666666",Phloem_companion_cell="#1b9e77") 
p=pheatmap(plot_data, annotation_col = meta ,cluster_rows = F,cluster_cols = F,show_colnames = F,gaps_col=bb,annotation_colors=ann_colors)
ggsave("heatmap.png",p,width=14)
ggsave("heatmap.pdf",p,width=14)
