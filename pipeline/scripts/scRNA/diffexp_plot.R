suppressPackageStartupMessages( library("dplyr") )
suppressPackageStartupMessages( library("ggplot2") )

### 1.2 

Sign67 <-  read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2018H2930S_HT2020-13912-h/SA20210604/diff_FC_1.2/clusters_6-vs-7-diff-padj-1e-04-FC-1.2_anno.xls",sep="\t",header=T,quote="",row.names=1)

all67 <- read.delim("/public/scRNA/works/tangxuan/backup_bendi/done/OE2018H2930S/houxu20200619/Diffexp/clusters_6-vs-7-all_diffexp_genes_anno.xls",sep="\t",header=T,quote="",row.names=1)

Sign37 <-  read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2018H2930S_HT2020-13912-h/SA20210604/diff_FC_1.2/clusters_3-vs-7-diff-padj-1e-04-FC-1.2_anno.xls",sep="\t",header=T,quote="",row.names=1)

all37 <- read.delim("/public/scRNA/works/tangxuan/backup_bendi/done/OE2018H2930S/houxu20200619/Diffexp/clusters_3-vs-7-all_diffexp_genes_anno.xls",sep="\t",header=T,quote="",row.names=1)

### 显著差异的 both， 各自 uniq 
both_sign <- intersect( rownames(Sign67),rownames(Sign37) )
sign_unique_67 <- setdiff( rownames(Sign67),rownames(Sign37) )
sign_unique_37 <- setdiff( rownames(Sign37),rownames(Sign67) )

### 同时都是非显著的差异基因
both_all <- intersect(rownames(all67),rownames(all37))
union_sign <-  union(rownames(Sign67),rownames(Sign37))
both_no_sign <- setdiff( both_all,union_sign )

### 数据集的总和 ###
union_sign <-  union(rownames(Sign67),rownames(Sign37))
all_plot_gene <- union(union_sign,both_no_sign) 

# all_diff <- rbind(all37,all67)

#all671 <- all67[all_plot_gene,] %>% select(log2FoldChange) %>% dplyr::rename("log2FoldChange_clusters_6-vs-7"="log2FoldChange")
all671 <- all67[all_plot_gene,] 
#all371 <- all37[all_plot_gene,] %>% select(log2FoldChange) %>% dplyr::rename("log2FoldChange_clusters_3-vs-7"="log2FoldChange")
all371 <- all37[all_plot_gene,]

plot <- as.data.frame(all_plot_gene) %>% mutate(type=1) %>%  tibble::column_to_rownames("all_plot_gene")
plot$`log2FoldChange_clusters_6_vs_7`=""
plot$`log2FoldChange_clusters_6_vs_7` = all671$log2FoldChange
plot$`log2FoldChange_clusters_3_vs_7`=""
plot$`log2FoldChange_clusters_3_vs_7` =  all371$log2FoldChange

### 
plot$Sign=""
plot[which(rownames(plot) %in% both_sign),"Sign"] ="Sign_different_both_comparions"
plot[which(rownames(plot) %in% sign_unique_67),"Sign"] ="Sign_different_in_clusters_6_vs_7"
plot[which(rownames(plot) %in% sign_unique_37),"Sign"] ="Sign_different_in_clusters_3_vs_7"
plot[which(rownames(plot) %in% both_no_sign),"Sign"] ="no_Sign_different_both_comparions"

plot_data <- plot[,2:4]

plot_data$Sign=factor(plot_data$Sign,levels=c("Sign_different_both_comparions","Sign_different_in_clusters_6_vs_7","Sign_different_in_clusters_3_vs_7","no_Sign_different_both_comparions"))

# pp = ggplot(plot_data) +
# 		 geom_point(aes(x=log2FoldChange_clusters_6_vs_7,
# 						y=log2FoldChange_clusters_3_vs_7,
# 						color=Sign),size=0.5)+
# 		 scale_color_manual(values=c("#3852A1","#44B549","#F02D19","grey"))+
# 		 theme_bw() +
# 		 theme(text = element_text(size=14),
# 				plot.title = element_text(hjust = 0.5,
# 										  size = 24),
# 				axis.title = element_text(size = 18),
# 				panel.grid = element_blank()) + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(0, 15))

pp = ggplot(plot_data) +
		 geom_point(aes(x=log2FoldChange_clusters_6_vs_7,
						y=log2FoldChange_clusters_3_vs_7,
						color=Sign))+
		 scale_color_manual(values=c("#3852A1","#44B549","#F02D19","grey"))+
		 theme_bw() +
		 theme(text = element_text(size=14),
				plot.title = element_text(hjust = 0.5,
										  size = 24),
				axis.title = element_text(size = 18),
				panel.grid = element_blank())  +scale_x_continuous(limits = c(-2, 4))
                
###  +scale_x_continuous(limits = c(-4, 6)) + scale_y_continuous(limits = c(-4, 6))

output_dir="./"
ggsave(file.path(output_dir,"plot_FC1.2.pdf"), width=11,plot= pp)
ggsave(file.path(output_dir,"plot_FC1.2.png"), width=11,plot= pp, dpi=1000)

data <- plot_data %>% tibble::rownames_to_column("gene")

write.table(data,file.path(output_dir,"plot_data_FC1.2.xls"),sep="\t",quote=F,row.names=F)



#### FC 1.1

Sign67 <-  read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2018H2930S_HT2020-13912-h/SA20210604/diff_FC_1.1/clusters_6-vs-7-diff-padj-1e-04-FC-1.1_anno.xls",sep="\t",header=T,quote="",row.names=1)

all67 <- read.delim("/public/scRNA/works/tangxuan/backup_bendi/done/OE2018H2930S/houxu20200619/Diffexp/clusters_6-vs-7-all_diffexp_genes_anno.xls",sep="\t",header=T,quote="",row.names=1)

Sign37 <-  read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2018H2930S_HT2020-13912-h/SA20210604/diff_FC_1.1/clusters_3-vs-7-diff-padj-1e-04-FC-1.1_anno.xls",sep="\t",header=T,quote="",row.names=1)

all37 <- read.delim("/public/scRNA/works/tangxuan/backup_bendi/done/OE2018H2930S/houxu20200619/Diffexp/clusters_3-vs-7-all_diffexp_genes_anno.xls",sep="\t",header=T,quote="",row.names=1)

### 显著差异的 both， 各自 uniq 
both_sign <- intersect( rownames(Sign67),rownames(Sign37) )
sign_unique_67 <- setdiff( rownames(Sign67),rownames(Sign37) )
sign_unique_37 <- setdiff( rownames(Sign37),rownames(Sign67) )

### 同时都是非显著的差异基因
both_all <- intersect(rownames(all67),rownames(all37))
union_sign <-  union(rownames(Sign67),rownames(Sign37))
both_no_sign <- setdiff( both_all,union_sign )

### 数据集的总和 ###
union_sign <-  union(rownames(Sign67),rownames(Sign37))
all_plot_gene <- union(union_sign,both_no_sign) 

# all_diff <- rbind(all37,all67)

#all671 <- all67[all_plot_gene,] %>% select(log2FoldChange) %>% dplyr::rename("log2FoldChange_clusters_6-vs-7"="log2FoldChange")
all671 <- all67[all_plot_gene,] 
#all371 <- all37[all_plot_gene,] %>% select(log2FoldChange) %>% dplyr::rename("log2FoldChange_clusters_3-vs-7"="log2FoldChange")
all371 <- all37[all_plot_gene,]

plot <- as.data.frame(all_plot_gene) %>% mutate(type=1) %>%  tibble::column_to_rownames("all_plot_gene")
plot$`log2FoldChange_clusters_6_vs_7`=""
plot$`log2FoldChange_clusters_6_vs_7` = all671$log2FoldChange
plot$`log2FoldChange_clusters_3_vs_7`=""
plot$`log2FoldChange_clusters_3_vs_7` =  all371$log2FoldChange

### 
plot$Sign=""
plot[which(rownames(plot) %in% both_sign),"Sign"] ="Sign_different_both_comparions"
plot[which(rownames(plot) %in% sign_unique_67),"Sign"] ="Sign_different_in_clusters_6_vs_7"
plot[which(rownames(plot) %in% sign_unique_37),"Sign"] ="Sign_different_in_clusters_3_vs_7"
plot[which(rownames(plot) %in% both_no_sign),"Sign"] ="no_Sign_different_both_comparions"

plot_data <- plot[,2:4]

plot_data$Sign=factor(plot_data$Sign,levels=c("Sign_different_both_comparions","Sign_different_in_clusters_6_vs_7","Sign_different_in_clusters_3_vs_7","no_Sign_different_both_comparions"))

# pp = ggplot(plot_data) +
# 		 geom_point(aes(x=log2FoldChange_clusters_6_vs_7,
# 						y=log2FoldChange_clusters_3_vs_7,
# 						color=Sign),size=0.5)+
# 		 scale_color_manual(values=c("#3852A1","#44B549","#F02D19","grey"))+
# 		 theme_bw() +
# 		 theme(text = element_text(size=14),
# 				plot.title = element_text(hjust = 0.5,
# 										  size = 24),
# 				axis.title = element_text(size = 18),
# 				panel.grid = element_blank()) + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(0, 15))

pp = ggplot(plot_data) +
		 geom_point(aes(x=log2FoldChange_clusters_6_vs_7,
						y=log2FoldChange_clusters_3_vs_7,
						color=Sign))+
		 scale_color_manual(values=c("#3852A1","#44B549","#F02D19","grey"))+
		 theme_bw() +
		 theme(text = element_text(size=14),
				plot.title = element_text(hjust = 0.5,
										  size = 24),
				axis.title = element_text(size = 18),
				panel.grid = element_blank())  +scale_x_continuous(limits = c(-2, 4))
                
##+scale_x_continuous(limits = c(-4, 6)) + scale_y_continuous(limits = c(-4, 6))

output_dir="./"
ggsave(file.path(output_dir,"plot_FC1.1.pdf"), width=11,plot= pp)
ggsave(file.path(output_dir,"plot_FC1.1.png"), width=11,plot= pp, dpi=1000)

data <- plot_data %>% tibble::rownames_to_column("gene")

write.table(data,file.path(output_dir,"plot_data_FC1.1.xls"),sep="\t",quote=F,row.names=F)
