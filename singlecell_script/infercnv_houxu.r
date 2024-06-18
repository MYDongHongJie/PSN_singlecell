library(RColorBrewer)
library(infercnv)
#需求是删除Erythrocytes细胞，重出infercnv结果

infercnv_obj2 = readRDS("/PERSONALBIO/work/singlecell/s04/Analysis/PN20230328006/SPA2024042901/inferCNV/infercnv_out/run.final.infercnv_obj")

data_ob = readRDS("/PERSONALBIO/work/singlecell/s04/Analysis/PN20230328006/SPA2024042901/infer_data_Pre.rds")
data_ob = subset(data_ob,Barcode %in% colnames(infercnv_obj@expr.data))
#拿出对应细胞类型的barcode
infercnv_obj = infercnv_obj2
infercnv_obj@tumor_subclusters$subcluster$Erythrocytes = NULL
infercnv_obj@observation_grouped_cell_indices$Erythrocytes=NULL

barcc = names(infercnv_obj2@tumor_subclusters$subclusters$Erythrocytes$Erythrocytes_s1)
exprdata = infercnv_obj@expr.data[,!(colnames(infercnv_obj@expr.data) %in% barcc)]
infercnv_obj@expr.data =exprdata

##重新对应expr.data的位置
metadata = data_ob@meta.data %>% select(celltype,Barcode) %>% filter(celltype!="Erythrocytes")
Tc_meta = metadata[metadata$celltype == "TCells",]
#infercnv_obj@reference_grouped_cell_indices["TCells"]  = NULL
infercnv_obj@reference_grouped_cell_indices[["TCells"]] = which(colnames(infercnv_obj@expr.data) %in% Tc_meta$Barcode)


for (num in c(1:14)){
	Temp_meta = metadata[metadata$celltype == paste0("EpithelialCells_",num),]
	infercnv_obj@observation_grouped_cell_indices[[paste0("EpithelialCells_",num)]] = 
    which(colnames(infercnv_obj@expr.data) %in% Temp_meta$Barcode)
}

plot_cnv(infercnv_obj, k_obs_groups = 1, cluster_by_groups = TRUE, 
            cluster_references = TRUE, plot_chr_scale = FALSE, 
            chr_lengths = NULL, out_dir = out_dir, 
            x.range = "auto", title = "inferCNV", 
            output_filename = "infercnv", output_format = "png", 
            write_expr_matrix = FALSE, png_res = 300, 
            useRaster = TRUE)


plot_cnv(infercnv_obj, k_obs_groups = 1, cluster_by_groups = TRUE, 
            cluster_references = TRUE, plot_chr_scale = FALSE, 
            chr_lengths = NULL, out_dir = out_dir, 
            x.range = "auto", title = "inferCNV", 
            output_filename = "infercnv", output_format = "pdf", 
            write_expr_matrix = FALSE, png_res = 300, 
            useRaster = TRUE)
saveRDS(infercnv_obj,"run.final.infercnv_obj")


#infercnv_obj = readRDS(file.path(opt$lbj,"infercnv_out/run.final.infercnv_obj"))
seurat_obj = data_ob
expr <- infercnv_obj@expr.data


infercnv::plot_cnv(infercnv_obj, 
                   out_dir = "test",
                   plot_chr_scale = T, 
                   output_filename = "better_plot",output_format = "pdf", 
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2)))

infercnv::plot_cnv(infercnv_obj, 
                   plot_chr_scale = T, 
                   out_dir = "test",
                   output_filename = "better_plot",output_format = "png", 
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2)))
 