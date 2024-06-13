#3. 出一张类似下图的比例图，横坐标是每个celltype中有克隆型细胞占总celltype细胞数的比例，纵坐标是对应celltype中，两组存在共享克隆型细胞的比例
#下载数据
obs://oe-scrna/works/scRNA/PROJECT/DOE202211093_Human/liuxuan/DOE202211093_Human/VDJ_Report/Clonotypes/merged_vdj_contig_annotation.xls

data = read.delim("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/3.Clonotype_dotplot/merged_vdj_contig_annotation.xls")
#先将表格筛选一下
data_ob = readRDS("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230629/PBMC_MAIT_sc.rds")
DefaultAssay(data_ob) = "RNA"
data_ob@meta.data$rawbc = rownames(data_ob@meta.data)
data_ob@meta.data$rawbc = gsub('-1_2', '-2',data_ob@meta.data$rawbc)
data_ob@meta.data$rawbc = gsub('-1_1', '-1',data_ob@meta.data$rawbc)

data$barcode = gsub('HD_TCR-', '',data$barcode)
data$barcode = gsub('sepsis_TCR-', '',data$barcode)
data$sampleid = gsub('HD_TCR', 'HD',data$sampleid)
data$sampleid = gsub('sepsis_TCR', 'sepsis',data$sampleid)

data[which(data$sampleid =="HD"),"barcode"] = paste(data[which(data$sampleid =="HD"),"barcode"],"-1",sep = "")
data[which(data$sampleid =="sepsis"),"barcode"] = paste(data[which(data$sampleid =="sepsis"),"barcode"],"-2",sep = "")

data_ob@meta.data$TCR= "No_TCR"

data_ob@meta.data$TCR = as.character(data_ob@meta.data$TCR)
for (i in intersect(data$barcode, data_ob@meta.data$rawbc)){
data_ob@meta.data[which(data_ob@meta.data$rawbc ==i),"TCR"]= data[which(data$barcode ==i),"clonotype_id"]
}

#横坐标是每个celltype中有克隆型细胞占总celltype细胞数的比例，纵坐标是对应celltype中，两组存在共享克隆型细胞的比例
df = data_ob@meta.data %>% select(celltype,orig.ident,TCR,rawbc) %>% rename("group" =orig.ident)

result_df <- df %>% group_by(celltype) %>% summarize(Clonal_cells = sum(TCR != "No_TCR"),Cell_number = n())

#计算两个group(HD和sepsis)中除了NO_TCR以外，相同TCR的个数，返回数据框
shared_TCR <- df %>% filter(TCR != "No_TCR") %>%  group_by(celltype) %>% summarize(num_shared_TCR = n())
shared_TCR  = shared_TCR %>% group_by(celltype) %>% summarize(num_shared_TCR = sum(num_shared_TCR))

#即为共享克隆型
TCR_df = df %>% filter(TCR != "No_TCR")

li_HD = list()
li_sepsis = list()
for (i in unique(TCR_df$celltype)){
li_HD[[i]] = TCR_df[which(TCR_df$celltype ==i & TCR_df$group =="HD"),"TCR"]
li_sepsis[[i]] = TCR_df[which(TCR_df$celltype ==i & TCR_df$group =="sepsis"),"TCR"]
}

clonotype_sum = list()
for (i in unique(TCR_df$celltype)){
clonotype_sum[[i]] = intersect(li_HD[[i]],li_sepsis[[i]])
}

#将list转换为dataframe格式
clonotype_sum <- data.frame( celltype = rep(names(clonotype_sum), lengths(clonotype_sum)), clonotype = unlist(clonotype_sum, use.names = FALSE) )

#按照celltype分组，计算clonotype_sum数据框中clonotype列与df数据框中对应celltype和TCR
data_df = data.frame()
n = 1
for (i in unique(clonotype_sum$celltype)){
sub_df = df[which(df$celltype ==i),]
sub_clonotype_sum = clonotype_sum[which(clonotype_sum$celltype ==i),]
data_df[n,1] = i
data_df[n,2] = sum(sub_df$TCR %in% sub_clonotype_sum$clonotype)
n = n+1
}
colnames(data_df)[1] = "celltype"
colnames(data_df)[2] = "Unique_clonotype_cell_numbers"

#把没有共享克隆型的细胞数量加回去
for ( i in setdiff(unique(df$celltype), unique(clonotype_sum$celltype))){
data_df[dim(data_df)[1]+1,1] = i
data_df[dim(data_df)[1],2] = 0
}

#把result_df和data_df合并到一起
data = data_df %>% left_join(result_df,by = "celltype")
data$Shared_clonal_cells_in_each_celltype = data$Unique_clonotype_cell_numbers / data$Cell_number
data$Clonal_cells_in_each_celltype = data$Clonal_cells / data$Cell_number

nlevel = length(unique(data$celltype))
CustomCol2 <- function(n){
  my_palette=c(
"#F57D06","#A05734","#B4794D","#C79C65","#D8C17D","#61438D","#B8A3C4","#D04545","#93BDD1","#5993B2","#33709F","#489544","#639F95","#9EC989")
  return(my_palette[n])
}
p = ggplot(data = data, mapping = aes(x =Clonal_cells_in_each_celltype*100, y = Shared_clonal_cells_in_each_celltype*100,colour= celltype)) + geom_point(size = 6)+ scale_colour_manual( values = CustomCol2(1:nlevel))+theme_classic()+  labs( x ="Clonal cells in each celltype(%)", y = "Shared clonal cells in each celltype (%)")+  theme( axis.title =element_text(size = 14),axis.text =element_text(size = 16) ,legend.text = element_text(size = 14),legend.title=element_text(size=16))+ theme(panel.background = element_rect( size = 3))

setwd("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/3.Clonotype_dotplot")
ggsave("Clonotype_dotplot.png",p,width = 10,height = 7)
ggsave("Clonotype_dotplot.pdf",p,width = 10,height = 7)
write.csv(data,file="data.csv",quote=F, row.names = F)
