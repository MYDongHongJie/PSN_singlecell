#路径"/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/2.CDR3"
#2. 每组CDR3长度的统计柱状图，例图如下，2张图
#思路：按照barcodes筛选表格，按照clonotype_aa列中的TRA和TRB拆分，TRA中有哪些AA，TRB中有哪些AA。根据AA的名字，匹配CDR3aa列，把每个unique AA的count数进行加和，然后除以总的counts数量获得频率。绘图
#下载数据
#obs://oe-scrna/works/scRNA/PROJECT/DOE202211093_Human/liuxuan/DOE202211093_Human/VDJ_Report/Clonotypes/HD_TCR.xls
#obs://oe-scrna/works/scRNA/PROJECT/DOE202211093_Human/liuxuan/DOE202211093_Human/VDJ_Report/Clonotypes/sepsis_TCR.xls


HD_table = read.delim("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/2.CDR3/HD_TCR.xls")


sepsis_table = read.delim("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/2.CDR3/sepsis_TCR.xls")


#先将表格筛选一下
data_ob = readRDS("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230629/PBMC_MAIT_sc.rds")
DefaultAssay(data_ob) = "RNA"
data_ob@meta.data$rawbc = rownames(data_ob@meta.data)
data_ob@meta.data$rawbc = gsub('-1_2', '-2',data_ob@meta.data$rawbc)
data_ob@meta.data$rawbc = gsub('-1_1', '-1',data_ob@meta.data$rawbc)


#HD是-1， sepsis是-2
HD_table$barcode = gsub('HD_TCR-', '',HD_table$barcode)
HD_table$barcode = paste(HD_table$barcode,"-1",sep = "")
sepsis_table$barcode = gsub('sepsis_TCR-', '',sepsis_table$barcode)
sepsis_table$barcode = paste(sepsis_table$barcode,"-2",sep = "")


HD_Barcodes = intersect(HD_table$barcode,data_ob@meta.data$rawbc)
HD_table= HD_table[HD_table$barcode %in%HD_Barcodes,]
setwd("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/2.CDR3")
write.csv(HD_table,file="HD_filtered_TCR.csv",quote=F, row.names = F)


sepsis_Barcodes = intersect(sepsis_table$barcode,data_ob@meta.data$rawbc)
sepsis_table = sepsis_table[sepsis_table$barcode %in%sepsis_Barcodes,]
setwd("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/2.CDR3")
write.csv(sepsis_table,file="sepsis_filtered_TCR.csv",quote=F, row.names = F)


clonotype_aa = HD_table$clonotype_aa
clonotype_aa  = append(clonotype_aa,sepsis_table$clonotype_aa)


#按照clonotype_aa列，对TRB进行拆分
trb_strings <- sapply(clonotype_aa , function(x) {
trb_substrs <- grep("TRB", unlist(strsplit(x, ";")), value = TRUE)
return(trb_substrs)
})


trb_strings = unlist(trb_strings)
names(trb_strings) = NULL
trb_strings = gsub('TRB:', '',trb_strings)


TRB = unique(trb_strings )


#按照clonotype_aa列，对TRA进行拆分
tra_strings <- sapply(clonotype_aa , function(x) {
tra_substrs <- grep("TRA", unlist(strsplit(x, ";")), value = TRUE)
return(tra_substrs)
})


tra_strings = unlist(tra_strings)
names(tra_strings) = NULL
tra_strings = gsub('TRA:', '',tra_strings)


TRA = unique(tra_strings)


#根据AA的名字，匹配CDR3aa列，把每个unique AA的count数进行加和，然后除以总的counts数量获得频率。绘图
HD_table$sampleid = "HD"
sepsis_table$sampleid = "sepsis"


HD_df = HD_table %>% select(count,CDR3aa,barcode,sampleid)
sepsis_df = sepsis_table %>% select(count,CDR3aa,barcode,sampleid)


df = rbind(HD_df,sepsis_df)
TRA = as.data.frame(TRA)
colnames(TRA)[1] = "CDR3aa"
TRA$Type = "TRA"


TRB = as.data.frame(TRB)
colnames(TRB)[1] = "CDR3aa"
TRB$Type = "TRB"


TRA_TRB = rbind(TRA,TRB)


df = df %>% left_join(TRA_TRB,by = "CDR3aa")


summarized_data <- df %>% group_by(sampleid, CDR3aa) %>% summarise(total_count = sum(count))


data = ungroup(summarized_data)
#将TRA和TRB的信息加回去
data= data %>% left_join(TRA_TRB,by = "CDR3aa")
data = data%>% group_by(sampleid,Type) %>%  mutate(frequency = total_count / sum(total_count))


data$CDR3_length = nchar(data$CDR3aa)


#针对TRA
Plot_data = data[which(data$Type =="TRA"),]
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}
nlevel = length(unique(Plot_data$sampleid))
p = ggplot(Plot_data ,aes(x=CDR3_length,y=100*frequency,fill=sampleid))+geom_bar(position="dodge",stat="identity")+scale_y_continuous(expand = c(0,0.001))+theme_classic()+  labs( x ="CDR3 length(aa)", y = "Frequency (%)")+  theme( axis.title =element_text(size = 14),axis.text =element_text(size = 16) ,legend.text = element_text(size = 14),legend.title=element_text(size=16)) + scale_fill_manual( values = CustomCol2(1:nlevel))+ scale_x_continuous(breaks = 1:max(Plot_data$CDR3_length))


setwd("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/2.CDR3")
write.csv(Plot_data,file="TRA_Plot_data.csv",quote=F, row.names = F)
ggsave("TRA.png",p)
ggsave("TRA.pdf",p)


#针对TRB
Plot_data = data[which(data$Type =="TRB"),]
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}
nlevel = length(unique(Plot_data$sampleid))
p = ggplot(Plot_data ,aes(x=CDR3_length,y=100*frequency,fill=sampleid))+geom_bar(position="dodge",stat="identity")+scale_y_continuous(expand = c(0,0.001))+theme_classic()+  labs( x ="CDR3 length(aa)", y = "Frequency (%)")+  theme( axis.title =element_text(size = 14),axis.text =element_text(size = 16) ,legend.text = element_text(size = 14),legend.title=element_text(size=16)) + scale_fill_manual( values = CustomCol2(1:nlevel))+ scale_x_continuous(breaks = 1:max(Plot_data$CDR3_length))




setwd("/public/scRNA_works/works/liuxuan/Project/conjoint/DOE202211093_Human/houxu20230808/2.CDR3")
ggsave("TRB.png",p)
ggsave("TRB.pdf",p)
write.csv(Plot_data,file="TRB_Plot_data.csv",quote=F, row.names = F)
