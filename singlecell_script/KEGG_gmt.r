library(dplyr)
library("optparse")
#library(jsonlite)

option_list <- list(
	make_option(c("-i", "--input"), help="Keggid"),
	make_option(c("-k", "--file"), help="kegg or go  file"),default=NULL,
	make_option(c("-n", "--name"), help="name of ouput file"),default=NULL,
	make_option(c("-t", "--type"), help="go or kegg"),default=NULL,
	make_option(c("-o", "--output"), help="output file"),default=NULL,
	make_option(c("-s", "--species"), help="species,only for go,we can choese （ssc,hsa,mmu,rno）",default=NULL)
)

opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
#读取文件

data = read.delim(opt$input,header = F)

##kegg
if (opt$type == 'kegg'){  
	if(is.null(opt$file)){
		print('if you select kegg,please input the kegg  file')
    stop()
	}else {
		 kegg_go_gmt = read.delim(opt$file)
	}
	temp = kegg_go_gmt[kegg_go_gmt$level3_pathway_id %in% data$V1,]
	gmt_process = temp %>% select(ko_name, level3_pathway_id, level3_pathway_name) %>% unique(.)
	gmt_process$name = paste0('(',gmt_process$level3_pathway_id,')',gmt_process$level3_pathway_name)
	gmtlist = list()
	for (Go in unique(gmt_process$name)){
		gmtlist[[Go]] = gmt_process$ko_name[gmt_process$name == Go]
	}
	filtered_list <- lapply(gmtlist, function(x) x[!grepl("\\s{1,}", x)])#删除含有空格的基因
}else if (opt$type == 'go') {
	library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(org.Ss.eg.db)
	if (is.null(opt$species)){
    print('Please input the species')
		stop()
  }else{
    species=opt$species
  }
  if( species=="ssc" ||species=="hsa" | species=="mmu" | species=="rno"){
		species.org=switch(species,
									"ssc"=org.Ss.eg.db,
									"hsa"=org.Hs.eg.db,
									"mmu"=org.Mm.eg.db,
									"rno"=org.Rn.eg.db,
					)
  }
	kegg_go_gmt$pathways = paste0(paste0("(",kegg_go_gmt$V1,")"),kegg_go_gmt$V2)
	k = kegg_go_gmt[,1]

	data = AnnotationDbi::select(species.org,keys = k,
	columns = c("GO","SYMBOL"),keytype="GO")

	gmt_process = data %>% dplyr::select(GO, SYMBOL) 
	gmt_process = gmt_process %>% unique()
	gmtlist = list()
	
	for (Go in unique(gmt_process$GO)){
		gmtlist[[Go]] = gmt_process$SYMBOL[gmt_process$GO == Go]
	}

	##排序且去除未匹配到的数据
	kegg_go_gmt <- kegg_go_gmt[!is.na(match(kegg_go_gmt$V1,names(gmtlist))), ]
	kegg_go_gmt <- kegg_go_gmt[order(match(kegg_go_gmt$V1, names(gmtlist))), ]	#对kegg_go_gmt进行排序，确保和list的names保持一致		
	##对list的names进行重新命名
	names(gmtlist) <- kegg_go_gmt$pathways[match(names(gmtlist), kegg_go_gmt$V1)]
	filtered_list = gmtlist
}
 

if(!is.null(opt$output)){
	path = './'
}else{
	path = opt$output
	dir.create(path,recursive = T)
}

filename = opt$name
file = file.path(path,filename)
output_gmt<-function(geneset,file){
sink(file) #将输出结果重定向到file
lapply(names(geneset),function(i){
#按照列名，将要连接的变量转化为向量型，用制表符连接成一个字符串
cat(paste(c(i,'tmp',geneset[[i]]),collapse='\t')) 
cat('\n') #输出后新起一行
})
sink() #结束重定向
}
output_gmt(filtered_list,file)
