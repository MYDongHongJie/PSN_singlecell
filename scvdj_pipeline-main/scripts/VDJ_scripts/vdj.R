suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2"))
#suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

option_list = list(
    make_option(c("-i", "--input"), type="character",
          help="The input exprssion matrix in RDS format."),
    make_option(c("-c", "--contifile"),type="character",
           help="The vdj merge anno file."),
    make_option( c("--reduct","-r"), type ="character",default="tsne",
                help="reduction used to visualize."),
    make_option( c("--topn"), type="integer", default = 10,
                 help = "[OPTIONAL] the number of top clonotype"),
    make_option( c("--type"), type="character", default = "TCR",
                 help = "[OPTIONAL] the choose TCR or BCR"),
    make_option( c("--groupby", "-g"), type = "character", default = "clusters",
                help = "[OPTIONAL]The variable in the metadata used to split the graph by the variable levels to
                        comparing the gene expression difference in different levels."),
    make_option( c("--splitby", "-y"), type = "character", default = "clusters",
                help = "[OPTIONAL]The variable in the metadata used to split the graph by the variable levels to
                        comparing the gene expression difference in different levels."),
    make_option( c("--analysiscontent","-a"), type= "character",default="vlnplot,featureplot",
                 help = "the visulization methods for the marker genes of each cell cluster.
                 he methods can be topn_clontype_id,clontype_id"),
    make_option( c("--num","-n"), type= "integer",default= 3,
                 help = "the visulization methods for the marker genes of each cell cluster.
                 he methods can be topn_clontype_id,clontype_id"),
    make_option( c("--outdir","-o"),type="character", default = "./",
        help="the output directory of  results." )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if ( is.null(opt$outdir) ){
    output_dir = getwd()
} else {
    if ( file.exists(opt$outdir) ){
        output_dir = opt$outdir
    } else {
        output_dir = opt$outdir
        dir.create(output_dir, recursive= T )
    }
}

if ( is.null( opt$topn )){
    topn = 10
} else {
    topn = opt$topn
}

if ( is.null(opt$input) ){
    stop("the seurat object is NOT AVAILABLE!")
}else{
    seurat_ob = OESingleCell::ReadX(input =opt$input, informat = 'h5seurat',verbose = F)
    if ( seurat_ob@version < 3){
        seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
    }
}

if (dim(seurat_ob)[2] < 500){
  pointsize=1.5
}else{
  pointsize=0.5
}

if ( is.null(opt$analysiscontent) ){
    print("NO analysis content choosed,the default method vlnplot and featureplot will be used!")
    analysiscontent = c("topn_clonotype_id","clonotype_id","clonotype_id_num")
}else if( opt$analysiscontent  == "all" ){
    analysiscontent = c("topn_clonotype_id","clonotype_id","clonotype_id_num")
}else{
    analysiscontent = unlist(strsplit(opt$analysiscontent ,","))
}

if ( is.null( opt$splitby)){
    facetby = "clusters"
}else{
    facetby = unlist( strsplit( opt$splitby, ",", perl =T ) )
}

if ('rawbc' %in% colnames(seurat_ob@meta.data)) {
        Barcode_content = 'rawbc'
}else{
        Barcode_content = 'orig.ident'
}   

clone=read.delim(opt$contifile,sep='\t')
clone$barcode <- gsub(paste0("_",opt$type,"-"),"-", clone$barcode )
rownames(clone)=clone$barcode

seurat_ob@meta.data$clonotype = paste0("non-",opt$type)
seurat_ob@meta.data[intersect(rownames(clone), rownames(seurat_ob@meta.data)),]$clonotype = opt$type

seurat_ob@meta.data[,"clonotype"]=factor(seurat_ob@meta.data[,"clonotype"],levels=c(opt$type,paste0("non-",opt$type)))
ggtsne = DimPlot(object = seurat_ob, reduction = opt$reduct,pt.size = pointsize,group.by="clonotype") + ggplot2::coord_fixed(ratio = 1) +
              ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::scale_color_manual( values =  c("red","lightgrey"))
ggplot2::ggsave(file.path(output_dir,paste0("clonotype.pdf")), ggtsne, width=7.5)
system(glue::glue("convert  -verbose -density 500 -trim {output_dir}/clonotype.pdf  -quality 100  -flatten  {output_dir}/clonotype.png"))
        #ggplot2::ggsave(file.path(output_dir1,paste0("top_",topn,"_clonotype_id",".png")), ggtsne, width = 7.5, dpi=500) 


if ("topn_clonotype_id" %in% analysiscontent){
    if ( file.exists("topn_clonotype_id") ){
        output_dir1 = file.path(output_dir,"topn_clonotype_id") 
    } else {
        output_dir1 = file.path(output_dir,"topn_clonotype_id") 
        dir.create(output_dir1, recursive= T )
    }
    if (!("topn_clonotype_id" %in% colnames(seurat_ob@meta.data))){
        #挑选topn
        clon1 = clone  %>%  dplyr::filter(clonotype_id %in% paste0("clonotype",c(1:topn)) ) 
        #单细胞转录组映射
        seurat_ob@meta.data$topn_clonotype_id = "Unselected"
        seurat_ob@meta.data[intersect(rownames(clon1), rownames(seurat_ob@meta.data)),]$topn_clonotype_id = as.character(clon1[intersect(rownames(clon1), rownames(seurat_ob@meta.data)),"clonotype_id"])
        seurat_ob@meta.data$topn_clonotype_id = factor(as.character(seurat_ob@meta.data$topn_clonotype_id)  , levels=c( paste0("clonotype",sort(as.numeric(gsub( "clonotype","", unique(seurat_ob@meta.data$topn_clonotype_id)[-which(unique(seurat_ob@meta.data$topn_clonotype_id)=="Unselected")])) ) ) ,"Unselected") ) 
       #simplified_meta = seurat_ob@meta.data %>%
       #                 dplyr::rename( "Barcode" = Barcode_content) %>%
       #                 dplyr::select( Barcode, sampleid, clusters,group,!!facetby,topn_clonotype_id)
       # write.table(simplified_meta, quote = F,sep =",",row.names = F,
       #      file.path(output_dir1,paste0("top_",topn,"_clonotype_id",".metadata.csv",collapse = "")))
    }
    nlevel = length(unique(seurat_ob@meta.data[,"topn_clonotype_id"]))
    #绘图      
    ggtsne = DimPlot(object = seurat_ob, reduction = opt$reduct,pt.size = pointsize,group.by="topn_clonotype_id") + ggplot2::coord_fixed(ratio = 1) + 
              ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::scale_color_manual( values =  c(CustomCol2(1:(nlevel-1) ),"lightgrey"))
        ggplot2::ggsave(file.path(output_dir1,paste0("top_",topn,"_clonotype_id",".pdf")), ggtsne, width=7.5)
        system(glue::glue("convert  -verbose -density 500 -trim {output_dir1}/top_{topn}_clonotype_id.pdf  -quality 100  -flatten  {output_dir1}/top_{topn}_clonotype_id.png"))
        #ggplot2::ggsave(file.path(output_dir1,paste0("top_",topn,"_clonotype_id",".png")), ggtsne, width = 7.5, dpi=500) 
        
    for (i in facetby){ 
        if (length(unique(seurat_ob@meta.data[,i]))>1) {
            nrow = ceiling(length(unique(seurat_ob@meta.data[,i]))/2)
            ggtsne = DimPlot(object = seurat_ob, reduction = opt$reduct,pt.size = pointsize,group.by="topn_clonotype_id",split.by=i,ncol = 2) +
                     ggplot2::coord_fixed(ratio = 1) +
                     ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
                     ggplot2::scale_color_manual( values =  c(CustomCol2(1:(nlevel-1) ),"lightgrey"))
            ggplot2::ggsave(file.path(output_dir1,paste0("top_",topn,"_splitby_",i,"_clonotype_id",".pdf")), ggtsne, width=14,height=5*nrow)
            ggplot2::ggsave(file.path(output_dir1,paste0("top_",topn,"_splitby_",i,"_clonotype_id",".png")), ggtsne, width = 14,height=5*nrow, dpi=1000)   
        }
    }
}

if ("clonotype_id" %in% analysiscontent){
    seurat_ob@meta.data$clonotype_id = "Unselected"
    seurat_ob@meta.data[intersect(rownames(clone), rownames(seurat_ob@meta.data)),]$clonotype_id = as.character(clone[intersect(rownames(clone), rownames(seurat_ob@meta.data)),"clonotype_id"])

    simplified_meta = seurat_ob@meta.data %>%
                        dplyr::rename( "Barcode" = Barcode_content) %>%
                        dplyr::select( Barcode, clonotype_id,sampleid, clusters,group,!!facetby)
    #venn
    if(length(unique(seurat_ob@meta.data[,i]))>1){
        if ( file.exists("venn") ){
            output_dir2 = file.path(output_dir,"venn")
        } else {
            output_dir2 = file.path(output_dir,"venn")
            dir.create(output_dir2, recursive= T )
        }
        for (i in facetby){
            if ( file.exists(file.path(output_dir,"venn",i) )){
                output_dir_venn = file.path(output_dir,"venn",i)
            } else {
                output_dir_venn = file.path(output_dir,"venn",i)
                dir.create(output_dir_venn, recursive= T )
            }
         
            for ( j in unique(seurat_ob@meta.data[,i]) ) {
                barcode=seurat_ob@meta.data[which(seurat_ob@meta.data[,i]==j),Barcode_content]
                subset <- subset(simplified_meta, simplified_meta[,"Barcode"] %in% barcode ) %>% dplyr::filter( clonotype_id != "Unselected" ) %>% dplyr::distinct(clonotype_id, .keep_all=T)
            
                write.table(subset, file.path(output_dir_venn,paste0(i,"_",j,"_clonotype_id.xls")), sep="\t", row.names=F,quote=F)
            }
            print("venn...")
            file <- paste(i,"_",unique(seurat_ob@meta.data[,i]),"_clonotype_id.xls",sep="",collapse=" ")
            files <- paste(output_dir_venn,paste0(i,"_",unique(seurat_ob@meta.data[,i]),"_clonotype_id.xls"),sep="/",collapse=" ")
            print(files)
            names=paste(unique(seurat_ob@meta.data[,i]),sep="",collapse=" ")
            print(names)
            ncol=paste(rep(2,times=length(unique(seurat_ob@meta.data[,i]))),sep="",collapse=" ")
            print(ncol)
        #system(glue::glue("echo -e 'module purge && source /home/lipeng/miniconda3/bin/activate base && perl  /public/scRNA_works/works/lipeng/script/venn_graph.pl -i '{file}' -c '{ncol}' -o ./ -l '{names}' -f sets' >{output_dir_venn}/w.sh"))
            system(glue::glue("module purge && /home/lipeng/miniconda3/bin/perl  scripts/VDJ_scripts/venn_graph.pl -i {files} -c {ncol} -o {output_dir_venn} -l {names} -f sets "))
        }
    }

   #柱形图
    if ( file.exists("stat") ){
        output_dir3 = file.path(output_dir,"stat")
    } else {
        output_dir3 = file.path(output_dir,"stat")
        dir.create(output_dir3, recursive= T )
    }

    subset <- simplified_meta %>% dplyr::filter( clonotype_id != "Unselected" )

    groupby=opt$groupby  #group
    bartype="dodge"
 
    for (facetbyx in facetby) {
        RemoveDuplicates <- subset[,c(groupby, facetbyx,"clonotype_id")] %>%  dplyr::group_by( .dots= c(groupby, facetbyx)) %>% dplyr::distinct(clonotype_id, .keep_all=T) 
        DATA <- as.data.frame(RemoveDuplicates[,c(groupby, facetbyx)]) %>% table() %>% as.data.frame()  %>% dplyr::rename(clonotype_unique_number=Freq)

        write.table(DATA, file.path(output_dir3,paste0(groupby,"_", facetbyx,"clonotype_id_freq_info.xls")), sep="\t", row.names=F,quote=F)
        nlevel=length(unique(DATA[,groupby]))

        sum_all <- ggplot(DATA,aes_string(x=facetbyx,y="clonotype_unique_number",fill=groupby) )+
                            geom_bar(stat="identity", position = bartype ) +
                            labs(x=" ",y=paste0("clonotype_unique_number"), title="") +
                            theme(plot.title = element_text(hjust = 0.5),
                            panel.grid.major =element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust = 1),
                            axis.line = element_line(color = "black"))+
                            scale_y_continuous(expand = c(0,0)) + ## //这个可以去掉与X轴间隙
                            #scale_x_discrete(expand = c(0.5,0)) +
                            scale_fill_manual(values=CustomCol2(1:nlevel)) +
                            guides(fill=guide_legend(title=groupby))

            ggsave(file.path(output_dir3,paste0("groupby_",groupby,"_splitby_",facetbyx,"_summary_clonotype_id_freq_plot.pdf",collapse="")),plot=sum_all,width=10,height=5)
            ggsave(file.path(output_dir3,paste0("groupby_",groupby,"_splitby_",facetbyx,"_summary_clonotype_id_freq_plot.png",collapse="")),dpi=1000, plot = sum_all,width=10,height=5)  
         }
     #write.table(simplified_meta, quote = F,sep =",",row.names = F,
     #        file.path(output_dir3,paste0("clonotype_id",".metadata.csv",collapse = "")))
       
}

if ("clonotype_id_num" %in% analysiscontent){
    if (!("clonotype_id" %in% colnames(seurat_ob@meta.data))){
        seurat_ob@meta.data[intersect(rownames(clon), rownames(seurat_ob@meta.data)),]$clonotype_id = as.character(clone[intersect(rownames(clone), rownames(seurat_ob@meta.data)),"clonotype_id"])    
    }
    if ( file.exists("clonotype_id_num") ){
        output_dir4 = file.path(output_dir,"clonotype_id_num")
    } else {
        output_dir4 = file.path(output_dir,"clonotype_id_num")
        dir.create(output_dir4, recursive= T )
    }
    num <- as.data.frame(table(seurat_ob@meta.data$clonotype_id))
    seurat_ob@meta.data$clonotype_id_num <- paste0(">=",opt$num)
    seurat_ob@meta.data[which(seurat_ob@meta.data$clonotype_id=='Unselected'),'clonotype_id_num']="Unselected"
    for (i in seq(1,opt$num-1)){
        seurat_ob@meta.data[which(seurat_ob@meta.data$clonotype_id %in% num[which(num$Freq==i),'Var1'] ),'clonotype_id_num'] = i
    }
    seurat_ob@meta.data$clonotype_id_num <- factor(seurat_ob@meta.data$clonotype_id_num,levels=c("Unselected",seq(1,opt$num-1),paste0(">=",opt$num)))
    ggtsne = DimPlot(object = seurat_ob, reduction = opt$reduct, pt.size = 1, group.by="clonotype_id_num") + ggplot2::coord_fixed(ratio = 1) +
        ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_colour_manual( values = c("lightgrey",CustomCol2(1:opt$num)))
    ggplot2::ggsave(file.path(output_dir4,paste0("clonotype_id_num",".pdf")), ggtsne,width=7.5)
    system(glue::glue("convert  -verbose -density 500 -trim {output_dir4}/clonotype_id_num.pdf  -quality 100  -flatten  {output_dir4}/clonotype_id_num.png"))
    #ggplot2::ggsave(file.path(output_dir4,paste0("clonotype_id_num",".png")), ggtsne,width=7.5)

    #simplified_meta = seurat_ob@meta.data %>%
    #                    dplyr::rename( "Barcode" = Barcode_content) %>%
    #                    dplyr::select(Barcode,sampleid,clusters,group,!!facetby,clonotype_id_num)

    #write.table(simplified_meta, quote = F,sep =",",row.names = F,file.path(output_dir4,paste0("clonotype_id_num.metadata.csv")))

}

info=colnames(seurat_ob@meta.data)[which( colnames(seurat_ob@meta.data) %in% c("clonotype_id","topn_clonotype_id","clonotype_id_num"))]
simplified_meta = seurat_ob@meta.data %>% rownames_to_column("barcode") %>%
                        dplyr::select(barcode,Barcode_content,sampleid, clusters,group,clonotype,!!info)
#去掉两张表重复的信息
clone=clone[,-which(colnames(clone) %in% c("sampleid","clonotype_id"))]
data=left_join(simplified_meta,clone,by="barcode")

write.table(data, quote = F,sep ="\t",row.names = F,file.path(output_dir,paste0("clonotype.metadata.xls")))
#保存rds
#OESingleCell::SaveX(seurat_ob, output = output_dir,
#                      outformat = "h5seurat", prefix = paste0("clonotype_id",Sys.Date()), update = FALSE )
#saveRDSMC(seurat_ob,file.path(output_dir,paste0("clonotype_id",Sys.Date(),".rds")))
