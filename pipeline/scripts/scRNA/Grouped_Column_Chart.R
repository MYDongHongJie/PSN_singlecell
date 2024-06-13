suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("OESingleCell"))

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
    make_option( c("--RDS", "-v"), type = "character", default = "TRUE",
                 help = "the seurat object saved as R object in RDS format."),
    make_option( c("--outdir","-o"),type="character", default = "./",
                help="the outdir directory of results.", metavar="character"),
    make_option( c("--groupby", "-g"), type = "character", default = NULL,
                help = "[OPTIONAL]The grouppinig variable in the metadata for,separate the cells to visulize marker genes."),
    make_option( c("--topby", "-c"), type = "character", default = NULL,
                 help="the column used to pick top n marker gene to visulize.The
                 option can be one of the column in the input marker genes table."),
    make_option( c("--splitby", "-y"), type = "character", default = NULL,
                help = "[OPTIONAL]The variable in the metadata used to split the graph by the variable levels to
                        comparing the gene expression difference in different levels."));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if ( is.null(opt$RDS) ){
    stop("the seurat object is NOT AVAILABLE!")
}else{
    seurat_ob = readRDSMC(opt$RDS)
    }

if ( is.null(opt$outdir) ){
    output_dir = getwd()
} else {
    if ( file.exists(opt$outdir) ){
        output_dir = opt$outdir
    } else {
        output_dir = opt$outdir
        dir.create(output_dir,recursive=T)
    }
}

if ( is.null( opt$groupby ) ){
    print( "NO groupping variable AVAILABLE for cell groupping! The default cell clusters id will be used!")
    groupby = "clusters"
}else{
    groupby = opt$groupby
}

if ( is.null( opt$topby ) ){
    print( "NO groupping variable AVAILABLE for cell groupping! The default cell group id will be used!")
    topby = "group"
}else{
    topby = opt$topby
}

if ( is.null( opt$splitby ) ){
    data = as.data.frame(seurat_ob@meta.data[,c(opt$groupby,opt$topby)]) %>% 
                    group_by(.dots = c(opt$groupby,opt$topby)) %>% 
                    dplyr::summarize(cell_number = n()) %>% 
                    mutate(freq = round((cell_number / sum(cell_number)) * 100,2))
    write.table(data,file.path(output_dir,'grouped_data1.xls'), sep = '\t',col.names=T, row.names =F)

    pp = ggplot(data, aes(x = 0, y = freq , fill = get(opt$topby))) + 
           labs(x = "", y = "", title = "")+coord_polar(theta = "y")+
           geom_bar(stat = "identity" , position = "stack",color = "white" , width = 0.1 ) + 
           theme_bw()+scale_fill_manual(values=CustomCol2(1:length(unique(seurat_ob@meta.data[[opt$topby]]))))+
           geom_text(stat="identity",color = "black",aes(x = 0.06,y = freq ,label = scales::percent(freq/100)), size=1, position=position_stack(vjust = 0.5))+
           facet_wrap(~get(opt$groupby),ncol = 2)+theme_void()+theme(legend.title = element_blank())

           ggsave(file.path(output_dir,'Grouped_Column_Chart.png'),pp)
           ggsave(file.path(output_dir,'Grouped_Column_Chart.pdf'),pp)
}else{
    splitby=opt$splitby
    DATA = as.data.frame(seurat_ob@meta.data[,c(opt$groupby,opt$topby,opt$splitby)]) %>% 
                    group_by(.dots = c(opt$groupby,opt$topby,opt$splitby)) %>% 
                    dplyr::summarize(cell_number = n()) %>% 
                    mutate(freq = (cell_number / sum(cell_number)) * 100)

    write.table(DATA,file.path(output_dir,'grouped_data.xls'), sep = '\t',col.names=T, row.names =F)

    p = ggplot(DATA , aes(x = get(opt$topby) , weight = freq ,fill = get(opt$splitby) ))+ 
                coord_flip()+geom_bar(position = "stack")+
                facet_wrap(~get(opt$groupby),ncol = 1,strip.position="top")+
                scale_fill_manual(values=CustomCol2(1:length(unique(seurat_ob@meta.data[[opt$splitby]]))))+ theme_minimal()+
                scale_y_continuous(expand = c(0,0))+ theme(axis.title = element_blank(),legend.title = element_blank())

                ggsave(file.path(output_dir,'Grouped_Column_Chart_plot.png'), p,height = 15)
                ggsave(file.path(output_dir,'Grouped_Column_Chart_plot.pdf'), p,height = 15)
}