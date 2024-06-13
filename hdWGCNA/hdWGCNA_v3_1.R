#module purge && source /home/liuxuan/miniconda3/bin/activate hdWGCNA_env

#Loading Packages
if (!requireNamespace('argparse', quietly = TRUE)) {
  stop("Please install argparse to parse the command line paramters!")
}
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(hdWGCNA))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(SeuratDisk, lib.loc = "/home/dongjiaoyang/miniconda3/envs/OESingleCell/lib/R/library"))
suppressPackageStartupMessages(library(UCell))
suppressPackageStartupMessages(library(qlcMatrix))
suppressPackageStartupMessages(library(corrplot, lib.loc = "/home/dongjiaoyang/miniconda3/envs/OESingleCell/lib/R/library"))
#加载PlotDMEsLollipop用到的ggforestplot包
suppressPackageStartupMessages(library(ggforestplot))


# ======================= COMMAND LINE PARAMETERS SETTING =======================
# ======================= GLOBAL command line parameters setting=================
parser = ArgumentParser(description = "hdWGCNA_pipeline",
                        usage = "%(prog)s [global options]" )
parser$add_argument("-i", "--input", type = "character",
             help = "The input exprssion matrix in several possible format.")
parser$add_argument("-f", "--informat", type = "character", default = "h5seurat",
             help = "The format of data object, the possible choices can be:h5seurat,(seurat)rds.[default: %(default)s]")
parser$add_argument("-o", "--output", type = "character", default = "./",
             help = "the output directory")
parser$add_argument("-j", "--ncores", type = "integer", default = 8,
             help = "the number of CPUs used to improve the performace.[default: %(default)s]")
parser$add_argument( "-g", "--groupby", type = "character", default = "clusters",
             help = "specify the columns in seurat_obj@meta.data to group by.  It should be a character vector of Seurat metadata column names representing groups for which metacells will be computed.")
parser$add_argument( "--ident", type = "character", default = "clusters",
             help = "set the Idents of the metacell seurat object. Tips: ident.group must be in group.by")
parser$add_argument( "--ident2use","-q", type = "character", default = NULL,
             help = "The name of the group of interest in the group.by column. 需要截取的分组变量，如有多个可用,分隔")
parser$add_argument( "--networkType","-n", type = "character", default = "signed",
             help = "you can also use unsigned or signed hybrid")
parser$add_argument( "--reduct1",type = "character",default = "pca",
             help = "the (primary) reduction methods used, will be used in metacell construction.")
parser$add_argument( "--reduct2",type = "character",default = "umap",
             help = "the results of computed reduction methods used as input for secondary reduction,choice can be tsne or umap.")      
parser$add_argument( "--soft_power",type = "character",default = NULL,
             help = "[optional] the soft power used for network construction. Automatically selected by default.") 
parser$add_argument( "--topn_hub",type = "integer",default = 10,
             help = "the number of hub genes by kME to select for each module.") 
parser$add_argument( "--Gene_scoring",type = "character",default = "UCell",
             help = " compute gene scoring for the topn hub genes by kME for each module. The possible choices can be: UCells or Seurat.") 
parser$add_argument( "--MetaCells",type = "character",default = "TRUE",
             help = "A logical determining if we use the metacells (TRUE) or the full expression matrix (FALSE)") 
parser$add_argument("--traits",type = "character",default = "nCount_RNA,nFeature_RNA",
             help = "List of traits to correlate")
parser$add_argument("--rds",type = "character",default = NULL,
             help = "hdWGCNA_object.rds")
parser$add_argument("--TOMFiles",type = "character",default = NULL,
             help = "hdWGCNA_TOM.rds file path")
subparsers = parser$add_subparsers(help = "subcommands:")
# Downstream Custome visualization ========
# create the parsers for subcommand create
sub_Customize = subparsers$add_parser("Customize",
     help = "Downstream Custome visualization")
sub_Customize$add_argument("-m", "--module", type = "character", default = NULL,
     help = "Plot customized module in Vlnplot")
sub_Customize$add_argument("--n_hubs", type = "integer", default = 10,
     help = "The number of hub genes to plot for each module.")
sub_Customize$add_argument("-c","--contrast", type = "character",default = NULL,
          help = "[Optional]levels of a factor used to compare with for final differenetial results. The format is Factor:interesting_levle:reference_level.")
opt = parser$parse_args()
args<-commandArgs(TRUE)
#================== GLOBAL PARAMETERS PARSING ===================================
# setting the output directory
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir,recursive = T)
    }
}
output_dir = normalizePath(output_dir )

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = opt$ncores)

if (is.null(opt$rds)){
  informat = opt$informat
  if (informat == "h5seurat" ){
      seurat_obj = SeuratDisk::LoadH5Seurat(opt$input, assays = "RNA",verbose = F )
  }else if (informat == "rds") {
      seurat_obj <- readRDS(opt$input)
  }else{
      stop("Inputfile should in rds or h5seurat format.")
  }

  #需要设置感兴趣的groups
  if ( is.null(opt$ident2use ) ){
      print("NO cell identity column name AVAILABLE! All cells in group.by column will be used as default.")
  }else{
      ident2use = unlist(strsplit( opt$ident2use,",",perl = T))
      Seurat::Idents(object = seurat_obj)<- opt$groupby 
      seurat_obj= subset(seurat_obj, idents = ident2use )
      if( is.factor(Seurat::Idents(object = seurat_obj))){
        seurat_obj@meta.data[[opt$groupby ]] = droplevels(seurat_obj@meta.data[[opt$groupby ]])
      }
      group_name = unique(seurat_obj@meta.data[[opt$groupby ]])
      if (all(ident2use %in% group_name)){
        ident <- ident2use
      } else {
        not_in_group <- ident2use[!ident2use %in% group_name]
        stop(paste("The following group of interest are not in group.by column:", not_in_group))
      }
  }

  if ( !opt$reduct1 %in% names(seurat_obj@reductions) ) stop(paste0("Primary reduction ",opt$reduct1," is not present in the imput file." ))
  if ( !opt$reduct2 %in% names(seurat_obj@reductions) ) stop(paste0("Secondary reduction ",opt$reduct2," is not present in the imput file." ))


  #初始化配置
  start_time <- Sys.time() # 记录初始时间
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "hdWGCNA" # the name of the hdWGCNA experiment
  )
  end_time <- Sys.time() # 记录终止时间
  print(paste0("初始化配置完成，共耗时: ",end_time - start_time))
  #Time difference of 2.125303 mins

  #根据细胞数设置k值的大小, k should between 25 and 75
  if(dim(seurat_obj)[2]<=50000){
      k = 25
  }else if (dim(seurat_obj)[2]>50000 && dim(seurat_obj)[2]<=80000) {
      k = 55
  }else if (dim(seurat_obj)[2]>120000) {
      k = 75
  }


  #metacells构建
  # construct metacells  in each group
  #ident.group设置使用到的分群seurat_obj@misc$hdWGCNA$wgcna_metacell_obj@active.ident
  #Default slot is slot='counts'
  start_time <- Sys.time() # 记录初始时间
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("sampleid",opt$groupby), # specify the columns in seurat_obj@meta.data to group by
    reduction = opt$reduct1, # select the dimensionality reduction to perform KNN on
    k = k, # nearest-neighbors parameter
    max_shared = 10, # maximum number of shared cells between two metacells
    ident.group = opt$ident # set the Idents of the metacell seurat object
  )

  # normalize metacell expression matrix:
  seurat_obj <- NormalizeMetacells(seurat_obj)

  end_time <- Sys.time() # 记录终止时间
  print(paste0("Meta细胞构建完成,共耗时: ",end_time - start_time,sep = ""))

  #共表达网络分析
  group_name = unique(seurat_obj@misc$hdWGCNA$wgcna_metacell_obj[[opt$groupby]][[opt$groupby]])
  
  use_metacells = as.logical(opt$MetaCells) 
  start_time <- Sys.time() # 记录初始时间
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = group_name, # the name of the group of interest in the group.by column
    group.by= opt$groupby , # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'RNA', # using RNA assay
    slot = 'data',
    use_metacells = use_metacells
  )
  end_time <- Sys.time() # 记录终止时间
  print(paste0("SetDatExpr运行完成, 共耗时: ",end_time - start_time,sep = ""))

  #选择软阈值
  use_metacells = as.logical(opt$MetaCells) 
  start_time <- Sys.time() # 记录初始时间
  # Test different soft powers:
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = opt$networkType,
    use_metacells = use_metacells
  )

  # plot the results:
  # assemble with patchwork
  pdf(paste0(output_dir,"/PlotSoftPowers.pdf",sep = ""))
  plot_list <- PlotSoftPowers(seurat_obj)
  p = wrap_plots(plot_list, ncol=2)
  plot(p)
  dev.off()

  end_time <- Sys.time() # 记录终止时间
  print(paste0("软阈值选择完成，共耗时：",end_time - start_time," ,结果见PlotSoftPowers.pdf文件"))
  #[1] "软阈值选择"
  #Time difference of 42.52228 secs
  #save.image(paste0(output_dir,"/3.SoftPowers.RData",sep = ""))


  start_time <- Sys.time() # 记录初始时间
  # construct co-expression network:
  seurat_obj <- ConstructNetwork(
    seurat_obj, soft_power=opt$soft_power,
    setDatExpr=FALSE,
    tom_outdir = output_dir,
    wgcna_name = "hdWGCNA",
    networkType = opt$networkType
  )

  pdf(paste0(output_dir,"/PlotDendrogram.pdf",sep = ""))
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
  dev.off()



  #4. 模块特征基因与连通性
  #4.1 计算模块特征基因
  start_time <- Sys.time() # 记录初始时间
  # need to run ScaleData first or else harmony throws an error:
  seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

  # compute all MEs in the full single-cell dataset,默认进行harmony去批次
  seurat_obj <- ModuleEigengenes(seurat_obj,group.by.vars="sampleid")

  # harmonized module eigengenes:
  hMEs <- GetMEs(seurat_obj)
  hMEs$rawbc = rownames(hMEs)
  hMEs = hMEs%>% select(rawbc, everything())
  write.table(hMEs, quote = F,sep ="\t",row.names = F,paste0(output_dir,"/hMEs.xls",sep = ""))

  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)
  MEs$rawbc = rownames(hMEs)
  MEs = MEs%>% select(rawbc, everything())
  write.table(MEs, quote = F,sep ="\t",row.names = F,paste0(output_dir,"/MEs.xls",sep = ""))

  end_time <- Sys.time() # 记录终止时间
  print(paste0("模块特征基因，共耗时：",end_time - start_time))

  #计算连通性
  #需要加载qlcMatrix R包
  start_time <- Sys.time() # 记录初始时间
  # compute eigengene-based connectivity (kME):
  #使用的 group.by参数与 SetDatExpr 函数中使用的 group.by 和 group_name 参数要保持一致，以保证在计算 kME 时仅考虑相关的细胞群。
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = opt$groupby, group_name = group_name
  )

  # plot genes ranked by kME for each module
  modules = GetModules(seurat_obj) %>% subset(module != 'grey')
  mods <- levels(modules$module); mods <- mods[mods != 'grey']

  pdf(paste0(output_dir,"/PlotKMEs.pdf",sep = ""),width = 3*length(mods),height =  0.7*length(mods))
  p = PlotKMEs(seurat_obj, ncol=5,wgcna_name = "hdWGCNA")
  plot(p)
  dev.off()

  end_time <- Sys.time() # 记录终止时间
  print(paste0("计算连通性, 共耗时：",end_time - start_time))

  #获取基因的模块分布表
  modules <- GetModules(seurat_obj,wgcna_name = "hdWGCNA")
  write.table(modules, quote = F,sep ="\t",row.names = F,paste0(output_dir,"/gene_modules_allocation.xls",sep = ""))

  #Extract the top N hub genes for a given set of modules.
  modules_topn <- GetHubGenes(seurat_obj, n_hubs = opt$topn_hub, wgcna_name = "hdWGCNA")
  write.table(modules_topn, quote = F,sep ="\t",row.names = F,paste0(output_dir,"/modules_topn_hubgenes.xls",sep = ""))

  #计算hub-gene打分
  # compute gene scoring for the top 25 hub genes by kME for each module
  # with Seurat method
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = opt$topn_hub ,
    method= opt$Gene_scoring
  )

  saveRDS(seurat_obj, file=paste0(output_dir,"/hdWGCNA_object.rds",sep = ""))
} else  seurat_obj = readRDS(opt$rds)

###############################################################
###############################################################
###########################下游可视化###########################

if ( file.exists(paste0(output_dir,"/Plot",sep = "")) ){
     allfile = dir(paste0(output_dir,"/Plot",sep = ""))
     target <- grep("*", allfile)
     file.remove(paste0(output_dir,"/Plot/",allfile[target],sep = ""))
}else{
    dir.create(paste0(output_dir,"/Plot",sep = ""),recursive = T)
}

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
seurat_obj,
wgcna_name = "hdWGCNA",
features='hMEs', # plot the hMEs
order=TRUE, # order so the points with highest hMEs are on top
reduction = opt$reduct2
)


pdf(paste0(output_dir,"/Plot/1.ModuleFeaturePlot.pdf",sep = ""))
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)
dev.off()   

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  wgcna_name = "hdWGCNA",
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = ifelse(opt$Gene_scoring == "UCell" ,TRUE, FALSE),# depending on Seurat vs UCell for gene scoring
  reduction= opt$reduct2
)

pdf(paste0(output_dir,"/Plot/2.ModuleFeaturePlot_score.pdf",sep = ""))
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)
dev.off()

pdf(paste0(output_dir,"/Plot/3.Correlation.pdf",sep = ""))
# plot module correlagram
colfunc <- grDevices::colorRampPalette(c('#004A97',  'white', '#005E00'))
col = colfunc(200)
ModuleCorrelogram(seurat_obj,col = col)
dev.off()

# get hMEs from seurat object 
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = opt$groupby)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
coord_flip() +
RotatedAxis() +
scale_color_gradient2(high='#005E00', mid='white', low='#004A97')

pdf(paste0(output_dir,"/Plot/4.Dotplot.pdf",sep = ""))
# plot output
p
dev.off()

####网络分析
#由于修改过TOM文件的默认位置，需要修改seurat_obj里对应的文件路径
if (is.null(opt$rds)){
seurat_obj@misc[["hdWGCNA"]]$wgcna_net$TOMFiles =  paste0(output_dir,"/hdWGCNA_TOM.rda",sep = "")
}else{ 
  if (!file.exists(opt$TOMFiles)) stop(paste0("TOM file ",opt$TOMFiles," not found. Please update path to TOM file." ))
  seurat_obj@misc[["hdWGCNA"]]$wgcna_net$TOMFiles = opt$TOMFiles 
}

ModuleNetworkPlot(seurat_obj,outdir = paste0(output_dir,"/Plot/ModuleNetworks"),wgcna_name = "hdWGCNA" )

# hubgene network
pdf(paste0(output_dir,"/Plot/5.HubGeneNetworkPlot.pdf"))
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

#使用Umap同时可视化共表达中的所有基因
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 , wgcna_name = "hdWGCNA" # min distance between points in UMAP space
)

#可视化
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
pdf(paste0(output_dir,"/Plot/6.RunModuleUMAP.pdf"))
p = ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
  color=umap_df$color, # color each point by WGCNA module
  size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +umap_theme()
print(p)
dev.off()

pdf(paste0(output_dir,"/Plot/7.ModuleUMAPPlot.pdf"))
ModuleUMAPPlot(
seurat_obj,
edge.alpha=0.25,
sample_edges=TRUE,
edge_prop=0.1, # proportion of edges to sample (20% here)
label_hubs=2 ,# how many hub genes to plot per module?
keep_grey_edges=FALSE
)
dev.off()

#性状相关性
#Warnings: 在进行相关性分析前，确保traits为Factor或者numeric格式
traits = unlist(strsplit(opt$traits, ",", perl =T))

if ("group" %in% traits){
  traits[which(traits =="group")] = "Group"
  colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) =="group")]  ="Group"
}

seurat_obj <- ModuleTraitCorrelation(
seurat_obj,
traits = traits,
group.by= opt$groupby
)
# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(seurat_obj)
print("names(mt_cor)")
names(mt_cor)

print("names(mt_cor$cor)")
names(mt_cor$cor)

pdf(paste0(output_dir,"/Plot/8.PlotModuleTraitCorrelation.pdf"))
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 4,
  text_digits = 2,
  text_color = 'white',
  high_color = '#FDE725',
  mid_color = '#20908C',
  low_color = '#440154',
  plot_max = 0.2,
  combine=TRUE
)
dev.off()

#一对多的DME分析
DMEs_all <- FindAllDMEs(
  seurat_obj,
  group.by = opt$groupby,
  wgcna_name = 'hdWGCNA'
)

#由于DMEs_all的分组列默认为“group”,为防止与seurat中的group列混淆，出表格时改为groupby作为列名
DMEs_all_df = DMEs_all
colnames(DMEs_all_df)[ which(colnames(DMEs_all_df) =="group")] = opt$groupby
write.table(DMEs_all_df , quote = F,sep ="\t",row.names = F, paste0(output_dir,"/Plot/9.DMEs_all.xls",sep = ""))

pdf(paste0(output_dir,"/Plot/9.FindAllDMEs.pdf"))
p <- PlotDMEsVolcano(
  seurat_obj,
  DMEs_all,
  wgcna_name = 'hdWGCNA',
  plot_labels=TRUE,
  show_cutoff=FALSE
)+ theme_bw() + NoLegend() + facet_wrap(~group, ncol=3)
plot(p)
dev.off()

file.copy("/public/scRNA_works/works/liuxuan/Test/hdWGCNA_Script/hdWGCNA结果说明.doc",output_dir)

setwd(paste0(output_dir,"/Plot/"))
system("for i in  `ls *.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")

setwd(paste0(output_dir,"/Plot/ModuleNetworks"))
system("for i in  `ls *.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")

setwd(output_dir)
system("for i in  `ls *.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")

if ( "Customize" %in% args ){
  if ( file.exists(paste0(output_dir,"/Customize",sep = "")) ){
    allfile = dir(paste0(output_dir,"/Customize",sep = ""))
    target <- grep("*", allfile)
    file.remove(paste0(output_dir,"/Customize/",allfile[target],sep = ""))
  }else{
    dir.create(paste0(output_dir,"/Customize",sep = ""),recursive = T)
  }

  if ( is.null(opt$module) ){
    print("No module specifed, use all modules for ploting")
    modules = GetModules(seurat_obj) %>% subset(module != 'grey')
    mods <- levels(modules$module)
    module = levels(modules$module)[levels(modules$module) !="grey"]
    all_module = TRUE
  
    pdf( paste0(output_dir,"/Customize/1.HubGeneNetworkPlot.pdf",sep = ""))
    HubGeneNetworkPlot(
      seurat_obj,
      n_hubs = opt$n_hubs, n_other=10,
      edge_prop = 0.75,
      mods = "all"
    )
    dev.off()
  }else{
  modules = GetModules(seurat_obj) %>% subset(module != 'grey')
  mods <- levels(modules$module)
  all_module = FALSE

  module = unlist(strsplit(opt$module, ",", perl =T))
  
  if (all(module %in% mods)){
    print(paste0("将对 ", module," 模块进行个性化绘图。"))
  }else{
    stop(paste0(setdiff(module, mods)," 模块不在module list中，请核对。"))
  }
  # hubgene network
  pdf( paste0(output_dir,"/Customize/1.HubGeneNetworkPlot.pdf",sep = ""))
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = opt$n_hubs, n_other=10,
    edge_prop = 0.75,
    mods = module
  )
  dev.off()
  }

  if ( is.null(opt$contrast ) ){ 
    print("未提供比较组和实验组信息，无法执行两组比较的DME分析")
  }else{
    contrast = opt$contrast
    contrasts = unlist( strsplit(contrast,":",perl = T) )
    all_levels = as.vector(unique(seurat_obj@meta.data[,contrasts[1]]))

    if ( contrasts[2] == "all" & contrasts[3] != "all" ){
      all_levels = all_levels[-which(all_levels==contrasts[3])]
      all_comparisions = paste(contrasts[1],all_levels,contrasts[3],sep = ":")
    }else if( contrasts[2] == "all" & contrasts[3] == "all" ){
      all_comparisions = lapply(all_levels,
                      function(x) paste(contrasts[1],x,paste0(all_levels[-which(all_levels==x)]),sep = ":"))
      all_comparisions = unlist(all_comparisions)
      
      #去重
      all_comparisions <- unique(lapply(all_comparisions, function(x) {
        y <- unlist(strsplit(x, split = ":"))
        paste(sort(y), collapse = ":")
      }))
      all_comparisions = unlist(all_comparisions)
      remove_last_semicolon <- function(x) {
        gsub(":([^:]+)$", "", x)
      }
      all_comparisions <- lapply(all_comparisions, remove_last_semicolon)
      all_comparisions  = paste(opt$groupby,unlist(all_comparisions),sep = ":")
    }else if ( contrasts[2] != "all" & contrasts[3] == "all" ){
      ref_levels = paste0(all_levels[-which(all_levels==contrasts[2])])
      all_comparisions = paste(contrasts[1],contrasts[2],ref_levels,sep = ":")
    }else{
      all_comparisions = contrast
    }

    for (i in all_comparisions){
    contrasts = unlist( strsplit(i,":",perl = T) )
    predicate_group1 = paste0(contrasts[1]," == '", contrasts[2],"'")
    group1 <- seurat_obj@meta.data %>% subset(eval(parse(text = predicate_group1))) %>% rownames
    predicate_group2 = paste0(contrasts[1]," == '", contrasts[3],"'")
    group2 <- seurat_obj@meta.data %>% subset(eval(parse(text = predicate_group2))) %>% rownames

    DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    wgcna_name='hdWGCNA'
    )
    compare_file_name = paste(unlist( strsplit(i,":",perl = T) ),collapse = "_")
    write.table(DMEs, quote = F,sep ="\t",row.names = F,paste0(output_dir,"/Customize/2.DME_",compare_file_name,".xls",sep = ""))

    pdf(paste0(output_dir,"/Customize/2.DME_",compare_file_name,"_Volcano.pdf",sep = ""))
    p = PlotDMEsVolcano(
      seurat_obj,
      DMEs,
      wgcna_name='hdWGCNA'
    )+theme_bw()+ NoLegend()
    plot(p)
    dev.off()
    }
  }
setwd(paste0(output_dir,"/Customize/"))
print("Convert pdf to png...")
system("for i in  `ls *.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")
}