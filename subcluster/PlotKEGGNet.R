#!/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/scanpy/bin/Rscript
#PlotKEGGNet.R
#version: 0.1
#pipeline: single cell VDJ, single cell v3
#author: CaoWei
#update: 2024/09/13

suppressMessages(library("argparse"))
suppressMessages(library(pathview))


kegg.dir<-'/PERSONALBIO/work/singlecell/s00/software/4.Annotation/KEGG_fig/'

parser <- ArgumentParser(description='Process names and file path')

parser$add_argument('-m','--allmarkers',required=T,type='character',
                    help='Input all markers file')
parser$add_argument('-o','--output',required=T,type='character',
                    help='Figure and table output directory')
parser$add_argument('-s','--species',required=T,type='character',
                    help='species ')
parser$add_argument('-n','--topn',type='character',
                    help='The top number of KEGG pathways')
parser$add_argument('-t','--method',required=T,type='character',
                    help='If it is marker analysis, this option is marker. If it is differential analysis, this option is A-vs-B ')
args = parser$parse_args()
topnum = args$topn
outdir <- normalizePath(args$output)
if(args$method=='marker'){
	markers<-read.table(args$allmarkers,sep="\t",header=T,row.names = 1)
}else{
	markers<-read.table(args$allmarkers,sep="\t",header=T)
}


species <- args$species
species.org=switch(species,
				"hsa"={suppressMessages(library(org.Hs.eg.db));org.Hs.eg.db}, #human
				"mmu"={suppressMessages(library(org.Mm.eg.db));org.Mm.eg.db}, #mouse
				"rno"={suppressMessages(library(org.Rn.eg.db));org.Rn.eg.db}, #rat
				"dme"={suppressMessages(library(org.Dm.eg.db));org.Dm.eg.db}, #fruit_fly
				"dre"={suppressMessages(library(org.Dr.eg.db));org.Dr.eg.db}, #zebrefish
				"ath"={suppressMessages(library(org.At.tair.db));org.At.tair.db}, #Arabidopsis
				"sce"={suppressMessages(library(org.Sc.sgd.db));org.Sc.sgd.db}, #yeast
				"cel"={suppressMessages(library(org.Ce.eg.db));org.Ce.eg.db}, #C.elegans
				"bta"={suppressMessages(library(org.Bt.eg.db));org.Bt.eg.db}, #Bovine
				"mcc"={suppressMessages(library(org.Mmu.eg.db));org.Mmu.eg.db}, #monkey
				"cfa"={suppressMessages(library(org.Cf.eg.db));org.Cf.eg.db}, #dog
				"ssc"={suppressMessages(library(org.Ss.eg.db));org.Ss.eg.db}, #pig
				"gga"={suppressMessages(library(org.Gg.eg.db));org.Gg.eg.db}, #chicken
				"xla"={suppressMessages(library(org.Xl.eg.db));org.Xl.eg.db}, #frog
				"ptr"={suppressMessages(library(org.Pt.eg.db));org.Pt.eg.db}, #chimpanzee
				"aga"={suppressMessages(library(org.Ag.eg.db));org.Ag.eg.db}, #mosquito
				#org.EcK12.eg.db nocode Ecoli_strain_K12
				#org.EcSakai.eg.db nocode Ecoli_strain_Sakai
				#org.Pf.plasmo.db pfa/pfd/pfh Malaria
				#org.Mxanthus.db nocode Myxococcus xanthus
                stop(paste('Species',species,"is not support now")))

if(!all(c('SYMBOL','ENTREZID') %in% keytypes(species.org))){stop('gene symbol and entrezid not in database')}

clust_count<-table(markers$cluster)
rawdir=getwd()
for(clust_num in names(clust_count)){
  if(clust_count[[clust_num]]<=1){next}
  print(paste("Drawing cluster",clust_num))
  cluster_markers = subset(markers,cluster==clust_num & p_val < 0.05)
  rownames(cluster_markers) <- cluster_markers$gene
  matchid_data<-select(species.org,columns='ENTREZID',keytype='SYMBOL',keys=cluster_markers$gene)
  remove_row <- which(is.na(matchid_data$ENTREZID)|duplicated(matchid_data$ENTREZID))
  draw_df <- setNames(cluster_markers[-remove_row,'avg_log2FC'], matchid_data[-remove_row,'ENTREZID'])
  
  #Due to stupid pathview function cannot specify output folder and output file name, can only generate figure in current dir
	if(args$method=="marker"){
		pro_wk = file.path(normalizePath(outdir),"Each_celltype_marker",paste0('cluster_',clust_num),"enrichment/KEGG_NET")
	}else{
		pro_wk = file.path(normalizePath(outdir),args$method,paste0('cluster',clust_num),"enrichment/KEGG_NET")
	}
	dir.create(pro_wk)
  setwd(pro_wk)
  df <- read.table('../all/KEGG_enrichment.xls',header=T,sep='\t')
	if (toupper(topnum) != 'ALL'){
		df<-df[order(df$Pvalue),]
  	df<-df[1:min(nrow(df),as.numeric(topnum)),]
	}


  
  for(pathway_id in df$PathwayID){
    id=substr(pathway_id,4,8)
    draw_df<-draw_df/max(abs(draw_df))
    pathview(gene.data = draw_df, pathway.id = id, species = species,out.suffix = "KEGG_net", kegg.native = T,kegg.dir=paste0(kegg.dir,species),same.layer = F)   
  }
	setwd(rawdir)
}







