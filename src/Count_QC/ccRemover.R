#!/usr/bin/env Rscript

#############################################
#Author: guangye gong 
#Creat Time: 2018-7-25 14:45:00
#############################################


library("optparse") 
option_list = list(
	make_option(c("-e", "--expression"), type="character", default=NULL, 
          help="expression matrix file name", metavar="character"),
	make_option(c("-s", "--species"), type="character", default=NULL,
	  help="The species which the gene names are from. Either human or mouse.", metavar="character"),
	make_option(c("-t", "--nametype"), type="character", default=NULL, 
          help="The type of gene name considered either, Ensembl gene IDS (ensembl), offical gene symbols (symbol), Entrez gene IDS (entrez), or Unigene IDS (unigene).", metavar="character"),
	make_option(c("-c", "--cutoff"),type="double", default=3,
          help="The significance cutoff for identifying sources of variation related to the cellcycle.The default value is 3, which roughly corresponds to a p-value of 0.01.", metavar="double"),
	make_option(c("-m", "--maxit"),type="double", default=4,
          help="The maximum number of iterations for the algorithm. The default value is 4.", metavar="double"),
	make_option(c("-b", "--nboot"),type="double", default=200,
          help="The number of bootstrap repititions to be carried out on each iteration to determine the significance of each component. The default value is 200.", metavar="double"),
	make_option(c("-n", "--ntop"),type="double", default=10,
          help="The number of components considered tested at each iteration as cell-cycle effects. The default value is 10.", metavar="double")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$expression) | is.null(opt$species) | is.null(opt$nametype)){
  print_help(opt_parser)
  stop("expression matrix file ,species, type of gene name must be supplied", call.=FALSE)
}

library(ccRemover)
expression<-read.delim(opt$expression,stringsAsFactors = F,header = T ,row.names = 1)
mean_gene_exp <- rowMeans(expression)
# Center data and select small sample for example
expression_data_cen <- t(scale(t(expression), center=TRUE, scale=FALSE))
# Extract gene names
gene_names <- rownames(expression_data_cen)
# Determine which genes are annotated to the cell-cycle
cell_cycle_gene_indices <- gene_indexer(gene_names,species = opt$species, name_type = opt$nametype)
# Create "if_cc" vector
if_cc <- rep(FALSE,nrow(expression_data_cen))
if_cc[cell_cycle_gene_indices] <- TRUE
# Move data into list
dat <- list(x=expression_data_cen, if_cc=if_cc)
# Run ccRemover
xhat <- ccRemover(dat, cutoff = opt$cutoff, max_it = opt$maxit, nboot = opt$nboot, ntop = opt$ntop)
xhat <- xhat + mean_gene_exp
xhat= cbind(rownames(xhat),xhat)
colnames(xhat)[1]=opt$nametype
out<-unlist(strsplit(opt$expression,split='/'))
out<-out[length(out)]
write.table(xhat,paste("ccRemover_",out,sep=""),quote=F,row.names=F,sep="\t",na="")
