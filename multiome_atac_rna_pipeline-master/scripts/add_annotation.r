suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))


option_list = list(
    make_option( c("--input", "-i" ), type = "character",
                 help = "The input exprssion matrix in RDS format."),
    make_option( c("--species", "-e" ), type = "character", 
	             help = "e.g. Human, Mouse or Other."),
    make_option( c("--anno","-a" ), type = "character", default = NULL,
                 help = "[OPTIONAL]Annotation file submitted when species is Other.")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if ( is.null(opt$species) ){
    stop("Please provide the species information with parameter \"-e\" !")
}else{
    species = opt$species
}

marker = read.delim(opt$input,sep="\t",stringsAsFactors=F)

if(species == "Human" ){
	anno = read.delim("/data/database/cellranger-refdata/refdata-cellranger-GRCh38-3.0.0/annotation/gene_annotation.xls",stringsAsFactors=F,sep='\t',quote = "")
	if( colnames(marker[1]) == "gene" ){
		res = left_join(marker, anno, by = c("gene"="id") )
	}
	else{
		res = left_join(marker, anno, by = c("GeneID"="id") )
	}
	res[is.na(res)] = "--" 
	write.table(res, paste0(strsplit(opt$input, "\\.xls")[[1]][1], "_anno.xls", collapse=""),sep='\t',quote=F,row.names=F)
} else if(species == "Mouse" ){
	anno = read.delim("/data/database/cellranger-refdata/refdata-cellranger-mm10-3.0.0/annotation/gene_annotation.xls",stringsAsFactors=F,sep='\t',quote = "")
	if( colnames(marker[1]) == "gene" ){
		res = left_join(marker, anno, by = c("gene"="id") )
	}
	else{
		res = left_join(marker, anno, by = c("GeneID"="id") )
	}
	res[is.na(res)] = "--" 
	write.table(res, paste0(strsplit(opt$input, "\\.xls")[[1]][1], "_anno.xls", collapse=""),sep='\t',quote=F,row.names=F)
} else {
    if ( is.null(opt$anno) ){
    stop("choosing species other than 'Human' and 'Mouse', please provide the annotation file with parameter \"-a\" !")
    }
    else{
        anno = read.delim(opt$anno,stringsAsFactors=F,sep='\t',quote = "")
        if( colnames(marker[1]) == "gene" ){
            res = left_join(marker, anno, by = c("gene"="id") )
        }
        else{
            res = left_join(marker, anno, by = c("GeneID"="id") )
        }
        res[is.na(res)] = "--"
        write.table(res, paste0(strsplit(opt$input, "\\.xls")[[1]][1], "_anno.xls", collapse=""),sep='\t',quote=F,row.names=F)
    }
}

file.remove(opt$input)

