######################################
# 自定义表输出函数
######################################
write_table <- function(data, outtablename){
    data <- data
    outtablename <- outtablename
    write.table(x=data, file=outtablename, quote = F, col.names = T, row.names = F, sep=",")
}

######################################
# 删除字符串前后的空格，将字符串中间的空格转换为下划线
######################################

black2underline <- function(x){
    x <- x
    
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub(" ", "_", x)
    
    return(x)
}

######################################
# 自定义绘制表达量密度函数
######################################
plot_density <- function(data, control, case, gene_or_SJ){
  data <- data
  data$tag <- ""
  data[data$cell.group == "cell.group.g1",]$tag <- control
  data[data$cell.group == "cell.group.g2",]$tag <- case
  
  x <- data$pct.cells.expr
  z <- data$tag
  maintitle <- ""
  ytitle <- "Density"
  xtitle <- paste0(gene_or_SJ, " Expressed (% Cells)")
  xmin <- 0 ; xmax <- max(x) ; xinterval <- 10
  legendtitle <- ""
  
  d <- density(x)
  max_height_location <- round(d$x[which.max(d$y)])
  
  # Plot
  plot <- ggplot() +
    geom_density(data, mapping=aes(x=x, color=z), alpha=0.5) +
    scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
    labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
    geom_vline(xintercept = max_height_location,  colour="#990000", linetype="dashed") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_blank(),
          plot.title=element_text(hjust = 0.5, size=15),
          plot.subtitle=element_text(hjust = 0.5, size=15),
          axis.line.y.left = element_line(color="black"),
          axis.line.x = element_line(color="black"),
          axis.title=element_text(size=12),
          axis.text=element_text(size=12),
          axis.text.x=element_text(size=8, colour="black"),
          axis.text.y=element_text(size=8, colour="black"),
          legend.title=element_text(size=8),
          legend.text=element_text(size=8)
    )
  out <- list(number=max_height_location, figure=plot, data=data)

  return(out)
  
}

######################################
# 判断降维类型;  数据准备时和绘制散点图时需要判断降维类型
######################################

which_dimension_type <- function(dimension_type){
    dimension_type <- dimension_type

    if(dimension_type %in% c("umap", "tsne")){
    NULL
    }else{
    message("您没有正确设置降维方法，目前只支持两种降维方法，分别为 umap 和 tsne")
    quit(save = "no", status = 1, runLast = TRUE)
    }
}

######################################
# 将变量以及变量类型打印出来，方便检查
######################################

check_paramaters <- function(par) {
  message(paste("变量", deparse(substitute(par)), "的类型是", class(par), "，值是", par))
}

######################################
# 创建文件夹
######################################

create_folder <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

######################################
# 提取位置匹配外显子的 junction， 返回 junction 列表
######################################

extract_juntions_by_annotate <- function(gtf, junctions, genename_or_geneid){

    gtf <- gtf 
    junctions <- junctions
    genename_or_geneid <- genename_or_geneid

    #######################
    # message("根据gtf文件注释junction对应的基因")
    #######################

    ## 解析gtf文件
    gtf <- gtf[which(gtf$V3=="exon"), ]
    if (genename_or_geneid == "gene_name"){
        # Retrieve gene names
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("gene_name", x, value=TRUE))
        . <- gsub("gene_name", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)

        gtf$gene_short_name <- .   
    }

    if (genename_or_geneid == "gene_id"){
        # Retrieve gene names        
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("gene_id", x, value=TRUE))
        . <- gsub("gene_id", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)

        gtf$gene_short_name <- .
    }

    gtf$V4 <- gtf$V4 - 1
    gtf$V5 <- gtf$V5 + 1

    # logical <- any(grepl("chr", gtf$V1, fixed=TRUE))

    # Collapse start position
    # Keep unique entries
    gtf.small <- gtf[, c("V1", "V5", "gene_short_name")]

    # if(logical==FALSE) {
    #     gtf.small$chr.pos <- paste("chr", gtf$V1, ":", gtf$V5, sep="")
    # } else if(logical==TRUE) {
        gtf.small$chr.pos <- paste(gtf$V1, ":", gtf$V5, sep="")
    # }
        
    gtf.small <- gtf.small[,c("chr.pos", "gene_short_name")]
    gtf.small <- unique(gtf.small)

    # Collapse
    gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
    . <- by(gtf.small[,"gene_short_name"], gtf.small[,"chr.pos"], function(x) {paste(x, collapse="|")})
    gtf.small.collapsed <- data.frame("chr"=as.character(names(.)), "gene_short_name"=as.character(.), stringsAsFactors=FALSE)

    # Save as new object
    gtf.small.collapsed.start <- gtf.small.collapsed
    names(gtf.small.collapsed.start) <- paste(names(gtf.small.collapsed.start), ".start", sep="")

    # Collapse end position
    # Keep unique entries
    gtf.small <- gtf[, c("V1", "V4", "gene_short_name")]

    # if(logical==FALSE) {
    #     gtf.small$chr.pos <- paste("chr", gtf$V1, ":", gtf$V4, sep="")
    # } else if(logical==TRUE) {
        gtf.small$chr.pos <- paste(gtf$V1, ":", gtf$V4, sep="")
    # }

    gtf.small <- gtf.small[,c("chr.pos", "gene_short_name")]
    gtf.small <- unique(gtf.small)

    # Collapse
    gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
    . <- by(gtf.small[,"gene_short_name"], gtf.small[,"chr.pos"], function(x) {paste(x, collapse="|")})
    gtf.small.collapsed <- data.frame("chr"=as.character(names(.)), "gene_short_name"=as.character(.), stringsAsFactors=FALSE)

    # Save as new object
    gtf.small.collapsed.end <- gtf.small.collapsed
    names(gtf.small.collapsed.end) <- paste(names(gtf.small.collapsed.end), ".end", sep="")

    ################

    df <- as.data.frame(junctions)
    . <- strsplit(junctions, split=":", fixed=TRUE)
    df$chr <- sapply(., function(x) {as.character(x[1])})
    df$start <- sapply(., function(x) {as.character(x[2])})
    df$end <- sapply(., function(x) {as.character(x[3])})
    df$chr.start <- paste(df$chr, ":", df$start, sep="")
    df$chr.end <- paste(df$chr, ":", df$end, sep="")
    df <- join(df, gtf.small.collapsed.start, by="chr.start", type="left")
    df <- join(df, gtf.small.collapsed.end, by="chr.end", type="left")
    df <- na.omit(df)

    keep_junctions <- df$junctions

    return(keep_junctions)

}

######################################
# 按照差异类型提取细胞
######################################

subset_cell_according_diff_type <- function(df, diff_type, g1, g2, ...){
  
  df <- df    # sample.metadata
  diff_type <- diff_type
  g1 <- g1
  g2 <- g2
  
  result_list <- list()
  
  if (!diff_type %in% colnames(df)) {
    message("您没有正确设置差异分组的类型，请检查")
    message(paste0("差异分组类型可以从以下类型中选择：", colnames(df)))
    quit(save = "no", status = 1, runLast = TRUE)
  }
  
  
  df <- df[df[,diff_type] %in% c(g1, g2), ]
  df$group_for_plot <- df[,diff_type]
  cells.g1 <- df[df[,diff_type] == g1, "cell.id"]
  cells.g2 <- df[df[,diff_type] == g2, "cell.id"]
  
  if (diff_type == "clusters"){
    names(g1) <- paste0("Cluster", g1)
    names(g2) <- paste0("Cluster", g2)
    df$group_for_plot <- paste0("Cluster", df$group_for_plot)
  } else {
    names(g1) <- g1
    names(g2) <- g2
  }
  
  
  result_list[["data"]] <- df
  result_list[["cells.g1"]] <- cells.g1
  result_list[["cells.g2"]] <- cells.g2
  result_list[["g1"]] <- g1
  result_list[["g2"]] <- g2
  
  # 默认情况下返回数据集和对应的细胞ID，若需要其他类型信息，在函数的后面添加相应的名称即可按照比较组返回各自的信息
  if (length(list(...)) > 0){
    for(p in list(...)){
      
      if (! p %in% colnames(df)) {
        message(paste0("您没有正确设置标签 ", p ,"，请检查"))
        quit(save = "no", status = 1, runLast = TRUE)
      }
      
      result_list[[paste0(p, '.g1')]] <- df[df[,diff_type] == g1, p] 
      result_list[[paste0(p, '.g2')]] <- df[df[,diff_type] == g2, p]

    }
  }
  
  return(result_list)
}

######################################
# 解析GTF文件
######################################

parsing_gtf <- function(gtf, gene.metadata, gene_or_exon){
    gtf <- gtf
    gene_or_exon <- gene_or_exon
    gene.metadata <- gene.metadata
    #############

    # Subset gene entries
    if (!is.null(gene_or_exon)){
        gtf <- gtf[which(gtf$V3==gene_or_exon), ]
    }
    

    # gene_id
    . <- strsplit(gtf$V9, split=";")
    . <- sapply(., function(x) grep("gene_id", x, value=TRUE))
    . <- gsub("gene_id", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)
    # cellranger在处理gtf文件时去掉了gene_id的版本号，因此所有10X的项目都要把版本号去掉
    # M20项目不对gene_id做处理
    if(!all(grepl(pattern = "\\.", x = gene.metadata$gene_id))){
        . <- str_split_fixed(., pattern = "\\.", n=2)[,1]
    }
    
    gtf$gene_id <- .

    # gene_type
    . <- strsplit(gtf$V9, split=";")
    . <- sapply(., function(x) grep("gene_biotype", x, value=TRUE))
    . <- gsub("gene_biotype", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)
        
    if(.[1] == "character(0)") {
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("gene_type", x, value=TRUE))
        . <- gsub("gene_type", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)
    }
    gtf$gene_type <- .

    if (!is.null(gene_or_exon)){
        return(gtf)
    }

    # exon_id
    . <- strsplit(gtf$V9, split=";")
    . <- sapply(., function(x) grep("exon_id", x, value=TRUE))
    . <- gsub("exon_id", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)
  
    gtf$exon_id <- .
  
    # transcript_id
    . <- strsplit(gtf$V9, split=";")
    . <- sapply(., function(x) grep("transcript_id", x, value=TRUE))
    . <- gsub("transcript_id", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)
    gtf$transcript_id <- .

    # transcript_biotype
    . <- strsplit(gtf$V9, split=";")
    . <- sapply(., function(x) grep("transcript_biotype", x, value=TRUE))
    . <- gsub("transcript_biotype", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)
  
    if(length(unique(.))==1 & unique(.)[1] == "character(0)") {
    
        . <- strsplit(gtf$V9, split=";")
        . <- sapply(., function(x) grep("transcript_type", x, value=TRUE))
        . <- gsub("transcript_type", "", .)
        . <- gsub(" ", "", .)
        . <- gsub("\"", "", .)
    
        gtf$transcript_biotype <- .
    } else {
        gtf$transcript_biotype <- .
    }

    return(gtf)

}

######################################
# 创建MARVEL对象
######################################

CreateMarvelObject.10x.ZCY <- function(gene.norm.matrix=NULL,
                                   gene.norm.pheno=NULL,
                                   gene.norm.feature=NULL,
                                   gene.count.matrix=NULL,
                                   gene.count.pheno=NULL,
                                   gene.count.feature=NULL,
                                   sj.count.matrix=NULL,
                                   sj.count.pheno=NULL,
                                   sj.count.feature=NULL,
                                   pca=NULL,
                                   gtf=NULL
                                   ) {
        
    # Annotate row and columns for matrices
        # Track progress
        message("Annotating rows and columns for matrices...")
        
        # Normalised gene expression matrix
        colnames(gene.norm.matrix) <- gene.norm.pheno$cell.id
        rownames(gene.norm.matrix) <- gene.norm.feature$gene_short_name
        
        # Raw (counts) gene expression matrix
        colnames(gene.count.matrix) <- gene.count.pheno$cell.id
        rownames(gene.count.matrix) <- gene.count.feature$gene_short_name
        
        # Raw (counts) SJ matrix
        colnames(sj.count.matrix) <- sj.count.pheno$cell.id
        rownames(sj.count.matrix) <- sj.count.feature$coord.intron
        
    # Create s3 object
    s3 <- list()
    class(s3) <- "Marvel"
    
    # Fill slots
    message("Filling in slots for MARVEL object...")
    
    s3$sample.metadata <- gene.norm.pheno
    s3$gene.metadata <- gene.norm.feature
    s3$gene.norm.matrix <- gene.norm.matrix
    s3$gene.count.matrix <- gene.count.matrix
    s3$sj.count.matrix <- sj.count.matrix
    s3$pca <- pca
    s3$gtf <- gtf
    
    # Returen final object
    message("MARVEL object created")
    return(s3)
        
}

######################################
# 给基因做注释
######################################

AnnotateGenes.10x.ZCY <- function(MarvelObject, gene_or_exon = "gene") {
    
    # Define arguments
    MarvelObject <- MarvelObject
    gene_or_exon <- gene_or_exon
    gtf <- MarvelObject$gtf
    gene.metadata <- MarvelObject$gene.metadata
    df.gene.norm <- MarvelObject$gene.norm
    
    # Example arguments
    # MarvelObject <- marvel
    # gtf <- MarvelObject$gtf
    # gene.metadata <- MarvelObject$gene.metadata
    # df.gene.norm <- MarvelObject$gene.norm
    
    ###############################################################
    
    gtf <- parsing_gtf(gtf=gtf, gene.metadata=gene.metadata, gene_or_exon=gene_or_exon)
  
    # Subset overlapping genes
    # overlap <- intersect(gene.metadata$gene_short_name, unique(gtf$gene_short_name))
    overlap <- intersect(unique(gtf$gene_id), gene.metadata$gene_id) 
    message(paste(nrow(gene.metadata), " genes found in gene metadata", sep=""))
    message(paste(length(overlap), " genes overlapped with GTF and will be subset-ed", sep=""))
    
    gene.metadata <- gene.metadata[which(gene.metadata$gene_id %in% overlap), , drop=FALSE]
    df.gene.norm <- df.gene.norm[gene.metadata$gene_short_name, ]
    
    # Annotate gene metadata
    # gene.metadata <- join(gene.metadata, gtf, by="gene_short_name", type="left", match="first")
    sum(is.na(gene.metadata$gene_type))
    
    message("Gene type annotated")
    
    # Check alignment
    table(rownames(df.gene.norm)==gene.metadata$gene_short_name)
    
    ###############################################################
    
    # Update slots
    MarvelObject$gene.metadata <- gene.metadata
    MarvelObject$gene.norm.matrix <- df.gene.norm
    
    message("Gene metadata and normalised gene expression matrix updated")
        
    # Return final object
    return(MarvelObject)
        
}

######################################
# 给junction做注释
######################################

AnnotateSJ.10x.ZCY <- function(MarvelObject, gene_or_exon = "exon") {
    
    # Define arguments
    MarvelObject <- MarvelObject
    gene_or_exon <- gene_or_exon
    gtf <- MarvelObject$gtf
    df.sj.count <- MarvelObject$sj.count.matrix
    gene.metadata <- MarvelObject$gene.metadata
    # Example arguments
    #MarvelObject <- marvel
    #gtf <- MarvelObject$gtf
    #df.sj.count <- MarvelObject$sj.count.matrix
    
    #########################################################
    
    # Create SJ metadata
    message("Creating splice junction metadata...")
    
    df <- data.frame("coord.intron"=rownames(df.sj.count), stringsAsFactors=FALSE)
    
    . <- strsplit(rownames(df.sj.count), split=":", fixed=TRUE)
    df$chr <- sapply(., function(x) {as.character(x[1])})
    df$start <- sapply(., function(x) {as.character(x[2])})
    df$end <- sapply(., function(x) {as.character(x[3])})
    
    # Parse GTF
        # Track progress
    message("Parsing GTF...")
    gtf <- parsing_gtf(gtf=gtf, gene.metadata=gene.metadata, gene_or_exon=gene_or_exon)
    gtf <- gtf[, c("V1", "V4", "V5", "gene_id",  "gene_type")]

    gtf <- plyr::join(gtf, gene.metadata)
    # Subset exon entries
    # gtf <- gtf[which(gtf$V3=="exon"), ]

    # Convert exon coordinates to intron (to match STAR output)
    gtf$V4 <- gtf$V4 - 1
    gtf$V5 <- gtf$V5 + 1

    # Check for chr prefix
    # logical <- any(grepl("chr", gtf$V1, fixed=TRUE))
            
    # Collapse start position
    # Track progress
    message("Matching gene names with SJ start coordinates in GTF...")
                    
    # Keep unique entries
    gtf.small <- gtf[, c("V1", "V5", "gene_short_name")]
        
        gtf.small$chr.pos <- paste(gtf$V1, ":", gtf$V5, sep="")
            
        gtf.small <- gtf.small[,c("chr.pos", "gene_short_name")]
        gtf.small <- unique(gtf.small)
        
        # Collapse
        gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
        . <- by(gtf.small[,"gene_short_name"], gtf.small[,"chr.pos"], function(x) {paste(x, collapse="|")})
        gtf.small.collapsed <- data.frame("chr"=as.character(names(.)), "gene_short_name"=as.character(.), stringsAsFactors=FALSE)
        
        # Save as new object
        gtf.small.collapsed.start <- gtf.small.collapsed
        names(gtf.small.collapsed.start) <- paste(names(gtf.small.collapsed.start), ".start", sep="")
        
    # Collapse end position
        # Track progress
        message("Matching gene names with SJ end coordinates in GTF...")
        
        # Keep unique entries
        gtf.small <- gtf[, c("V1", "V4", "gene_short_name")]
        
        gtf.small$chr.pos <- paste(gtf$V1, ":", gtf$V4, sep="")

        gtf.small <- gtf.small[,c("chr.pos", "gene_short_name")]
        gtf.small <- unique(gtf.small)
        
        # Collapse
        gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
        . <- by(gtf.small[,"gene_short_name"], gtf.small[,"chr.pos"], function(x) {paste(x, collapse="|")})
        gtf.small.collapsed <- data.frame("chr"=as.character(names(.)), "gene_short_name"=as.character(.), stringsAsFactors=FALSE)
        
        # Save as new object
        gtf.small.collapsed.end <- gtf.small.collapsed
        names(gtf.small.collapsed.end) <- paste(names(gtf.small.collapsed.end), ".end", sep="")
        
    # Annotate splice junctions
        # Track progress
        message("Annotating splice junctions...")
        
        # Retrieve SJ start, end coordinates
        df$chr.start <- paste(df$chr, ":", df$start, sep="")
        df$chr.end <- paste(df$chr, ":", df$end, sep="")
        #df$chr.start.end <- paste(df$chr, ":", df$start, ":", df$end, sep="")

        # Indicate id for re-ordering later
        df$id <- c(1:nrow(df))

        # Annotate SJ start,end with gene names
        #dim(gtf.small.collapsed.start) ; dim(gtf.small.collapsed.end) ; dim(df)
        df <- join(df, gtf.small.collapsed.start, by="chr.start", type="left")
        df <- join(df, gtf.small.collapsed.end, by="chr.end", type="left")
        #dim(df) ; sum(is.na(df$gene_short_name.start)) ; sum(is.na(df$gene_short_name.end))

        # Merge gene annotation
            # Create new column
            df$sj.type <- NA
            
            # Single-end record: Same gene at SJ start & end
            index.l <- !grepl("|", df$gene_short_name.start, fixed=TRUE) & !grepl("|", df$gene_short_name.end, fixed=TRUE) & !is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & df$gene_short_name.start==df$gene_short_name.end
            index <- which(index.l==TRUE)
            
            if(length(index) >= 1) {
                df$sj.type[index] <- "start_known.single.gene|end_known.single.gene|same"
            }
            
            # Single-end record: Different gene at SJ start & end
            index.l <- is.na(df$sj.type) & !grepl("|", df$gene_short_name.start, fixed=TRUE) & !grepl("|", df$gene_short_name.end, fixed=TRUE) & !is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & !df$gene_short_name.start==df$gene_short_name.end
            index <- which(index.l==TRUE)
            
            if(length(index) >= 1) {
                
                df$sj.type[index] <- "start_known.single.gene|end_known.single.gene|different"
                
            }
            
            # No record: Both sj
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & is.na(df$gene_short_name.end)
            index <- which(index.l==TRUE)
            
            if(length(index) >= 1) {
                
                df$sj.type[index] <- "start_unknown.gene|end_unknown.gene"
                
            }
                        
            # No record (start SJ) + end single gene
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & !grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_unknown.gene|end_known.single.gene"
                 
            }
            
            # No record (start SJ) + end multi gene
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_unknown.gene|end_known.multi.gene"
                 
            }
            
            # No record (end SJ) + start single gene
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & !grepl("|", df$gene_short_name.start, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.single.gene|end_unknown.gene"
                 
            }
            
            # No record (end SJ) + start multi gene
            index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.multi.gene|end_unknown.gene"
                 
            }
            
            # Multiple genes start and end
            index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, fixed=TRUE) & grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.multi.gene|end_known.multi.gene"
                 
            }
            
            # Multiple genes start + single gene end
            index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, fixed=TRUE) & !grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                 
                df$sj.type[index] <- "start_known.multi.gene|start_known.single.gene"
                 
            }

            # Multiple genes end + single gene start
            index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & !is.na(df$gene_short_name.start) & !grepl("|", df$gene_short_name.start, fixed=TRUE) & grepl("|", df$gene_short_name.end, fixed=TRUE)
            index <- which(index.l==TRUE)
             
            if(length(index) >= 1) {
                df$sj.type[index] <- "start_known.single.gene|end_known.multi.gene"
            }

            # Check for missing sj annotations
            if(sum(is.na(df$sj.type)) == 0) {
                message("All SJ successfully annotated")
            } else {
                message("Some SJ NOT successfully annotated")
            }

        # Formet input for MARVEL
        df <- df[,c("coord.intron", "gene_short_name.start", "gene_short_name.end", "sj.type")]
    
    #########################################################
    
    # Update slot
    MarvelObject$sj.metadata <- df
            
    # Return final object
    return(MarvelObject)
        
}

######################################
# 根绘制基因在百分之多少的细胞中表达密度图，用于自动确定最低基因表达量
######################################

PlotPctExprCells.Genes.10x.ZCY <- function(MarvelObject, diff_type, g1, g2, min.pct.cells=1) {

    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    diff_type <- diff_type
    g1 <- g1  # 样本1
    g2 <- g2  # 样本2
    min.pct.cells <- min.pct.cells  # 至少在百分之多少的细胞中表达的基因才保留
    
    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #cell.group.g1 <- cell.ids.1
    #cell.group.g2 <- cell.ids.2
    #min.pct.cells <- 10
    
    ################################################################
    
    sample.metadata.list <- subset_cell_according_diff_type(df=sample.metadata, 
                                                        diff_type=diff_type,
                                                        g1=g1,
                                                        g2=g2
                                                        )
    sample.metadata <- sample.metadata.list[["data"]]
    cell.group.g1 <- sample.metadata.list[["cells.g1"]]
    cell.group.g2 <- sample.metadata.list[["cells.g2"]]

    message("统计基因在细胞中的表达情况")
    # Compute num. of cells in which gene is expressed: Group 1
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g1]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .
    
    # Compute num. of cells in which gene is expressed: Group 2
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g2]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g2),  
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
        
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Censor non-expressing genes
    message(paste0("至少在 ", min.pct.cells , "% 细胞中表达的基因才被保留，用于后续密度图绘制"))
    results <- results[which(results$pct.cells.expr > min.pct.cells), ]  # 至少在百分之多少的细胞中表达的基因才保留 —— 筛基因

    # Density plot
        # Definitions
        data <- results
        x <- data$pct.cells.expr
        z <- data$cell.group
        maintitle <- ""
        ytitle <- "Density"
        xtitle <- "Gene Expressed (% Cells)"
        xmin <- 0 ; xmax <- max(x) ; xinterval <- 10
        legendtitle <- ""
       
        # Plot
        plot <- ggplot() +
            geom_density(data, mapping=aes(x=x, color=z), alpha=0.5) +
            scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
            labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=12),
                axis.text.x=element_text(size=8, colour="black"),
                axis.text.y=element_text(size=8, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                    
    # Save into new slot
    MarvelObject$pct.cells.expr$Gene$Plot <- plot
    MarvelObject$pct.cells.expr$Gene$Data <- results
    
    # Return final object
    return(MarvelObject)

}

######################################
# 根绘制junction在百分之多少的细胞中表达密度图，用于自动确定最低基因表达量
######################################

PlotPctExprCells.SJ.10x.ZCY <- function(MarvelObject, diff_type, g1, g2, min.pct.cells.genes=10, min.pct.cells.sj=10, downsample=FALSE, downsample.pct.sj=10, seed=1) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix   # 基因表达矩阵
    df.sj.count <- MarvelObject$sj.count.matrix     # SJ矩阵
    diff_type <- diff_type
    g1 <- g1  # 样本1
    g2 <- g2  # 样本2
    min.pct.cells.genes <- min.pct.cells.genes  # 上一步密度分布图得到的阈值
    min.pct.cells.sj <- min.pct.cells.sj
    downsample <- downsample
    downsample.pct.sj <- downsample.pct.sj
    seed <- seed
    
    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #sj.metadata <- MarvelObject$sj.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #cell.group.g1 <- cell.ids.1
    #cell.group.g2 <- cell.ids.2
    #min.pct.cells.genes <- 10
    #min.pct.cells.sj <- 5
    #downsample <- TRUE
    #downsample.pct.sj <- 10
    #seed <- 1
    
    ################################################################
    
    sample.metadata.list <- subset_cell_according_diff_type(df=sample.metadata, 
                                                        diff_type=diff_type,
                                                        g1=g1,
                                                        g2=g2
                                                        )
    sample.metadata <- sample.metadata.list[["data"]]
    cell.group.g1 <- sample.metadata.list[["cells.g1"]]
    cell.group.g2 <- sample.metadata.list[["cells.g2"]]

    # 像函数 PlotPctExprCells.Genes.10x() 一样，再按照样本重新统计一下各基因的表达量，并筛选至少在5%以上细胞中表达的基因
    # Compute num. of cells in which gene is expressed: Group 1
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g1]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .
    
    # Compute num. of cells in which gene is expressed: Group 2
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g2]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g2),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
        
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Censor non-expressing genes
    message(paste0("至少在 ", min.pct.cells.genes , "% 细胞中表达的基因才被保留，用于后续密度图绘制"))
    results <- results[which(results$pct.cells.expr > min.pct.cells.genes), ]
    
    # Save as new object
    # 基因表达表达，内含有样本信息、基因名字、共有多少个细胞、有多少个细胞表达了对应基因、表达基因占细胞总量的百分比
    # 每个junction都对应基因，这里保留主要表达的基因，并按照基因名字提取junction
    results.genes <- results

    ####################### DOWN-SAMPLE SJ ###########################
    
    # SJ 矩阵输出的为原始矩阵，未对细胞和UMI进行筛选
    # 上一步分析发现，在11%的细胞中表达的基因最多，因此就
    # 在 SJ 矩阵中随机抽取 11% 的SJ 【相当于抽基因，并没有对细胞进行筛选】
    # SJ 矩阵的行变少了
    # 随机抽取过程肯定会对结果带来影响

    if(downsample==TRUE) {
        
        # Set seed
        set.seed(seed)
        
        # Find no. of SJ to downsample
        size <- round(nrow(df.sj.count) * (downsample.pct.sj / 100), digits=0)
        coord.introns <- sample(rownames(df.sj.count), size=size, replace=FALSE)
        
        # Subset
        sj.metadata <- sj.metadata[which(sj.metadata$coord.intron %in% coord.introns), ]
        df.sj.count <- df.sj.count[coord.introns, ]
        
    }
    
    ################################################################
    
    # Compute num. of cells in which SJ is expressed: Group 1
        # Subset cells
        # 从随机抽取后的小矩阵中抽取 样本1 的矩阵，再按照至少在百分之5的细胞中表达，再次筛选 junction
        df.sj.count.small <- df.sj.count[, cell.group.g1]
        
        # Subset expressed genes
        # gene 矩阵和 SJ 矩阵都有基因名称，用基因名称对 SJ 矩阵进行过滤，只保留 gene 矩阵和 SJ 矩阵共有的内容
        gene_short_names <- results.genes[which(results.genes$cell.group=="cell.group.g1"), "gene_short_name"]
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        # 统计每个junction在细胞中的表达百分比
        . <- apply(df.sj.count.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "coord.intron"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .

    # Compute num. of cells in which SJ is expressed: Group 2
    # 以相同的思路处理样本2
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.g2]
        
        # Subset expressed genes
        gene_short_names <- results.genes[which(results.genes$cell.group=="cell.group.g2"), "gene_short_name"]
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "coord.intron"=names(.),
                        "n.cells.total"=length(cell.group.g2),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
            
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Censor non-expressing genes
    results <- results[which(results$pct.cells.expr > min.pct.cells.sj), ]
    
    # Save as new object
    results.sj <- results
        
    ################################################################
    
    # Density plot
        # Definitions
        data <- results.sj
        x <- data$pct.cells.expr
        z <- data$cell.group
        maintitle <- ""
        ytitle <- "Density"
        xtitle <- "SJ Expressed (% Cells)"
        xmin <- 0 ; xmax <- max(x) ; xinterval <- 10
        legendtitle <- ""
       
        # Plot
        plot <- ggplot() +
            geom_density(data, mapping=aes(x=x, color=z), alpha=0.5) +
            scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
            labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=12),
                axis.text.x=element_text(size=8, colour="black"),
                axis.text.y=element_text(size=8, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                    
    # Save into new slot
    MarvelObject$pct.cells.expr$SJ$Plot <- plot
    MarvelObject$pct.cells.expr$SJ$Data <- results
    
    # Return final object
    return(MarvelObject)
}

######################################
# 绘制 marvel pipeline 中的火山图；在源码的基础上进行修改
######################################

PlotDEValues.SJ.10x.ZCY <- function(MarvelObject, diff_type, g1, g2, min.pct.cells.gene , min.pct.cells.sj,  pval=0.05, log2fc=NULL, delta=NULL, min.gene.norm=NULL, anno=FALSE, anno.coord.intron=NULL, label.size=2) {

    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    diff_type <- diff_type
    g1 <- g1    # control
    g2 <- g2    # case
    min.pct.cells.gene <- min.pct.cells.gene
    min.pct.cells.sj <- min.pct.cells.sj
    df <- MarvelObject$DE$SJ$Table
    pval <- pval
    log2fc <- log2fc
    delta <- delta
    anno <- anno
    anno.coord.intron <- anno.coord.intron
    label.size <- label.size
    min.gene.norm <- min.gene.norm
    
    # Define arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$SJ$Table
    #pval <- 0.05
    #log2fc <- NULL
    #delta <- 5
    #anno <- FALSE
    #anno.coord.intron <- df$coord.intron[c(1:10)]
    #label.size <- 2
    #min.gene.norm <- 1
    
    ############################################################
    
    sample.metadata.list <- subset_cell_according_diff_type(df=sample.metadata, 
                                                        diff_type=diff_type,
                                                        g1=g1,
                                                        g2=g2
                                                        )
    g1 <- sample.metadata.list[['g1']]
    g2 <- sample.metadata.list[['g2']]

    # Subset/plot expressed genes
    df <- df[which(df$mean.expr.gene.norm.g1.g2 > min.gene.norm), ] # mean.expr.gene.norm.g1.g2 两个样本平均log(表达量)
    df <- df[which(df$pct.cells.expr.gene.g1 > min.pct.cells.gene), ]
    df <- df[which(df$pct.cells.expr.gene.g2 > min.pct.cells.gene), ]
    df <- df[which(df$pct.cells.expr.sj.g1 > min.pct.cells.sj), ]
    df <- df[which(df$pct.cells.expr.sj.g2 > min.pct.cells.sj), ]

    # Indicate sig events and direction
    # 有两种筛选方法，可以按照PSI差的绝对值进行筛选，也可以按照PSI差异倍数进行筛选

    # 只对SJ做评估
    if(!is.null(log2fc)) {
        # 按照差异可变剪切倍数筛选
        df$sig <- NA
        df$sig[which(df$pval < pval & df$log2fc > log2fc)] <- "Up"
        df$sig[which(df$pval < pval & df$log2fc < (log2fc*-1))] <- "Down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("Up", "Down", "n.s."))
        table(df$sig)
        
    } else if(!is.null(delta)){
        # 按照差异可变剪切绝对值筛选
        df$sig <- NA
        df$sig[which(df$pval < pval & df$delta > delta)] <- "Up"
        df$sig[which(df$pval < pval & df$delta < (delta*-1))] <- "Down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("Up", "Down", "n.s."))
        table(df$sig)
    }
    
    # Indicate color scheme
    sig.up <- which(df$sig=="Up")
    sig.down <- which(df$sig=="Down")
    if(length(sig.up) != 0 & length(sig.down) != 0) {
        col.breaks <- c("red", "blue", "gray")
    } else if(length(sig.up) != 0 & length(sig.down) == 0) {
        col.breaks <- c("red", "gray")
    } else if(length(sig.up) == 0 & length(sig.down) != 0) {
        col.breaks <- c("blue", "gray")
    } else if(length(sig.up) == 0 & length(sig.down) == 0) {
        col.breaks <- "gray"
    }
    
    # Create labels
    # 添加注释
    # 给junction添加基因名，需要注释的junction需要老师自己确定，并放在文件中
    if(anno==TRUE) {
       df$label <- NA
       df$label[which(df$coord.intron %in% anno.coord.intron)] <- df$gene_short_name[which(df$coord.intron %in% anno.coord.intron)]
    } else {
        df$label <- NA
    }

    # Compute mean norm gene expression
    # 计算两个样本的平均表达量
    numerator <- (df$mean.expr.gene.norm.g1 * df$n.cells.expr.gene.norm.g1) + (df$mean.expr.gene.norm.g2 * df$n.cells.expr.gene.norm.g2)
    denominator <- df$n.cells.expr.gene.norm.g1 + df$n.cells.expr.gene.norm.g2
    df$mean.expr.gene.norm.g1.g2 <- numerator/denominator
    
    ############################################################
    
    # Plot
       # Definition
       data <- df
       x <- data$mean.expr.gene.norm.g1.g2  # 默认图片的横坐标是基因在两个样本中的平均表达量
       x_lim_min = min(range(x))
       x_lim_max = max(range(x))
       label_location_x <- x_lim_min + (x_lim_max - x_lim_min)/15

       message(paste0("x轴范围从 ", x_lim_min, " 到 ", x_lim_max, "标签位置在 ", label_location_x ))

       z <- data$sig
       label <- data$label
       maintitle <- ""
       xtitle <- "Mean log2(Norm. GEX + 1)"
       
       if(!is.null(log2fc)) {
           y <- data$log2fc
           ytitle <- "log2fc"
       } else if(!is.null(delta)){
           y <- data$delta
           y[y > 100] <- 100    # 绘制火山图时，将大于100的点改为100
           ytitle <- "delta PSI"
       }

        g1 <- names(g1)
        g2 <- names(g2)

        range.y <- c(range(y), range(y)*-1)
        y_lim_max <- range.y[which.max(range.y)] 
        label_location_y.g2 <- y_lim_max
        y_lim_max <- ifelse(y_lim_max > 0, y_lim_max * 1.15, y_lim_max * 0.85)
        label_location_y.g2 <-  label_location_y.g2 + abs(y_lim_max - label_location_y.g2)/2

        y_lim_min <- range.y[which.min(range.y)] 
        label_location_y.g1 <- y_lim_min
        y_lim_min <- ifelse(y_lim_min > 0, y_lim_min * 0.85, y_lim_min * 1.15)
        label_location_y.g1 <- label_location_y.g1 - abs(label_location_y.g1-y_lim_min)

        message(paste0("y轴范围从 ", y_lim_min, " 到 ", y_lim_max ))
        message(paste0("标签的位置为 ", label_location_y.g2, "和", label_location_y.g1))
        message(paste0("标签分别为 ", g2, ' 和 ', g1))

       # Plot
       plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z), shape=20, alpha = 0.75, size=0.2) +
                ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                scale_colour_manual(values=col.breaks) +
                # coord_cartesian(ylim(c(y_lim_min, y_lim_max))) + 
                scale_x_continuous(n.breaks = 5, limits=c(x_lim_min, x_lim_max)) +
                scale_y_continuous(n.breaks = 7, limits=c(y_lim_min, y_lim_max)) +  # 使用 scale_y_continuous() 函数时直接使用 ylim 控制范围无效
                coord_cartesian(clip = 'off') + # 防止 annotate() 添加的文字溢出绘图区时显示不全
                # annotate("text", x = label_location_x, y = label_location_y.g2, label = g2, colour = "black") + 
                # annotate("text", x = label_location_x, y = label_location_y.g1, label = g1, colour = "black") + 
                annotate("text", x = -Inf, y = label_location_y.g1, label = g1, colour = "black", hjust = -.1) + 
                annotate("text", x = -Inf, y = label_location_y.g2, label = g2, colour = "black", hjust = -.1) + 
                labs(title=maintitle, x=xtitle, y=ytitle) +
                theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border=element_blank(),
                    # plot.title=element_text(hjust = 0.5, size=15),
                    # plot.subtitle=element_text(hjust = 0.5, size=15),
                    axis.line.y.left = element_line(color="black"),
                    axis.line.x = element_line(color="black"),
                    axis.title=element_text(size=12),
                    axis.text=element_text(size=12, colour="black"),
                    axis.text.x=element_text(size=10, colour="black"),
                    axis.text.y=element_text(size=10, colour="black"),
                    legend.position="none",
                    legend.title=element_text(size=8),
                    legend.text=element_text(size=8)
                    ) 


    # Save to new slot
    MarvelObject$DE$SJ$VolcanoPlot$SJ$Plot <- plot
    MarvelObject$DE$SJ$VolcanoPlot$SJ$Data <- df
    
    # Return final object
    return(MarvelObject)
}

######################################
# 图片展示细胞样本来源
######################################

PlotValues.PCA.CellGroup.10x.ZCY <- function(MarvelObject, g1, g2, dimension_type, diff_type, legendtitle, alpha=0.75, point.size=1.0, point.stroke=0.1, point.colors=NULL, point.size.legend=2) {

    # Example aruguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$sample.metadata
    g1 <- g1    # 样本1名称（control），用于绘图
    g2 <- g2   # 样本2名称（case），用于绘图
    dimension_type <- dimension_type
    diff_type <- diff_type
    legendtitle <- legendtitle
    alpha <- alpha
    point.size <- point.size
    point.stroke <- point.stroke
    point.colors <- point.colors
    point.size.legend <- point.size.legend
    
    # Example aruguments
    # MarvelObject <- marvel
    # g1 <- diff_control
    # g2 <- diff_case
    # df <- MarvelObject$sample.metadata
    # legendtitle <- ""
    # alpha <- 0.75
    # point.size <- 0.2
    # point.stroke <- 0
    # point.colors <- NULL
    # point.size.legend <- 2
    # dimension_type <- "tsne"
    
    ##########################################################################
    ########################### SUBSET CELLS #################################
    ##########################################################################
    
    # Annotate coordinate table
    . <- subset_cell_according_diff_type(df=df, diff_type=diff_type, g1=g1, g2=g2)
    df <- .[["data"]]
    
    # Set factor levels
    df$group_for_plot <- factor(df$group_for_plot, levels=c(unique(df$group_for_plot)))
    
    # Definition
    data <- df
    x <- data$x
    y <- data$y
    z <- data$group_for_plot
    
    if(dimension_type=="umap"){
        xtitle <- "UMAP-1"
        ytitle <- "UMAP-2"
    } else if(dimension_type=="tsne"){
        xtitle <- "tSNE-1"
        ytitle <- "tSNE-2"
    }
    
    # Define color scheme
    if(is.null(point.colors[1])) {
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
                            
      n = length(levels(z))
      point.colors = gg_color_hue(n)
    }

    # No. of columns for legend
    if(length(unique(z)) > 20) {
        ncol.legend <- 3
    } else if(length(unique(z)) > 10) {
        ncol.legend <- 2
    } else if(length(unique(z)) <= 10) {
        ncol.legend <- 1
    }

    # Define margins before removing cells
    xmin <- min(data$x); xmax <- max(data$x)
    ymin <- min(data$y); ymax <- max(data$y)
    
    # Plot (with legends)
    plot <- ggplot() +
          geom_point(data, mapping=aes(x=x, y=y, fill=z), color="black", pch=21, size=point.size, alpha=alpha, stroke=point.stroke) +
          scale_fill_manual(values=point.colors) +
          scale_x_continuous(limits=c(xmin, xmax)) +
          scale_y_continuous(limits=c(ymin, ymax)) +
          labs(title=NULL, x=xtitle, y=ytitle, fill=legendtitle) +
          ggtitle("Dimensionality reduction") + 
          theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              plot.title = element_text(size=15, hjust=0.5),
              axis.line=element_line(colour = "black"),
              axis.title=element_text(size=12),
              axis.text.x=element_text(size=10, colour="black"),
              axis.text.y=element_text(size=10, colour="black"),
              legend.title=element_text(size=8),
              legend.text=element_text(size=8)
              ) +
           guides(fill = guide_legend(override.aes=list(size=point.size.legend, alpha=alpha, stroke=point.stroke), ncol=ncol.legend))

    # Save into new slot
    MarvelObject$adhocPlot$PCA$CellGroup <- plot

    # Return final object
    return(MarvelObject)
}

######################################
# 基因表达量映射到降维图形中
######################################

PlotValues.PCA.Gene.10x.ZCY <- function(MarvelObject, g1, g2, gene_short_name, dimension_type, diff_type, log2.transform=TRUE, cell.ids=NULL, point.size=0.1, color.gradient=c("grey90","blue","red")) {

    # Example aruguments
    MarvelObject <- MarvelObject
    g1 <- g1
    g2 <- g2
    gene_short_name <- gene_short_name
    dimension_type <- dimension_type
    diff_type <- diff_type
    df.coord <- MarvelObject$sample.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    log2.transform <- log2.transform
    cell.ids <- cell.ids
    point.size <- point.size
    color.gradient <- color.gradient
    
    # Example aruguments
    # MarvelObject <- marvel
    # g1 <- diff_control
    # g2 <- diff_case
    # df.coord <- MarvelObject$pca
    # df.gene.norm <- MarvelObject$gene.norm.matrix
    # cell.ids <- NULL
    # gene_short_name <- "ATP5F1C"
    # point.size <- 0.1
    # color.gradient <- c("grey","cyan","green","yellow","red")
    # dimension_type <- "tsne"
    # diff_type <- "sample"
    # log2.transform <- TRUE
    
    ##########################################################################
    ###################### ANNOTATE EXPRESSION ###############################
    ##########################################################################
    
    # Subset relevant SJ
    df.gene.norm <- df.gene.norm[gene_short_name, ]
    
    # Tabulate counts
    df.gene.norm  <- data.frame("cell.id"=names(df.gene.norm),
                                "expr.gene.norm"=as.numeric(df.gene.norm),
                                stringsAsFactors=FALSE
                                )
    
    # Annotate coordinate table
    df.coord <- join(df.coord, df.gene.norm, by="cell.id", type="left")

    df.coord <- subset_cell_according_diff_type(df=df.coord, diff_type=diff_type, g1=g1, g2=g2)[["data"]]

    ##########################################################################
    ####################### SET THRESHOLD FOR COUNTS #########################
    ##########################################################################
    
    # Define margins before removing cells
    xmin <- min(df.coord$x); xmax <- max(df.coord$x)
    ymin <- min(df.coord$y); ymax <- max(df.coord$y)
    
    ##########################################################################
    ################################ PLOT ####################################
    ##########################################################################
    
    # Reorder by expression
    # 排序对点的展示效果有作用，相互重叠的很多，让极端值更多的展示出来
    df.coord <- df.coord[order(df.coord$expr.gene.norm), ]
        
    # Definitions
    data <- df.coord
    x <- data$x
    y <- data$y
    
    # 原始的norm矩阵没有经过log2转换，这里只能设置为true
    if(log2.transform==TRUE) {
        z <- log2(data$expr.gene.norm + 1)
    }
    
    if(dimension_type=="umap"){
        xtitle <- "UMAP-1"
        ytitle <- "UMAP-2"
    } else if(dimension_type=="tsne"){
        xtitle <- "tSNE-1"
        ytitle <- "tSNE-2"
    }

    legendtitle  <- "log2(expr)"
    
    # Plot
    plot <- ggplot() +
        geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size) +
        scale_color_gradientn(colors=color.gradient) +
        scale_x_continuous(limits=c(xmin, xmax)) +
        scale_y_continuous(limits=c(ymin, ymax)) +
        labs(title=NULL, x=xtitle, y=ytitle, color=legendtitle) +
        ggtitle(gene_short_name) +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            plot.title=element_text(hjust = 0.5, size=15),
            plot.subtitle=element_text(hjust = 0.5, size=15),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_text(size=12),
            axis.text.x=element_text(size=10, colour="black"),
            axis.text.y=element_text(size=10, colour="black"),
            legend.title=element_text(size=8),
            legend.text=element_text(size=8)
            )

    ##########################################################################

    # Save into new slot
    MarvelObject$adhocPlot$PCA$Gene <- plot

    # Return final object
    return(MarvelObject)

}

######################################
# PSI 量映射到降维图形中
######################################

PlotValues.PCA.PSI.10x.ZCY <- function(MarvelObject, g1, g2, dimension_type, diff_type, psi_or_expr="psi", cell.ids=NULL, coord.intron, min.gene.count=NULL, point.size=0.1, log2.transform=TRUE, color.gradient=c("grey90","blue","red")) {

    # Example aruguments
    MarvelObject <- MarvelObject
    g1 <- g1
    g2 <- g2
    dimension_type <- dimension_type
    diff_type <- diff_type
    psi_or_expr <- psi_or_expr
    df.coord <- MarvelObject$sample.metadata
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.ids <- cell.ids
    coord.intron <- coord.intron
    min.gene.count <- min.gene.count
    point.size <- point.size
    log2.transform <- log2.transform
    color.gradient <- color.gradient
    

    # Example aruguments
    # MarvelObject <- marvel
    # g1 <- "1"
    # g2 <- "2"
    # dimension_type <- "tsne"
    # diff_type <- "cluster"
    # df.coord <- MarvelObject$sample.metadata
    # sj.metadata <- MarvelObject$sj.metadata
    # df.gene.count <- MarvelObject$gene.count.matrix
    # df.sj.count <- MarvelObject$sj.count.matrix
    # coord.intron <- "chr6:30671275:30672634"
    # min.gene.count <- 0.01
    # point.size <- 0.1
    # log2.transform <- TRUE
    # color.gradient <- c("grey90","blue","red")
    
    ##########################################################################
    ####################### ANNOTATE SJ COUNTS ###############################
    ##########################################################################
    
    # Subset relevant SJ
    df.sj.count <- df.sj.count[coord.intron, ]
    
    # Tabulate counts
    df.sj.count  <- data.frame("cell.id"=names(df.sj.count),
                               "sj.count"=as.numeric(df.sj.count),
                               stringsAsFactors=FALSE
                               )
    
    # Annotate coordinate table
    df.coord <- join(df.coord, df.sj.count, by="cell.id", type="left")

    df.coord <- subset_cell_according_diff_type(df=df.coord, diff_type=diff_type, g1=g1, g2=g2)[["data"]]

    ##########################################################################
    ####################### ANNOTATE SJ COUNTS ###############################
    ##########################################################################
    
    # Retrieve gene
    # 找到 junction 对应的 gene
    gene_short_name <- sj.metadata[which(sj.metadata$coord.intron==coord.intron), "gene_short_name.start"]
    
    # Subset relevant gene
    # 提取对应 gene 表达矩阵
    df.gene.count <- df.gene.count[gene_short_name, ]
    
    # Tabulate counts
    df.gene.count  <- data.frame("cell.id"=names(df.gene.count),
                                 "gene.count"=as.numeric(df.gene.count),
                                  stringsAsFactors=FALSE
                                  )
    
    # Annotate coordinate table
    df.coord <- join(df.coord, df.gene.count, by="cell.id", type="left")
    
    ##########################################################################
    ############################# COMPUTE PSI ################################
    ##########################################################################
    
    # 表格中记录的PSI值为cluster水平的，并是细胞水平的，这里在细胞水平重新计算
    df.coord$psi <- df.coord$sj.count / df.coord$gene.count * 100

    # 对应基因的表达量为0时，psi就是NA，绘图时不显示，需要修改为0
    df.coord[which(is.na(df.coord$psi)),]$psi <- 0

    ##########################################################################
    ####################### SET THRESHOLD FOR COUNTS #########################
    ##########################################################################
    
    # Define margins before removing cells
    xmin <- min(df.coord$x); xmax <- max(df.coord$x)
    ymin <- min(df.coord$y); ymax <- max(df.coord$y)
    
    # Censor lowly expressed gene
    if( !is.null(min.gene.count )){
        df.coord$psi[which(df.coord$gene.count < min.gene.count)] <- NA
        df.coord <- df.coord[!is.na(df.coord$psi), ]
    }

    ##########################################################################
    ############################## SUBSET CELLS ##############################
    ##########################################################################
    
    # if(!is.null(cell.ids[1])) {
    #     df.coord <- df.coord[which(df.coord$cell.id %in% cell.ids), ]
    # }
    
    ##########################################################################
    ################################ PLOT ####################################
    ##########################################################################
    
    # Cap max PSI
    df.coord$psi[which(df.coord$psi > 100)] <- 100
    
    # Reorder by expression
    if (psi_or_expr=="psi"){
        df.coord <- df.coord[order(df.coord$psi), ]
    } else if(psi_or_expr=="expr"){
        df.coord <- df.coord[order(df.coord$sj.count), ]
    }
    
    # Definitions
    data <- df.coord
    x <- data$x
    y <- data$y
    
    if(log2.transform==TRUE & psi_or_expr=="psi") {
        z <- log2(data$psi + 1)
        legendtitle <- "log2(PSI)"
    } else if (log2.transform==TRUE & psi_or_expr=="expr"){
        z <- log2(data$sj.count + 1)
        legendtitle <- "log2(expr)"
    } else if(log2.transform==FALSE & psi_or_expr=="psi"){
        z <- data$psi
        legendtitle <- "PSI"
    } else if (log2.transform==FALSE & psi_or_expr=="expr"){
        z <- data$sj.count
        legendtitle <- "expr"
    }


    if(dimension_type=="umap"){
        xtitle <- "UMAP-1"
        ytitle <- "UMAP-2"
    } else if(dimension_type=="tsne"){
        xtitle <- "tSNE-1"
        ytitle <- "tSNE-2"
    }
    

    # Plot
    plot <- ggplot() +
        geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size) +
        scale_color_gradientn(colors=color.gradient) +
        scale_x_continuous(limits=c(xmin, xmax)) +
        scale_y_continuous(limits=c(ymin, ymax)) +
        labs(title=NULL, x=xtitle, y=ytitle, color=legendtitle) +
        ggtitle(coord.intron) +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "white", fill = "white"),
            panel.border=element_blank(),
            plot.title=element_text(hjust = 0.5, size=15),
            plot.subtitle=element_text(hjust = 0.5, size=15),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_text(size=12),
            axis.text.x=element_text(size=10, colour="black"),
            axis.text.y=element_text(size=10, colour="black"),
            legend.title=element_text(size=8),
            legend.text=element_text(size=8),
            legend.background = element_rect(colour = "white", fill = "white")
            )

    ##########################################################################

    # Save into new slot
    MarvelObject$adhocPlot$PCA$PSI <- plot

    # Return final object
    return(MarvelObject)

}

######################################
# 统计junction表达量，并做显著性分析
######################################

CompareValues.SJ.10x.ZCY <- function(MarvelObject, g1, g2, diff_type, min.pct.cells.genes=0, min.pct.cells.sj=0, min.gene.norm=-100, seed=1, n.iterations=100, coord.introns=NULL, downsample=FALSE, show.progress=TRUE) {
    
    # 先按照基因表达量对基因进行简单筛选，然后按照基因名从SJ矩阵中提取SJ，
    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    g1 <- g1  # control
    g2 <- g2  # case
    diff_type <- diff_type
    min.pct.cells.genes <- min.pct.cells.genes
    min.pct.cells.sj <- min.pct.cells.sj
    seed <- seed
    n.iterations <- n.iterations
    downsample <- downsample
    show.progress <- show.progress
    min.gene.norm <- min.gene.norm
    
    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #sj.metadata <- MarvelObject$sj.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.gene.count <- MarvelObject$gene.count.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #coord.introns.custom <- NULL
    #cell.group.g1 <- cell.ids.g1
    #cell.group.g2 <- cell.ids.g2
    #min.pct.cells.genes <- 0.01
    #min.pct.cells.sj <- 0.01
    #seed <- 1
    #n.iterations <- 100
    #downsample <- FALSE
    #min.gene.norm <- 0.01
    #show.progress <- FALSE
    
    #################################################################
    ################### SUBSET SPECIFC SJs ##########################
    #################################################################
    

    
    sample.metadata.list <- subset_cell_according_diff_type(df=sample.metadata, 
                                                        diff_type=diff_type,
                                                        g1=g1,
                                                        g2=g2
                                                        )
    sample.metadata <- sample.metadata.list[["data"]]
    cell.group.g1 <- sample.metadata.list[["cells.g1"]]
    cell.group.g2 <- sample.metadata.list[["cells.g2"]]


    # Report progress
    message(paste(length(cell.group.g1), " cells from ", g1," and ", length(cell.group.g2), " cells from ", g2," included", sep=""))
    
    # Subset matrices
    df.gene.norm <- df.gene.norm[, c(cell.group.g1, cell.group.g2)]
    df.gene.count <- df.gene.count[, c(cell.group.g1, cell.group.g2)]
    df.sj.count <- df.sj.count[, c(cell.group.g1, cell.group.g2)]
    
    # Check alignment
    table(colnames(df.gene.norm)==colnames(df.gene.count))
    table(colnames(df.gene.count)==colnames(df.sj.count))

    #################################################################
    ################## SUBSET EXPRESSED GENES (1) ###################
    #################################################################
    
    # Compute num. of cells in which gene is expressed: Group 1
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g1]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .
    
    # Compute num. of cells in which gene is expressed: Group 2
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g2]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g2),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
        
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Subset expressed genes
    results <- results[which(results$pct.cells.expr > min.pct.cells.genes), ]
    
    # Subset genes expressing in both cell groups
    gene_short.names.1 <- results[which(results$cell.group=="cell.group.g1"), "gene_short_name"]
    gene_short.names.2 <- results[which(results$cell.group=="cell.group.g2"), "gene_short_name"]
    gene_short_names <- intersect(gene_short.names.1, gene_short.names.2)
    
    # Report progress
    message(paste(length(gene_short.names.1), " genes expressed in cell group 1", sep=""))
    message(paste(length(gene_short.names.2), " genes expressed in cell group 2", sep=""))
    message(paste(length(gene_short_names), " genes expressed in BOTH cell group and retained", sep=""))
    
    #################################################################
    ################## SUBSET EXPRESSED GENES (2) ###################
    #################################################################
    
    # Subset expressed genes from part 1
    df.gene.norm.small <- df.gene.norm[gene_short_names, ]
    
    # Compute combined average
    . <- apply(df.gene.norm.small, 1, function(x) {mean(log2(x + 1))})
    mean.combined.df <- data.frame("gene_short_name"=names(.),
                                   "mean.expr.gene.norm.g1.g2"=as.numeric(.),
                                   stringsAsFactors=FALSE
                                   )
    
    # Subset expressed genes

    index <- which(mean.combined.df$mean.expr.gene.norm.g1.g2 > min.gene.norm)
    mean.combined.df <- mean.combined.df[index, ]
    gene_short_names <- mean.combined.df$gene_short_name
    
    # Report progress
    message(paste(length(gene_short_names), " genes with mean log2(expression + 1) > ", min.gene.norm, " retained", sep=""))
    
    #################################################################
    ###################### SUBSET EXPRESSED SJ ######################
    #################################################################
    
    # Compute num. of cells in which SJ is expressed: Group 1
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.g1]
        
        # Subset expressed genes
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "coord.intron"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .

    # Compute num. of cells in which SJ is expressed: Group 2
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.g2]
        
        # Subset expressed genes
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "coord.intron"=names(.),
                        "n.cells.total"=length(cell.group.g2),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
            
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Subset expressed SJ
    results <- results[which(results$pct.cells.expr > min.pct.cells.sj), ]
    
    # Subset SJ expressing in either cell groups
    coord.introns.1 <- results[which(results$cell.group=="cell.group.g1"), "coord.intron"]
    coord.introns.2 <- results[which(results$cell.group=="cell.group.g2"), "coord.intron"]
    coord.introns <- unique(c(coord.introns.1, coord.introns.2))
    
    # Report progress
    message(paste(length(coord.introns.1), " SJ expressed in cell group 1", sep=""))
    message(paste(length(coord.introns.2), " SJ expressed in cell group 2", sep=""))
    message(paste(length(coord.introns), " SJ expressed in EITHER cell groups and retained", sep=""))
    
    # Report final numbers
    n.sj <- length(coord.introns)
    n.genes <- length(unique(sj.metadata[which(sj.metadata$coord.intron %in% coord.introns), "gene_short_name.start"]))
    message(paste("Total of ", n.sj, " SJ from ", n.genes, " genes included for DE analysis", sep=""))
        
    #################################################################
    ############## SUBSET EXPRESSED GENES, SJ, CELLS ################
    #################################################################
    
    # Subset cells, expressed genes, SJs
        # Metadata
        sj.metadata <- sj.metadata[which(sj.metadata$coord.intron %in% coord.introns), ]
    
        # Group 1
        df.sj.count.g1 <- df.sj.count[sj.metadata$coord.intron, cell.group.g1]
        df.gene.count.g1 <- df.gene.count[unique(sj.metadata$gene_short_name.start), cell.group.g1]
        
        # Group 2
        df.sj.count.g2 <- df.sj.count[sj.metadata$coord.intron, cell.group.g2]
        df.gene.count.g2 <- df.gene.count[unique(sj.metadata$gene_short_name.start), cell.group.g2]
        
        # Check alignment
        table(row.names(df.sj.count.g1)==row.names(df.sj.count.g2))
        table(row.names(df.gene.count.g1)==row.names(df.gene.count.g2))

    #################################################################
    ######################### COMPUTE PSI ###########################
    #################################################################
    
    # Compute PSI: Group 1
        # Report progress
        message("Computing PSI for cell group 1...")
        
        # Compute cell group size
        n.cells.total <- ncol(df.sj.count.g1)
        
        # Tabulate SJ metrices
            # Compute % expressed SJ
            n.cells.expr.sj <- apply(df.sj.count.g1, 1, function(x) {sum(x != 0)})
            pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
            
            # Compute total SJ counts
            sj.count.total <- apply(df.sj.count.g1, 1, function(x) {sum(x)})
            
            # Save into data frame
            results.sj <- data.frame("coord.intron"=names(sj.count.total),
                                     "n.cells.total"=n.cells.total,
                                     "n.cells.expr.sj"=n.cells.expr.sj,
                                     "pct.cells.expr.sj"=pct.cells.expr.sj,
                                     "sj.count.total"=sj.count.total,
                                     stringsAsFactors=FALSE
                                     )
            row.names(results.sj) <- NULL
            
            # Annotate gene
            results.sj <- join(results.sj, sj.metadata[, c("coord.intron", "gene_short_name.start")], by="coord.intron", type="left")
            
        # Tabulate gene metrices
            # Compute % expressed genes
            n.cells.expr.gene <- apply(df.gene.count.g1, 1, function(x) {sum(x != 0)})
            pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
            
            # Compute total gene counts
            gene.count.total <- apply(df.gene.count.g1, 1, function(x) {sum(x)})
        
            # Tabulate results
            results.gene <- data.frame("gene_short_name.start"=names(gene.count.total),
                                       "n.cells.expr.gene"=n.cells.expr.gene,
                                       "pct.cells.expr.gene"=pct.cells.expr.gene,
                                       "gene.count.total"=gene.count.total,
                                       stringsAsFactors=FALSE
                                       )
                                       
        # Merge
        results <- join(results.sj, results.gene, by="gene_short_name.start", type="left")
        
        # Compute PSI
        results$psi <- round(results$sj.count.total / results$gene.count.total * 100, digits=2)
        
        # Reorder columns
        cols <- c("coord.intron", "gene_short_name.start", "n.cells.total",
                  "n.cells.expr.sj", "pct.cells.expr.sj",
                  "n.cells.expr.gene", "pct.cells.expr.gene",
                  "sj.count.total", "gene.count.total", "psi"
                  )
        
        results <- results[, cols]
        
        names(results)[which(names(results)=="gene_short_name.start")] <- "gene_short_name"
        
        # Indicate cell group
        names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))], ".g1", sep="")
     
         # Save as new object
         results.g1 <- results
     
    # Compute PSI: Group 2
        # Report progress
        message("Computing PSI for cell group 2...")
        
        # Compute cell group size
        n.cells.total <- ncol(df.sj.count.g2)
        
        # Tabulate SJ metrices
            # Compute % expressed SJ
            n.cells.expr.sj <- apply(df.sj.count.g2, 1, function(x) {sum(x != 0)})
            pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
            
            # Compute total SJ counts
            sj.count.total <- apply(df.sj.count.g2, 1, function(x) {sum(x)})
            
            # Save into data frame
            results.sj <- data.frame("coord.intron"=names(sj.count.total),
                                     "n.cells.total"=n.cells.total,
                                     "n.cells.expr.sj"=n.cells.expr.sj,
                                     "pct.cells.expr.sj"=pct.cells.expr.sj,
                                     "sj.count.total"=sj.count.total,
                                     stringsAsFactors=FALSE
                                     )
            row.names(results.sj) <- NULL
            
            # Annotate gene
            results.sj <- join(results.sj, sj.metadata[, c("coord.intron", "gene_short_name.start")], by="coord.intron", type="left")
            
        # Tabulate gene metrices
            # Compute % expressed genes
            n.cells.expr.gene <- apply(df.gene.count.g2, 1, function(x) {sum(x != 0)})
            pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
            
            # Compute total gene counts
            gene.count.total <- apply(df.gene.count.g2, 1, function(x) {sum(x)})
        
            # Tabulate results
            results.gene <- data.frame("gene_short_name.start"=names(gene.count.total),
                                       "n.cells.expr.gene"=n.cells.expr.gene,
                                       "pct.cells.expr.gene"=pct.cells.expr.gene,
                                       "gene.count.total"=gene.count.total,
                                       stringsAsFactors=FALSE
                                       )
                                       
        # Merge
        results <- join(results.sj, results.gene, by="gene_short_name.start", type="left")
        
        # Compute PSI
        results$psi <- round(results$sj.count.total / results$gene.count.total * 100, digits=2)
        
        # Reorder columns
        cols <- c("coord.intron", "gene_short_name.start", "n.cells.total",
                  "n.cells.expr.sj", "pct.cells.expr.sj",
                  "n.cells.expr.gene", "pct.cells.expr.gene",
                  "sj.count.total", "gene.count.total", "psi"
                  )
        
        results <- results[, cols]
        
        names(results)[which(names(results)=="gene_short_name.start")] <- "gene_short_name"
        
        # Indicate cell group
        names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))], ".g2", sep="")
        
        # Save as new object
        results.g2 <- results
        
    # Merge group 1, 2
    table(results.g1$coord.intron==results.g2$coord.intron)
    table(results.g1$gene_short_name==results.g2$gene_short_name)
    
    
    index.l <- table(results.g1$coord.intron==results.g2$coord.intron)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
    
        results.g2$coord.intron <- NULL
        results.g2$gene_short_name <- NULL
        results <- cbind.data.frame(results.g1, results.g2)
        
    } else {
        message("Error in merging tables from Group 1 and 2")
    }
    
    # Compute log2fc, delta
    results$log2fc <- log2( (results$psi.g2 + 1) / (results$psi.g1 + 1) )
    results$delta <- results$psi.g2 - results$psi.g1
    
    # Save as new object
    results.obs <- results
    
    #################################################################
    ##################### COMPUTE P-VALUES ##########################
    #################################################################
    
    # Report progress
    message("Creating null distributions...")
     
    # Set random num. generator
    set.seed(seed)
    
    if(show.progress==TRUE) {
        pb <- txtProgressBar(1, n.iterations, style=3)
    }
    
    .list.results.perm <- list()
    
    for(i in 1:n.iterations) {
        
        # Create null distribution
            # Shuffle cell ids
            cell.ids.shuffled <- sample(colnames(df.sj.count), size=ncol(df.sj.count), replace=FALSE)
            
            # Shuffle SJ matrix
            df.sj.count.shuffled <- df.sj.count
            colnames(df.sj.count.shuffled) <- cell.ids.shuffled
            
            # Match column in gene matrix
            df.gene.count.shuffled <- df.gene.count
            colnames(df.gene.count.shuffled) <- cell.ids.shuffled
            
            # Check alignment
            table(colnames(df.sj.count.shuffled)==colnames(df.gene.count.shuffled))
            
            # Report progress
            index.l <- table(colnames(df.sj.count.shuffled)==colnames(df.gene.count.shuffled))
            index.true <- length(which(names(index.l)==TRUE))
            index.false <- length(which(names(index.l)==FALSE))
             
            if(index.true==1 & index.false==0) {
                #message(paste("Iteration ", 1, " ...", sep=""))
            } else {
                return(message(paste("Error in iteration ", 1, " ...", sep="")))
            }
        
        # Subset cells, expressed genes, SJs
            # Group 1
            df.sj.count.g1 <- df.sj.count.shuffled[results.obs$coord.intron, cell.group.g1]
            df.gene.count.g1 <- df.gene.count.shuffled[unique(results.obs$gene_short_name), cell.group.g1]
            
            # Group 2
            df.sj.count.g2 <- df.sj.count.shuffled[results.obs$coord.intron, cell.group.g2]
            df.gene.count.g2 <- df.gene.count.shuffled[unique(results.obs$gene_short_name), cell.group.g2]
            
            # Check alignment
            table(row.names(df.sj.count.g1)==row.names(df.sj.count.g2))
            table(row.names(df.gene.count.g1)==row.names(df.gene.count.g2))
            
        # Compute PSI: Group 1
            # Report progress
            #message("Computing PSI for cell group 1...")
            
            # Compute cell group size
            #n.cells.total <- ncol(df.sj.count.g1)
            
            # Tabulate SJ metrices
                # Compute % expressed SJ
                #n.cells.expr.sj <- apply(df.sj.count.g1, 1, function(x) {sum(x != 0)})
                #pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
                
                # Compute total SJ counts
                sj.count.total <- apply(df.sj.count.g1, 1, function(x) {sum(x)})
                
                # Save into data frame
                results.sj <- data.frame("coord.intron"=names(sj.count.total),
                                         #"n.cells.total"=n.cells.total,
                                         #"n.cells.expr.sj"=n.cells.expr.sj,
                                         #"pct.cells.expr.sj"=pct.cells.expr.sj,
                                         "sj.count.total"=sj.count.total,
                                         stringsAsFactors=FALSE
                                         )
                row.names(results.sj) <- NULL
                
                # Annotate gene
                results.sj <- join(results.sj, sj.metadata[, c("coord.intron", "gene_short_name.start")], by="coord.intron", type="left")
                
            # Tabulate gene metrices
                # Compute % expressed genes
                #n.cells.expr.gene <- apply(df.gene.count.g1, 1, function(x) {sum(x != 0)})
                #pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
                
                # Compute total gene counts
                gene.count.total <- apply(df.gene.count.g1, 1, function(x) {sum(x)})
            
                # Tabulate results
                results.gene <- data.frame("gene_short_name.start"=names(gene.count.total),
                                           #"n.cells.expr.gene"=n.cells.expr.gene,
                                           #"pct.cells.expr.gene"=pct.cells.expr.gene,
                                           "gene.count.total"=gene.count.total,
                                           stringsAsFactors=FALSE
                                           )
                                           
            # Merge
            results <- join(results.sj, results.gene, by="gene_short_name.start", type="left")
            
            # Compute PSI
            results$psi <- round(results$sj.count.total / results$gene.count.total * 100, digits=2)
            
            # Reorder columns
            #cols <- c("coord.intron", "gene_short_name.start", "n.cells.total",
                      #"n.cells.expr.sj", "pct.cells.expr.sj",
                      #"n.cells.expr.gene", "pct.cells.expr.gene",
                      #"sj.count.total", "gene.count.total", "psi"
                      #)
                                  
            #results <- results[, cols]
            
            names(results)[which(names(results)=="gene_short_name.start")] <- "gene_short_name"
            
            # Indicate cell group
            names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))], ".g1", sep="")
         
            # Save as new object
            results.g1 <- results
            #results.g1 <- results[, c("coord.intron", "psi.g1")]

        # Compute PSI: Group 2
            # Report progress
            #message("Computing PSI for cell group 2...")
            
            # Compute cell group size
            #n.cells.total <- ncol(df.sj.count.g2)
            
            # Tabulate SJ metrices
                # Compute % expressed SJ
                #n.cells.expr.sj <- apply(df.sj.count.g2, 1, function(x) {sum(x != 0)})
                #pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
                
                # Compute total SJ counts
                sj.count.total <- apply(df.sj.count.g2, 1, function(x) {sum(x)})
                
                # Save into data frame
                results.sj <- data.frame("coord.intron"=names(sj.count.total),
                                         #"n.cells.total"=n.cells.total,
                                         #"n.cells.expr.sj"=n.cells.expr.sj,
                                         #"pct.cells.expr.sj"=pct.cells.expr.sj,
                                         "sj.count.total"=sj.count.total,
                                         stringsAsFactors=FALSE
                                         )
                row.names(results.sj) <- NULL
                
                # Annotate gene
                results.sj <- join(results.sj, sj.metadata[, c("coord.intron", "gene_short_name.start")], by="coord.intron", type="left")
                
            # Tabulate gene metrices
                # Compute % expressed genes
                #n.cells.expr.gene <- apply(df.gene.count.g2, 1, function(x) {sum(x != 0)})
                #pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
                
                # Compute total gene counts
                gene.count.total <- apply(df.gene.count.g2, 1, function(x) {sum(x)})
            
                # Tabulate results
                results.gene <- data.frame("gene_short_name.start"=names(gene.count.total),
                                           #"n.cells.expr.gene"=n.cells.expr.gene,
                                           #"pct.cells.expr.gene"=pct.cells.expr.gene,
                                           "gene.count.total"=gene.count.total,
                                           stringsAsFactors=FALSE
                                           )
                                           
            # Merge
            results <- join(results.sj, results.gene, by="gene_short_name.start", type="left")
            
            # Compute PSI
            results$psi <- round(results$sj.count.total / results$gene.count.total * 100, digits=2)
            
            # Reorder columns
            #cols <- c("coord.intron", "gene_short_name.start", "n.cells.total",
                      #"n.cells.expr.sj", "pct.cells.expr.sj",
                      #"n.cells.expr.gene", "pct.cells.expr.gene",
                      #"sj.count.total", "gene.count.total", "psi"
                      #)
            
            #results <- results[, cols]
            
            names(results)[which(names(results)=="gene_short_name.start")] <- "gene_short_name"
            
            # Indicate cell group
            names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))], ".g2", sep="")
            
            # Save as new object
            results.g2 <- results
            #results.g2 <- results[, c("coord.intron", "psi.g2")]
            
        # Merge
        results <- join(results.g1, results.g2, by="coord.intron", type="left")
        
        # Compute permutated delta
        results$delta.perm <- results$psi.g2 - results$psi.g1
        
        # Save into list
        row.names(results) <- results$coord.intron
        results <- results[, "delta.perm", drop=FALSE]
        
        .list.results.perm[[i]] <- results
             
        if(show.progress==TRUE) {
            # Track progress
            setTxtProgressBar(pb, i)
        }
    }
        
    results.perm <- do.call(cbind.data.frame, .list.results.perm)
    
    # Compute p-value
        # Report progress
        message("Computing P values...")
        
        # Annotate delta observed
        index.l <- table(results.obs$coord.intron==row.names(results.perm))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
            results.perm <- cbind.data.frame(results.perm, results.obs[,"delta",drop=FALSE])
        } else {
            return(message("Error in consolidating observed and permutated results"))
        }
        
        # Compute indices
        index.delta.perm <- grep("perm", names(results.perm))
        index.delta.obs <- which(names(results.perm)=="delta")

        # Compute p-value
        pval <- apply(results.perm[, ], 1, function(x) {
                
                    #x <- as.numeric(results.perm[2, ]) # Test
                                
                    # Retrieve delta observed, perm
                    delta.obs <- x[index.delta.obs]
                    delta.perm <- x[index.delta.perm]
                    
                    # Compute pval (1-sided)
                    #if(delta.obs < 0) {
                        
                        #pval <- sum(delta.perm < delta.obs) / n.iterations
                        
                    #} else if(delta.obs >= 0) {
                        
                        #pval <- sum(delta.perm > delta.obs) / n.iterations
                        
                    #}
                    
                    # Compute pval (2-sided)
                    pval <- sum(abs(delta.perm) > abs(delta.obs)) / n.iterations
                        
                    return(pval)
                })
    
        results.obs$pval <- pval
        
        # Check results
        #results.small <- results.obs[which(results.obs$pval < 0.05 & results.obs$delta > 0), ]
        #results.small <- results.obs[which(results.obs$pval < 0.05 & results.obs$delta < 0), ]
        #table(results.obs$pval < 0.05)
        
    # Re-order by pval
    results.obs <- results.obs[order(results.obs$pval), ]
    
    ################################################################

    # Compute mean norm gene expression
    #. <- apply(df.gene.norm, 1, function(x) {mean(log2(x + 1))})
    #. <- data.frame("gene_short_name"=names(.),
                    #"mean.expr.gene.norm.g1.g2"=as.numeric(.),
                    #stringsAsFactors=FALSE
                    #)
    
    # Annotate
    results.obs <- join(results.obs, mean.combined.df, by="gene_short_name", type="left")
    
    ################################################################
 
    # Save into new slot
    MarvelObject$sample.metadata <- sample.metadata
    MarvelObject$DE$SJ$Table <- results.obs
    MarvelObject$DE$SJ$cell.group.g1 <- cell.group.g1
    MarvelObject$DE$SJ$cell.group.g2 <- cell.group.g2
    
    # Return final object
    return(MarvelObject)
            
}

######################################
# 统计gene表达量，并做显著性分析
######################################

CompareValues.Genes.10x.ZCY <- function(MarvelObject, log2.transform=TRUE, show.progress=TRUE, method="wilcox", mast.method="bayesglm", mast.ebayes=TRUE) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    de <- MarvelObject$DE$SJ$Table
    df.gene.norm <- MarvelObject$gene.norm.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.group.g1 <- MarvelObject$DE$SJ$cell.group.g1
    cell.group.g2 <- MarvelObject$DE$SJ$cell.group.g2
    
    # Example arguments
    #MarvelObject <- marvel
    #de <- MarvelObject$DE$SJ$Table
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #cell.group.g1 <- MarvelObject$DE$SJ$cell.group.g1
    #cell.group.g2 <- MarvelObject$DE$SJ$cell.group.g2
    #log2.transform <- TRUE
    #method <- "wilcox"
    #mast.method <- "bayesglm"
    #mast.ebayes <- TRUE
    
    #################################################################
    
    # Subset genes of SJ included for DE
    # 从差异可变剪切的结果中提取基因名字，首先要保证我们分析的基因是差异可变剪切，再去考虑是不是差异表达
    gene_short_names <- unique(de$gene_short_name)
    
    # Subset cell group 1
    df.gene.norm.g1 <- df.gene.norm[gene_short_names, cell.group.g1]
    
    # Subset cell group 2
    df.gene.norm.g2 <- df.gene.norm[gene_short_names, cell.group.g2]
    
    # log2 transform values
    if(log2.transform==TRUE | method=="mast") {
        df.gene.norm.g1 <- log2(df.gene.norm.g1 + 1)
        df.gene.norm.g2 <- log2(df.gene.norm.g2 + 1)
    }
    
    # Tabulate expression: Group 1
        # Compute group size
        n.cells.total <- length(cell.group.g1)
        
        # Compute % expr cells
        n.cells.expr.gene.norm <- apply(df.gene.norm.g1, 1, function(x) {sum(x != 0)})
        pct.cells.expr.gene.norm <- round(n.cells.expr.gene.norm/n.cells.total * 100, digits=2)
        
        # Compute mean expr
        mean.expr.gene.norm <- apply(df.gene.norm.g1, 1, function(x) {mean(x)})
        
        # Tabulate results
        results <- data.frame("gene_short_name"=gene_short_names,
                              "n.cells.total.norm"=n.cells.total,
                              "n.cells.expr.gene.norm"=n.cells.expr.gene.norm,
                              "pct.cells.expr.gene.norm"=pct.cells.expr.gene.norm,
                              mean.expr.gene.norm=mean.expr.gene.norm,
                              stringsAsFactors=FALSE
                              )
        row.names(results) <- NULL
        
        # Indicate group
        names(results)[which(names(results) != "gene_short_name")] <- paste(names(results)[which(names(results) != "gene_short_name")], ".g1", sep="")
        
        # Save as new object
        results.g1 <- results
        
    # Tabulate expression: Group 2
        # Compute group size
        n.cells.total <- length(cell.group.g2)
        
        # Compute % expr cells
        n.cells.expr.gene.norm <- apply(df.gene.norm.g2, 1, function(x) {sum(x != 0)})
        pct.cells.expr.gene.norm <- round(n.cells.expr.gene.norm/n.cells.total * 100, digits=2)
        
        # Compute mean expr
        mean.expr.gene.norm <- apply(df.gene.norm.g2, 1, function(x) {mean(x)})
        
        # Tabulate results
        results <- data.frame("gene_short_name"=gene_short_names,
                              "n.cells.total.norm"=n.cells.total,
                              "n.cells.expr.gene.norm"=n.cells.expr.gene.norm,
                              "pct.cells.expr.gene.norm"=pct.cells.expr.gene.norm,
                              mean.expr.gene.norm=mean.expr.gene.norm,
                              stringsAsFactors=FALSE
                              )
        row.names(results) <- NULL
        
        # Indicate group
        names(results)[which(names(results) != "gene_short_name")] <- paste(names(results)[which(names(results) != "gene_short_name")], ".g2", sep="")
        
        # Save as new object
        results.g2 <- results
    
    # Merge group 1, 2
    index.l <- table(results.g1$gene_short_name==results.g2$gene_short_name)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
        results.g2$gene_short_name <- NULL
        results <- cbind.data.frame(results.g1, results.g2)
    } else {
        return(message("Error in merging tables from Group 1 and 2"))
    }
    
    #################################################################
    
    if(method=="wilcox"){
        
        # Compute log2fc
        results$log2FoldChange.gene <- results$mean.expr.gene.norm.g2 - results$mean.expr.gene.norm.g1
        
        # Compute p-values
        message("Performing Wilcox rank sum test...")
        
        if(show.progress==TRUE) {
            pb <- txtProgressBar(1, length(gene_short_names), style=3)
        }

        pval <- NULL
        
        # 每个基因单独检验
        for(i in 1:length(gene_short_names)) {
            
            # Retrieve values
            values.g1 <- as.numeric(df.gene.norm.g1[gene_short_names[i], ])
            values.g2 <- as.numeric(df.gene.norm.g2[gene_short_names[i], ])
            
            # Wilcox
            pval[i] <- wilcox.test(values.g1, values.g2)$p.value
            
            if(show.progress==TRUE) {
                # Track progress
                setTxtProgressBar(pb, i)
            }
        }
        
        results$p-value.gene <- pval
        
    } else if(method=="mast"){
        
        # Prepare cdata
            # Indicate group 1,2
            cdata <- data.frame("cell.id"=c(cell.group.g1, cell.group.g2))
            cdata$condition <- ifelse(cdata$cell.id %in% cell.group.g1, "g1", "g2")
            cdata$condition <- factor(cdata$condition, levels=c("g1", "g2"))
            
            # Retrieve expression values
            df.exp.master <- cbind(df.gene.norm.g1, df.gene.norm.g2)
            df.exp.master <- df.exp.master[,cdata$cell.id]
            
            # Compute n expressed genes
            . <- apply(df.exp.master, 2, function(x) {sum(x != 0)})
            . <- data.frame("cell.id"=names(.),
                            "n.genes"=as.numeric(.),
                            stringsAsFactors=FALSE
                            )
                            
            # Compute gene detection rate
            .$cngeneson <- scale(.$n.genes)
            cdata <- join(cdata, ., by="cell.id", type="left")
            
            # Format for MAST
            names(cdata)[which(names(cdata)=="cell.id")] <- "wellKey"
            row.names(cdata) <- cdata$wellKey
            
        # Prepare fdata
        fdata <- data.frame("primerid"=row.names(df.exp.master),
                            stringsAsFactors=FALSE
                            )
        row.names(fdata) <- fdata$primerid
        
        # Prepare SingleCellAssay object
        sca <- MAST::FromMatrix(as.matrix(df.exp.master), cdata, fdata)
        
        # Build model
        zlmCond <- MAST::zlm(~condition + cngeneson, sca, method=mast.method, ebayes=mast.ebayes, silent=TRUE)

        # Only test the condition coefficient
        summaryCond <- summary(zlmCond, doLRT='conditiong2')

        # Retrieve DE table
        summaryDt <- summaryCond$datatable
        #fcHurdle <- merge(summaryDt[contrast=='conditiong2' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditiong2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
        #fcHurdle <- as.data.frame(fcHurdle)
        
        index <- which(summaryDt$contrast=='conditiong2' & summaryDt$component=='H')
        summaryDt.small.1 <- summaryDt[index, c("primerid", "Pr(>Chisq)")]
        
        index <- which(summaryDt$contrast=='conditiong2' & summaryDt$component=='logFC')
        summaryDt.small.2 <- summaryDt[index, c("primerid", "coef", "ci.hi", "ci.lo")]
        
        fcHurdle <- join(summaryDt.small.1, summaryDt.small.2, by="primerid")
        
        fcHurdle <- as.data.frame(fcHurdle)
        
        # Subset relevant columns
        results. <- fcHurdle
        results. <- results.[,c("primerid", "Pr(>Chisq)", "coef")]
        names(results.) <- c("gene_short_name", "p-value.gene", "log2FoldChange.gene")
        results. <- results.[,c("gene_short_name", "log2FoldChange.gene", "p-value.gene")]
        
        # Annotate original result table
        results <- join(results, results., by="gene_short_name", type="left")
        
    }
    
    #################################################################

    # Adjust for multiple testing
    results$q-value.gene <- p.adjust(results$p-value.gene, method="fdr", n=length(results$p-value.gene))
    
    ################################################################
 
    # Update SJ DE table
    de <- join(de, results, by="gene_short_name", type="left")
    
    ################################################################
    # Save into new slot
    MarvelObject$DE$SJ$Table <- de
    
    # Return final object
    return(MarvelObject)
            
}

######################################
# 给差异表达的junction和基因定性
######################################

IsoSwitch.10x.ZCY <- function(MarvelObject, pval.sj=0.05, log2fc.sj=NULL, delta.sj=NULL, min.gene.norm=-100, min.pct.cells.sj=0.01, min.pct.cells.gene=0.01, pval.gene=0.05, log2fc.gene=NULL, pval.adj.gene = NULL) {

    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$DE$SJ$Table
    pval.sj <- pval.sj
    log2fc.sj <- log2fc.sj
    delta.sj <- delta.sj
    pval.adj.gene <- pval.adj.gene
    pval.gene <- pval.gene
    log2fc.gene <- log2fc.gene
    min.gene.norm <- min.gene.norm
    min.pct.cells.gene <-  min.pct.cells.gene
    min.pct.cells.sj <- min.pct.cells.sj
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$SJ$Table
    #pval.sj <- 0.05
    #log2fc.sj <- NULL
    #delta.sj <- 5
    # pval.adj.gene <- NULL
    #pval.gene <- 0.05
    #log2fc.gene <- 0.5
    #min.gene.norm <- -100
    
    ######################### CLASSIFY RELATIONSHIP ###########################
    
    # 首先提取显著上下调的SJ，并将其提取出来
    # Indicate sig events and direction: SJ


    # Subset/plot expressed genes
    df <- df[which(df$mean.expr.gene.norm.g1.g2 > min.gene.norm), ] # mean.expr.gene.norm.g1.g2 两个样本平均log(表达量)
    df <- df[which(df$pct.cells.expr.gene.g1 > min.pct.cells.gene), ]
    df <- df[which(df$pct.cells.expr.gene.g2 > min.pct.cells.gene), ]
    df <- df[which(df$pct.cells.expr.sj.g1 > min.pct.cells.sj), ]
    df <- df[which(df$pct.cells.expr.sj.g2 > min.pct.cells.sj), ]


    if(!is.null(log2fc.sj)) {
        
        df$sig <- NA
        df$sig[which(df$pval < pval.sj & df$log2fc > log2fc.sj & df$mean.expr.gene.norm.g1.g2)] <- "Up"
        df$sig[which(df$pval < pval.sj & df$log2fc < (log2fc.sj*-1) & df$mean.expr.gene.norm.g1.g2)] <- "Down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("Up", "Down", "n.s."))
        table(df$sig)
        
    } else if(!is.null(delta.sj)){
        
        df$sig <- NA
        df$sig[which(df$pval < pval.sj & df$delta > delta.sj & df$mean.expr.gene.norm.g1.g2)] <- "Up"
        df$sig[which(df$pval < pval.sj & df$delta < (delta.sj*-1) & df$mean.expr.gene.norm.g1.g2)] <- "Down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("Up", "Down", "n.s."))
        table(df$sig)
        
    }
    
    names(df)[which(names(df)=="sig")] <- "sig.sj"
    
    # Subset sig SJ
    df <- df[which(df$sig.sj %in% c("Up", "Down")), ]
    
    # 在SJ显著差异的前提下，再去看基因的差异【基因可以差异，也可以不差异】
    # Indicate sig events and direction: Gene
    if( !is.null(pval.adj.gene) ){
        df$sig <- NA
        df$sig[which(df$pval.adj.gene.norm < pval.adj.gene & df$log2FoldChange.gene > log2fc.gene)] <- "Up"
        df$sig[which(df$pval.adj.gene.norm < pval.adj.gene & df$log2FoldChange.gene < (log2fc.gene*-1))] <- "Down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("Up", "Down", "n.s."))
        table(df$sig)
    } else if( !is.null(pval.gene) ){
        df$sig <- NA
        df$sig[which(df$pval.gene.norm < pval.gene & df$log2FoldChange.gene > log2fc.gene)] <- "Up"
        df$sig[which(df$pval.gene.norm < pval.gene & df$log2FoldChange.gene < (log2fc.gene*-1))] <- "Down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("Up", "Down", "n.s."))
        table(df$sig)      
    }

    
    names(df)[which(names(df)=="sig")] <- "sig.gene"
    
    # Classify SJ-gene relationship
    table(df$sig.sj, df$sig.gene)
    
    df$cor <- NA
    
    df$cor[which(df$sig.sj=="Up" & df$sig.gene=="Up")] <- "Coordinated"
    df$cor[which(df$sig.sj=="Down" & df$sig.gene=="Down")] <- "Coordinated"
    df$cor[which(df$sig.sj=="Down" & df$sig.gene=="Up")] <- "Opposing"
    df$cor[which(df$sig.sj=="Up" & df$sig.gene=="Down")] <- "Opposing"
    df$cor[which(df$sig.sj=="Up" & df$sig.gene=="n.s.")] <- "Iso-Switch"
    df$cor[which(df$sig.sj=="Down" & df$sig.gene=="n.s.")] <- "Iso-Switch"
    
    sum(is.na(df$cor))
    
    # Find mixed relationship
    gene_short_names <- unique(df$gene_short_name)

    ###################
    corr <- NULL
    for(i in 1:length(gene_short_names)) {
        # Subset gene
        df.small <- df[which(df$gene_short_name == gene_short_names[i]), ]
        # Stratify cor
        test.unique <- unique(df.small$cor)
        if(length(test.unique) != 1) { 
            corr[i] <- "Complex"
        } else {
            corr[i] <- as.vector(test.unique)
        }
    }
    
    # Update SJ-gene cor column
    results <- data.frame("gene_short_name"=gene_short_names, "cor.complete"=corr, stringsAsFactors=FALSE)
    df <- join(df, results, by="gene_short_name", type="left")
    ###################

    # 只要junction差异上下调的方向相反，就定义为 Complex
    corr <- NULL
    corr <- hash::hash()
    for (gene in gene_short_names){
        df.small <- df[which(df$gene_short_name == gene), ]
        test.unique <- unique(df.small$sig.sj)
        if(length(test.unique) != 1){
            corr[[gene]] <- gene
        }
    }

    for (gene in hash::keys(corr)) {
        df[which(df$gene_short_name == gene),]$cor.complete <- "Complex"
    }

    corr <- NULL

    ############################ PLOT PROPORTION ###########################
    
    # Set factor levels
    levels <- intersect(c("Coordinated", "Opposing", "Iso-Switch", "Complex"), unique(df$cor.complete))

    # Tabulate freq (for genes, not splicing)
    #. <- as.data.frame(table(df$cor.complete), stringsAsFactors=FALSE)
    df.small <- unique(df[, c("gene_short_name", "cor.complete")])
    . <- as.data.frame(table(df.small$cor.complete), stringsAsFactors=FALSE)
    names(.) <- c("sj.gene.cor", "freq")
    .$pct <- .$freq / sum(.$freq) * 100
    
    # Set factor levels
    .$sj.gene.cor <- factor(.$sj.gene.cor, levels=levels)
    . <- .[order(.$sj.gene.cor), ]
    
    # Compute statistics for plot
    .$fraction <- .$freq / sum(.$freq)
    .$ymax <- cumsum(.$fraction)
    .$ymin = c(0, .$ymax[-length(.$ymax)])
    
    # Definitions
    data <- .
    xmax <- nrow(data) + 1
    xmin <- nrow(data)
    ymax <- data$ymax
    ymin <- data$ymin
    z <- data$sj.gene.cor
    maintitle <- ""
    xtitle <- ""
    ytitle <- ""
    legendtitle <- "Gene-SJ Relationship"
    
    # Plot
    plot <- ggplot() +
        geom_rect(data=data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=z), color="black") +
        coord_polar(theta="y") +
        xlim(c(2, 4)) +
        #scale_fill_manual(values=colors) +
        labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            plot.title=element_text(hjust = 0.5, size=15),
            plot.subtitle=element_text(hjust = 0.5, size=15),
            axis.line = element_blank(),
            axis.ticks=element_blank(),
            axis.text=element_blank(),
            legend.title=element_text(size=9),
            legend.text=element_text(size=9)
            )
            
    #########################################################################
    
    # Save into new slot
    # df <- select(df, -c("cor"))
    MarvelObject$SJ.Gene.Cor$Data <- df
    MarvelObject$SJ.Gene.Cor$Proportion$Plot <- plot
    MarvelObject$SJ.Gene.Cor$Proportion$Table <- data[,c("sj.gene.cor", "freq", "pct")]
    
    # Return final object
    return(MarvelObject)

}

######################################
# 筛选含有可变剪切的候选基因
######################################

#【后续】还需要添加最低表达量参数

Get_Top_DE_SJ_Gene <- function(MarvelObject, DE_type, top_DE_junctions, log2fc.sj=NULL, delta=NULL){
    
    MarvelObject <- MarvelObject
    DE_type <- DE_type
    top_DE_junctions <- top_DE_junctions
    log2fc.sj <- log2fc.sj
    delta <- delta

    #########################

    filter_duplicate_when_get_list <- function(x, top_DE_junctions){
        # 只提取前 min(length(unique(x)), top_DE_junctions) 个基因
        x <- x 
        top_DE_junctions <- top_DE_junctions
        my_vector <- x
        my_hash <- hash()
        i <- 0
        for (v in my_vector) {
            if (v %in% keys(my_hash)){
                next()
            } else {
                i = i + 1
                my_hash[[v]] <- i
            }
            if(length(keys(my_hash)) == top_DE_junctions){
                break()
            }
        }
        return(keys(my_hash))
    }


    top_gene_list <- function(df.small){
        df.small <- df.small
        # 需要预先判读方向，如果只有低表达或者只有高表达，只返回单边基因
        if(length(df.small$log2fc) == 0){
            return(c())
        } else {
            if(max(df.small$log2fc) < 0){
                down_gene_list <- filter_duplicate_when_get_list(x=rev(df.small$gene_short_name), top_DE_junctions=top_DE_junctions)
                return(down_gene_list)
            } else if(min(df.small$log2fc) > 0){
                up_gene_list <- filter_duplicate_when_get_list(x=df.small$gene_short_name, top_DE_junctions=top_DE_junctions)
                return(up_gene_list)
            } else {
                df.small.negtivate <- df.small[df.small$log2fc < 0,]
                down_gene_list <- filter_duplicate_when_get_list(x=rev(df.small.negtivate$gene_short_name), top_DE_junctions=top_DE_junctions)
                df.small.positive <- df.small[df.small$log2fc > 0,]
                up_gene_list <- filter_duplicate_when_get_list(x=df.small.positive$gene_short_name, top_DE_junctions=top_DE_junctions)
                gene_list <- c(down_gene_list, up_gene_list)
                return(gene_list)
            }
        }
    }

    # 将一个基因含有多个junction的提取出来，暂时不绘制只有一个junction的基因
    df <- MarvelObject$SJ.Gene.Cor$Data
    df.small <- df[df$cor.complete == DE_type,]
    multi_sj_one_gene <- names(table(df.small$gene_short_name)[(table(df.small$gene_short_name) > 1)])
    df.small <- df.small[df.small$gene_short_name %in% multi_sj_one_gene,]

    if (!is.null(log2fc.sj) & is.null(delta)) {

        df.small <- df.small[order(df.small$log2fc, decreasing = T), ]
        gene_list <- top_gene_list(df.small=df.small)
        return(gene_list)

    } else if(is.null(log2fc.sj) & !is.null(delta)){

        df.small <- df.small[order(df.small$delta, decreasing = T), ]
        gene_list <- top_gene_list(df.small=df.small)
        return(gene_list)

    } else {

        message("您没有在函数 Get_Top_DE_SJ_Gene 的使用中正确设置可变剪切的筛选标准")
        quit(save = "no", status = 1, runLast = TRUE)

    }
}


get_candidate_genes <- function(MarvelObject, top_DE_junctions, delta=NULL, log2fc=NULL){

  MarvelObject <- MarvelObject
  candidate_table <- MarvelObject$SJ.Gene.Cor$Data
  top_DE_junctions <- top_DE_junctions
  delta <- delta
  log2fc <- log2fc

  ####
  DE_genes_hash <- hash()
  DE_types <- unique(candidate_table$cor.complete)

  for( i in c(1:length(DE_types))){

    DE_type=DE_types[i]

    if(DE_type == "Complex"){

      .set(DE_genes_hash, keys = DE_type, values = unique(candidate_table[candidate_table$cor.complete == "Complex",]$gene_short_name))
    
    } else {
        
        if (is.null(log2fc) && !is.null(delta)){
            top_SJ_gene_list <- Get_Top_DE_SJ_Gene(MarvelObject=MarvelObject, 
                                                DE_type=DE_type, 
                                                top_DE_junctions=top_DE_junctions, 
                                                delta=delta
                                                )
        }

        if (!is.null(log2fc) && is.null(delta)){
            top_SJ_gene_list <- Get_Top_DE_SJ_Gene(MarvelObject=MarvelObject, 
                                                DE_type=DE_type, 
                                                top_DE_junctions=top_DE_junctions, 
                                                log2fc=log2fc
                                                )
        }

      .set(DE_genes_hash, keys = DE_type, values = top_SJ_gene_list )

    }
  }

  return(DE_genes_hash)
}



######################################
# 绘制目的基因gene表达量点图
######################################

adhocGene.TabulateExpression.Gene.10x.ZCY <- function(MarvelObject, gene_short_name, diff_type, g1, g2, log2.transform=TRUE, min.pct.cells=10, downsample=FALSE, seed=1) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    df.gene.norm <- MarvelObject$gene.norm.matrix
    sample.metadata <- MarvelObject$sample.metadata
    gene.metadata <- MarvelObject$gene.metadata
    gene_short_name <- gene_short_name
    diff_type <- diff_type
    g1 <- g1
    g2 <- g2
    min.pct.cells <- min.pct.cells
    downsample <- downsample
    seed <- seed
    log2.transform <- log2.transform

    # Example arguments
    #MarvelObject <- marvel
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #gene_short_name <- "TPM1"
    #min.pct.cells <- 10
    #downsample <- TRUE
    #seed <- 1
    
    ##########################################################################
    
    sample.metadata.list <- subset_cell_according_diff_type(df=sample.metadata, 
                                                        diff_type=diff_type,
                                                        g1=g1,
                                                        g2=g2
                                                        )
    sample.metadata <- sample.metadata.list[["data"]]
    cell.group.g1 <- sample.metadata.list[["cells.g1"]]
    cell.group.g2 <- sample.metadata.list[["cells.g2"]]
    g1 <- sample.metadata.list[["g1"]]
    g2 <- sample.metadata.list[["g2"]]

    cell.group.list <- list("control" = cell.group.g1, "case" = cell.group.g2)

    names(cell.group.list) <- c(names(g1), names(g2))

    # Downsample
    .list <- list()
    
    set.seed(seed)
    
    if(downsample==TRUE){
        
        # Find lowest common demoninator
        n.cells <- min(sapply(cell.group.list, length))
        
        # Track progress
        message(paste("Downsampling to ", n.cells, " cells per group", sep=""))
        
        # Downsample
        for(i in 1:length(cell.group.list)) {
            
            cell.ids <- cell.group.list[[i]]
            
            cell.ids <- sample(cell.ids, size=n.cells, replace=FALSE)
            
            .list[[i]] <- cell.ids
            
        }
        
        names(.list) <- names(cell.group.list)
        cell.group.list <- .list
            
    }

    # Subset relevant gene
    df.gene.norm <- df.gene.norm[gene_short_name, ]
    df.gene.norm <- data.frame("cell.id"=names(df.gene.norm),
                               "exp.gene.norm"=as.numeric(df.gene.norm),
                               stringsAsFactors=FALSE
                               )
    row.names(df.gene.norm) <- NULL
    
    # Tranform values
    if(log2.transform==TRUE) {
        
        df.gene.norm$exp.gene.norm <- log2(df.gene.norm$exp.gene.norm + 1)
    
    }
    
    # Annotate cell groups
        # Create ref table
        .list <- list()
        
        for(i in 1:length(cell.group.list)) {
            
            . <- data.frame("cell.id"=cell.group.list[[i]],
                            "group"=names(cell.group.list)[i]
                            )
                            
            .list[[i]] <- .
            
        }
        
        ref <- do.call(rbind.data.frame, .list)
        
        # Annotate
        df.gene.norm <- join(df.gene.norm, ref, by="cell.id", type="left")
        df.gene.norm <- df.gene.norm[!is.na(df.gene.norm$group), ]
        
        # Set factor levels
        df.gene.norm$group <- factor(df.gene.norm$group, levels=names(cell.group.list))
        df.gene.norm <- df.gene.norm[order(df.gene.norm$group), ]
        
    # Compute mean expression, % cells expressed
    groups <- levels(df.gene.norm$group)
    
    mean.expr <- NULL
    pct.cells.expr <- NULL
    
    for(i in 1:length(groups)) {
        
        # Define cell group
        group <- groups[i]
        
        # Subset
        df.gene.norm.small <- df.gene.norm[which(df.gene.norm$group==group), ]
        
        # Compute mean
        mean.expr[i] <- mean(df.gene.norm.small$exp.gene.norm)
        
        # Compute % cells expressing gene
        pct.cells.expr[i] <- sum(df.gene.norm.small$exp.gene.norm != 0) / nrow(df.gene.norm.small) * 100
        
    }
    
    results <- data.frame("cell.group"=groups,
                          "mean.expr"=mean.expr,
                          "pct.cells.expr"=pct.cells.expr,
                          stringsAsFactors=FALSE
                          )
    
    # Censor lowly expressing cell types
    results$pct.cells.expr[which(results$pct.cells.expr < min.pct.cells)] <- NA
     
    # Set factor levels
    results$cell.group <- factor(results$cell.group, levels=names(cell.group.list))
    results <- results[order(results$cell.group), ]
    results$gene_name <- gene_short_name
    results$gene_id <- gene.metadata[gene.metadata$gene_short_name == gene_short_name,]$gene_id
    
    # Bubble plot
        # Definition
        data <- results
        x <- data$gene_name
        y <- data$cell.group
        z1 <- data$pct.cells.expr
        z2 <- data$mean.expr
        maintitle <- gene_short_name
        xtitle <- ""
        ytitle <- ""
        legendtitle.size <- "% cells"
        legendtitle.color <- "mean[log2(expr)]"
        labels.y <- data$cell.group
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, size=z1, color=z2)) +
            scale_color_gradientn(colors=c("gray","cyan","green","yellow","red")) +
            scale_y_discrete(limits = rev(levels(y))) +
            labs(x=xtitle, y=ytitle, size=legendtitle.size, color=legendtitle.color) +
            ggtitle("Gene") +
            theme(panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    plot.title = element_text(size=12, hjust=0.5),
                    axis.line.x=element_line(colour = "black"),
                    axis.line.y=element_line(colour = "black"),
                    axis.text.x=element_text(size=10, colour="black", angle = 45, hjust = 1),
                    axis.text.y=element_text(size=10, colour="black"),
                    legend.title=element_text(size=10),
                    legend.text=element_text(size=8),
                    legend.key=element_blank()
                    )  +
                scale_size_area(breaks=c(0.1,1,25,50,75,100), limits=c(0, 100)) +
                guides(size = guide_legend(order = 2)) 
 
    ##################################################################
 
    # Save into new slots
    MarvelObject$adhocGene$Expression$Gene$Table <- results
    MarvelObject$adhocGene$Expression$Gene$Plot <- plot
    MarvelObject$adhocGene$cell.group.list <- cell.group.list
    MarvelObject$adhocGene$gene_short_name <- gene_short_name
    
    # Return final object
    return(MarvelObject)
            
}

######################################
# 绘制目的基因junction表达量点图
######################################

adhocGene.TabulateExpression.PSI.10x.ZCY <- function(MarvelObject, diff_exp_junction_only=TRUE, min.pct.cells=NULL ) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    diff_exp_junction_only <- diff_exp_junction_only
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.group.list <- MarvelObject$adhocGene$cell.group.list
    gene_short_name <- MarvelObject$adhocGene$gene_short_name
    min.pct.cells <- min.pct.cells
    
    # # Example arguments
    # MarvelObject <- marvel
    # sj.metadata <- MarvelObject$sj.metadata
    # df.gene.count <- MarvelObject$gene.count.matrix
    # df.sj.count <- MarvelObject$sj.count.matrix
    # cell.group.list <- MarvelObject$adhocGene$cell.group.list
    # gene_short_name <- MarvelObject$adhocGene$gene_short_name
    # min.pct.cells <- 0.01
    # target_junction_in_tabulate <- TRUE
    
    ##########################################################################
    ########################### TABULATE SJ COUNTS ###########################
    ##########################################################################
    
    # Define SJ
    # 按照基因名从原始数据中提取junction
    # 注意这里的问题，我们通过统计junction的差异表达找到了候选的junction，是否有必要仅绘制候选junction的图？
    if (diff_exp_junction_only == TRUE){
        DE_table <- MarvelObject$SJ.Gene.Cor$Data
        coord.introns <- DE_table[DE_table$gene_short_name == gene_short_name,]$coord.intron
    } else {
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start==gene_short_name), "coord.intron"]
    }
    
    
    # Subset relevant SJ
    df.sj.count <- df.sj.count[coord.introns, ]
    
    # Compute counts, % expr for each cell group
    .list <- list()
    
    for(i in 1:length(cell.group.list)) {
        
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.list[[i]]]
        
        # Compute counts
        sj.count.total <- apply(df.sj.count.small, 1, function(x) {sum(x)})
        
        # % cells expressing SJ
        n.cells.total <- ncol(df.sj.count.small)
        n.cells.expr.sj <- apply(df.sj.count.small, 1, function(x) {sum(x != 0)})
        pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
        
        # Save into data frame
        results <- data.frame("cell.group"=names(cell.group.list)[i],
                              "coord.intron"=row.names(df.sj.count.small),
                              "n.cells.total"=n.cells.total,
                              "n.cells.expr.sj"=n.cells.expr.sj,
                              "pct.cells.expr.sj"=pct.cells.expr.sj,
                              "sj.count.total"=sj.count.total,
                              stringsAsFactors=FALSE
                              )
                              
        row.names(results) <- NULL
        
        # Save into list
        .list[[i]] <- results
        
    }
    
    results.sj.count <- do.call(rbind.data.frame, .list)
    
    ##########################################################################
    ########################## TABULATE GENE COUNTS ##########################
    ##########################################################################
        
    # Subset relevant gene
    df.gene.count <- df.gene.count[gene_short_name, , drop=FALSE]
    
    # Compute counts, % expr for each cell group
    .list <- list()
    
    for(i in 1:length(cell.group.list)) {
        
        # Subset cells
        df.gene.count.small <- df.gene.count[, cell.group.list[[i]], drop=FALSE]
        
        # Compute counts
        gene.count.total <- apply(df.gene.count.small, 1, function(x) {sum(x)})
        
        # % cells expressing SJ
        n.cells.total <- ncol(df.gene.count.small)
        n.cells.expr.gene <- apply(df.gene.count.small, 1, function(x) {sum(x != 0)})
        pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
        
        # Save into data frame
        results <- data.frame("cell.group"=names(cell.group.list)[i],
                              "gene_short_name"=row.names(df.gene.count.small),
                              "n.cells.total"=n.cells.total,
                              "n.cells.expr.gene"=n.cells.expr.gene,
                              "pct.cells.expr.gene"=pct.cells.expr.gene,
                              "gene.count.total"=gene.count.total,
                              stringsAsFactors=FALSE
                              )
                              
        row.names(results) <- NULL
        
        # Save into list
        .list[[i]] <- results
        
    }
    
    results.gene.count <- do.call(rbind.data.frame, .list)
    
    ##########################################################################
    ################################# MERGE ##################################
    ##########################################################################
    
    # Annotate gene count
    results <- join(results.sj.count, results.gene.count[,c("cell.group", "gene.count.total")], by="cell.group", type="left")
    
    # Compute PSI
    results$psi <- round(results$sj.count.total/results$gene.count.total * 100, digits=2)
    
    # Remove lowly expressed SJ
    results$psi[which(results$pct.cells.expr.sj < min.pct.cells)] <- NA
    results <- results[!is.na(results$psi), ]
    
    # Re-order by % expressing cells
        # Reorder by % expressing cells
        #. <- aggregate(results$pct.cells.expr.sj, list(results$coord.intron), function(x) {sum(!is.na(x))})
        . <- aggregate(results$n.cells.expr.sj, list(results$coord.intron), function(x) {sum(x)})
        . <- .[order(.[,2], decreasing=TRUE), ]
        results$coord.intron <- factor(results$coord.intron, levels=.[,1])
                        
        # Reorder by cell group
        results$cell.group <- factor(results$cell.group, levels=names(cell.group.list))
        results <- results[order(results$cell.group, results$coord.intron), ]
        
        # Annotate column no.
        coord.introns <- levels(results$coord.intron)
        . <- data.frame("coord.intron"=coord.introns,
                        "figure.column"=paste("SJ-", c(1:length(coord.introns)), sep=""),
                        stringsAsFactors=FALSE
                        )
        results <- join(results, ., by="coord.intron", type="left")
        
        # Reorder columns
        cols.1 <- c("cell.group", "figure.column", "coord.intron")
        cols.2 <- setdiff(names(results), cols.1)
        results <- results[, c(cols.1, cols.2)]
            
    # Bubble plot
        # Definition
        data <- results
        x <- data$coord.intron
        y <- data$cell.group
        z1 <- data$pct.cells.expr.sj
        z2 <- data$psi
        maintitle <- gene_short_name
        xtitle <- ""
        ytitle <- ""
        legendtitle.size <- "% cells"
        legendtitle.color <- "PSI"
        labels.y <- data$cell.group
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, size=z1, color=z2)) +
            scale_color_gradientn(colors=c("gray","cyan","green","yellow","red")) +
            scale_y_discrete(limits = rev(levels(y))) +
            labs( x=xtitle, y=ytitle, size=legendtitle.size, color=legendtitle.color) +
            ggtitle("Junctions") +
            theme(panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_rect(colour = "white", fill = "white"),
                    plot.title = element_text(size=12, hjust=0.5),
                    axis.line.x=element_line(colour = "black"),
                    axis.line.y=element_line(colour = "black"),
                    axis.text.x=element_text(size = 10, colour="black", angle = 45, hjust = 1),
                    axis.text.y=element_text(size=10, colour="black"),
                    legend.title=element_text(size=10),
                    legend.text=element_text(size=8),
                    legend.key=element_blank(),
                    legend.background = element_rect(colour = "white", fill = "white")
            ) + 
            scale_size_area(breaks=c(0.1,1,25,50,75,100), limits=c(0, 100)) +
            guides(size = guide_legend(order = 2)) 

    ##################################################################
    
    # Save into new slots
    MarvelObject$adhocGene$Expression$PSI$Table <- results
    MarvelObject$adhocGene$Expression$PSI$Plot <- plot
    
    # Return final object
    return(MarvelObject)
            
}

######################################
# 绘制目的基因和junction结构图
######################################

adhocGene.PlotSJPosition.10x.ZCY <- function(MarvelObject, gene_short_name, diff_exp_junction_only,  coord.intron.ext=50, rescale_introns=FALSE, show.protein.coding.only=TRUE) {
  
  # Define arguments
  MarvelObject <- MarvelObject
  gene_short_name <- gene_short_name
  diff_exp_junction_only <- diff_exp_junction_only
  sj.metadata <- MarvelObject$sj.metadata
  gtf <- MarvelObject$gtf
  gene.metadata <- MarvelObject$gene.metadata
  coord.intron.ext <- coord.intron.ext
  rescale_introns <- rescale_introns
  show.protein.coding.only <- show.protein.coding.only
  
  
    #   # Example arguments
    #   MarvelObject <- marvel
    #   gene_short_name <- gene
    #   sj.metadata <- MarvelObject$sj.metadata
    #   gtf <- MarvelObject$gtf
    #   coord.intron.ext <- 50
    #   rescale_introns <- TRUE
    #   show.protein.coding.only <- TRUE
  
  #################################################################
  #################### RETRIEVE TRANSCRIPTS #######################
  #################################################################
  
  if (diff_exp_junction_only == TRUE){
    DE_table <- MarvelObject$SJ.Gene.Cor$Data
    coord.intron <- DE_table[DE_table$gene_short_name == gene_short_name,]$coord.intron
  } else {
    coord.intron <- sj.metadata[which(sj.metadata$gene_short_name.start==gene_short_name), "coord.intron"]
  }
  

  message(paste0("基因 ", gene_short_name, " 具有 ", length(coord.intron), " 个junction，分别为：", paste(coord.intron, collapse = ", ") ))
  
  
  # Subset (by approximation) releavnt gene
  # 用gene_id提取信息，这样每个基因的内容都是唯一的
  gene_id <- gene.metadata[gene.metadata$gene_short_name == gene_short_name,]$gene_id

  gtf <- gtf[grep(gene_id, gtf$V9), ]
  gtf <- parsing_gtf(gtf=gtf, gene.metadata=gene.metadata, gene_or_exon=NULL)

  gtf <- gtf[which(gtf$gene_id==gene_id), ]
  
  # Report progress
  message("Retrieving transcripts from GTF file...")
  
  
  # Annotate transcripts with ORF status
  anno <- unique(gtf[,c("transcript_id", "transcript_biotype")])
  anno <- anno[grep("ENST|ENSMUST", anno$transcript_id), ]
  
  # Report progress
  message(paste(nrow(anno), " transcripts identified", " 转录本是从gtf文件中提取的", sep=""))
  
  #################################################################
  ############# PREPARE WIGGLEPLOTR INPUT: EXON FILE ##############
  #################################################################
  
  # Define transcript ids
  transcript_ids <- anno$transcript_id
  grange.exon.list <- list()
  
  for(i in 1:length(transcript_ids)) {
    
    # Retrieve exons of relevant transcript
    gtf.small <- gtf[which(gtf$transcript_id==transcript_ids[i]), ]
    gtf.small <- gtf.small[which(gtf.small$V3=="exon"), ]
    # Create GRange object
    grange <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(gtf.small$V1),
                                     ranges=IRanges::IRanges(gtf.small$V4,
                                                             width=(gtf.small$V5-gtf.small$V4)+1
                                     ),
                                     strand=gtf.small$V7[1],
                                     exon_id=gtf.small$exon_id,
                                     exon_name=gtf.small$exon_id,
                                     exon_rank=c(1:length(gtf.small$exon_id))
    )
    # Save into list
    grange.exon.list[[i]] <- grange
  }
  
  # Convert list to GenomicRanges::GRanges list
  grange.exon.list <- GenomicRanges::GRangesList(grange.exon.list)
  names(grange.exon.list) <- transcript_ids
  
  #################################################################
  ########### PREPARE WIGGLEPLOTR INPUT: EXON FILE + SJ ###########
  #################################################################
  
  junction_table <- data.frame()
  n <- 0
  
  for(j in 1:length(coord.intron)){
    for(i in 1:length(grange.exon.list)) {
      # Retrieve transcript
      grange <- grange.exon.list[[i]]
      exon <- as.data.frame(grange)
      exon$start <- as.numeric(exon$start)
      exon$end <- as.numeric(exon$end)
      chr.gene <- as.character(unique(exon$seqnames))
      transcript_id.sj <- names(grange.exon.list)[i]
      
      # Retrieve SJ chr, start, end
      . <- strsplit(coord.intron[j], split=":", fixed=TRUE)[[1]]
      chr.sj <- .[1]
      start.sj <- as.numeric(.[2])
      end.sj <- as.numeric(.[3])
      
      if (chr.gene != chr.sj){
        MarvelObject$adhocGene$SJPosition$Plot.1 <- NULL
        MarvelObject$adhocGene$SJPosition$Plot.2 <- NULL
        MarvelObject$adhocGene$SJPosition$metadata <- NULL
        MarvelObject$adhocGene$SJPosition$junction_table <- NULL
        MarvelObject$adhocGene$SJPosition$structure_table <- NULL
        MarvelObject$adhocGene$SJPosition$coord.intron <- NULL
  
        return(MarvelObject)
      }


      # 找到junction的起始位置和终止位置
      if (any(start.sj-1 == exon$end) == TRUE &  any(end.sj+1 == exon$start) == TRUE){
        if (abs(which(start.sj-1 == exon$end)-which(end.sj+1 == exon$start)) == 1 ){
          n = n + 1
            exon.sj.start <- exon$end[which(exon$end==start.sj-1)]
            exon.sj.end <- exon$start[which(exon$start==end.sj+1)]
            left_exon <- exon[which(start.sj-1 == exon$end),]$exon_rank
            right_exon <- exon[which(end.sj+1 == exon$start),]$exon_rank
            left_exon_id <- exon[which(start.sj-1 == exon$end),]$exon_id
            right_exon_id <- exon[which(end.sj+1 == exon$start),]$exon_id
            strand <- exon[which(start.sj-1 == exon$end),]$strand
            
            junction_table[n, "transcript_id"] <- transcript_id.sj
            junction_table[n, "chr"] <- chr.sj
            junction_table[n, "start"] <- start.sj
            junction_table[n, "end"] <- end.sj
            junction_table[n, "left_exon"] <- left_exon
            junction_table[n, "right_exon"] <- right_exon
            junction_table[n, "strand"] <- strand
            junction_table[n, "left_exon_id"] <- left_exon_id
            junction_table[n, "right_exon_id"] <- right_exon_id

        } else {
            message(paste0("转录本 ",transcript_id.sj," 没有可变剪切 ", coord.intron[j] ))
        }

      } else {
        message(paste0("转录本 ",transcript_id.sj," 没有可变剪切 ", coord.intron[j] ))
      }
    }
  }

  junction_table$start <- as.integer(junction_table$start)
  junction_table$end <- as.integer(junction_table$end)
  
  #################################################################
  ################ PREPARE WIGGLEPLOTR INPUT: CDS #################
  #################################################################
  
  # Define transcript ids
  transcript_ids <- anno$transcript_id
  grange.cds.list <- list()
  transcript.ids <- NULL
  
  for(i in 1:length(transcript_ids)) {
    
    # Retrieve exons of relevant transcript
    gtf.small <- gtf[which(gtf$transcript_id==transcript_ids[i]), ]
    gtf.small <- gtf.small[which(gtf.small$V3=="CDS"), ]
    
    if(nrow(gtf.small) != 0) {
      
      # Create GRange object
      grange <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(gtf.small$V1),
                                       ranges=IRanges::IRanges(gtf.small$V4,
                                                               width=(gtf.small$V5-gtf.small$V4)+1
                                       ),
                                       strand=gtf.small$V7[1],
                                       exon_id=gtf.small$exon_id,
                                       exon_name=gtf.small$exon_id,
                                       exon_rank=c(1:length(gtf.small$exon_id))
      )
      
      # Save into list
      grange.cds.list[[i]] <- grange
      transcript.ids[i] <- transcript_ids[i]
      
    } else {
      # Save into list
      grange.cds.list[[i]] <- FALSE
      transcript.ids[i] <- FALSE
      
    }
  }
  
  # Convert list to GenomicRanges::GRanges list
  index.keep <- which(transcript.ids != FALSE)
  grange.cds.list <- grange.cds.list[index.keep]
  transcript.ids <- transcript.ids[index.keep]
  
  if(length(grange.cds.list) != 0) {
    
    grange.cds.list <- GenomicRanges::GRangesList(grange.cds.list)
    names(grange.cds.list) <- transcript.ids
    
  }
  
  
  #################################################################
  ############# PREPARE WIGGLEPLOTR INPUT: METADATA ###############
  #################################################################
  
  # Retrieve transcript ids
  metadata <- data.frame("transcript_id"=names(grange.exon.list), stringsAsFactors=FALSE)
  
  # Annotate gene id, name
  #metadata$gene_id <- gene_id
  metadata$gene_short_name <- gene_short_name
  
  # Indicate strand
  strand <- gtf$V7[1]
  
  if(strand=="+") {
    metadata$strand <- 1
  } else {
    metadata$strand <- -1
  }
  
  # Annotate transcript type
  # Retrieve transcript type

  metadata <- plyr::join(metadata, anno, type="left")
  metadata$transcript_id.biotype <- paste(metadata$transcript_id, " (", metadata$transcript_biotype, ")", sep="")
  
  # Update GRange exon list
  . <- data.frame("transcript_id"=names(grange.exon.list), stringsAsFactors=FALSE)
  . <- join(., metadata[,c("transcript_id", "transcript_id.biotype")], by="transcript_id", type="left")
  names(grange.exon.list) <- .$transcript_id.biotype
  
  # Updateexon-SJ table
  junction_table$transcript_id.biotype <- metadata[match(junction_table$transcript_id, metadata$transcript_id),]$transcript_id.biotype
  junction_table$transcript_id <- NULL

  # Update GRange CDS list
  if (length(grange.cds.list) > 0){
    
    # 当候选基因为非编码蛋白的基因时，就没有CDS，执行 if 条件内的命令就会报错
    . <- data.frame("transcript_id"=names(grange.cds.list), stringsAsFactors=FALSE)
    . <- join(., metadata[,c("transcript_id", "transcript_id.biotype")], by="transcript_id", type="left")
    names(grange.cds.list) <- .$transcript_id.biotype

  }

  # Update metadata
  metadata$transcript_id <- metadata$transcript_id.biotype
  # Remove intermediate column
  metadata$transcript_id.biotype <- NULL
  metadata$transcript_biotype <- NULL
  
  #################################################################
  ################# FILTER FOR SPECIFIC BIOTYPE ###################
  #################################################################
  
  if(show.protein.coding.only==TRUE) {
    transcript_ids <- metadata[grep("protein_coding", metadata$transcript_id, fixed=TRUE), "transcript_id"]
    junction_protein_coding_only <- junction_table[grep("protein_coding", junction_table$transcript_id, fixed=TRUE), ]
    
    if(length(transcript_ids) != 0 & dim(junction_protein_coding_only)[1] != 0) {
      metadata <- metadata[grep("protein_coding", metadata$transcript_id, fixed=TRUE), ]
      grange.exon.list <- grange.exon.list[metadata$transcript_id]
      overlap <- intersect(names(grange.cds.list), transcript_ids)
      grange.cds.list <- grange.cds.list[overlap]
      junction_table <- junction_table[grep("protein_coding", junction_table$transcript_id, fixed=TRUE), ]
    } else {
      message(paste0("基因 ", gene_short_name, " 的junction并不位于编码蛋白质的转录本上，绘制所有转录本" ))
      # MarvelObject$adhocGene$SJPosition$metadata <- metadata
      # return(MarvelObject)
    }
  }
  
  #################################################################
  ######################## WIGGLEPLOTR ############################
  #################################################################
  
  if ( length(grange.cds.list) > 0 ){
    plot <- wiggleplotr::plotTranscripts(exons=grange.exon.list,
                                        cdss=grange.cds.list,
                                        transcript_annotations=metadata,
                                        rescale_introns=rescale_introns,
                                        new_intron_length=50,
                                        flanking_length=c(50,50),
                                        connect_exons=TRUE,
                                        transcript_label=TRUE,
                                        region_coords = NULL
    ) + scale_fill_manual(values = c("gray", "black"), name = "box") +
        theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        )
  } else {

    plot <- wiggleplotr::plotTranscripts(exons=grange.exon.list,
                                        cdss=NULL,
                                        transcript_annotations=metadata,
                                        rescale_introns=rescale_introns,
                                        new_intron_length=50,
                                        flanking_length=c(50,50),
                                        connect_exons=TRUE,
                                        transcript_label=TRUE,
                                        region_coords = NULL
    ) + scale_fill_manual(values = c("black", "gray"), name = "box") +  
        theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        )
  }
  
  plot_data <- plot$data
  
  if (rescale_introns == TRUE){
    junction_table$rescale_start = ""
    junction_table$rescale_end = ""
    transcript_id.sj <- junction_table$transcript_id
    for (id in transcript_id.sj) {
      left_exon <- junction_table[junction_table$transcript_id == id,]$left_exon
      right_exon <- junction_table[junction_table$transcript_id == id,]$right_exon
      junction_table[junction_table$transcript_id == id,]$rescale_start <- plot_data[plot_data$transcript_id == id,][left_exon,]$end
      junction_table[junction_table$transcript_id == id,]$rescale_end <- plot_data[plot_data$transcript_id == id,][right_exon,]$start
    }
    if("" %in% junction_table$rescale_start | "" %in% junction_table$rescale_end){
      message("Error: junction位置rescale失败")
      quit(save = "no", status = 1, runLast = TRUE)
    }

    junction_table$rescale_start <- as.integer(junction_table$rescale_start)
    junction_table$rescale_end <- as.integer(junction_table$rescale_end)
    junction_table$row_start <- junction_table$start
    junction_table$row_end <- junction_table$end
    junction_table$start <- junction_table$rescale_start
    junction_table$end <- junction_table$rescale_end
    junction_table$rescale_end <- NULL
    junction_table$rescale_start <- NULL
  }
  
  
  junction_table$transcript_rank <- plot_data[match(junction_table$transcript_id.biotype, plot_data$transcript_id),]$transcript_rank

  if ("row_start" %in% colnames(junction_table)){
    junction_table$color <- paste0(junction_table$chr, ":", junction_table$row_start, ":", junction_table$row_end)
    junction_table$width <- junction_table$row_end - junction_table$row_start + 1
    junction_table$biotype <- "junction"
    output_junction_table <- subset(junction_table, select = c(transcript_id.biotype, chr, row_start, row_end, width, strand, left_exon_id, left_exon, right_exon_id, right_exon, biotype ))
    names(output_junction_table) <- c("transcript_name", "chr", "start", "end", "width", "strand", "left_exon_id","left_exon_rank", "right_exon_id", "right_exon_rank", "biotype")
  } else {
    junction_table$color <- paste0(junction_table$chr, ":", junction_table$start, ":", junction_table$end)
    junction_table$width <- junction_table$end - junction_table$start + 1
    junction_table$biotype <- "junction"
    output_junction_table <- subset(junction_table, select = c(transcript_id.biotype, chr, start, end, width, strand, left_exon_id, left_exon, right_exon_id, right_exon, biotype ))
    names(output_junction_table) <- c("transcript_name", "chr", "start", "end", "width", "strand", "left_exon_id","left_exon_rank", "right_exon_id", "right_exon_rank", "biotype")
  }
  
  output_junction_table$gene_name <- gene_short_name
  output_junction_table$gene_id <- gene_id

  plot.2 <- plot + guides(color="none", fill="none") + new_scale_color() +  geom_segment(data=junction_table, linewidth = 2,
                aes(x=start, y=transcript_rank, xend=end, yend=transcript_rank, color = color )) +
                #scale_color_manual(values = c("red", "gold", "green", "blue",  "#94A684", "#FFEEF4", "#98D8AA", "#D67BFF"), name = "segment")
                scale_color_manual(values = rainbow(length(unique(junction_table$color))), name = "Junctions") +
                theme(legend.position="bottom") +
                guides(color = guide_legend( ncol=4, byrow=T ) )
  
  # Save into new slots
  transcript_table <- as.data.frame(grange.exon.list)
  transcript_table$group <- NULL
  transcript_table$exon_name <- NULL
  names(transcript_table)[names(transcript_table) == "group_name"] <- "transcript_name"
  names(transcript_table)[names(transcript_table) == "seqnames"] <- "chr"
  transcript_table$biotype <- "exon"

  if(length(grange.cds.list) > 0){
    cds_table <- as.data.frame(grange.cds.list)
    cds_table$group <- NULL
    cds_table$exon_name <- NULL
    names(cds_table)[names(cds_table) == "group_name"] <- "transcript_name"
    names(cds_table)[names(cds_table) == "seqnames"] <- "chr"
    cds_table$biotype <- "cds"
    structure_table <- rbind(transcript_table, cds_table)
  } else {
    structure_table <- transcript_table
  }

  structure_table$gene_name <- gene_short_name
  structure_table$gene_id <- gene_id

  MarvelObject$adhocGene$SJPosition$Plot.1 <- plot
  MarvelObject$adhocGene$SJPosition$Plot.2 <- plot.2
  MarvelObject$adhocGene$SJPosition$metadata <- metadata
  MarvelObject$adhocGene$SJPosition$exonfile <- grange.exon.list
  MarvelObject$adhocGene$SJPosition$cdsfile <- grange.cds.list
  MarvelObject$adhocGene$SJPosition$junction_table <- output_junction_table
  MarvelObject$adhocGene$SJPosition$structure_table <- structure_table
  MarvelObject$adhocGene$SJPosition$coord.intron <- coord.intron
  
  # Return final object
  return(MarvelObject)
  
}

