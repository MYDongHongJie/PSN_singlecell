library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(org.Ss.eg.db)
  library(ggplot2)
  library(dplyr)
  library(rjson)
  library(DOSE)
  library(jsonlite)
  library(stringr)
	library(RColorBrewer)
  source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/color/color.R")
  colors = colorls$"NPG"
GeneTranslate <- function(genelist,translate){
    transdf <- read.table(translate,sep="\t")
    transgene <- transdf[genelist %in% transdf$V1,V2]
}

enrichment<-function(species,outDir,geneList){
    library(clusterProfiler)
    library(org.Ss.eg.db)
    library(ggplot2)
    library(dplyr)
    library(rjson)
    library(DOSE)
    library(jsonlite)
    # Go enrichment
    if( species=="ssc" ||species=="hsa" | species=="mmu" | species=="rno"){
        species.org=switch(species,
                     "ssc"=org.Ss.eg.db,
                     "hsa"=org.Hs.eg.db,
                     "mmu"=org.Mm.eg.db,
                     "rno"=org.Rn.eg.db,
    		      )
    }
    if(exists("species.org")){
        processGOdb(geneList,species.org,outDir)
    }else{
        processGOtxt(geneList,species,outDir)
    }

   # KEGG enrichment
    processKEGGtxt(geneList,species,outDir)
}

PlotGo <- function(out_df,outDir){
  if( !"Term" %in% colnames(out_df)){
		if(!"Category" %in% colnames(out_df)){
			colnames(out_df) <- c("GO","Term","List","Total","Pvalue","adjustPvalue")
		}else{
			colnames(out_df) <- c("GO","Term","Category","List","Total","Pvalue","adjustPvalue")
		}
    #out_df$Category <- "None"
  }
    ##删除Category为NA的数据
    out_df<-out_df[!is.na(out_df$Category),]
    out_df$Term=substring(out_df$Term,1,50)
		
    out_df<-out_df[!duplicated(out_df[,"Term"]),]
    out_df$rich=out_df$List / out_df$Total
    out_df$Number=out_df$List
    out_df<-out_df[order(out_df[,"rich"]),]
    
    out_df$Term=factor(out_df$Term,levels=unique(as.character(out_df$Term)))
    #print(out_df)
		#Term
	
		text_size=9
		nTermSum= sum(nchar(as.character(out_df$Term))>30)


    q<-qplot(rich,Term,data=out_df,colour= Pvalue,size=Number,main="GO Enrichment")+
      scale_colour_gradient(low="red",high="green")+
      theme(panel.background=element_rect(fill="white",color="black"),
			panel.grid.major=element_line(color="grey80",linetype="dotted"),
			axis.text.y=element_text(colour="black",size=text_size)
			,aspect.ratio = 6/4
			)+
      scale_y_discrete(labels=function(x) str_wrap(x, width=30))
    ggsave(paste(outDir,"GO.richfactor.pdf",sep="/"),width = 6+nTermSum*0.06, height = (6+nTermSum*0.06)*4/6)
    ggsave(paste(outDir,"GO.richfactor.png",sep="/"),width = 6+nTermSum*0.06, height = (6+nTermSum*0.06)*4/6)

    out_df=out_df[order(out_df$Category,out_df$List),]
		
    out_df$Category <- factor(out_df$Category, levels = unique(out_df$Category))
    out_df$Term <- factor(out_df$Term, levels = out_df$Term)
    
		
		p<- ggplot(data = out_df) +
      geom_col(aes(x = Term, y = List, fill = Category),width =0.8) +
      scale_fill_manual(values=rev(c("#FFEDA0", "#FEB24C", "#F03b01")))+
      #scale_color_brewer(type="seq",palette="Dark2")+
      
      scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
      ylab("Counts")+coord_flip()+scale_y_continuous(expand=c(0, 0))+ theme_bw()+theme(panel.grid.major=element_line(colour=NA),
            panel.grid.minor = element_blank(),
						panel.border = element_blank(),
        panel.grid = element_blank())+theme( strip.text.y = element_text(angle = 0),axis.text.x=element_text(colour="black"),
			axis.text.y=element_text(colour="black",size=text_size),
			aspect.ratio = 3/4) 
    ggsave(plot=p,filename=paste(outDir,"GO_enrichment_pvalue_barplot.pdf",sep="/"), width = 8+nTermSum*0.06, height =(8+nTermSum*0.06)/2)
    ggsave(plot=p,filename=paste(outDir,"GO_enrichment_pvalue_barplot.png",sep="/"), width = 8+nTermSum*0.06, height = (8+nTermSum*0.06)/2)

}

prepareMapping<-function(pathwayInfo,tag="term2gene",header=T){
      mat = pathwayInfo
      mat=switch(tag,
                 term2gene=mat[,c(2,1)],
                 term2name=mat[,c(2,3)])
      return(mat)
    }
processGOtxt <- function(geneList,species,outDir){
    json_dic<-read_json("/PERSONALBIO/work/singlecell/s00/software/script/1.source/ref.json")
	if('desc' %in% names(json_dic[[species]])){
		desc_file <- json_dic[[species]][['desc']]
		pathwayInfo_GO = read.csv(desc_file, colClasses = "character")
		pathwayInfo_GO <- pathwayInfo_GO[pathwayInfo_GO$database=='GO',]
	}else{
		pathwayInfo_GO_file <- json_dic[[species]][['go_desc']]
		pathwayInfo_GO = read.delim(pathwayInfo_GO_file, header=F, sep="\t",check.names=F,quote="",colClasses = "character")
	}

    term2gene=prepareMapping(pathwayInfo=pathwayInfo_GO,tag="term2gene")
    term2name=prepareMapping(pathwayInfo=pathwayInfo_GO,tag="term2name")

    go_enrichment<-enricher(geneList,
                            TERM2GENE = term2gene,
                            TERM2NAME = term2name,
                            pAdjustMethod = "BH",
                            minGSSize=1 ,maxGSSize=100000,
                            qvalueCutoff = 1, pvalueCutoff=1)
    if(!is.null(go_enrichment)){
        res=go_enrichment@result
        res$Total=apply(res,1,function(x){getBGnumber(x[4])})
        df=data.frame( PathwayID=res$ID,
                     Pathway=res$Description,
                     List=res$Count,
                     Total=res$Total,
                     Pvalue=res$pvalue,
                     adjustPvalue=res$p.adjust,
                     Gene=res$geneID
        )
				Category_list = read.delim("/PERSONALBIO/work/singlecell/s04/Test/go_hierarchy.txt",header =T)
				df$Category <- Category_list$ONTOLOGY[match(df$PathwayID, Category_list$GO)]
				df =df[,c("PathwayID","Pathway","Category","List","Total","Pvalue","adjustPvalue","Gene")]
        write.table(df,file=paste(outDir,"/","GO_enrichment.xls",sep=""),col.names=T,row.names=F,quote=F,sep='\t')
        if(!'Category' %in% colnames(df)){
					if(nrow(df) < 20){out_df <- df[seq(1,nrow(df),1),]}else{out_df <- df[seq(1,20,1),]} 
				}else{
					#分组取每个组的top5
					out_df <- df %>%  filter(Pvalue < 0.05) %>%  group_by(Category) %>% arrange(Category,  desc(List), Pvalue) %>% slice_head(n=5) %>%  ungroup()
          out_df = as.data.frame(out_df)
				}
        PlotGo(out_df,outDir)
    }
}

processGOdb<- function(geneList,species.org,outDir){
  en <- enrichGO(gene       = geneList,
                 OrgDb         = species.org,
                 keyType       = 'SYMBOL',
                 ont           = "ALL",
                 minGSSize     = 1 ,
                 maxGSSize     = 100000 ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
    if(!is.null(en)){
        res=en@result
        res$Total=apply(res,1,function(x){getBGnumber(x[5])})
        df=data.frame( Category=res$ONTOLOGY,
                   GO=res$ID,
                   Term=res$Description,
                   List=res$Count,
                   Total=res$Total,
                   Pvalue=res$pvalue,
                   adjustPvalue=res$p.adjust,
                   Gene=res$geneID
        )
        df<-df[order(df$Pvalue),]
        write.table(df,file=paste(outDir,"/","GO_enrichment.xls",sep=""),col.names=T,row.names=F,quote=F,sep='\t')
    
        if(!'Category' %in% colnames(df)){
					if(nrow(df) < 20){out_df <- df[seq(1,nrow(df),1),]}else{out_df <- df[seq(1,20,1),]} 
				}else{
					#分组取每个组的top5
					out_df <- df %>%  filter(Pvalue < 0.05) %>%  group_by(Category) %>% arrange(Category,  desc(List), Pvalue) %>% slice_head(n=5) %>%  ungroup()
          out_df = as.data.frame(out_df)
				}
        PlotGo(out_df,outDir)
    }
}



PlotKegg <- function(out_df,outDir){
    if(!("Category" %in% colnames(out_df))){
        out_df$Category <- "None"
    }
    #out_df$Pathway=factor(out_df$Pathway,levels=unique(out_df$Pathway))
    out_df$Pathway=substring(out_df$Pathway,1,50)
    out_df<-out_df[!duplicated(out_df[,"Pathway"]),]
    out_df$rich=out_df$List / out_df$Total
    out_df$Number=out_df$List
    out_df<-out_df[order(out_df[,"rich"]),]
    out_df$Pathway=factor(out_df$Pathway,levels=unique(as.character(out_df$Pathway)))

		
		text_size=8
		nTermSum = sum(nchar(as.character(out_df$Pathway))>30)

    q<-qplot(rich,Pathway,data=out_df,colour= Pvalue,size=Number,main="KEGG Enrichment")+
      scale_colour_gradient(low="red",high="green")+
      theme(panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey80",linetype="dotted"),
			axis.text.y=element_text(colour="black",size=text_size),
			aspect.ratio = 6/4)+
      scale_y_discrete(labels=function(x) str_wrap(x, width=30))
    ggsave(plot=q,filename=paste(outDir,"KEGG.richfactor.pdf",sep="/"),width = 6+nTermSum*0.3, height = (6+nTermSum*0.3)*4/6)
    ggsave(plot=q,filename=paste(outDir,"KEGG.richfactor.png",sep="/"),width = 6+nTermSum*0.3, height = (6+nTermSum*0.3)*4/6)

    out_df=out_df[order(out_df$List),]
    out_df$Category <- factor(out_df$Category, levels = unique(out_df$Category))
    out_df$Pathway <- factor(out_df$Pathway, levels = out_df$Pathway)

		custom_colors <- colorRampPalette(rev(brewer.pal(3, "YlOrRd")))(299)
    p<- ggplot(data = out_df) +
      geom_col(aes(x = Pathway, y = List, fill = Pvalue),width =0.8) +
       
      scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
      ylab("Counts") +coord_flip()+scale_y_continuous(expand=c(0, 0))+ theme_bw()+theme(panel.grid.major=element_line(colour=NA),
            panel.grid.minor = element_blank(),
						panel.border = element_blank(),
        panel.grid = element_blank())+theme(strip.text.y = element_text(angle = 0),axis.text.x=element_text(colour="black"),
			axis.text.y=element_text(colour="black",size=text_size),
			aspect.ratio = 3/4) + scale_fill_gradientn(colors = custom_colors)
    ggsave(plot = p,filename =paste(outDir,"KEGG_enrichment_pvalue_barplot.pdf",sep="/"), width = 8+nTermSum*0.4, height = (8+nTermSum*0.4)/2)
    ggsave(plot = p,filename =paste(outDir,"KEGG_enrichment_pvalue_barplot.png",sep="/"), width = 8+nTermSum*0.4, height = (8+nTermSum*0.4)/2)
}

processKEGGtxt <- function(geneList,species,outDir){
    json_dic<-read_json("/PERSONALBIO/work/singlecell/s00/software/script/1.source/ref.json")
	if('desc' %in% names(json_dic[[species]])){
		desc_file <- json_dic[[species]][['desc']]
		pathwayInfo_KEGG = read.csv(desc_file, colClasses = "character")
		pathwayInfo_KEGG <- pathwayInfo_KEGG[pathwayInfo_KEGG$database=='KEGG',]
	}else{
		pathwayInfo_KEGG_file <- json_dic[[species]][['kegg_desc']]

		pathwayInfo_KEGG = read.delim(pathwayInfo_KEGG_file, header=F, sep="\t",check.names=F,quote="",colClasses = "character")
	}
    
    term2gene=prepareMapping(pathwayInfo=pathwayInfo_KEGG,tag="term2gene")
    term2name=prepareMapping(pathwayInfo=pathwayInfo_KEGG,tag="term2name")
    kegg_enrichment<-enricher(geneList,
                              TERM2GENE = term2gene,
                              TERM2NAME = term2name,
                              pAdjustMethod = "BH",
                              minGSSize=1 ,maxGSSize=100000,
                              qvalueCutoff = 1, pvalueCutoff=1)

    if(!is.null(kegg_enrichment)){ #判断是否为空
        res=kegg_enrichment@result
        res$Total=apply(res,1,function(x){getBGnumber(x[4])})
        df=data.frame( PathwayID=res$ID,
                     Pathway=res$Description,
                     List=res$Count,
                     Total=res$Total,
                     Pvalue=res$pvalue,
                     adjustPvalue=res$p.adjust,
                     Gene=res$geneID
        )
		if('desc' %in% names(json_dic[[species]])){
			df$Category <- pathwayInfo_KEGG$category[match(df$PathwayID, pathwayInfo_KEGG$term_id)]
		}
        write.table(df,file=paste(outDir,"/","KEGG_enrichment.xls",sep=""),col.names=T,row.names=F,quote=F,sep='\t')
        if(nrow(df) < 20){out_df <- df[seq(1,nrow(df),1),]}else{out_df <- df[seq(1,20,1),]}
        PlotKegg(out_df,outDir) 
    }
}

processKEGGdb<- function(geneList,species,outDir){
    species.kegg=switch(species,
                      "ssc"="/PERSONALBIO/work/singlecell/s00/software/script/1.source/KEGG.db/ssc/",
                      "cgib"="/PERSONALBIO/work/singlecell/s00/software/script/1.source/KEGG.db/cgib/",
                      "dosa"="/PERSONALBIO/work/singlecell/s00/software/script/1.source/KEGG.db/dosa/",
    )
    library(KEGG.db,lib.loc=species.kegg)
    en <- enrichKEGG(gene     = geneList,
                   organism      = species,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   minGSSize     = 1 ,
                   maxGSSize     = 100000 ,
                   qvalueCutoff  = 1,
                   use_internal_data =T)
    if(species=="ssc"){
        en <- DOSE::setReadable(en, OrgDb='org.Ss.eg.db',keyType='ENTREZID')
    }
    kegg_level=read.table("/PERSONALBIO/work/singlecell/s00/software/script/1.source/KEGG.db/pathway_level",header=F,sep="\t",stringsAsFactors=F,colClasses="character",row.names=1)
    if(!is.null(en)){
        res=en@result
        res$Total=apply(res,1,function(x){getBGnumber(x[4])})
        df=data.frame( Category=kegg_level[gsub(species,"",res$ID),1] ,
                   PathwayID=res$ID,
                   Pathway=res$Description,
                   List=res$Count,
                   Total=res$Total,
                   Pvalue=res$pvalue,
                   adjustPvalue=res$p.adjust,
                   Gene=res$geneID
        )
        write.table(df,file=paste(outDir,"/","KEGG_enrichment.xls",sep=""),col.names=T,row.names=F,quote=F,sep='\t')
        if(nrow(df) < 20){out_df <- df[seq(1,nrow(df),1),]}else{out_df <- df[seq(1,20,1),]}
        PlotKegg(out_df,outDir)
    }
}

processREACTOME <- function(geneList,species,outDir){
  library(ReactomePA)
  species.reactome=switch(species,
                          'hsa'='human',
                          'mmu'='mouse',
                          'rno'="rat",
  )
  en <- enrichPathway(gene       = geneList,
                      organism      = species.reactome,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      minGSSize     = 1 ,
                      maxGSSize     = 100000 ,
                      qvalueCutoff  = 1)
  if(species=="hsa"){
    en <- DOSE::setReadable(en, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  }
  if(species=="mmu"){
    en <- DOSE::setReadable(en, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
  }
  #   if(length(unique(en@result$p.adjust<0.05))>1){
  if(!is.null(en)){
    res=en@result
    res$Total=apply(res,1,function(x){getBGnumber(x[4])})
    df=data.frame( PathwayID=res$ID,
                   Pathway=res$Description,
                   List=res$Count,
                   Total=res$Total,
                   Pvalue=res$pvalue,
                   adjustPvalue=res$p.adjust,
                   Gene=res$geneID
    )
    write.table(df,file=paste(outDir,"/","Reactome_enrichment.xls",sep=""),col.names=T,row.names=F,quote=F,sep='\t')
    
    #Plot
    if(nrow(df) < 20){out_df <- df[seq(1,nrow(df),1),]}else{out_df <- df[seq(1,20,1),]}
    out_df$Pathway=factor(out_df$Pathway,levels=unique(out_df$Pathway))
    out_df$Pathway=substring(out_df$Pathway,1,50)
    out_df<-out_df[!duplicated(out_df[,"Pathway"]),]
    out_df$rich=out_df$List / out_df$Total
    out_df$Number=out_df$List
    q<-qplot(rich,Pathway,data=out_df,colour=Pvalue,size=Number,main="Reactome Enrichment")+
      scale_colour_gradient(low="red",high="green",limits=c(0,1))+
      theme(panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey80",linetype="dotted"))
    ggsave(paste(outDir,"Reactome.richfactor.pdf",sep="/"),width = 10, height = 10)
    ggsave(paste(outDir,"Reactome.richfactor.png",sep="/"),width = 10, height = 10)
    
    
    p<- ggplot(data = out_df) +
      geom_col(aes(y = Pathway, x = -log10(as.numeric(Pvalue))),width =0.8) +
      theme( strip.text.y = element_text(angle = 0),axis.text.x=element_text(size=10,angle=80,hjust=1)) +
      xlab("-log10(P-value)")
    ggsave(paste(outDir,"Reactome_enrichment_pvalue_barplot.pdf",sep="/"), width = 10, height = 10)
    ggsave(paste(outDir,"Reactome_enrichment_pvalue_barplot.png",sep="/"), width = 10, height = 10)
    
  }}

richfactor<-function(en){
  total <- apply(en,1,function(x){getBGnumber(x[4])})
  return(en$Count/total)
}
getBGnumber <- function(ratio,split="/"){
  list<-strsplit(ratio, split = split)[[1]]
  return(as.numeric(list[1]))
}
