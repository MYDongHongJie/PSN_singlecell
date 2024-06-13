### example
docstring <- "example1:\\n\\n\\
 sctool -i singlecell_object.clustering_resolution0.4.rds \\\\\\n\\
        -f rds \\\\\\n\\
        --assay SCT \\\\\\n\\
        --dataslot data \\\\\\n\\
        -o result  \\\\\\n\\
        score  --gmt  celltype.gmt  \\\\\\n\\
               --min.size 1  \\\\\\n\\
               --plot featureplot  \\\\\\n\\
               -g clusters \\\\\\n\\
               --ccolors spectral"
# calculate the score of gene modules
sub_score <- subparsers$add_parser("score",
                                   help = "calculate the gene module score",
                                   description = docstring,
                                   formatter_class = 'argparse.RawTextHelpFormatter',
                                   #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                   argument_default = "True")
# sub_score$add_argument("--name", type = "character", default = NULL,
#                        help = paste0("the  name of internal gene set. ",
#                                      "Choices can be:Macosko_mouse, Macosko_human, Jeny2020.[default: %(default)s] "))
sub_score$add_argument("--gmt", type = "character", default = NULL,
                       help = "the gene sets in gmt format for calculation.[default: %(default)s] ")
sub_score$add_argument("--min.size", type = "integer", default = 5,
                       help = "the minimal number of genes in each gene set.[default: %(default)s]")
sub_score$add_argument("--plot", type = "character", default = NULL,
                       #choices=c("featureplot","vlnplot","NULL"),
                       help = paste0("methods used to visualize the gene module score of each cell. Choices:NULL,vlnplot,",
                                     "featureplot. [default: %(default)s]  no plot"))
sub_score$add_argument("-g", "--groupby", type = "character", default = "clusters",
                       help = paste0("[OPTIONAL]The grouppinig variable in the metadata for separate the cells to visulize",
                                     "marker genes.[default: %(default)s] "))
sub_score$add_argument("--reduct", type = "character", default = "umap",
                       help = paste0("[OPTIONAL]name(s) of reductions[default: %(default)s].NULL:all reductions,NA:only ",
                                     "global reductions,FALSE:no reductions.  Note: Only reductions associated with an assay",
                                     " loaded in assays or marked as global will be loaded."))
sub_score$add_argument("-s", "--pointsize", type = "double", default = 0.01,
                       help = "[OPTIONAL]the point size in the plot.")
sub_score$add_argument("-r", "--spointsize", type = "double", default = 0.8,
                       help = "[OPTIONAL]the point size in the plot.")
sub_score$add_argument("-a", "--HE", type = "character", default = "T", help = "切片背景")
sub_score$add_argument("-t", "--textsize", type = "double", default = 20,
                       help = "[OPTIONAL]the text size in the plot.")
sub_score$add_argument("--dodge", type = "character", default = "FALSE",
                       help = paste0("[OPTIONAL]visualize the feature between the contrast groups separately for each ",
                                     "level in variable specified by --groupby.[default: %(default)s] "))
sub_score$add_argument("-y", "--splitby", type = "character", default =  NULL,
                       help = paste0("[OPTIONAL]the variable in the metadata used to split the graph by the variable ",
                                     "levels to comparing the module difference in different levels.[default: %(default)s] "))
sub_score$add_argument("--ccolors", type = "character", default = "spectral",
                       help = paste0("[OPTIONAL]the name of customized continious color palatte, recommandations:spectral,",
                                     " solar_extra,flame_light or color scale used to map the continious expression value ",
                                     "for feature plot or dotplot,format:low_color,high_color.[default: %(default)s]"))
sub_score$add_argument("--palette", type = "character", default = "customecol2",
                      help = paste0("the discrete color schema mapped to the cell annotations specified by --groupby.[default: %(default)s]",
                                    " Choices:customecol2:50,blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,",
                                    " tableau20:20,Buen:17,UKBB:18,TF1:17,paired:12"))
sub_score$add_argument("--min.cutoff", type = "character", default = "q10",
                       help = "min.cutoff value in featureplot show ")
sub_score$add_argument("--max.cutoff", type = "character", default = "q90",
                       help = "max.cutoff value in featureplot show  ")
sub_score$add_argument(
    "--crop",
    default = 'FALSE',
    help = "whether to crop in spatialplot, data from cytassist project should be 'TRUE'.[default: %(default)s]"
)
sub_score$add_argument(
        "--pvalue",
        type = "character",
        default = NULL,
        help = "[OPTIONAL]use like GROUP1:GROUP2+GROUP2:GROUP3 to add pvalue on vlnplot."
)

# =============== Subcmd: score, gene module score calculation =========
args <- commandArgs(TRUE)
if ("score" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "score") {
    futile.logger::flog.info("step1:read the specified assay and data slot in data object into memory") #=================================
    suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                    informat = opt$informat,
                                                    assays = assays,
                                                    data.use = dataslots,  # data slot is enough
                                                    verbose = F))
    Seurat::DefaultAssay(data_ob) <- assays
    futile.logger::flog.info("step2:import gmt information by use gmt = GSEABase::getGmt(con=opt$gset)") #================================
    gset_list <- GSEABase::geneIds(GSEABase::getGmt(con = opt$gmt))
    gene_list <- unlist(gset_list, use.names = FALSE)
    match <- Seurat::CaseMatch(search = as.vector(unlist(gset_list, use.names = FALSE)), match = rownames(data_ob@assays$SCT))
    ## if have unmatch genes,save 
    unmatch <- gene_list[!gene_list %in% names(match)]
    if (length(unmatch) > 0) {
      futile.logger::flog.warn("There are some gene name unmatch gene in query object,please check:", unmatch, capture = TRUE)
      unmatch %>% as.data.frame %>% readr::write_tsv(file = glue::glue("{output_dir}/unuse_gene.list.xls"), col_names = FALSE)
      #quit()
    }
    for (name in names(gset_list)){
        match <- Seurat::CaseMatch(search = gset_list[[name]], match = rownames(data_ob@assays$SCT))
        gset_list[[name]]=as.vector(match)
    }
    futile.logger::flog.info(glue::glue("step3:filter gene sets by {opt$min.size}")) #==================================
    gset_list <- gset_list[lengths(gset_list) > opt$min.size]
    futile.logger::flog.info("step4:using seurat::Addmodulescore function to get module score") #=======================
    data_ob <- Seurat::AddModuleScore(data_ob, features = gset_list, name = "module")
    cellmeta <- OESingleCell::colData(data_ob)
    cellmeta <- cellmeta %>%  dplyr::rename_at(dplyr::vars(dplyr::starts_with("module")), ~names(gset_list))
    data_ob <- Seurat::AddMetaData(data_ob, metadata = cellmeta)

    if (opt$groupby %in% c("clusters", "sampleid", "group")  ){
          gset_score <- cellmeta %>%
                tibble::rownames_to_column(var = "barcode") %>% 
                dplyr::select(barcode, clusters, sampleid, group, names(gset_list) )
    }else{
          gset_score <- cellmeta %>%
                tibble::rownames_to_column(var = "barcode") %>% 
                dplyr::select(barcode, clusters, sampleid, group, !!opt$groupby, names(gset_list) )
    }
    futile.logger::flog.info("output gene_module_scores.xls")
    readr::write_tsv(gset_score, file = glue::glue("{output_dir}/gene_module_scores.xls"))
    futile.logger::flog.info(glue::glue("step5:ploting: {opt$plot} ")) #================================================
    if (!is.null(opt$plot)) {
      vismethods <- unlist(strsplit(opt$plot, ",")) %>% as.vector()
      # violin plot or featureplot for each gene module in specified group
      ccolors <- unlist(strsplit(opt$ccolors, ","))
      if (length(ccolors) == 1) { # this means to use customized discrete color schema in OESingleCell package
        # info <- glue::glue("the customized continious colors are: \n {paste(names(OESingleCell::continuous_palette),
        # lengths(OESingleCell::continuous_palette), sep =': ', collapse = '\n')}.")
        # warning("NO specified color schema found!")
        # warning(info)
        # stop()
        continious_colors <- OESingleCell::SelectColors(palette = ccolors, is.discrete = FALSE)
      }else { # this means to use the customized color schema from command line specified by user
        continious_colors <- ccolors
      }

      lapply(vismethods, function(method) {
        if (method == "featureplot") {
            lapply(names(gset_list), function(x) {
           # parallel::mclapply(names(gset_list),function(x) {
            pdf(NULL)
            ##for featureplot 目前版本seurat中FeaturePlot split.by使用时存在bug,无法统一分图标尺,需要人为修改
            min.feature.value <- Seurat::SetQuantile(opt$min.cutoff, data_ob@meta.data[, x])
            max.feature.value <- Seurat::SetQuantile(opt$max.cutoff, data_ob@meta.data[, x])
            ###merge
            ggf1 <- Seurat::FeaturePlot(data_ob,
                                        reduction = opt$reduct,
                                        features = x,
                                        #split.by = opt$splitby,
                                        keep.scale = "all",
                                        min.cutoff = min.feature.value,
                                        max.cutoff = max.feature.value,
                                        label=FALSE,
                                        #cols =continious_colors,
                                        #blend=TRUE,
                                        combine = TRUE,
                                        order = T,
                                        pt.size = opt$pointsize) &
              ggplot2::scale_color_gradientn(colours = continious_colors,
                                             limits = c(min.feature.value, max.feature.value))
            ggf1 <- ggf1 + theme(legend.position = "top",legend.key.width=unit(1,'cm'),legend.direction ="horizontal",legend.title=element_blank(),legend.text=element_text(size=10))+ggtitle("")
            ##split
            if(!is.null(opt$splitby)){
              spitby_list<-sort(unique(cellmeta[[opt$splitby]]))
              print(spitby_list[stringr::str_order(spitby_list)])
              data_ob@meta.data[[opt$splitby]] <- factor(data_ob@meta.data[[opt$splitby]], levels = spitby_list[stringr::str_order(spitby_list)])
              ggf2 <- lapply(spitby_list[stringr::str_order(spitby_list)],function(y) {
                            desired_cells <- subset(cellmeta, eval(parse(text = glue::glue('{opt$splitby} =="{y}"'))))
                            data_ob.split <- subset(data_ob,cells=rownames(desired_cells))
                            ggf <- Seurat::FeaturePlot(data_ob.split,
                                                reduction = opt$reduct,
                                                features = x,
                                                keep.scale = "all",
                                                min.cutoff = min.feature.value,
                                                max.cutoff = max.feature.value,
                                                label=FALSE,
                                                #cols =continious_colors,
                                                #blend=TRUE,
                                                combine = TRUE,
                                                order = T,
                                                pt.size = opt$pointsize) & ggtitle(y) & theme(plot.title = element_text(size = 15,hjust = 0.5,face="plain"))
                            return(ggf)
                        })

              ggf2<- patchwork::wrap_plots(ggf2, ncol = length(unique(cellmeta[[opt$splitby]])),guides = "collect") &
                            ggplot2::scale_color_gradientn(colours = continious_colors,
                                                       limits = c(min.feature.value, max.feature.value)) &
                    theme(legend.position = "right",legend.title=element_blank(),legend.text=element_text(size=10))  &
                    theme(legend.position = "none")

            }
            ##for spatial featureplot
            if (!is.null(Seurat::Images(data_ob))) {
              crop <- as.logical(opt$crop)
              if ( opt$HE == FALSE){
                gsp <- OESingleCell::SpatialPlot(data_ob ,
                                      features = x,
                                      alpha = 1,
                                      pt.size.factor = opt$spointsize,
                                      stroke = 0,
                                      min.cutoff = min.feature.value,
                                      max.cutoff = max.feature.value,
                                      combine = FALSE,
                                      HE = F,
                                      cols = ccolors,
                                      ncol = length(Images(seurat_ob)),
                                      crop = crop
                )}else{
                  gsp <- OESingleCell::SpatialPlot(data_ob ,
                                                features = x,
                                                combine = FALSE,
                                                alpha = 1,
                                                min.cutoff = min.feature.value,
                                                max.cutoff = max.feature.value,
                                                cols = ccolors,
                                                stroke = 0,
                                                HE = T,
                                                pt.size.factor = opt$spointsize,
                                                crop = crop
                )
                }
              gsp<- patchwork::wrap_plots(gsp, ncol = length(Seurat::Images(data_ob)),guides = "collect")  &
                    # theme(legend.direction ="horizontal")
                     theme(legend.position = "none")
                    # ggplot2::scale_colour_continuous(limits = c(min.feature.value, max.feature.value))&
                    #theme(legend.position = "right",legend.title=element_blank(),legend.text=element_text(size=10))
            }

            splitby_len <-ifelse(!is.null(opt$splitby),length(unique(cellmeta[[opt$splitby]])),1)
            image_len <- ifelse(!is.null(Seurat::Images(data_ob)),length(Seurat::Images(data_ob)),1)

            if( !is.null(Seurat::Images(data_ob))){
              max_colnum <- max(splitby_len, image_len)
              layout <- glue::glue("A{strrep('#',max_colnum-1)}\n{strrep('B',splitby_len)}{strrep('#',max_colnum-splitby_len)}\n{strrep('C',image_len)}{strrep('#',max_colnum-image_len)}")
              if(!is.null(opt$splitby)){
                  ggf_merge <- patchwork::wrap_plots(A = ggf1,
                                                     B = ggf2,
                                                     C = gsp,
                                                     heights = c(1, 1, 1.3),
                                                     guides = "collect",
                                                     design = layout)
              }else {
                  ggf_merge <- patchwork::wrap_plots(A = ggf1,
                                                     C = gsp,
                                                     heights = c(1,0, 1.3),
                                                     guides = "collect",
                                                     design = layout)
              }
              ggf_merge <- ggf_merge + patchwork::plot_annotation(title = glue::glue('Addmodulescore:{x}'),
                                                                  theme = theme(legend.direction = "horizontal",
                                                                                legend.position = "top",
                                                                                legend.key.width = unit(1, 'cm'),
                                                                                plot.title = element_text(size = opt$textsize,
                                                                                                          hjust = 0.5,
                                                                                                          face = "bold")))
                           #theme(text = element_text('mono'))
              width <- 3.1 * max_colnum
              height <- 3.1 *3+1
            }
            else if( !is.null(opt$splitby)){
              layout <- glue::glue("A{strrep('#',splitby_len-1)}\n{strrep('B',splitby_len)}")
              ggf_merge <- patchwork::wrap_plots(A = ggf1,
                                                 B = ggf2,
                                                 heights = c(1, 1),
                                                 design = layout,
                                                 guides="collect")+
                           patchwork::plot_annotation(title = glue::glue('Addmodulescore:{x}'),
                                                      theme = theme(legend.direction ="horizontal",
                                                                    legend.position = "top",
                                                                    legend.key.width=unit(1,'cm'),
                                                                    plot.title = element_text(size = opt$textsize,
                                                                                              hjust = 0.5,
                                                                                              face="bold")))
              width <- 3.1 * splitby_len
              height <- 3.1 *2+1
            }
            else{
              ggf_merge <- ggf1 + patchwork::plot_annotation(title = glue::glue('Addmodulescore:{x}'),
                                                             theme = theme(legend.direction ="horizontal",
                                                                    legend.position = "top",
                                                                    legend.key.width=unit(1,'cm'),
                                                                    plot.title = element_text(size = opt$textsize,
                                                                                              hjust = 0.5,
                                                                                              face="bold")))

              width <- 3.1
              height <- 3.1+1
            }
            #ggm<-patchwork::wrap_plots(plotlist = ggmodules, ncol = 1)
            OESingleCell::save_ggplots(plot = ggf_merge,
                                       filename = glue::glue("{output_dir}/gene_module_for_{x}_featureplot"),
                                       width = width,
                                       height =height,
                                       dpi = 300,
                                       limitsize =F)
            } )
           # } , mc.cores = 8)
        }
        if (!is.null(opt$pvalue)) {
            print("Hello,开始展示P值")

            all_comparisions <- list()
            contrasts_list <- unlist(strsplit(opt$pvalue, "\\+", perl = T))
            #print(contrasts_list)
            # contrasts_list = unlist(strsplit(pvalue, "\\+", perl = T)) ##delete
            for (contrast in contrasts_list) {
                contrast <- paste0(as.character(opt$groupby), ":", contrast)
                contrasts <- unlist(strsplit(contrast, ":", perl = T))
                assay_metadata <- data_ob@meta.data
                #print(contrasts[1])
                #print(contrasts[2])
                #print(contrasts[3])
                all_levels <- as.vector(unique(assay_metadata[, opt$groupby]))
                print(all_levels)
                if (contrasts[2] == "all" & contrasts[3] != "all") {
                    all_levels <- all_levels[-which(all_levels == contrasts[3])] #delete the reference level
                    all_comparisions <- paste(all_levels, contrasts[3], sep = ":")
                    break
                }else if (contrasts[2] == "all" & contrasts[3] == "all") {
                    all_comparisions <- "all"
                    break
                }else if (contrasts[2] != "all" & contrasts[3] == "all") {
                    all_levels <- all_levels[-which(all_levels == contrasts[2])]
                    all_comparisions <- paste(contrasts[2], all_levels, sep = ":")
                    break
                }else {
                    if (!contrasts[2] %in% all_levels | !contrasts[3] %in% all_levels) {
                        print(paste0(contrasts[2], ":", contrasts[3], "所选分组中细胞数为0,请检查分组比较信息。已跳过该分组。"))
                    }else if (table(assay_metadata[, opt$groupby])[contrasts[2]] <= 1 | table(assay_metadata[, opt$groupby])[contrasts[3]] <= 1) {
                        print(paste0(contrasts[2], ":", contrasts[3], "所选分组中细胞数小于2,请检查分组比较信息。已跳过该分组。"))
                    }else {
                        all_comparisions <- c(all_comparisions, paste0(contrasts[2], ":", contrasts[3]))
                    }
                }
            }
            print(paste0("确定比较分组：",all_comparisions))

            #找分组信息
            my_comparisons <- sort(unique((assay_metadata[, opt$groupby])))
            comp <- list()
            ##若为all:all,则通过循环生成各个分组两两结合比对的list
            `%!in%` <- Negate(`%in%`)
            if ("all" %!in% all_comparisions) {
                for (i in all_comparisions) {
                  list <- strsplit(i, ":", perl = TRUE)
                  comp <- c(comp, list)
                }
            }else {
                for (a in 1:(length(my_comparisons) - 1)) {
                  for (b in 1:(length(my_comparisons) - a)) {
                    list <- c(as.character(my_comparisons[a]), as.character(my_comparisons[a + b]))
                    comp <- c(comp, list(list))
                  }
                }  
            }
            print(paste0("最终画图比较分组信息：",comp))
        }
        if (method == "vlnplot"){
            print("vlnplot")
            plot_data <- data_ob@meta.data %>% dplyr::select(c(opt$groupby, names(gset_list)))
            plot_data[[opt$groupby]] <- as.factor(plot_data[[opt$groupby]])
            gs <- lapply(gsub("-", ".", names(gset_list)), function(gset) {
                plot_value <- plot_data %>% dplyr::select(c(opt$groupby, tidyselect::all_of(gset)))
                colnames(plot_value) <- c("group", "value")
                #print(head(plot_data))
                #print(gset)
                Pvalue = opt$pvalue
                if (is.null(Pvalue)) {
                  print("不进行P值展示")
                  gs <- ggstatsplot::ggbetweenstats(
                          data = plot_value,
                          x = group,
                          y = value,
                          plot.type = "boxviolin",
                          results.subtitle = FALSE,
                          messages = FALSE,
                          pairwise.comparisons = FALSE,
                          mean.label.size = 0, # size of the label for mean
                          centrality.plotting = F,
                          ylab = gset
                  ) +
                          scale_color_manual(values = OESingleCell::SelectColors(unique(plot_data[[opt$groupby]]))) +
                          theme(axis.text.x = element_text(size = 8, colour = "black"),
                                axis.text.y = element_text(size = 8, colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                axis.line = element_line(colour = "black"))
                }else{
                    print("进行P值展示")
                    gs <- ggstatsplot::ggbetweenstats(
                        data = plot_value,
                        x = group,
                        y = value,
                        plot.type = "boxviolin",
                        results.subtitle = FALSE,
                        messages = FALSE,
                        pairwise.comparisons = FALSE,
                        mean.label.size = 0, # size of the label for mean
                        centrality.plotting = F,
                        ylab = gset
                ) +
                        scale_color_manual(values = OESingleCell::SelectColors(unique(plot_data[[opt$groupby]]))) +
                        theme(axis.text.x = element_text(size = 8, colour = "black"),
                              axis.text.y = element_text(size = 8, colour = "black"),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              axis.line = element_line(colour = "black"))
                    gs <- gs + ggsignif::geom_signif(comparisons = comp, #指定比较对象
                                                     test = "wilcox.test", #指定检验方法
                                                     size = 0.4, #指定标记中线条的尺寸
                                                     textsize = 2.6, #指定标记中文字部分的大小
                                                     vjust = -0.05, #指定标记中文字部分与横线之间的距离指定标记中文字部分与横线之间的距离
                                                     step_increase = 0.1, #指定每根线条距离高低
                                                     #tip_length = c(0.2, 0.2), #指定短竖线的长度
                                                     map_signif_level = T) #F显示具体数值，T显示显著性
                }
            })
          ggb <- do.call(ggpubr::ggarrange,
                         c(gs, list(ncol = 1, align = "v", legend = "none")))
          OESingleCell::save_ggplots(plot = ggb,
                                     file.path(output_dir, "gene_module_violin_plot"),
                                     limitsize =F,
                                     width = 1 * length(unique(cellmeta[[opt$groupby]])),
                                     # 如果all:all，那么height判断
                                     height = ifelse(is.null(opt$pvalue),length(names(gset_list)) * 3,length(names(gset_list)) * 5 ))
        }
      })
    }
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}


#top10marker转gmt
# marker=read.table("top10_markers_for_each_cluster.xls",
#                   sep="\t",quote="",header=1)%>%.[,c('gene',"cluster")]%>%unstack()
#
# lapply(1:length(marker),function(x){
#     print(x)
#     vom=c(names(marker)[x]," ",marker[[x]])%>%as.data.frame() %>% t()
#     write.table( vom,file = "top10_markers_for_each_cluster.gmt",
#                 sep = "\t",row.names = F,col.names = F,quote = F, append = TRUE)
#
#
#
# })
