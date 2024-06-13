## 1.函数定义============================================================================================================
data_concater2 <- function(x) {
    x <- names(sort(table(x), decreasing = T))
    paste(x[1])
}

data_concater <- function(x) {
    x <- levels(factor(x))
    paste(x, collapse = "+")
}

###  对compass分数进行预处理==================
# 1）log10(c+1)
# 2）过滤掉变化不大的反应（min_range=1e-3）
# 3）取相对量(相对于整体所有反应的最小通量)
#' @return   celltype's plot color for this celltype-celltype interaction  matrix
#' @param df a matrix  compass score for all cells
#' @param min.range min cutoff for reaction, default=1e-3
rxn.consistency <- function(df, min.range = 1e-3) {
    df <- -log(df + 1)
    df <- df[apply(df, MARGIN = 1, FUN = max, na.rm = T) - apply(df, MARGIN = 1, FUN = min, na.rm = T) >= min.range,]
    df <- df - min(apply(df, MARGIN = 1, FUN = min, na.rm = T))
    return(df)
}

## cohens_d效应量统计函数=====================
#' @param x compass score matrix A
#' @param y compass score matrix  B
#' @return cohens_d matrix
cohens_d <- function(x, y) {
    pooled_std <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))
    return((mean(x) - mean(y)) / pooled_std)
}

## wilcox.test统计检验======================
#' @param  df  compass score matrix
#' @param groupA  group A barcodes
#' @param groupB  group B barcodes
#' @return wilcox test results
my.wilcox <- function(df, groupA, groupB) {
    results <- matrix(ncol = 3, nrow = nrow(df))
    for (i in 1:nrow(df)) {
        groupA.val <- as.numeric(df[i, groupA])
        groupB.val <- as.numeric(df[i, groupB])
        obj <- try(wilcox.test(groupA.val, groupB.val, alternative = "two.sided"))
        if (!is(obj, "try-error")) {
            results[i, 1] <- obj$statistic
            results[i, 2] <- obj$p.value
            results[i, 3] <- cohens_d(groupA.val, groupB.val)
        }
    }
    results <- data.frame(results)
    rownames(results) <- rownames(df)
    colnames(results) <- c("wilcox_stat", "wilcox_pval", "cohens_d")
    results$adjusted_pval <- p.adjust(results$wilcox_pval, method = "BH")
    return(results)
}

### 氨基酸代谢通路汇总======================================
amino_acid_metab <- c("Alanine and aspartate metabolism",
                      "Arginine and Proline Metabolism",
                      "beta-Alanine metabolism",
                      "Cysteine Metabolism",
                      "D-alanine metabolism",
                      "Folate metabolism",
                      "Glutamate metabolism",
                      "Glycine, serine, alanine and threonine metabolism",
                      "Histidine metabolism",
                      "Lysine metabolism",
                      "Methionine and cysteine metabolism",
                      "Taurine and hypotaurine metabolism",
                      "Tryptophan metabolism",
                      "Tyrosine metabolism",
                      "Urea cycle",
                      "Valine, leucine, and isoleucine metabolism")
### subparsers
docstring <- " example1:\\n\\n\\
sctool -i  sc.rds -f rds --assay RNA -o result compass -x ~new_clusters --pvalue 0.05 --species mus_musculus --labels input_data/rxn.names.csv --metacell False --group.by 'new_clusters,sampleid'  --contrast new_clusters:all:1"
sub_compass <- subparsers$add_parser("compass", help = "differential test between specified groups for compass")
# sub_diffexp$add_argument("--modules",
#                          type = "character",
#                          choices=c("prepare","downstream" ), help = "选择模块:1)准备matrix矩阵，用于compass分析 2）进行下游分析，如PCA和差异筛选")
#sub_compass$add_argument("--downsample", type = "character",  default = NULL, help = "")
sub_compass$add_argument("--targetN", type = "integer", default = NULL, help = "the number of cell to keep after downsample.")
sub_compass$add_argument("--reaction_metadata", type = "character", default = "input_data/database_RECON2/reaction_metadata.csv", help = "Recon2数据库meta信息.")
sub_compass$add_argument("--group.by", type = "character", default = "clusters",
                         help = "The grouppinig variable in the metadata for groupping cells in the PCA plot.")
sub_compass$add_argument("-x", "--design", type = "character", default = NULL,
                         help = "The group design for cell clusters or samples to make differential expression analysis.[default: %(default)s]")
sub_compass$add_argument("-M", "--addition_metadata", type = "character", default = NULL,
                         help = "[Optional]additional metadata for each sample which includes sample id and additional sample groupping info.[default: %(default)s]")
sub_compass$add_argument("-c", "--contrast", type = "character", default = NULL,
                         help = "[Optional]levels of a factor used to compare with for final differenetial results. The format is Factor:interesting_levle:reference_level.[default: %(default)s]")
sub_compass$add_argument("-p", "--adjusted_pval", type = "double", default = 0.1,
                         help = "the adjusted_pval of the reaction differential.[default: %(default)s]", metavar = "P-value")
sub_compass$add_argument("--cohens_d", type = "double", default = 0.2,
                         help = "the cohens_d of the reaction differential .[default: %(default)s]")
sub_compass$add_argument("--species", type = "character", help = "e.g. homo_sapiens,mus_musculus")
#sub_compass$add_argument("--metacell", type = "character", help = "e.g. TRUE or FALSE")
sub_compass$add_argument("--labels", type = "character", help = "rxn labels in MAplot", default = "input_data/rxn.names.csv")
#Subcmd: compass analysis for single cell RNA ================================================================
args <- commandArgs(TRUE)
if ("compass" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "compass") {
        futile.logger::flog.info("step1:加载seurat对象，数据筛选处理") ##======================================================
        suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                        informat = opt$informat,
                                                        assays = assays[1],
                                                        data.use = dataslots,
                                                        verbose = F
        ))
        Seurat::DefaultAssay(data_ob) <- assays[1]
        if (!is.null(opt$predicate)) {
            predicate <- opt$predicate
            futile.logger::flog.info(glue::glue("提供了提取细胞参数--predicate {predicate}，将按照要求提取特定细胞"))
            df <- slot(data_ob, "meta.data")
            desired_cells <- subset(df, eval(parse(text = opt$predicate)))
            data_ob <- subset(data_ob, cells = rownames(desired_cells))
        }
        if (!is.null(opt$targetN)) {
            targetN <- opt$targetN
            futile.logger::flog.info(glue::glue("针对数据进行降采样,目标细胞数目{targetN}"))
            data_ob <- OESingleCell::DownSample(data_ob,
                                                opt$targetN,
                                                slot = dataslots,
                                                use_PCs = TRUE,
                                                PCs = 100,
                                                k = "auto",
                                                seed = 4563,
                                                replace = FALSE,
                                                alpha = 0.1,
                                                max_iter = 200,
                                                verbose = F) }
        futile.logger::flog.info("step2:运行compass软件，每个细胞所需运行资源约为0.5 cpu.h") #===================================
        if (!file.exists(glue::glue("{output_dir}/1.compass_output/reactions.tsv"))) {
            expression <- Seurat::GetAssayData(object = data_ob, assay = assays[1], slot = "data") %>%
                    as.data.frame %>%
                    tibble::rownames_to_column("Gene Symbols")
            readr::write_tsv(expression, glue::glue("{output_dir}/expression.tsv"))
            metadata <- data_ob@meta.data %>%
                    tibble::rownames_to_column("barcode") %>%
                    dplyr::select("barcode", "clusters", "group")
            readr::write_tsv(metadata, glue::glue("{output_dir}/metadata.tsv"))
            # if (opt$metacell) {
            #   cmd <- glue::glue('compass ',
            #                     ' --data {output_dir}/expression.tsv ',
            #                     ' --num-processes {opt$ncores}',
            #                     ' --microcluster-size 10 ',
            #                     ' --output-dir {output_dir}/1.compass_output',
            #                     ' --species {opt$species}')
            #   futile.logger::flog.info(cmd)
            #   system(cmd)
            # }else {
            cmd <- glue::glue('module load OESingleCell/dev && compass ',
                              ' --data {output_dir}/expression.tsv ',
                              ' --num-processes {opt$ncores}',
                              ' --output-dir {output_dir}/1.compass_output',
                              ' --species {opt$species}', " && module unload OESingleCell/dev")
            futile.logger::flog.info(cmd)
            system(cmd)
            #}
        }else {
            futile.logger::flog.info(glue::glue("经检测,compass运行结果已存在,跳过：{output_dir}/1.compass_output/reactions.tsv "))
        }
        futile.logger::flog.info("step3:读入compass分数、分组信息以及反应标注信息") #=============================================
        reaction.penalties <- data.table::fread(glue::glue("{output_dir}/1.compass_output/reactions.tsv"), na.strings = c("NA", "N/A")) %>% tibble::column_to_rownames("V1")
        reaction.meta <- data.table::fread(opt$reaction_metadata, na.strings = c("NA", "N/A"))
        rxn.labels.names <- readr::read_csv(opt$labels)
        #rxn.labels.names<-read_csv("input_data/rxn.names.csv")
        #get rid duplicates
        rxn.labels.names <- unique(rxn.labels.names)
        #only keep first occurance of remaining duplicates
        rxn.labels.names <- rxn.labels.names[match(unique(rxn.labels.names$rxn), rxn.labels.names$rxn),]

        futile.logger::flog.info("step4::对compass分数进行预处理") #==========================================================
        reaction.consistencies <- rxn.consistency(reaction.penalties)

        futile.logger::flog.info("step5:分别利用正反应及负反应，进行主成分分析, 需要注意的是该过程中利用的是核心反应") #===============
        ## Filter core reactions
        core.rxns <- reaction.meta %>%
                dplyr::filter(confidence %in% c(0, 4) & !is.na(EC_number))
        core.reaction.consistencies <- subset(reaction.consistencies,
                                              rownames(reaction.consistencies) %in% paste0(core.rxns$reaction_no_direction, "_pos") |
                                                      rownames(reaction.consistencies) %in% paste0(core.rxns$reaction_no_direction, "_neg"))
        ## run PCA
        lapply(c("_pos", "_neg", "_"), function(pos_neg) {
            df2 <- core.reaction.consistencies %>%
                    t() %>%
                    as.data.frame() %>%
                    dplyr::select(matches(pos_neg))
            if (pos_neg == "_") { pos_neg <- "all" }
            res.pca <- FactoMineR::PCA(df2, graph = FALSE)
            var <- factoextra::get_pca_var(res.pca) #查看样本的主成分分析结果
            plot_screplot <- factoextra::fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 60)) #绘制主成分碎石图，查看每一个主成分能在多大程度上代表原来的特征
            OESingleCell::save_ggplots(glue::glue("{output_dir}/2.pca/screeplot_for{pos_neg}_reaction"),
                                       plot = plot_screplot,
                                       width = 10,
                                       height = 8)
            for (group_by in  stringr::str_split(opt$group.by, ",")[[1]]) {
                metadata <- data_ob@meta.data %>%
                        tibble::rownames_to_column("barcode") %>%
                        dplyr::select(group_by)
                plot_pca <- factoextra::fviz_pca_ind(res.pca,
                                                     mean.point = F,
                                                     label = "none",
                                                     addEllipses = TRUE,
                                                     habillage = metadata[, group_by] %>% factor(),
                                                     palette = OESingleCell::SelectColors(n = metadata[, group_by] %>% unique %>% length,
                                                                                          palette = "blindless")) +
                        theme_bw(base_size = 20) +
                        theme(plot.title = element_text(hjust = 0.5),
                              strip.background = element_rect(colour = "black", fill = "white", linetype = "solid"),
                              panel.border = element_rect(colour = "black", size = 1, fill = NA),
                              panel.background = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.key = element_rect(fill = "white", colour = NA),
                              legend.title = element_blank())
                OESingleCell::save_ggplots(
                        glue::glue("{output_dir}/2.pca/PCA_group.by_{group_by}_for{pos_neg}_reaction"),
                        plot = plot_pca,
                        width = 10,
                        height = 8)
            }
        })
        futile.logger::flog.info("stpe6:进行wilcoxn检验") #==================================================================
        assay_metadata <- OESingleCell::colData(data_ob)
        if (is.null(opt$contrast)) { #no contrast is provided
            factors_indesign <- strsplit(opt$design, "[~+ ]+", perl = T)
            last_factor_indesign <- factors_indesign[length(factors_indesign)]
            if (is.null(assay_metadata[, last_factor_indesign])) {
                stop("The factor in design formula does not exist in assay metadata.")
            }
            variable_levels <- levels(assay_metadata[, last_factor_indesign])
            contrast <- paste(last_factor_indesign, variable_levels[length(variable_levels)], variable_levels[1], sep = ":")
        }else {
            contrast <- opt$contrast
        }
        contrasts <- unlist(strsplit(contrast, ":", perl = T))
        all_levels <- as.vector(unique(assay_metadata[, contrasts[1]]))
        if (contrasts[2] == "all" & contrasts[3] != "all") {
            all_levels <- all_levels[-which(all_levels == contrasts[3])] #delete the reference level
            all_comparisions <- paste(contrasts[1], all_levels, contrasts[3], sep = ":")
        }else if (contrasts[2] == "all" & contrasts[3] == "all") {
            all_comparisions <- lapply(all_levels,
                                       function(x) paste(contrasts[1], x, paste0(all_levels[-which(all_levels == x)], collapse = ","), sep = ":"))
            all_comparisions <- unlist(all_comparisions)
        }else if (contrasts[2] != "all" & contrasts[3] == "all") {
            #delete the interested level in the  reference level
            ref_levels <- paste0(all_levels[-which(all_levels == contrasts[2])], collapse = ",")
            all_comparisions <- paste(contrasts[1], contrasts[2], ref_levels, sep = ":")
        }else {
            all_comparisions <- contrast
        }
        future.apply::future_lapply(all_comparisions, function(j) {
            futile.logger::flog.info("1.运行wilcox.test") #==========================
            group <- stringr::str_split(j, ":")[[1]][1]
            case <- stringr::str_split(j, ":")[[1]][2]
            control <- stringr::str_split(j, ":")[[1]][3]
            group1 <- assay_metadata %>%
                    dplyr::filter(!!rlang::sym(group) == case) %>%
                    rownames
            group2 <- assay_metadata %>%
                    dplyr::filter(!!rlang::sym(group) == control) %>%
                    rownames
            wilcox.results <- my.wilcox(reaction.consistencies,
                                        group1,
                                        group2)
            futile.logger::flog.info("2.添加各反应的meta信息") #======================
            wilcox.results$rxn <- rownames(wilcox.results)
            for (i in 1:nrow(wilcox.results)) {
                r <- wilcox.results$rxn[i]
                if (r %in% reaction.meta$reaction_no_direction) {
                    wilcox.results$metadata_r_id[i] <- r
                }else if (substr(r, 1, nchar(r) - 4) %in% reaction.meta$reaction_no_direction) {
                    wilcox.results$metadata_r_id[i] <- substr(r, 1, nchar(r) - 4)
                }

                W <- dplyr::left_join(wilcox.results, reaction.meta, by = c("metadata_r_id" = "reaction_no_direction"))
                rownames(W) <- W$rxn
            }

            diff_dir <- glue::glue("{output_dir}/3.diff_compass/{group}_{case}-vs-{control}/volcano/")
            if (!file.exists(diff_dir)) { dir.create(diff_dir, recursive = TRUE) }
            write.csv(W, file = glue::glue("{output_dir}/3.diff_compass/{group}_{case}-vs-{control}/{group}_{case}-vs-{control}-all-reactions.csv"))
            futile.logger::flog.info("3:对反应进行过滤") #=============================
            #  Next we join the metadata to the reactions in a new dataframe W, so that we can filter out non-core reactions.
            #  More specifically, we remove reactions with confidence other than 0 or 4(4：most confident; 0：unassigned confidence) and filter
            #  out reactions in the citric acid cycle subsystem which are outside of the mitochondria. W<-W[(W$confidence %in% c(0,4)) , ]
            W <- W %>%
                    dplyr::mutate(subsystem = as.character(subsystem)) %>%
                    dplyr::filter(confidence %in% c(0, 4) & !is.na(EC_number)) %>%
                    dplyr::mutate(subsystem = ifelse(stringr::str_detect(subsystem, "Citric acid cycle") & !stringr::str_detect(formula, "\\[m\\]"),
                                                     'other',
                                                     subsystem))
            write.csv(W, file = glue::glue("{output_dir}/3.diff_compass/{group}_{case}-vs-{control}/{group}_{case}-vs-{control}-all-reactions-filtered.csv"))
            futile.logger::flog.info("4.保存结果、绘制热图、绘制火山图") #=================================================================
            #===================================保存结果
            futile.logger::flog.info("----保存差异筛选结果")
            data <- W %>% dplyr::mutate(subsystem.select = ifelse(subsystem %in% amino_acid_metab, "Amino Acid Metabolism", subsystem))
            case <- stringr::str_split(j, ":")[[1]][2]
            control <- stringr::str_split(j, ":")[[1]][3]
            data$logP <- -log10(data$adjusted_pval)
            data <- data[order(data$adjusted_pval),]
            for (i in 1:nrow(data)) { ##逐行处理
                if ( # data$rxn[i] %in% rxn.names$rxn &
                        !is.na(data$adjusted_pval[i]) &
                                data$adjusted_pval[i] < opt$adjusted_pval &
                                abs(data$cohens_d[i]) > opt$cohens_d) {
                    #select significant values to label 显著性差异padjust值小于0.1,并且cohen's效应量值大于0.2
                    #if(!(rxn.names$label[rxn.names$rxn==data$rxn[i]] %in% data$label)){ #only label each rxn once
                    if (!(data$reaction_name[i] %in% data$label)) {
                        #  data$label[i]<-rxn.names$label[rxn.names$rxn==data$rxn[i]]
                        data$label[i] <- data$reaction_name[i]
                    }else {
                        data$label[i] <- "" }
                }else {
                    data$label[i] <- ""
                }
            }
            padj <- opt$adjusted_pval
            cohens_d <- opt$cohens_d
            write.csv(data, file = glue::glue("{output_dir}/3.diff_compass/{group}_{case}-vs-{control}/{group}_{case}-vs-{control}-diff-padj-{padj}-cohens_d-{cohens_d}.csv"))
            #===================================生成火山图=====================================
            futile.logger::flog.info("----绘制热图")
            metadata <- data_ob@meta.data %>% tibble::rownames_to_column("barcode") %>% dplyr::select("barcode", "clusters", "group")
            select.systems <- data[data$subsystem %in% names(table(data$subsystem))[table(data$subsystem) > 5],] %>%
                    dplyr::filter(adjusted_pval < opt$adjusted_pval & abs(cohens_d) > opt$cohens_d) %>%
                    .$subsystem %>%
                    unique()
            data <- data %>%
                    dplyr::mutate(group = ifelse(cohens_d > 0, case, control)) %>%
                    dplyr::mutate(col = dplyr::case_when(logP > 1 & cohens_d > 0 ~ "red",
                                                         logP > 1 & cohens_d < 0 ~ "blue",
                                                         TRUE ~ "grey")) %>%
                    dplyr::mutate(col = as.factor(col))

            colors <- c("lightblue", "grey", "pink")
            names(colors) <- c("blue", "grey", "red")
            for (name in select.systems) {
                r.name <- sub("/", " ", name)
                ## 数据整理
                select_rec <- data %>%
                        dplyr::filter(subsystem  == name) %>%
                        rownames
                select_rec_score <- reaction.consistencies %>%
                        dplyr::filter(rownames(reaction.consistencies) %in% select_rec) %>%
                        t() %>%
                        scale(center = TRUE, scale = TRUE) %>%
                        t()
                sc_group_name <- metadata %>%
                        dplyr::filter(as.character(barcode) %in% colnames(select_rec_score)) %>%
                        dplyr::arrange(!!!syms(c(group)))%>%
                        dplyr::pull(!!sym(group)) %>%
                        as.character()
                score_range <- select_rec_score %>%
                        as.data.frame %>%
                        dplyr::summarise_all(max) %>%
                        quantile %>%
                        dplyr::select(`75%`) %>%
                        dplyr::mutate(`75%` = `75%` / 10) %>%
                        ceiling() * 10
                col_fun <- circlize::colorRamp2(c(-(score_range$`75%`), 0, score_range$`75%`), c("navyblue", "white", "red"))
                ##设置列注释
                clusters.colors <- structure(names = as.character(unique(sc_group_name)),
                                             OESingleCell::SelectColors(1:length(unique(sc_group_name)), palette = "ditto" ))
                colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(clusters = as.factor(sc_group_name)),
                                                            col = list(clusters = clusters.colors),
                                                            show_annotation_name = T,
                                                            annotation_legend_param = list(clusters = list(ncol = 1, direction = "horizontal"))
                )
                ##绘图
                pdf(NULL)
                plot_heatmap <- ggplotify::as.ggplot(grid::grid.grabExpr(
                        ComplexHeatmap::draw(
                                ComplexHeatmap::Heatmap(as.matrix(select_rec_score),
                                                        col = col_fun,
                                                        cluster_rows = T,
                                                        cluster_columns = T,
                                                        show_column_names = F,
                                                        column_title_side = "top",
                                                        show_row_names = T,
                                                        row_names_gp = grid::gpar(fontsize = 10),
                                                        row_names_rot = 0,
                                                        row_names_side = "left",
                                                        heatmap_legend_param = list(
                                                            title = group,
                                                            # title_position = "topcenter", # 图例标题位置
                                                            # at = c(0, 1), # 图例范围
                                                            legend_direction = "vertical",
                                                            legend_height = unit(2, "cm") # 图例长度
                                                        ),
                                                        row_dend_width = unit(0.8, "mm"),
                                                        top_annotation = colAnn
                                ),
                                merge_legend = FALSE,
                                heatmap_legend_side = "right",
                                annotation_legend_side = "top"
                        )))
                OESingleCell::save_ggplots(plot = plot_heatmap,
                                           filename = glue::glue("{output_dir}/3.diff_compass/{group}_{case}-vs-{control}/heatmap/{r.name}"),
                                           width = 7,
                                           height = dim(select_rec_score)[1] * 0.15)
                ##=======================================================================ss======================================
                futile.logger::flog.info("----绘制火山图")
                plot_volcano <- ggplot(data[data$subsystem == name,], aes(x = cohens_d, y = logP)) +
                        geom_point(data = data[data$subsystem == name,],
                                   aes(x = cohens_d,
                                       y = logP,
                                       color = col),
                                   #color=colors[name][[name]],
                                   size = 1) +
                        scale_color_manual(values = colors[levels(data$col)]) +
                        geom_hline(yintercept = 1, linetype = "dashed", size = 0.4) +
                        geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
                        xlim(-1, 1) +
                        labs(title = name,
                             x = glue::glue(" Cohen's d ({case}-vs-{control})"),
                             y = expression(-log[10]("Wilcoxon-adjusted p"))) +
                        cowplot::theme_cowplot() +
                        # theme(legend.title=element_blank())+
                        # guides(color=guide_legend(override.aes=list(size=1)))+
                        #scale_color_brewer(palette="Set1")+
                        theme(plot.title = element_text(hjust = 0.5),
                              legend.position = "none",
                              panel.grid = element_blank(),
                              panel.background = element_blank()) +
                        # theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")) +
                        ggrepel::geom_text_repel(size = 4,
                                                 aes(label = label),
                                                 box.padding = 0.5,
                                                 segment.size = 0.2,
                                                 max.overlaps = Inf,
                                                 force = 50) #, ,box.padding = 0.005,subsystem
                OESingleCell::save_ggplots(filename = glue::glue("{output_dir}/3.diff_compass/{group}_{case}-vs-{control}/volcano/{r.name}"),
                                           plot = plot_volcano,
                                           width = 8,
                                           height = 7)
                #=============================================================================================================
                reaction.consistencies
            }

            futile.logger::flog.info("5: 绘制点图") #=============================
            data <- W %>%
                    dplyr::filter(!subsystem %in% c("Miscellaneous", "Unassigned", "Other")) %>%
                    dplyr::filter(!stringr::str_detect(subsystem, "Transport|Exchange"))
            data <- data[data$subsystem %in% names(table(data$subsystem))[table(data$subsystem) > 5],] # only keep subsystems with at least 5 cells
            data <- data %>% dplyr::mutate(sig.level = ifelse(adjusted_pval < opt$adjusted_pval,
                                                              glue::glue("BH-adjusted p<{opt$adjusted_pval}"),
                                                              "No Significant"),
                                           up.down = ifelse(cohens_d >= 0,
                                                            "up",
                                                            "down"))
            plot_dotplot <- ggplot(data = data,
                                   aes(x = cohens_d, y = reorder(subsystem, cohens_d, FUN = median))) +
                    geom_point(mapping = aes_string(color = "up.down", alpha = "sig.level")) +
                    geom_vline(xintercept = 0, linetype = "dotted") +
                    scale_alpha_manual(values = c(`No Significant` = 0.1, `BH-adjusted p<0.1` = 1)) +
                    scale_color_manual(values = c(up = "red", down = "blue")) +
                    theme_bw() +
                    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
                    labs(x = "Cohen's d", y = "", title = gsub("W.", "", j)) +
                    guides(color = F, alpha = guide_legend(title = ""))
            OESingleCell::save_ggplots(filename = glue::glue("{output_dir}/3.diff_compass/{group}_{case}-vs-{control}/dotplot"),
                                       plot = plot_dotplot,
                                       width = 10,
                                       height = 7)
        }, future.seed = 2020)
        #==============================================================================================================
        ## save session information
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}
