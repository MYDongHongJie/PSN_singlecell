#' Scissor: Single-Cell Identification of Subpopulations with bulk Sample phenOtype coRrelation
#'
#' \code{Scissor} is a novel approach that utilizes the phenotypes, such as disease stage, tumor metastasis, treatment response, and survival
#' outcomes, collected from bulk assays to identify the most highly phenotype-associated cell subpopulations from single-cell data.
#'
#' Scissor is a novel algorithm to identify cell subpopulations from single-cell data that are most highly associated with the given phenotypes.
#' The three data sources for Scissor inputs are a single-cell expression matrix, a bulk expression matrix, and a phenotype of interest.
#' The phenotype annotation of each bulk sample can be a continuous dependent variable, binary group indicator vector, or clinical survival data.
#' The key step of Scissor is to quantify the similarity between the single-cell and bulk samples by Pearson correlation for each pair of cells and bulk samples.
#' After this, Scissor optimizes a regression model on the correlation matrix with the sample phenotype.
#' The selection of the regression model depends on the type of the input phenotype, i.e., linear regression for continuous variables,
#' logistic regression for dichotomous variables, and Cox regression for clinical survival data.
#' Based on the signs of the estimated regression coefficients, the cells with non-zero coefficients can be indicated as
#' Scissor positive (Scissor+) cells and Scissor negative (Scissor-) cells, which are positively and negatively associated
#' with the phenotype of interest, respectively.
#'
#' @param bulk_dataset Bulk expression matrix of related disease. Each row represents a gene and each column represents a sample.
#' @param sc_dataset Single-cell RNA-seq expression matrix of related disease. Each row represents a gene and each column represents a sample.
#' A Seurat object that contains the preprocessed data and constructed network is preferred. Otherwise, a cell-cell similarity network is
#' constructed based on the input matrix.
#' @param phenotype Phenotype annotation of each bulk sample. It can be a continuous dependent variable,
#' binary group indicator vector, or clinical survival data:
#'   \itemize{
#'   \item Continuous dependent variable. Should be a quantitative vector for \code{family = gaussian}.
#'   \item Binary group indicator vector. Should be either a 0-1 encoded vector or a factor with two levels for \code{family = binomial}.
#'   \item Clinical survival data. Should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable,
#'   with '1' indicating event (e.g.recurrence of cancer or death), and '0' indicating right censored.
#'   The function \code{Surv()} in package survival produces such a matrix.
#'   }
#' @param tag Names for each phenotypic group. Used for linear and logistic regressions only.
#' @param alpha Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector.
#' If \code{alpha = NULL}, a default searching vector is used. The range of alpha is in \code{[0,1]}. A larger alpha lays more emphasis on the l1 norm.
#' @param cutoff Cutoff for the percentage of the Scissor selected cells in total cells. This parameter is used to restrict the number of the
#' Scissor selected cells. A cutoff less than \code{50\%} (default \code{20\%}) is recommended depending on the input data.
#' @param family Response type for the regression model. It depends on the type of the given phenotype and
#' can be \code{family = gaussian} for linear regression, \code{family = binomial} for classification, or \code{family = cox} for Cox regression.
#' @param Save_file File name for saving the preprocessed regression inputs into a RData.
#' @param Load_file File name for loading the preprocessed regression inputs. It can help to tune the model parameter \code{alpha}.
#' Please see Scissor Tutorial for more details.
#'
#' @return This function returns a list with the following components:
#'   \item{para}{A list contains the final model parameters.}
#'   \item{Coefs}{The regression coefficient for each cell.}
#'   \item{Scissor_pos}{The cell IDs of Scissor+ cells.}
#'   \item{Scissor_neg}{The cell IDs of Scissor- cells.}
#'
#' @references Duanchen Sun and Zheng Xia (2021): Phenotype-guided subpopulation identification from single-cell sequencing data. Nature Biotechnology.
#' @import Seurat Matrix preprocessCore
#' @export
Scissor <- function(bulk_dataset,
                    sc_dataset,
                    # assay= "RNA",
                    phenotype,
                    tag = NULL,
                    alpha = NULL,
                    cutoff = 0.2,
                    family = c("gaussian","binomial","cox"),
                    Save_file = "Scissor_inputs.RData",
                    Load_file = NULL){
    library(Seurat)
    library(Scissor)
    library(Matrix)
    library(preprocessCore)
    if (is.null(Load_file)){
        common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
        if (length(common) == 0) {
            stop("There is no common genes between the given single-cell and bulk samples.")
        }
        if (class(sc_dataset) == "Seurat"){
          as_matrix <- function(mat){
            tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
            row_pos <- mat@i+1
            col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
            val <- mat@x

            for (i in seq_along(val)){
              tmp[row_pos[i],col_pos[i]] <- val[i]
            }
            row.names(tmp) <- mat@Dimnames[[1]]
            colnames(tmp) <- mat@Dimnames[[2]]
            return(tmp) }

          sc_exprs <- as_matrix(sc_dataset@assays$RNA@data)
          #sc_exprs <- as_matrix(GetAssayData(object = sc_dataset, assay=opt$assay, slot = "data"))
          if(is.null(sc_dataset@graphs$RNA_snn)){
            sc_dataset <- Seurat::FindVariableFeatures(sc_dataset, selection.method = "vst", verbose = F)
            sc_dataset <- Seurat::ScaleData(sc_dataset, verbose = F)
            sc_dataset <- Seurat::RunPCA(sc_dataset, features = VariableFeatures(sc_dataset), verbose = F)
            sc_dataset <- Seurat::FindNeighbors(sc_dataset, dims = 1:10, verbose = F)
          }
          network  <- as_matrix(sc_dataset@graphs$RNA_snn)
          #network  <- as_matrix(sc_dataset@graphs[[glue::glue("{assay}_snn")]])
        }else{
            sc_exprs <- as.matrix(sc_dataset)
            Seurat_tmp <- Seurat::CreateSeuratObject(sc_dataset)
            Seurat_tmp <- Seurat::FindVariableFeatures(Seurat_tmp, selection.method = "vst", verbose = F)
            Seurat_tmp <- Seurat::ScaleData(Seurat_tmp, verbose = F)
            Seurat_tmp <- Seurat::RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), verbose = F)
            Seurat_tmp <- Seurat::FindNeighbors(Seurat_tmp, dims = 1:10, verbose = F)
            network  <- as.matrix(Seurat_tmp@graphs$RNA_snn)

        }

        diag(network) <- 0
        network[which(network != 0)] <- 1

        dataset0 <- as.matrix(cbind(bulk_dataset[common,], sc_exprs[common,]))     # Dataset before quantile normalization.
        dataset1 <- preprocessCore::normalize.quantiles(dataset0)                           # Dataset after  quantile normalization.
        rownames(dataset1) <- rownames(dataset0)
        colnames(dataset1) <- colnames(dataset0)

        Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
        Expression_cell <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
        X <- cor(Expression_bulk, Expression_cell)

        quality_check <- quantile(X)
        print("|**************************************************|")
        print("Performing quality-check for the correlations")
        print("The five-number summary of correlations:")
        print(quality_check)
        print("|**************************************************|")
        if (quality_check[3] < 0.01){
            warning("The median correlation between the single-cell and bulk samples is relatively low.")
        }
        if (family == "binomial"){
            Y <- as.numeric(phenotype$status)
            z <- table(Y)
            if (length(z) != length(tag)){
                stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
            }else{
                print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
                print("Perform logistic regression on the given phenotypes:")
            }
        }
        if (family == "gaussian"){
            print("OK")
            Y <- as.numeric(phenotype$status)
            z <- table(Y)
            print(Y)
            print(tag)
            if (length(z) != length(tag)){
                stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
            }else{
                tmp <- paste(z, tag)
                print(paste0("Current phenotype contains ", paste(tmp[1:(length(z)-1)], collapse = ", "), ", and ", tmp[length(z)], " samples."))
                print("Perform linear regression on the given phenotypes:")
            }
        }
        if (family == "cox"){
            Y <- as.matrix(phenotype)
            if (ncol(Y) != 2){
                stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
            }else{
                print("Perform cox regression on the given clinical outcomes:")
            }
        }
        save(X, Y, network, Expression_bulk, Expression_cell, file = Save_file)
    }else{
        load(Load_file)
    }
    print(alpha)
    #for (i in 1:length(alpha)){
    #infos <- parallel::mclapply(1:length(alpha), function(i){
    set.seed(123)
    fit0 <- Scissor::APML1(X, Y, family = family, penalty = "Net", alpha = alpha , Omega = network, nlambda = 100, nfolds = min(10,nrow(X)))
    fit1 <- Scissor::APML1(X, Y, family = family, penalty = "Net", alpha = alpha , Omega = network, lambda = fit0$lambda.min)
    if (family == "binomial"){
        Coefs <- as.numeric(fit1$Beta[2:(ncol(X)+1)])
    }else{
        Coefs <- as.numeric(fit1$Beta)
    }
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2)) / ncol(X)
    print(sprintf("alpha = %s", alpha))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", formatC(percentage*100, format = 'f', digits = 3)))

    # if (percentage < cutoff){
    #     break
    # }
    cat("\n")
    infos <-list(para=list(alpha = alpha,
                           lambda = fit0$lambda.min,
                           family = family),
                 Coefs = Coefs,
                 Scissor_pos = Cell1,
                 Scissor_neg = Cell2)

    # },mc.cores =4)
    print("|**************************************************|")
    return(infos)
    # return(list(para = list(alpha = alpha_param[i], lambda = fit0$lambda.min, family = family),
    #             Coefs = Coefs,
    #             Scissor_pos = Cell1,
    #             Scissor_neg = Cell2))
}

docstring<- " example1(without refmarker):\\n\\n\\
  sctool  -i  query.rds   -f rds  -o results  -d rds   --assay RNA   scissor  --bulkexp demo.fpkm.txt  --phenotype  demo.phenotype.txt  --tag phenotype1 --family gaussian  --alpha  0.07  --cutoff 0.2  "

sub_scissor<- subparsers$add_parser("scissor",
                                       description = docstring,
                                       formatter_class= 'argparse.RawTextHelpFormatter' ,
                                       #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                       argument_default = "True" ,
                                       help="Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data")
sub_scissor$add_argument("--bulkexp",
                           type = "character",
                           help = "[REQUIRED] Bulk expression matrix of related disease. The gene name colnames is 'gene_id'.And each row represents a gene and each column represents a sample." )
sub_scissor$add_argument("--phenotype",
                           type = "character",
                           help = paste0("Phenotype annotation of each bulk sample. It can be a continuous dependent variable, binary group indicator vector, or clinical survival data:",
                                         "1) Continuous dependent variable. Should be a quantitative vector for family = gaussian.",
                                         "2) Binary group indicator vector. Should be either a 0-1 encoded vector or a factor with two levels for family = binomial.",
                                         "3) Clinical survival data. Should be a two-column matrix with columns named 'time' and 'status'.",
                                         "The latter is a binary variable, with '1' indicating event (e.g.recurrence of cancer or death), and '0' indicating right censored. ",
                                         "The function Surv() in package survival produces such a matrix." ))
sub_scissor$add_argument("--tag",
                           type = "character",
                           help = "the selected phenotype's colname in phenotype data ")
sub_scissor$add_argument("--family",
                           type = "character",
                           help = "gaussian,binomial,cox")
sub_scissor$add_argument("--alpha",
                         default=0.05,
                         help =paste0("Parameter used to balance the effect of the l1 norm and the network-based penalties.",
                                         "It can be a number or a searching vector. If alpha = NULL, a default searching vector is used. ",
                                         "The range of alpha is in [0,1]. A larger alpha lays more emphasis on the l1 norm."))
sub_scissor$add_argument("--cutoff",
                         type = "double",
                         default = 0.2,
                         help =paste0("Cutoff for the percentage of the Scissor selected cells in total cells. ",
                                         "This parameter is used to restrict the number of the Scissor selected cells. ",
                                         "A cutoff less than 0.5 (default 0.2) is recommended depending on the input data[default: %(default)s]"))
sub_scissor$add_argument("--test",
                         default = TRUE,
                         help = " test or not [default: %(default)s]")
sub_scissor$add_argument("--palette",
                          type = "character",
                          default = "grey,lightblue,orange",
                          help = paste0(
                              "the discrete color schema mapped to the cell annotations specified by --groupby.[default: %(default)s]",
                              " Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,",
                              " tableau20:20,Buen:17,UKBB:18,TF1:17,paired:12"))
sub_scissor$add_argument("--groupby",
                          type = "character",
                          default = "clusters",
                          help = paste0(
                              "[OPTIONAL]The grouppinig variable in the metadata for groupping cells in the dimplot.",
                              "If not specified, %(default)s will be used as default."))
# ================ Subcmd: st_deconv, deconvolute the cell type composition for spatial transcriptomics or bulk transcriptomics using scRNA-seq data ========
args <- commandArgs(TRUE)
if ( "scissor"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name == "scissor"){
    #===================================================================================================================
    futile.logger::flog.info("step1:read the specified assay and data slot in data object into memory")
    suppressMessages(sc_dataset <- OESingleCell::ReadX(input = opt$input,
                                                       informat = opt$informat,
                                                       assays = assays[1],
                                                       data.use = dataslots,
                                                       verbose = F))
    Seurat::DefaultAssay(sc_dataset) <- assays[1]
     # get the subset of cells used for visualization if necessay subset cells on the expression matrix using logical exprssion on features
     if (!is.null(opt$predicate)) {
            #示例： predicate sampleid %in% c(\"A\",\"B\")
            futile.logger::flog.info("get the subset of cells used for visualization if necessay subset cells on the expression matrix using logical exprssion on features")
            df <- OESingleCell::colData(sc_dataset)
            desired_cells <- subset(df, eval(parse(text = opt$predicate)))
            sc_dataset <- sc_dataset[, rownames(desired_cells)]
        }
    #===================================================================================================================
    futile.logger::flog.info("step2:load bulk RNAseq express and phenotype data")#=================================
    bulk_dataset<- data.table::fread(opt$bulkexp,data.table = F) %>% tibble::column_to_rownames("gene_id")
    phenotype_all <-  data.table::fread(opt$phenotype,data.table = F)
    if(all(colnames(bulk_dataset) == phenotype_all$sample)){
      futile.logger::flog.info("表型数据与bulk表达矩阵中的列名顺序相同")
    }else{
      futile.logger::flog.info("表型数据与bulk表达矩阵中的列名顺序不相同")
      bulk_dataset<-bulk_dataset%>%dplyr::select(phenotype_all$sample)
    }
    if(opt$family=="gaussian" | opt$family=="binomial"){
      phenotype<- phenotype_all %>% dplyr::select(opt$tag)%>% dplyr::rename(status=opt$tag)
    }else if(opt$family=="cox"){
      phenotype <- phenotype_all %>% dplyr::select("time","status")
    }
    #====================================================================================================================
    futile.logger::flog.info("step3:run scissor")#=====================
    #if (is.null(opt$alpha)){
        #alpha_param <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
    #     alpha_param <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    #}else{
   # alpha_param <- stringr::str_split(opt$alpha,",") %>% unlist%>%as.numeric()
    #}
    alpha_param<- opt$alpha %>% as.numeric()
    if(! file.exists(glue::glue("{output_dir}/Scissor_{opt$tag}_infos.RData"))){
         # Load_file<- ifelse(file.exists(glue::glue("{output_dir}/Scissor_{opt$tag}.RData")),
         #                    glue::glue("{output_dir}/Scissor_{opt$tag}.RData"),
         #                    NULL)
        infos <- Scissor( bulk_dataset,
                          sc_dataset,
                          #assay=opt$assay,
                          phenotype,
                          tag =  names(table(phenotype$status)),
                          alpha = alpha_param,
                          cutoff= opt$cutoff,
                          family = opt$family,
                          Save_file = glue::glue("{output_dir}/Scissor_{opt$tag}.RData"))#,
                          #Load_file = Load_file)
        if (!file.exists(  glue::glue("{output_dir}/{alpha_param}/"))) {
            dir.create(glue::glue("{output_dir}/{alpha_param}/"), recursive = T)
          }
        save( infos, file = glue::glue("{output_dir}/{alpha_param}/Scissor_{opt$tag}_infos.RData"))
    }else{
      futile.logger::flog.info(glue::glue("Skipping, loading existing files:{output_dir}/{alpha_param}/Scissor_{opt$tag}.RData" ))
      load(glue::glue("{output_dir}/{alpha_param}/Scissor_{opt$tag}_infos.RData"))
    }

    futile.logger::flog.info("step4:result visulization")#==================================================
      futile.logger::flog.info(glue::glue("visulization for {alpha_param}"))
      Scissor_select <- rep("Background cells", ncol(sc_dataset))
      names(Scissor_select) <- colnames(sc_dataset)
      Scissor_select[infos$Scissor_pos] <- "Scissor+ cell"
      Scissor_select[infos$Scissor_neg] <- "Scissor- cell"
      sc_dataset <- Seurat::AddMetaData(sc_dataset, metadata = Scissor_select, col.name = glue::glue("scissor_{opt$tag}_{alpha_param}"))
      ### dimplot
      groupby<-stringr::str_split(opt$groupby,",")%>%unlist%>% as.character
      plot2<-list()
      reordered_cell_count_by_scissor <- table(sc_dataset@meta.data[[glue::glue("scissor_{opt$tag}_{alpha_param}")]])
      cell_count_labels <- glue::glue("{names(reordered_cell_count_by_scissor)}-{reordered_cell_count_by_scissor} cells")
      plot1 <- Seurat::DimPlot(sc_dataset,
                      reduction = opt$reduction,
                      group.by = glue::glue("scissor_{opt$tag}_{alpha_param}"),
                      pt.size = 0.1)+
              ggplot2::scale_colour_manual(values = stringr::str_split(opt$palette,",")%>%unlist%>%as.character  ,
                                           breaks = names(reordered_cell_count_by_scissor),
                                           labels = cell_count_labels)
      for(j in  groupby){
        plot2[[j]] <- Seurat::DimPlot(sc_dataset,
                                       reduction = opt$reduction,
                                       group.by = j,
                                       cols=OESingleCell::SelectColors(levels(sc_dataset@meta.data[[j]]), palette ="customecol2") ,
                                       pt.size = 0.1)
      }
      plot<-patchwork::wrap_plots(plot2  ,nrow=1)
      OESingleCell::save_ggplots( plot =patchwork::wrap_plots(plot ,plot1, nrow=1,widths = c(length(groupby),1)),
                                  file=glue::glue("{output_dir}/{alpha_param}/dimplot_for_scissor_{opt$tag}_{alpha_param}"),
                                  limitsize =F,
                                  width = length(groupby)*6+6,
                                  height =5 )
      ## barplot
      p1<-list()
      p2<-list()
      for(j in  1:length(groupby)){
          plot_data<-sc_dataset@meta.data%>%
                 tibble::rownames_to_column("barcodes")%>%
                 dplyr::select(groupby[j], dplyr::starts_with("scissor"))%>%
                 dplyr::rename(scissor_status=glue::glue("scissor_{opt$tag}_{alpha_param}"))%>%
                 dplyr::group_by(!!!rlang::syms(groupby[j]),  scissor_status) %>%
                 dplyr::summarize(count = dplyr::n()) %>%
                 dplyr::mutate(per= round(prop.table(count) * 100,2))
         # print(head(plot_data))
          p1[[j]] <- ggplot() + geom_bar(aes(y = count, x =!! rlang::sym(groupby[j]), fill =scissor_status),
                                         data = plot_data,
                                         stat="identity")+
                scale_fill_manual(name = 'scissor_status',  values =stringr::str_split(opt$palette,",")%>%unlist%>%as.character   ) +
                scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
                coord_cartesian(clip = 'off') +
                theme_bw() +
                    labs(title=glue::glue("{groupby[j]}"))+
                theme(
                  legend.position = 'right',
                  plot.title = element_text(hjust = 0.5),
                  text = element_text(size = 15),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.title = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust=1),
                  plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
                )
          p2[[j]] <- ggplot() + geom_bar(aes(y = per, x = !! rlang::sym(groupby[j]), fill = scissor_status),
                                         data = plot_data,
                                         stat="identity")+
                scale_fill_manual(name = 'scissor_status',  values =stringr::str_split(opt$palette,",")%>%unlist%>%as.character  ) +
                scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
                coord_cartesian(clip = 'off') +
                theme_bw() +
                labs(title=glue::glue("{groupby[j]}"))+
                theme(
                  legend.position = 'right',
                  plot.title = element_text(hjust = 0.5),
                  text = element_text(size = 15),
                  axis.text.x = element_text(angle = 45, hjust=1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.title = element_blank(),
                  plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
                )

        OESingleCell::save_ggplots(  plot = patchwork::wrap_plots(p1[[j]] ,p2[[j]],nrow=1),
                                     glue::glue("{output_dir}/{alpha_param}/barplot_for_scissor_{opt$tag}_{groupby[j]}_{alpha_param}"),
                                     limitsize =F,
                                     width = 16,
                                     height =7 )
      }
    #},mc.cores =4)
    ## 保存meta数据至predicate_Scissor_metadata.tsv"
    sc_dataset@meta.data %>%
      tibble::rownames_to_column("barcode")%>%
      #dplyr::select("barcode",stringr::str_split(opt$groupby,","), dplyr::starts_with("scissor"))%>%
      readr::write_tsv(glue::glue("{output_dir}/{alpha_param}/predicate_Scissor_metadata.tsv"))
  ##====================================================================================================================
    futile.logger::flog.info("step5: Reliability significance test")#=============================================
    if(opt$test){
      result_reliability_test<-list()
      resutl_evaluate_summary<-list()
      load(glue::glue("{output_dir}/Scissor_{opt$tag}.RData"))
      print("1")
      ## 1.显著性检验（Reliability significance test）10x 交叉验证法
      numbers <- length(infos$Scissor_pos) + length(infos$Scissor_neg)
      result_reliability_test<-Scissor::reliability.test(X,
                                                         Y,
                                                         network,
                                                         alpha =alpha_param,
                                                         family = opt$family,
                                                         cell_num = numbers,
                                                         n = 10,
                                                         nfold = 10)

      ## 2.Cell level evaluations 注意load文件路径不能为绝对路径，需要为相对路径
      print("2")
      resutl_evaluate_summary<- Scissor::evaluate.cell(glue::glue("{opt$output}/Scissor_{opt$tag}.RData"),
                                                       infos,
                                                       FDR = 0.05,
                                                       bootstrap_n = 100)
      resutl_evaluate_summary %>%
        tibble::rownames_to_column("barcode")%>%
        readr::write_tsv(glue::glue("{output_dir}/{alpha_param}/{opt$tag}_cell_level_{alpha_param}.tsv"))

      save(result_reliability_test,
           resutl_evaluate_summary,
           file = glue::glue("{output_dir}/{alpha_param}/Scissor_{opt$tag}_test.RData"))
    }

    ##====================================================================================================================
    ## save session information
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit() }
}
