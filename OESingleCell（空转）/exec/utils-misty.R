
#' Extracts data from an specific assay from a Seurat object and aligns the IDs to the geometry
#' and aligns the IDs to the geometry
#' @param visium.data Seurat object data of each sample.
#' @param assay matrix of counts.
#' @param geometry matrix of spatial locations for each spot.
#' @export
extract_seurat_data <- function(visium.data,
                                assay,
                                geometry) {
  print(assay)
  data <- SeuratObject::GetAssayData(visium.data, assay = assay) %>%
    as.matrix() %>%
    t() %>%
    tibble::as_tibble(rownames = NA)

  return(data %>% dplyr::slice(match(rownames(.), rownames(geometry))))
}

# Filters data to contain only features of interest
filter_data_features <- function(data,
                                 features) {
  if (is.null(features)) features <- colnames(data)

  return(data %>%
           tibble::rownames_to_column() %>%
           dplyr::select(rowname, tidyselect::all_of(features)) %>%
           dplyr::rename_with(make.names) %>%
           tibble::column_to_rownames())
}

#' Builds views depending on the paramaters defined
#' @param data Seurat object of each sample.
#' @param view.type views' interaction type (intra, juxta or para).
#' @param view.param views' distances threshold to intraview spots.
#' @param view.name views' names.
#' @param spot.ids all spots' id of intraview.
#' @param geometry matrix of spatial locations for each spot.
#' @export
create_default_views <- function(data,
                                 view.type,
                                 view.param,
                                 view.name,
                                 spot.ids,
                                 geometry) {

  mistyR::clear_cache()
  view.data.init <- mistyR::create_initial_view(data)
  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }
  if (view.type == "intra") {
    data.red <- view.data.init[["intraview"]]$data %>%
      tibble::rownames_to_column() %>%
      dplyr::filter(rowname %in% spot.ids) %>%
      dplyr::select(-rowname)
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>%
      mistyR::add_paraview(geometry, l = view.param)
    data.ix <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      dplyr::mutate(rowname = rownames(data)) %>%
      dplyr::filter(rowname %in% spot.ids) %>%
      dplyr::select(-rowname)
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>%
      mistyR::add_juxtaview(
        positions = geometry,
        neighbor.thr = view.param)
    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      dplyr::mutate(rowname = rownames(data)) %>%
      dplyr::filter(rowname %in% spot.ids) %>%
      dplyr::select(-rowname)
  }

  misty.view <- mistyR::create_view(
    paste0(view.name),
    data.red)

  return(misty.view)
}

#=======================================================================================================================
docstring <- " example1:\\n\\n\\
sctool  -i st.rds -f rds  -o results/misty misty   --misclist spotlight_results  --target celltype  --predictor celltype  --para 10 \\n\\n\\
example2: \\n\\n\\
sctool  -i st.rds -f rds --assay SCT  --dataslot data -o results/misty misty  --organism Human --misclist spotlight_results --target progeny --predictor celltype  --para 10 "
sub_misty <- subparsers$add_parser(
  "misty",
  description = docstring,
  formatter_class = "argparse.RawTextHelpFormatter",
  # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
  argument_default = "True",
  help = "Using misty to refer Spaital interactions between markers"
)
sub_misty$add_argument(
  "--organism",
  type = "character",
  default = "Human",
  help = "choose organism to perform progeny(Human or Mouse).[default: %(default)s]."
)
sub_misty$add_argument(
  "--misclist",
  type = "character",
  default = NULL,
  help = "the matrix name of the celltype annotation or other functional matrix in the misclist.[default: %(default)s]."
)
sub_misty$add_argument(
  "--target",
  type = "character",
  help = "choose which matrix as the target feature for the intraview (progeny or celltype)."
)
sub_misty$add_argument(
  "--predictor",
  type = "character",
  help = "choose which matrix as the predictor feature for the intraview (progeny or celltype)."
)
sub_misty$add_argument(
  "--remove_ct",
  type = "character",
  default = NULL,
  help = "filter out features which had no variance (noninformative).[default: %(default)s]."
)
sub_misty$add_argument(
  "--para",
  type = "integer",
  default = 15,
  help = "set the effective radius of spots in paraview.[default: %(default)s]."
)
#=======================================================================================================================
args <- commandArgs(TRUE)
if ("misty" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "misty") {
    # load data=========================================================================================================
    futile.logger::flog.info("step1: load query object data")
    suppressMessages(data_ob <- OESingleCell::ReadX(
            input = opt$input,
            informat = opt$informat,
            assays = assays[1],
            data.use = dataslots,
            verbose = F
    ))
    # creat annotation assay and get useful features====================================================================
    assay <- list()
    for (i in c(opt$predictor, opt$target)) {
      if (i == "celltype") {
        if (length(opt$misclist) > 0) {
          futile.logger::flog.info("step2: creat a misclist assay")
          anno <- data_ob@misc[[opt$misclist]]
          rownames(anno) <- anno[, 1]
          anno <- anno[, -1] %>% as.matrix() %>% t()
          data_ob[[opt$misclist]] <- SeuratObject::CreateAssayObject(data = anno)
          useful_features <- rownames(SeuratObject::GetAssayData(data_ob, assay = opt$misclist))
          if (!is.null(opt$remove_ct)) {
            useful_features <- useful_features[!useful_features %in% stringr::str_split(opt$remove_ct, ",")[[1]]]
          }
          tag <- "cell_types"
        }
        Seurat::DefaultAssay(data_ob) <- opt$misclist
      }
      if (i == "progeny") {
        library(progeny)
        data_ob <- progeny(data_ob,
                           assay_name = opt$assay,
                           scale = FALSE,
                           organism = glue::glue("{opt$organism}"),
                           top = 500,
                           perm = 1,
                           return_assay = TRUE)
        Seurat::DefaultAssay(data_ob) <- "progeny"
        useful_features <- rownames(data_ob)
        if (!is.null(opt$remove_ct)) {
          useful_features <- useful_features[!useful_features %in% stringr::str_split(opt$remove_ct, ",")[[1]]]
        }
        tag <- "pathway_activities"
      }

      if (opt$target == i) {
        assay[["target"]] <- Seurat::DefaultAssay(data_ob)
        annotation_target <- tag
        useful_features <- useful_features
      }else {
        assay[["predictor"]] <- Seurat::DefaultAssay(data_ob)
        annotation_predictor <- tag
        useful_features_ct <- useful_features
      }
    }
    #define parameters==================================================================================================
    if (identical(opt$predictor, opt$target)) {
      view_assays <- list("main" = assay[["target"]],
                          "juxta" = assay[["target"]],
                          "para" = assay[["target"]])
      view_features <- list("main" = useful_features,
                            "juxta" = useful_features,
                            "para" = useful_features)
      view_types <- list("main" = "intra",
                         "juxta" = "juxta",
                         "para" = "para")
      view_params <- list("main" = NULL,
                          "juxta" = 5,
                          "para" = opt$para)
      views <- c("intra", "juxta", "para")
    }else {
      view_assays <- list("main" = assay[["target"]],
                          "juxta" = assay[["target"]],
                          "para" = assay[["target"]],
                          "predictor_intra" = assay[["predictor"]],
                          "predictor_juxta" = assay[["predictor"]],
                          "predictor_para" = assay[["predictor"]])
      view_features <- list("main" = useful_features,
                            "juxta" = useful_features,
                            "para" = useful_features,
                            "predictor_intra" = useful_features_ct,
                            "predictor_juxta" = useful_features_ct,
                            "predictor_para" = useful_features_ct)
      view_types <- list("main" = "intra",
                         "juxta" = "juxta",
                         "para" = "para",
                         "predictor_intra" = "intra",
                         "predictor_juxta" = "juxta",
                         "predictor_para" = "para")
      view_params <- list("main" = NULL,
                          "juxta" = 5,
                          "para" = opt$para,
                          "predictor_intra" = NULL,
                          "predictor_juxta" = 5,
                          "predictor_para" = opt$para)
      rename_view <- c("main",
                       "juxta",
                       "para",
                       glue::glue("intra_{opt$predictor}"),
                       glue::glue("juxta_{opt$predictor}"),
                       glue::glue("para_{opt$predictor}"))
      names(view_assays) <- rename_view
      names(view_features) <- rename_view
      names(view_types) <- rename_view
      names(view_params) <- rename_view
      views <- c("intra", "juxta", "para",
                 glue::glue("intra_{opt$predictor}"),
                 glue::glue("juxta_{opt$predictor}"),
                 glue::glue("para_{opt$predictor}"))
    }

    #get subset of cells if necessary===================================================================================
    if (!is.null(opt$predicate)) {
       futile.logger::flog.info(glue::glue("get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
       #factors <- stringr::str_split(opt$predicate, ",")[[1]]
       df <- slot(data_ob, "meta.data")
       desired_cells <- subset(df, eval(parse(text = opt$predicate)))
       #desired_cells <- subset(df, df[[factors[1]]] %in% factors[2:length(factors)])
       data_ob <- subset(data_ob, cells = rownames(desired_cells))
    }
    if (!is.null(Seurat::Images(data_ob))) {
       unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique())]
       if (length(unuse_images) > 0) {
           data_ob@images[unuse_images] <- NULL
       }
    }

    # extract geometry and viewdata=================================================================================
    futile.logger::flog.info("step3: extract geometry and viewdata")
    sampleid <- unique(data_ob@meta.data$sampleid)
    #merge <- list()
    misty_result_out <- list()
    misty_result_outdir <- list()
    for (i in sampleid) {
      sub_data_ob <- subset(data_ob, subset = sampleid == i)
      geometry <- SeuratObject::GetTissueCoordinates(sub_data_ob,
                                                     cols = c("row", "col"),
                                                     image = i,
                                                     scale = NULL)
      view.data <- purrr::map(view_assays,
                              extract_seurat_data,
                              geometry = geometry,
                              visium.data = sub_data_ob)
      spot.ids <- rownames(view.data[[1]])
      view.data.filt <- purrr::map2(view.data,
                                    view_features,
                                    filter_data_features)
      view.names <- names(view.data.filt)
      # creat views and run misty===================================================================================
      futile.logger::flog.info("step4: creat initial view")
      views.main <- mistyR::create_initial_view(view.data.filt[[1]] %>%
                                                  tibble::rownames_to_column() %>%
                                                  dplyr::filter(rowname %in% spot.ids) %>%
                                                  dplyr::select(-rowname))

      futile.logger::flog.info("step5: creat default views")
      all.views <- purrr::pmap(list(view.data.filt[-1], view_types[-1], view_params[-1], view.names[-1]),
                               create_default_views,
                               spot.ids = spot.ids,
                               geometry = geometry)
      pline.views <- mistyR::add_views(views.main, unlist(all.views, recursive = FALSE))

      futile.logger::flog.info("step6: run misty")
      # misty_out <- glue::glue("{output_dir}/{i}_results")
      mistyR::run_misty(pline.views, glue::glue("{output_dir}/{i}_results"), cached = FALSE)

      misty_result_outdir[[i]] <- glue::glue("{output_dir}/{i}_results")

      #transform file format=============================================================================================
      fileData <- list()
      fileData[["coefficients.txt"]] <- read.table(glue::glue("{output_dir}/{i}_results/coefficients.txt"),header = TRUE)
      fileData[["performance.txt"]] <- read.table(glue::glue("{output_dir}/{i}_results/performance.txt"),header = TRUE)
      for (view in views) {
        for(feature in colnames(view.data.filt[[1]])) {
          fileData[[glue::glue("importances_{feature}_{view}.txt")]] <- read.table(
            glue::glue("{output_dir}/{i}_results/importances_{feature}_{view}.txt"),header = TRUE, sep = ",")
        }
      }
      for (j in names(fileData)){
                  temp <- fileData[[j]]
                  write.table(temp, glue::glue("{output_dir}/{i}_results/{j}.xls"),
                              sep = "\t",
                              quote = FALSE,
                              row.names = FALSE)
      }
      # collect results=================================================================================================
      futile.logger::flog.info("step7: collecting results")
      misty_result_out[[i]] <- mistyR::collect_results(glue::glue("{output_dir}/{i}_results"))
    }
    if (length(sampleid) > 1) {
      misty_result_out[["merge"]] <- mistyR::collect_results(misty_result_outdir)
      results_importance <- misty_result_out[["merge"]]$importances.aggregated$Importance
    }else {
      results_importance <- misty_result_out$importances.aggregated$Importance
    }
    max_importance <- max(results_importance[!is.na(results_importance)])
    min_importance <- min(results_importance[!is.na(results_importance)])
    print(glue::glue("max_importance:{max_importance}"))
    print(glue::glue("min_importance:{min_importance}"))

    for (i in names(misty_result_out)) {
      futile.logger::flog.info("step8: plotting")
      ### coefficients
      plot.coef1 <- mistyR::plot_improvement_stats(misty_result_out[[i]])
      plot.coef2 <- mistyR::plot_view_contributions(misty_result_out[[i]])
      OESingleCell::save_ggplots(
        glue::glue("{output_dir}/{i}_results/coefficients/improvement_stats_plot"),
        plot = plot.coef1,
        limitsize = FALSE,
        height = 6,
        dpi = 600)
      OESingleCell::save_ggplots(
        glue::glue("{output_dir}/{i}_results/coefficients/view_contributions_plot"),
        plot = plot.coef2,
        limitsize = FALSE,
        height = 6,
        dpi = 600)

      for (view  in  unique(misty_result_out[[i]]$importances.aggregated$view)) {

        plot.heatmap <- mistyR::plot_interaction_heatmap(misty_result_out[[i]],
                                                         view,
                                                         title = view,
                                                         max_imp = max_importance,
                                                         min_imp = min_importance,
                                                         cutoff = 0)
        OESingleCell::save_ggplots(
          glue::glue("{output_dir}/{i}_results/{view}/{view}_interaction_heatmap_plot"),
          plot = plot.heatmap,
          limitsize = FALSE,
          height = 6,
          dpi = 600
        )
        if (view %in% c("intra","juxta","para")) {
          pdf(file = glue::glue("{output_dir}/{i}_results//{view}/{view}_interaction_communities_plot.pdf"))
          mistyR::plot_interaction_communities(misty_result_out[[i]],
                                               view,
                                               title = view,
                                               cutoff = 0.5)
          dev.off()
          png(file = glue::glue("{output_dir}/{i}_results/{view}/{view}_interaction_communities_plot.png"))
          mistyR::plot_interaction_communities(misty_result_out[[i]],
                                               view,
                                               title = view,
                                               cutoff = 0.5)
          dev.off()
        }
        if(i!="merge"){system(glue::glue("mv {output_dir}/{i}_results/importances_*_{view}.txt.xls  {output_dir}/{i}_results/{view}/"))}
      }
      write.table(misty_result_out[[i]]$importances.aggregated,
                  glue::glue("{output_dir}/{i}_results/total_importances_{i}.xls"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
      if(i!="merge"){
        system(glue::glue("mv  {output_dir}/{i}_results/{{coefficients.txt.xls,performance.txt.xls}}  {output_dir}/{i}_results/coefficients"))
        system(glue::glue("rm {output_dir}/{i}_results/*.txt"))
      }
    }
    if (!file.exists(file.path(output_dir, "空间互作mistyR分析说明文档.docx"))) {
        file.copy("/public/dev_scRNA/oesinglecell3_test/document/空间互作mistyR分析说明文档.docx",
                  file.path(output_dir, "空间互作mistyR分析说明文档.docx")) }
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}
