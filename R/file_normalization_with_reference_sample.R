#' Normalize using reference sample
#'
#' @description Perform normalization using reference files. Takes advantage of
#' CytoNorm package.
#'
#' @param df Data frame containing following columns:
#' file_paths (the full path to the files to be normalized),
#' batch_label (batch label for each file), ref_ids (logical defining TRUE values
#' for reference sample).
#' @param markers_to_normalize Character vector, marker names to be normalized,
#' can be full marker name e.g. "CD45$" (only CD45 marker will be picked) or
#' "CD" (all markers containig "CD" will be used).
#' If NULL (default) all non-mass markers will be normalized.
#' @param transformList Transformation list to pass to the flowCore
#' transform function. Defult is set to NULL. Either transformList or
#' arcsine_transform needs to be defined.
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5, default is set to TRUE. Either
#' transformList or arcsine_transform needs to be defined.
#' @param nQ Numeric, as in CytoNorm, number of quantiles to use.
#' Default = 101, which results in quantiles for every percent of the data.
#' @param limit Numeric, as in CytoNorm, these values are modeled to map onto themselves
#' by the spline.
#' @param quantileValues Numeric, as in CytoNorm, If specified,
#' it should be a vector of length nQ with values between 0 and 1, giving the
#' percentages at which the quantiles should be computed. If NULL (default),
#' the quantiles will be evenly distributed, including 0 and 1.
#' @param goal Goal distribution.
#' Default "mean", can also be nQ numeric values or one of the batch labels.
#' @param to_plot Logical, if TRUE, a plot is generated (using the layout function)
#' showing all quantiles. Ic norm_with_cluster = TRUE, FlowSOM clustering quality
#' plots will be also generated. Default = FALSE.
#' @param norm_with_clustering Logical, if data should be normalized using
#' clustering algorithm, FlowSOM.Default set to FALSE.
#' @param seed Numeric, set to obtain reproducible results,
#' when norm_with_clustering set to TRUE. Default NULL.
#' @param nCells Numeric, the number of cells, to use for FlowSOM clustering.
#' This number is determined by total number of fcs files, as by default 1000 cells
#' is used per file. Only if norm_with_clustering set to TRUE.
#' @param xdim Numeric, parameter to pass to FlowSOM, width of the SOM grid.
#' Only if norm_with_clustering set to TRUE.
#' @param ydim Numeric, parameter to pass to FlowSOM, geight of the SOM grid.
#' Only if norm_with_clustering set to TRUE.
#' @param nClus Numeric, exact number of clusters for metaclustering.
#' Only if norm_with_clustering set to TRUE.
#' @param clustering_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45$" (only CD45 marker will be picked) or
#' "CD" (all markers containig "CD" will be used). Default (NULL), all the
#' markers defined in markers_to_normalize will be used.
#' Only if norm_with_clustering set to TRUE.
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#' be saved, default is set to working directory. If NULL, files will be saved
#' in file.path(getwd(), CytoNorm).
#' @param save_model Logical, if the model should be saved, if TRUE it will be
#' saved to our_dir. Default set to FALSE.
#'
#' @return model, describing the normalization function
#'
#' @export
#'
#' @examples
#' # Set input directory
#' gate_dir <- file.path(dir, "Gated")
#'
#' # Define reference samples
#' files_ref <- list.files(gate_dir,
#'                         pattern = "*_gated.fcs$",
#'                         full.names = TRUE,
#'                         recursive = T)
#'
#' df <- data.frame("file_paths" = files_ref,
#'                 "batch_labels" = stringr::str_match(files_ref, "day[0-9]*")[,1],
#'                 "ref_ids" = grepl("REF", files_ref))
#'
#'
#' model <- train_REF_model(df = df,
#'                         markers_to_normalize = c("CD", "HLA", "IgD",
#'                                                  "IL", "TN", "MCP", "MIP",
#'                                                  "Gran", "IFNa"),
#'                         arcsine_transform = TRUE,
#'                         nQ = 2,
#'                         limit = c(0,8),
#'                         quantileValues = c(0.05, 0.95),
#'                         goal = "mean",
#'                         norm_with_clustering = FALSE,
#'                         save_model = TRUE,
#'                         clustering_markers = c("CD", "HLA", "IgD"))
#'
train_REF_model <- function(df,
                            markers_to_normalize = NULL,
                            transformList = NULL,
                            arcsine_transform = TRUE,
                            nQ = 101,
                            limit = NULL,
                            quantileValues = NULL,
                            goal = "mean",
                            to_plot = TRUE,
                            norm_with_clustering = FALSE,
                            seed = NULL,
                            nCells = 10000,
                            xdim = 10,
                            ydim = 10,
                            nClus = 10,
                            clustering_markers = NULL,
                            out_dir = NULL,
                            save_model = FALSE){

  if(!is(df, "data.frame")){
    stop("df is not a data frame")
  }

  if(!all(colnames(df) == c("file_paths", "batch_labels", "ref_ids"))){
    stop("colnames in df does not match: file_paths, batch_ids, ref_ids
        please correct colnames")
  }

  if(!all(file.exists(df$file_paths))){
    id <- !file.exists(df$file_paths)
    print(df$file_paths[id])
    stop("above files have incorrect path")
  }


  flow_frame <- flowCore::read.FCS(df$file_paths[1])
  if (!is.null(markers_to_normalize)){

    matches <- paste(markers_to_normalize, collapse="|")

    m_to_keep <- names(grep(matches, FlowSOM::GetMarkers(flow_frame, flowCore::colnames(flow_frame)),
                            ignore.case = TRUE, value = TRUE))
    m_to_keep <- grep(pattern = "File", x = m_to_keep, ignore.case = TRUE,
                      invert = TRUE, value = TRUE)

  } else {
    m_to_keep <- grep(pattern = "Time|length|Center|Offset|Width|Residual|File|File_scattered",
                      x = flowCore::colnames(flow_frame),
                      ignore.case = TRUE, value = TRUE, invert = TRUE)
  }

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "CytoNormed")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  files_ref <- df$file_paths[df$ref_ids]
  message(paste("The following reference files were found"))
  print(basename(files_ref))

  labels_ref <- df$batch_labels[df$ref_ids]

  if(arcsine_transform & !is.null(transformList)){
    stop("define either arcsine_transform or transformList parameter")
  } else if(arcsine_transform==FALSE & is.null(transformList)){
    stop("define either arcsine_transform or transformList parameter")
  }

  if(arcsine_transform){
    trans <- flowCore::transformList(m_to_keep,
                                     CytoNorm::cytofTransform)

  } else {
    trans <- transformList
  }

  if(to_plot & !norm_with_clustering){

    png(file.path(out_dir, "REF_normalization_taining_model.png"),
        width = length(m_to_keep) * 300,
        height = (length(files_ref) * 2 + 1) * 300)
  }


  if(norm_with_clustering){
    print("Clustering the data using FlowSOM")

    if(!is.null(clustering_markers)){
      markers_all <- FlowSOM::GetMarkers(flow_frame, flowCore::colnames(flow_frame))
      clustering_pattern <- paste(clustering_markers, collapse = "|")
      clustering_pattern <- names(grep(clustering_pattern, markers_all, value = TRUE))
    } else {
      clustering_pattern <- m_to_keep
    }


    nC <- length(files_ref)*nCells
    model <- CytoNorm::CytoNorm.train(files = files_ref,
                                      labels = labels_ref,
                                      channels = m_to_keep,
                                      transformList = trans,
                                      plot = to_plot,
                                      seed = seed,
                                      normParams = list(nQ = nQ,
                                                        limit = limit,
                                                        quantileValues = quantileValues,
                                                        goal = goal),
                                      FlowSOM.params = list(nCells = nC,
                                                            xdim = xdim,
                                                            ydim = ydim,
                                                            nClus = nClus,
                                                            scale = FALSE,
                                                            colsToUse = clustering_pattern),
                                      outputDir = out_dir)

  } else {
    model <- CytoNorm::QuantileNorm.train(files = files_ref,
                                          labels = labels_ref,
                                          channels = m_to_keep,
                                          transformList = trans,
                                          nQ = nQ,
                                          limit = limit,
                                          quantileValues = quantileValues,
                                          goal = goal,
                                          plot = to_plot)
  }



  if(to_plot & !norm_with_clustering){
    dev.off()
  }
  print(out_dir)
  if (save_model){
    saveRDS(model, file.path(out_dir, "REF_normalization_model.RSD"))
  }
  return(model)
}

#' Normalize the data
#'
#' @param model Model of the batch effercts, as computed by train_REF_model.
#' @param df Data frame containing following columns:
#' file_paths (the full path to the files to be normalized),
#' batch_label (batch label for each file), ref_ids (logical defining TRUE values
#' for reference sample).
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5, default is set to TRUE. Either
#' transformList or arcsine_transform needs to be defined.
#' @param transformList Transformation list to pass to the flowCore
#' transform function. Defult is set to NULL. Either transformList or
#' arcsine_transform needs to be defined.
#' @param transformList.reverse Transformation list with the reverse function,
#' so the normalized files can be saved in the untransformed space
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#' be saved, default is set to working directory. If NULL, files will be saved
#' in file.path(getwd(), CytoNorm).
#' @param norm_with_clustering Logical, if teh model was built using
#' clustering algorithm, FlowSOM. Default is FALSE.
#'
#' @return normalized fcs files
#'
#' @export
#'
#' @examples
#' # Set input directory
#' gate_dir <- file.path(dir, "Gated")
#'
#' # Define reference samples
#' files_ref <- list.files(gate_dir,
#'                         pattern = "*_gated.fcs$",
#'                         full.names = TRUE,
#'                         recursive = T)
#'
#' df <- data.frame("file_paths" = files_ref,
#'                 "batch_labels" = stringr::str_match(files_ref, "day[0-9]*")[,1],
#'                 "ref_ids" = grepl("REF", files_ref))
#'
#'
#' model <- train_REF_model(df = df,
#'                         markers_to_normalize = c("CD", "HLA", "IgD",
#'                                                  "IL", "TN", "MCP", "MIP",
#'                                                  "Gran", "IFNa"),
#'                         arcsine_transform = TRUE,
#'                         nQ = 2,
#'                         limit = c(0,8),
#'                         quantileValues = c(0.05, 0.95),
#'                         goal = "mean",
#'                         norm_with_clustering = FALSE,
#'                         save_model = TRUE,
#'                         clustering_markers = c("CD", "HLA", "IgD"))
#'
#' # Normalize files
#' normalize_REF(model = model, df = df, arcsine_transform = TRUE,
#'               norm_with_clustering = FALSE)
#'
normalize_REF <- function(model,
                          df,
                          arcsine_transform = TRUE,
                          transformList = NULL,
                          transformList.reverse = NULL,
                          out_dir = NULL,
                          norm_with_clustering = FALSE){

  files <- df$file_paths

  if(!all(file.exists(files))){
    stop("the files does not exist, please specify the corect pathway")
  }

  labels <- df$batch_labels

  if(arcsine_transform & !is.null(transformList)){
    stop("define either arcsine_transform or transformList parameter")
  } else if(arcsine_transform==FALSE & is.null(transformList)){
    stop("define either arcsine_transform or transformList parameter")
  }

  if(arcsine_transform){
    flow_frame <- flowCore::read.FCS(df$file_paths[1])


    trans <- flowCore::transformList(flowCore::colnames(flow_frame),
                                     CytoNorm::cytofTransform)
    trans_rev <- flowCore::transformList(flowCore::colnames(flow_frame),
                                         CytoNorm::cytofTransform.reverse)

  } else {
    if(is.null(transformList)){
      stop("define transformList")
    }

    trans <- transformList

    if(is.null(transformList.reverse)){
      stop("define transformList")
    }
    trans_rev <- transformList.reverse
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "CytoNormed")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if(norm_with_clustering){
    CytoNorm::CytoNorm.normalize(model = model,
                                 files = files,
                                 labels = labels,
                                 transformList = trans,
                                 transformList.reverse = trans_rev,
                                 outputDir = out_dir)

  } else {
    CytoNorm::QuantileNorm.normalize(model = model,
                                     files = files,
                                     labels = labels,
                                     transformList = trans,
                                     transformList.reverse = trans_rev,
                                     outputDir = out_dir)
  }
}
