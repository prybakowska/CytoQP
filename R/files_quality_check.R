#' Check the quality of acquired files
#'
#' @description Wrapper function to perform sample quality scoring.
#' First, it clusters the data per each batch (if batch argument is defined) and
#' calculates the AOF scores and Quality scores per batch using AOF algorithm.
#' Next, based on Quality scores, it detects outliers across all the files,
#' regarding the batch.
#'
#' @param fcs_files Character, full path to fcs files.
#' @param file_batch_id Character vector with batch label for each fcs_file,
#' the order and the length needs to be the same as in fcs_files. if only batch
#' is processed can be prepared as e.g. file_batch_id <- rep("batch", length(files))
#' @param out_dir Character, pathway to where the plots should be saved,
#' default is set to NULL, which means that the following path will be created
#' file.path(getwd(), "Quality_Control").
#' @param phenotyping_markers Character vector, marker names to be used for
#' flowsom clustering including DNA marker Iridium and viability staining
#' if available. Can be full marker name e.g. "CD45" or pattern "CD" if
#' all CD-markers needs to be plotted. Default is set to NULL, thus all the mass
#' channels will be used.
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5. Default is set to TRUE.
#' If FALSE, the transform list to pass to the flowCore transform function must
#' be defined and pass as an additional argument to fsom_aof function.
#' @param sd Numeric, number of standard deviation allowed for file outlier
#' detection, default = 3.
#' @param nClus Numeric, as in FlowSOM, number of metaclusters to be obtained
#' @param ... Arguments to be passed to fsom_aof function for FlowSOM parameter
#' adjustment and plotting: xdim, ydim, transform_list, my_colors, seed, to_plot.
#'
#' @return Quality file scores, plots Quality AOF scores for all files and save .RDS and .csv Quality
#' scores for further analysis, files are saved in out_dir.
#'
#' @examples
#' # Set input directory
#' clean_dir <- file.path(dir, "Cleaned")
#'
#' # Define files for visualization
#' files <- list.files(clean_dir,
#'                    pattern = "_cleaned.fcs$",
#'                    full.names = TRUE)
#'
#' # Define batch_id for each file
#' file_batch_id <- stringr::str_match(basename(files),
#'                                    "(day[0-9]*).*.fcs")[,2]
#'
#' file_quality_check(fcs_files = files,
#'                   file_batch_id = file_batch_id,
#'                   phenotyping_markers = c("Ir","CD", "HLA", "IgD", "Pt"),
#'                   arcsine_transform = TRUE,
#'                   nClus = 10,
#'                   sd = 3)
#'
#' @importFrom flowCore transformList
#' @import ggplot2
#'
#' @export
file_quality_check <- function(fcs_files,
                               file_batch_id = NULL,
                               out_dir = NULL,
                               phenotyping_markers = NULL,
                               arcsine_transform = TRUE,
                               sd = 3,
                               nClus = 10,
                               ...){

  # Check parameters
  if(!is(fcs_files, "character") & !is(fcs_files, "list")) {
    stop("files must be a character vector or a list")
  }

  if(is(fcs_files, "list")){
    fcs_files <- unlist(fcs_files)
  }

  if(length(file_batch_id) != length(fcs_files)){
    stop("the lenght of the file_batch_id is not equal to the lenght of fcs_files")
  }
  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if (!is.null(file_batch_id)) {
    scores <- lapply(unique(file_batch_id), function(batch) {
      print(batch)

      files <- fcs_files[file_batch_id == batch]
      fsom <- fsom_aof(fcs_files = files,
                       phenotyping_markers = phenotyping_markers,
                       out_dir = out_dir,
                       arcsine_transform = arcsine_transform,
                       nClus = nClus,
                       batch = batch, ...)
      return(aof_scoring(fcs_files = files,
                         phenotyping_markers = phenotyping_markers,
                         fsom = fsom, out_dir = out_dir, batch = batch))
    })
    names(scores) <- unique(file_batch_id)
  }
  else {
    files <- fcs_files
    fsom <- fsom_aof(fcs_files = files,
                     phenotyping_markers = phenotyping_markers,
                     out_dir = out_dir, arcsine_transform = arcsine_transform,
                     nClus = nClus,
                     batch = NULL)

    scores <- aof_scoring(fcs_files = files,
                          phenotyping_markers = phenotyping_markers,
                          fsom = fsom, out_dir = out_dir, batch = NULL)
  }

  final_score <- file_outlier_detecion(scores = scores, out_dir = out_dir,
                                       sd = sd)
  return(final_score)
}


#' Prepares FlowSOM
#'
#' @description Builds FlowSOM object
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param phenotyping_markers Character vector, marker names to be used for
#' clustering, can be full marker name e.g. "CD45" or "CD" if all CD-markers
#' needs to be plotted.
#' @param nCells Numeric, the total number of cells, to use for FlowSOM clustering.
#' This number is determined by total number of fcs files, as a default 10000 cells
#' is used per file
#' @param xdim Numeric, as in FlowSOM, parameter to pass to FlowSOM,
#' width of the SOM grid
#' @param ydim Numeric, as in FlowSOM, parameter to pass to FlowSOM,
#' height of the SOM grid
#' @param nClus Numeric, exact number of clusters for metaclustering
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#'  be saved, default is set to working directory.
#' @param batch Character, the name of the acquisition batch for each fcs file.
#' This name is passed to FlowSOM plot name, default is set to NULL
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5. Default set to TRUE.
#' @param seed numeric, set to obtain reproducible results, default 1.
#' @param transform_list Transformation list to pass to the flowCore
#' transform function see flowCore::transformList.
#' @param my_colors An array specifying colors to be used for the background
#' coloring of metaclusters in FlowSOM and t-SNE plot. Must have a length equal
#' to the nClus.If NULL (default) colors defined in FlowSOM package will be used.
#'
#' @param to_plot Logical, if FlowSOM tree and t-SNE map should be plotted,
#' default set to TRUE.
#'
#' @return fsom object
#'
#' @export
fsom_aof <- function(fcs_files,
                     phenotyping_markers,
                     nCells = length(fcs_files)*10000,
                     xdim = 10,
                     ydim = 10,
                     nClus = 10,
                     out_dir = NULL,
                     batch = NULL,
                     arcsine_transform = TRUE,
                     transform_list = NULL,
                     my_colors = NULL,
                     seed = 1,
                     to_plot = TRUE){


  if(!exists("phenotyping_channels")){
    o <- capture.output(ff_tmp <- flowCore::read.FCS(file.path(fcs_files[1])))
    markers <- FlowSOM::GetMarkers(ff_tmp, flowCore::colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers,
                                       collapse = ("|")), markers, value = TRUE)
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }

  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if(arcsine_transform){
    trans <- flowCore::transformList(names(phenotyping_channels),
                                     CytoNorm::cytofTransform)
  }
  else {
    if(is.null(transform_list)){
      stop("parameter transform_list must be defined")
    }
    trans <- transform_list
  }

  fsom <- CytoNorm::prepareFlowSOM(file = fcs_files,
                                   colsToUse = names(phenotyping_channels),
                                   seed = seed,
                                   nCells = nCells,
                                   transformList = trans,
                                   FlowSOM.params = list(xdim = xdim,
                                                         ydim = ydim,
                                                         nClus = nClus,
                                                         scale = FALSE))

  if(to_plot){
    if(is.null(my_colors)){
      backgroundColors <- NULL
    }
    else {

      if(max(as.numeric(fsom$metaclustering)) < length(my_colors)){
        warning("The number of colors is greater than the number of metaclusters only
             needed number of colors will be used")
        nmcl <- max(as.numeric(fsom$metaclustering))
        backgroundColors <- my_colors[1:nmcl]

      }

      if(max(as.numeric(fsom$metaclustering)) > length(my_colors)){
        warning("The number of colors is lower than the number of metaclusters
              default colors will be used")
        backgroundColors <- NULL

      }

      if(max(as.numeric(fsom$metaclustering)) == length(my_colors)){
        backgroundColors <- my_colors
      }
    }

    if(!is.null(batch)){
      filename <- paste0(batch, "_FlowSOM_clustering.pdf")
    }
    else {
      filename <- "FlowSOM_clustering.pdf"
    }

    fsomPlot <- FlowSOM::PlotStars(fsom = fsom,
                                   title = "FlowSOM clustering",
                                   backgroundValues = fsom$metaclustering,
                                   maxNodeSize = 3,
                                   backgroundColors = backgroundColors)
    fsomTsne <- FlowSOM::PlotDimRed(fsom = fsom, plotFile = NULL, seed = seed,
                                    cTotal = 20000,
                                    title = "tSNE visualization of FlowSOM metaclusters")

    figure <- ggpubr::ggarrange(fsomPlot, fsomTsne,
                                # labels = c("FlowSOM clustering", "tsne"),
                                ncol = 2, nrow = 1)

    ggplot2::ggsave(filename = filename, plot = figure, device = "pdf",
                    path = out_dir,
                    width =24, height = 10)
  }

  return(fsom)
}


#' Calculates AOF scores and scaled AOF scores
#'
#' @description  Calculates AOF (Average Overlap Frequency) scores using flowSOM
#' object and greedy algorithm from cytutils package.
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param phenotyping_markers Character vector, marker names to be used for
#' clustering, can be full marker name e.g. "CD45" or "CD" if all CD-markers
#' needs to be plotted.
#' @param fsom FlowSOM object as generated by fsom_aof
#' @param out_dir Character, pathway to where the AOF scores and plots should
#' be saved, default is set to file.path(getwd(), "Quality_Control")
#' @param batch Character, acquisition batch for each fcs file.
#' This argument is passed to AOF plot names, default is set to NULL.
#'
#' @return Returns data frame with the scaled AOF scores and heatmap plots
#' representing AOF scores and scaled AOF scores.
#'
#' @export
aof_scoring <- function(fcs_files,
                        phenotyping_markers,
                        fsom,
                        out_dir = NULL,
                        batch = NULL){

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  if(!exists("phenotyping_channels")){

    ff_tmp <- flowCore::read.FCS(file.path(fcs_files[1]))
    markers <- FlowSOM::GetMarkers(ff_tmp, flowCore::colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers,
                                       collapse = ("|")), markers, value = TRUE)

    if(length(grep("Ir", phenotyping_channels)) > 1){
      phenotyping_channels <- phenotyping_channels[-(grep("Ir",
                                                          phenotyping_channels)[2])]
    }
  }

  aof_scores <- lapply(fcs_files, function(file) {
    print(paste("calculating AOF", file))
    File_ID <- which(fcs_files == file)
    idx <- which(fsom$data[,"File"] == File_ID)
    fcs_data <- fsom$data[idx,]
    MC <- fsom$metaclustering[fsom$map$mapping[idx, 1]]

    aof_tmp <- cytutils::greedyCytometryAof(fcs_data = fcs_data,
                                            y = MC,
                                            channel_names = names(phenotyping_channels),
                                            width = 0.05,
                                            cofactor = 5,
                                            verbose = TRUE)
    return(aof_tmp$Aof)
  })

  aof_scores <- do.call("rbind", aof_scores)
  rownames(aof_scores) <- fcs_files
  colnames(aof_scores) <- names(phenotyping_channels)

  scaled_aof_score(aof_scores = aof_scores,
                   out_dir = out_dir,
                   aof_channels = phenotyping_channels,
                   batch = batch)
}


#' Calculates scaled aof scores and sample quality scores
#'
#' @description Calculates scaled AOF scores and sample quality score.
#' Additionally plots heatmaps for both raw  aof scores and scaled aof scores.
#' @param aof_scores Matrix, array, Aof scores obtained using function
#' cytutils::greedyCytometryAof.
#' @param out_dir Character, pathway to where the plots and scores should be
#' saved, if NULL file.path(getwd(), "Quality_Control) will be used.
#' @param aof_channels Character vector with the markers and their corresponding
#' channel names used for aof_scoring. Used only for plotting markers instead
#' of channels. If NULL (default) the colnames(aof_scores) will be used.
#' @param batch Character, acquisition batch for each fcs file.
#' Pass to FlowSOM plot name, default is set to NULL
#' be saved, default is set to working directory.
#'
#' @return Returns data frame with sample scores.Saved heatmaps and data frames
#' for AOF scores and AOF scaled scores.
#'
#' @export
scaled_aof_score <- function(aof_scores, out_dir = NULL, aof_channels = NULL,
                             batch = NULL){
  aof_scores_scaled <- scale(aof_scores)
  aof_scores_scaled <- pmax(aof_scores_scaled, 0)^2
  sample_scores <- apply(aof_scores_scaled, 1, sum, na.rm = TRUE)

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  df <- as.data.frame(sample_scores)

  list_scores <- list("aof_scores_per_marker" = aof_scores,
                      "scaled_AOF" = aof_scores_scaled)

  for (name in names(list_scores)) {

    if(!is.null(batch)){
      filename <- file.path(out_dir, paste0(batch, "_", name, ".pdf"))
      main <- paste0(batch, "_", name)
    } else {
      filename <- file.path(out_dir, paste0(name, ".pdf"))
      main <- name
    }

    if(is.null(aof_channels)){
      phenotyping_channels <- colnames(list_scores[[name]])
    }
    else if (all(colnames(aof_scores) %in% names(aof_channels))){

      phenotyping_channels <- aof_channels[colnames(aof_scores)]
    } else {
      phenotyping_channels <- colnames(list_scores[[name]])
      warning("aof_channels do not correspond to the channels selected for
              AOF scoring, colnames from aof_scores will be used for plotting")
    }

    pheatmap::pheatmap(list_scores[[name]],
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       color = colorRampPalette(
                         RColorBrewer::brewer.pal(n = 9,
                                                  name = "YlGnBu"))(100),
                       display_numbers = TRUE,
                       labels_col = phenotyping_channels,
                       labels_row = basename(rownames(list_scores[[name]])),
                       filename = filename,
                       main = main,
                       number_format = "%.1f",
                       fontsize_number = 8,
                       number_color = "black",
                       width = 10)
  }

  if(is.null(batch)){
    saveRDS(list_scores, file.path(out_dir, "AOF_scores_and_Scaled_AOF_scores.RDS"))
  }
  else {
    saveRDS(list_scores, file.path(out_dir,
                                   paste0(batch, "_AOF_scores_and_Scaled_AOF_scores.RDS")))
  }
  return(df)
}

#' Detects outliers based on sample quality scores
#'
#' @description Detects outlier files based on sample quality scores,
#' generates plot for outlier and .csv file which indicates which fcs files
#' could be discarded from further analysis.
#'
#' @param scores List of scaled scores per acquisition batch or data frame
#' of scaled scores, both generated by scaled_aof_score or aof_scoring function
#' @param out_dir Character, pathway to where the plot and .csv files with
#' quality scores should be saved, default is set to NULL, thus
#' file.path(getwd(), "Quality_Control" will be generated.
#' @param sd How many standard deviation should be use to detect outliers
#' default is set to 3.
#'
#' @return Save .RDS and .csv Quality scores for further analysis, and
#' plots and save .png for Quality AOF scores for all files.
#'
#' @export
file_outlier_detecion <- function(scores, out_dir = NULL, sd) {

  if(!inherits(scores, "data.frame") & !inherits(scores, "list")){
    stop("df scores are neither data frame nor list of the data frames")
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if(inherits(scores, "list")){
    df_scores <- do.call(rbind, scores)
  }
  else {
    df_scores <- scores
  }

  df_scores$file_names <- basename(rownames(df_scores))

  scores_median <- stats::median(df_scores$sample_scores)
  scores_MAD <- stats::mad(df_scores$sample_scores)

  df_scores$quality <- ifelse(df_scores$sample_scores >
                                (scores_median + sd * scores_MAD),"bad","good")

  bad_scores <- sum(df_scores$quality == "bad")

  colors <- c("bad" = "red", "good" = "darkblue", "threshold= " = "orange")

  max_score <- max(df_scores$sample_scores)
  max_pctgs <- max_score + (max_score * 0.1)

  p <- ggplot2::ggplot(df_scores, aes(x = file_names, y = sample_scores,
                                      color = quality)) +
    geom_point(size = 4) +
    scale_colour_manual(values = colors) +
    ylim(-0.5, max_pctgs) +
    annotate(geom="text", x = mean(as.numeric(as.factor(df_scores$file_names))),
             y= max_score - 0.05*max_score, label=paste("N bad = ", bad_scores),
             color="red", size = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text = element_text(size = 11),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          # panel.border = fill = "black",
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
    scale_x_discrete(breaks = df_scores$file_names[df_scores$quality == "bad"])

  if(scores_median + sd * scores_MAD <= max_pctgs){
    p + geom_hline(yintercept = scores_median + sd * scores_MAD,
                   linetype = "dashed", color = "darkgreen", size = 1)
  }

  ggplot2::ggsave(filename = "Quality_AOF_score.png", plot = p,
                  path = file.path(out_dir))

  saveRDS(df_scores, file.path(out_dir, "Quality_AOF_score.RDS"))
  write.csv(df_scores, file = file.path(out_dir, "Quality_AOF_score.csv"))
  return(df_scores)
}
