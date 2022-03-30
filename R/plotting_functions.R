#' gate_out_beads
#'
#' @description removes beads from the files that contain them e.g files before
#' bead normalization
#'
#' @param bead_channel character, the mass for bead channel that is exclusively used for
#' beads identification, no marker is present at this channel, default 140.
#' @param flow_frame Flow frame, if unstransformed arcsine_transform should be
#' kept as default, TRUE.
#'
#' @return flow frame with beads removed
.gate_out_beads <- function(bead_channel,
                            flow_frame){
  ch <- grep(pattern = bead_channel, x = flowCore::colnames(flow_frame), value = TRUE)
  ids <- flow_frame[,ch] > 0
  # calculate threshold
  th <- flowDensity::deGate(obj = flow_frame[ids,], channel = ch)
  #remove beads
  cells_to_remove <- flow_frame@exprs[, ch] < th
  flow_frame <- flow_frame[cells_to_remove,]
  return(flow_frame)
}

#' Plots quantiles for the markers
#'
#' @description Calculates quantiles (0.01, 0.25, 0.5, 0.75, 0.99) for
#' selected markers and plots them as diagnostic plots.
#'
#' @param files_before_norm Character vector or list with the paths of the raw files.
#' @param files_after_norm Character vector or list with the paths of the normalized files.
#' @param batch_pattern Character, batch pattern to be match in the fcs file name
#' @param uncommon_prefix Character vector or string, uncommon prefix in
#' the basename of the fcs files. The file names need to match, so uncommon prefix
#' needs to be removed. If NULL (default) prefix like
#' "Norm|_CC_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS|.fcs" will be removed.
#' Default is set to NULL.
#' @param bead_channel character, the mass for bead channel that is exclusively
#' used for beads identification (no marker is assign to this channel),
#' Default 140.
#' @param remove_beads Logical, if beads needs to be removed. This needs to be
#' set to TRUE if files contain beads e.g before beads normalization,
#' default is set to TRUE. For the visualization purpose the beads will be
#' removed using channel set in bead_channel.
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5.
#' @param markers_to_plot character vector, marker names to be plotted, can be
#' full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param manual_colors character, vector of the colors to be used,
#' the number of colors needs to be equal to the length of batch_pattern
#' @param out_dir Character, pathway to where the plots should be saved,
#' default is set to working directory.
#' @param transform_list Transformation list to pass to the flowCore
#' transform function, see flowCore::transformList, if different transformation
#' than arcsine is needed. Only if arcsine_transform is FALSE. If NULL and
#' arcsine_transform = FALSE no transformation will be applied.
#'
#' @return Save the pdf with plots to out_dir.
#'
#' @examples
#' # Define files for visualization
#' # Before normalization
#' raw_data_dir <- file.path(dir, "RawFiles")
#' files_b <- list.files(raw_data_dir,
#'                       pattern = ".FCS$",
#'                       ignore.case = T,
#'                       full.names = TRUE)
#'
#' # After normalization
#' bead_norm_dir <- file.path(dir, "BeadNorm")
#' files_a <- list.files(bead_norm_dir,
#'                       pattern = "_beadNorm.fcs$",
#'                       ignore.case = T,
#'                       full.names = TRUE)
#'
#' # Define batch id and sample id for each file
#' batch_pattern <- stringr::str_match(basename(files_b), "(?i).*(day[0-9]*).*.FCS")[,2]
#'
#' plot_marker_quantiles(files_after_norm = files_a,
#'                      files_before_norm = files_b,
#'                      batch_pattern = batch_pattern,
#'                      arcsine_transform = TRUE,
#'                      remove_beads = TRUE,
#'                      bead_channel = "140",
#'                      uncommon_prefix = "_beadNorm.fcs|.FCS",
#'                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
#'                                          "TGF", "GR", "IFNa"),
#'                      manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
#'                      out_dir = bead_norm_dir)
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#'
#' @export
plot_marker_quantiles <- function(files_before_norm,
                                  files_after_norm,
                                  batch_pattern = NULL,
                                  batch_labels = NULL,
                                  remove_beads = FALSE,
                                  bead_channel = "140",
                                  uncommon_prefix = NULL,
                                  arcsine_transform = TRUE,
                                  markers_to_plot = NULL,
                                  plot_name = "Marker_distribution_across_aliquots_and_batches",
                                  manual_colors = NULL,
                                  out_dir = NULL,
                                  transform_list = NULL){


  if(is(files_after_norm, "list")){
    files_after_norm <- unlist(files_after_norm)
  }

  if(is(files_before_norm, "list")){
    files_before_norm <- unlist(files_before_norm)
  }

  if(!(length(files_after_norm) == length(files_before_norm))){
    stop("files_before and after does not have the same length")
  }

  fcs_files <- c(files_after_norm, files_before_norm)
  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if(!is.null(uncommon_prefix)){
    test_match_order(x = basename(gsub(uncommon_prefix,".fcs",files_after_norm)),
                     y = basename(gsub(".FCS",".fcs", files_before_norm)))
  } else {
    test_match_order(x = basename(gsub("Norm_|_beadNorm","",files_after_norm)),
                     y = basename(gsub(".FCS",".fcs",files_before_norm)))

  }

  if(is.null(batch_labels) & is.null(batch_pattern)){
    stop("define batch_labels or batch_pattern")
  } else if (!(is.null(batch_labels)) & !(is.null(batch_pattern))){
    stop("both batch_labels and batch_pattern are defined, desellect one option by  setting to NULL")
  }


  tmp <- c(paste0(files_after_norm, "_YES"), paste0(files_before_norm, "_NO"))

  if(!is.null(batch_labels)){
    if(length(tmp) != length(rep(batch_labels, 2))){
      stop("The lenght of batch labels is not equal to the lenght of files")
    }
  }

  ff_tmp <- flowCore::read.FCS(file.path(fcs_files[1]))

  if (!is.null(markers_to_plot)){

    if(!is.character(markers_to_plot)){
      stop ("markers are not a character vector")
    }

    matches <- paste(markers_to_plot, collapse="|")

    norm_markers <- grep(matches,
                         FlowSOM::GetMarkers(ff_tmp, find_mass_ch(ff_tmp,
                                                                  value = TRUE)),
                         value = TRUE, ignore.case = FALSE)
  } else {
    norm_markers <- find_mass_ch(ff_tmp, value = TRUE)
    norm_markers <- FlowSOM::GetMarkers(ff_tmp, norm_markers)
  }


  quantile_values <-  c(0.01, 0.25, 0.5, 0.75, 0.99)
  quantiles <- expand.grid(File = tmp,
                           Marker = norm_markers,
                           Quantile = quantile_values,
                           Value = NA)

  if(!is.null(batch_pattern)){
    quantiles <- cbind(quantiles, "Batch" = stringr::str_match(
      basename(as.character(quantiles$File)), batch_pattern)[,1])
  }

  if(!is.null(batch_labels)){
    quantiles <- cbind(quantiles, "Batch" = rep(batch_labels, 2))
  }

  quantiles$Normalization <- gsub(".*.fcs_|.*.FCS_", "", quantiles$File)
  quantiles$File <- gsub("_YES|_NO", "", quantiles$File)

  for (file in fcs_files) {
    print(file)

    ff <- flowCore::read.FCS(file, transformation = FALSE)

    if(arcsine_transform){
      ff <- flowCore::transform(ff,
                                flowCore::transformList(grep("Di", flowCore::colnames(ff),
                                                             value = TRUE),
                                                        CytoNorm::cytofTransform))
    } else if (!is.null(transform_list)){
      ff <- flowCore::transform(ff, transform_list)
    } else {
      ff <- ff
    }

    norm <- quantiles$Normalization[(which(quantiles$File == file)[1])]

    if(norm == "NO" & remove_beads){
      ff <- .gate_out_beads(bead_channel = bead_channel, flow_frame = ff)
    }

    for (marker in names(norm_markers)) {
      quantiles_res <-stats::quantile(flowCore::exprs(ff)[, marker],
                                      quantile_values)
      for (i in seq_along(quantiles_res)) {
        quantiles <- quantiles %>%
          dplyr::mutate(Value = replace(Value,
                                        File == file &
                                          Marker == norm_markers[marker] &
                                          Quantile == quantile_values[i],
                                        quantiles_res[i]))

      }
    }
  }

  if(is.null(uncommon_prefix)){
    quantiles$Sample <- gsub(pattern = "Norm_", replacement = "",
                             ignore.case = TRUE,
                             x = gsub(
                               pattern = "_CC_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS|.fcs",
                               replacement = "", ignore.case = TRUE,
                               x = basename(as.character(quantiles$File))))
  }
  else {
    uncommon_prefix <- paste(uncommon_prefix, collapse = ("|"))
    quantiles$Sample <- gsub(pattern = "Norm_", replacement = "",
                             ignore.case = TRUE,
                             x = gsub(pattern = uncommon_prefix,
                                      replacement = "",
                                      ignore.case = TRUE,
                                      x =  basename(as.character(
                                        quantiles$File))))
  }

  ncols <- length(unique(quantiles$Batch))
  p <- quantiles %>% dplyr::filter(Normalization == "YES") %>%
    ggplot2::ggplot(aes(x = Sample,
                        y = Value,
                        color = Batch)) +
    geom_point(data = quantiles %>% dplyr::filter(Normalization == "NO"),
               aes(alpha = ifelse(Quantile == "0.5", 2, 0)), color = "grey31") +
    geom_line(data = quantiles %>% dplyr::filter(Normalization == "NO"),
              aes(alpha = ifelse(Quantile != "0.5", 2, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 1,
                                ifelse(Quantile == "0.5", 0, 2))),
              alpha = 0.5, color = "grey31") +
    ylab(label = levels(quantiles$Marker))+
    geom_point(aes(alpha = ifelse(Quantile == "0.5", 2, 0))) +
    geom_line(aes(alpha = ifelse(Quantile != "0.5", 2, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 1,
                                ifelse(Quantile == "0.5", 0, 2))),
              alpha = 0.5) +
    scale_size_identity() +
    scale_alpha_identity() +
    facet_wrap(~ Marker + Batch, ncol = ncols, scales = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")

  if (!is.null(manual_colors)){
    p <- p + ggplot2::scale_colour_manual(values = c(manual_colors))
  }

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- getwd()
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}


  ggplot2::ggsave(filename = paste0(plot_name, ".pdf"),
                  plot = p,
                  path = out_dir,
                  width = length(fcs_files)*0.25,
                  height = length(norm_markers)*4, limitsize = FALSE)
}


#' Plot 2D scatter plots
#'
#' @description Plots biaxial plots
#'
#' @param fcs_files Character, pathway to fcs files.
#' @param markers_to_plot Character, pattern of the markers to be plotted e.g.
#' "CD" (all CD markers will be plotted), "CD41$" (only CD41 will be plotted).
#' @param y_marker Character, The marker to be plotted on the y-axis.
#' @param out_dir Character, path where fill are saved, if NULL (default)
#' files are saved in getwd()
#' @param out_put_name The name of the file that is saved in out_dir, default is
#' "marker_plots.png".
#'
#' @return Save plots in out_dir.
#' @export
#'
#' @examples
#' plot_2D_scatter_plots(fcs_files = fcs_files,
#'                       markers_to_plot = c("CD", "HLA"),
#'                       y_marker = "CD45$",
#'                       out_dir = getwd(),
#'                       out_put_name = "marker_plots.png")
#'
#'
plot_2D_scatter_plots <- function(fcs_files,
                                  markers_to_plot,
                                  y_marker,
                                  out_dir = NULL,
                                  out_put_name = "marker_plots.png"){


  if(!all(file.exists(fcs_files))){
    stop("fcs files do not exist")
  }

  ff <- flowCore::read.FCS(fcs_files[1], transformation = FALSE)

  if (!is.null(markers_to_plot)){

    if(!is.character(markers_to_plot)){
      stop ("markers are not a character vector")
    }

    matches <- paste(markers_to_plot, collapse="|")

    norm_markers <- grep(matches,
                         FlowSOM::GetMarkers(ff, find_mass_ch(ff,
                                                              value = TRUE)),
                         value = TRUE, ignore.case = FALSE)

  } else {
    norm_markers <- find_mass_ch(ff, value = TRUE)
    norm_markers <- FlowSOM::GetMarkers(ff, norm_markers)
  }

  matches_y <- paste(y_marker, collapse="|")
  ym <- grep(matches_y,
             FlowSOM::GetMarkers(ff, find_mass_ch(ff,
                                                  value = TRUE)),
             value = TRUE, ignore.case = FALSE)

  n_plots <- length(norm_markers) - 1
  grDevices::png(file.path(out_dir, paste0(out_put_name)),
      width = n_plots * 300, height = length(fcs_files) * 300)
  graphics::layout(matrix(1:(length(fcs_files) * n_plots), ncol = n_plots, byrow = TRUE))

  for(file in fcs_files){
    ff <- flowCore::read.FCS(file, transformation = FALSE)

    cols_to_trans <- grep(pattern = "Di", x = flowCore::colnames(ff), value = TRUE)

    ff_t <- flowCore::transform(ff, flowCore::transformList(cols_to_trans,
                                                            CytoNorm::cytofTransform))

    for(m in names(norm_markers)){

      if(m != names(ym)){
        flowDensity::plotDens(obj = ff_t, channels = c(m,names(ym)),
                              main = basename(file))
        print(paste("plotting", basename(file), m))
      }
    }
  }
  dev.off()
}

