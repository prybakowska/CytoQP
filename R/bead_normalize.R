#' Creates baseline file for bead normalization
#'
#' @description Creates the reference flow frame for which mean beads
#' values will be computed and used during the normalization.
#'
#' @param fcs_files Character, path to fcs files to be normalized.
#' @param beads Character, as in CATALYST::normCytof, "dvs"
#' (for bead masses 140, 151, 153 ,165, 175)
#' or "beta" (for bead masses 139, 141, 159, 169, 175)
#' or a numeric vector of masses. Default is set to "dvs".
#' @param to_plot Logical, indicates if plots should be generated,
#' default set to FALSE
#' @param out_dir Character, pathway to where the plots should be saved,
#' only if argument to_plot = TRUE, default is set to working directory.
#' @param k The same as in CATALYST::normCytof, integer width of the
#' median window used for bead smoothing (affects visualizations only).
#' @param ncells number of cells to be aggregated per each file, defaults is
#' set to 25000 per file.
#' @param ... Additional arguments to pass to CATALYST::normCytof.
#'
#' @return Returns reference, aggregated flow frame.
#'
#' @examples
#' # set input directory (pathway to the files that are going to be normalized)
#' raw_data_dir <- file.path(dir, "RawFiles")
#'
#' # set a directory where bead-normalized fcs files and plots will be saved
#' bead_norm_dir <- file.path(dir, "BeadNorm")
#'
#' # define full pathway to the files that you want to normalize
#' files <- list.files(raw_data_dir,
#'                     pattern = ".FCS$",
#'                     full.names = TRUE)
#'
#' # create baseline file to which all the files will be normalized
#' set.seed(2)
#' ref_sample <- baseline_file(fcs_files = files,
#'                             beads = "dvs",
#'                             out_dir = bead_norm_dir)
#' @export
baseline_file <- function(fcs_files, beads = "dvs", to_plot = FALSE,
                          out_dir = getwd(), k = 80, ncells = 25000, ...){

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "BeadNorm")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  ff <- FlowSOM::AggregateFlowFrames(fileNames = fcs_files,
                                     cTotal = length(fcs_files)*ncells)

  dat <- CATALYST::prepData(ff)

  dat_norm <- CATALYST::normCytof(x = dat,
                                  beads = beads,
                                  remove_beads = TRUE,
                                  norm_to = NULL,
                                  k = k,
                                  plot = to_plot,
                                  verbose = FALSE,
                                  transform = FALSE,
                                  ...)
  ff_ref <- CATALYST::sce2fcs(dat_norm$beads)
  rm(ff)

  if (to_plot == TRUE){

    # plot and save diagnostic plots
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))

    p <- dat_norm$scatter
    ggplot2::ggsave(filename = file.path(plot_dir,"RefBeadGate.png"),
                    plot = p, limitsize = FALSE)

    p <- dat_norm$lines
    ggplot2::ggsave(filename = file.path(plot_dir,"RefBeadLines.png"),
                    plot = p, limitsize = FALSE)
  }
  return(ff_ref)
}


#' Bead-based normalization
#'
#' @description Performs bead-based normalization using beads spiked in
#' the sample. It is based on functions from CATALYST package.
#'
#' @param files Character vector or list with the paths of the raw files.
#' @param cores Number of cores to be used.Works only for not-Widows users.
#' @param markers_to_keep Character vector, marker names to be kept after
#' the normalization, can be full marker name e.g. "CD45" or "CD".
#' If NULL (default) all markers will be normalized and kept in flowframe.
#' Selection of the markers will reduce file volume and speedup the analysis.
#' Non-mass channels like Time, Event_length, Gaussian parameters and in addition
#' palladium barcoding channels are kept if non_mass_ch set to NULL.
#' @param non_mass_channel Character vector, non-mass channels to keep for
#' further analysis. Can be full channel name like Eu151Di or 151.
#' By default "Time" and "event_length" will be always kept in the flow frame.
#' @param beads Character, as in CATALYST::normCytof, "dvs"
#' (for bead masses 140, 151, 153 ,165, 175)
#' or "beta" (for bead masses 139, 141, 159, 169, 175)
#' or a numeric vector of masses. Default is set to "dvs".
#' @param norm_to_ref flow frame, created by baseline_file function to which
#' input data will be normalized, default is set to NULL.
#' @param to_plot Logical if to plot bead gate and bead normalization lines
#' for each file.Defaults is set to TRUE.
#' @param out_dir Character, pathway to where the bead normalized fcs files
#' and plots should be saved, for plots only if argument to_plot = TRUE,
#' default is set to file.path(getwd(), BeadNorm).
#' @param k The same as in CATALYST::normCytof, integer width of the
#' median window used for bead smoothing (affects visualizations only).
#' @param remove_beads Logical, as in CATALYST::normCytof if beads should be
#' removed from fcs files. Default set to TRUE. Note, should be set to FALSE if
#' none of the channels is beads-specific.
#' @param ... Additional arguments to pass to normCytof.
#'
#' @return Save bead-normalized fcs files and plots to out_dir.
#'
#' @examples
#' # set input directory (pathway to the files that are going to be normalized)
#' raw_data_dir <- file.path(dir, "RawFiles")
#'
#' # set a directory where bead-normalized fcs files and plots will be saved
#' bead_norm_dir <- file.path(dir, "BeadNorm")
#'
#' # define full pathway to the files that you want to normalize
#' files <- list.files(raw_data_dir,
#'                     pattern = ".FCS$",
#'                     full.names = TRUE)
#'
#' # create baseline file to which all the files will be normalized
#' set.seed(2)
#' ref_sample <- baseline_file(fcs_files = files,
#'                             beads = "dvs",
#'                             out_dir = bead_norm_dir)
#'
#' # Normalize files
#' bead_normalize(files, cores = 1,
#'                out_dir = bead_norm_dir,
#'                non_mass_channel = NULL,
#'                norm_to_ref = ref_sample,
#'                to_plot = TRUE,
#'                remove_beads = TRUE,
#'                k = 80,
#'                markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir",
#'                                  "Viability","IL", "IFNa",
#'                                    "TNF", "TGF", "MIP", "MCP", "Granz"))
#'
#' @export
bead_normalize <- function(files,
                           cores = 1,
                           markers_to_keep = NULL,
                           non_mass_channel = NULL,
                           beads = "dvs",
                           norm_to_ref = NULL,
                           remove_beads = TRUE,
                           to_plot = TRUE,
                           out_dir = NULL,
                           k = 80,
                           ...){
  # Check parameters
  if(!is(files, "character") & !is(files, "list")) {
    stop("files must be a character vector or a list")
  }

  if (!all(file.exists(files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if (any(!is(cores, "numeric") | cores < 1)){
    stop("cores must be a positive number")
  }

  # Analysis with a single core
  if (cores == 1) {
    lapply(files, function(x) {
      .save_bead_normalize(x,
                           markers_to_keep,
                           beads,
                           non_mass_channel,
                           norm_to_ref,
                           remove_beads,
                           to_plot,
                           out_dir,
                           k)})
  }

  # Parallelized analysis
  else {
    BiocParallel::bplapply(files, function(x) {
      .save_bead_normalize(x,
                           markers_to_keep,
                           non_mass_channel,
                           beads,
                           norm_to_ref,
                           remove_beads,
                           to_plot,
                           out_dir,
                           k)},
      BPPARAM = MulticoreParam(workers = cores))
  }

}


.save_bead_normalize <- function(file,
                                 markers_to_keep,
                                 beads,
                                 non_mass_channel,
                                 norm_to_ref,
                                 remove_beads,
                                 to_plot,
                                 out_dir,
                                 k,
                                 ...) {

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "BeadNorm")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}


  # read fcs file
  ff <- flowCore::read.FCS(file, transformation = FALSE,
                           truncate_max_range = FALSE)

  # bead normalize the files
  ff_norm <- .bead_normalize_ind(flow_frame = ff,
                                 out_dir = out_dir,
                                 non_mass_channel = non_mass_channel,
                                 norm_to_ref = norm_to_ref,
                                 remove_beads = remove_beads,
                                 to_plot = to_plot,
                                 k = k,
                                 beads = beads,
                                 markers_to_keep = markers_to_keep)

  # save normalized FCS files
  flowCore::write.FCS(ff_norm, filename = file.path(out_dir,
                                                    gsub(".FCS","_beadNorm.fcs",
                                                         basename(file),
                                                         ignore.case = TRUE)))
}


# Función interna para leer, normalizar y guardar
.bead_normalize_ind <- function(flow_frame,
                                markers_to_keep = NULL,
                                non_mass_channel = NULL,
                                beads = "dvs",
                                norm_to_ref = NULL,
                                remove_beads = TRUE,
                                to_plot = TRUE,
                                out_dir = getwd(),
                                k = 80,
                                ...){

  if (!is.null(markers_to_keep)){

    matches <- paste(markers_to_keep, collapse="|")

    m_to_keep <- grep(matches, FlowSOM::GetMarkers(flow_frame, flowCore::colnames(flow_frame)),
                      ignore.case = TRUE, value = FALSE)

    if(is.null(non_mass_channel)){
      non_mass_ch <- grep("Time|length|Ce140|151|153|165|175|Center|Offset|Width|
                        |Residual|Pd",
                          flowCore::colnames(flow_frame),
                          ignore.case = TRUE, value = FALSE)

    } else {
      matches_ch <- paste(c(non_mass_channel, "Time", "length"), collapse="|")
      non_mass_ch <- grep(matches_ch,
                          flowCore::colnames(flow_frame),
                          ignore.case = TRUE, value = FALSE)

    }


    channels_to_keep <- c(m_to_keep, non_mass_ch)
    channels_to_keep <- flowCore::colnames(flow_frame)[sort(unique(channels_to_keep))]

    flow_frame <- flow_frame[, channels_to_keep]
  }

  if (is.null(norm_to_ref)){
    warning("the reference file is not defined. Each file will be normalized
       to its own bead mean but not across all files")
  }

  dat <- CATALYST::prepData(flow_frame)

  # normalize the data and remove beads
  dat_norm <- CATALYST::normCytof(x = dat,
                                  beads = beads,
                                  remove_beads = remove_beads,
                                  norm_to = norm_to_ref,
                                  k = k,
                                  plot = TRUE,
                                  transform = FALSE,
                                  ...)

  # convert back to .fcs files and save
  f <- CATALYST::sce2fcs(dat_norm$data)

  filename <- basename(flow_frame@description$FILENAME)

  if (to_plot == TRUE){

    # plot and save diagnostic plots
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))

    p <- dat_norm$scatter
    ggplot2::ggsave(filename = file.path(plot_dir,
                                         gsub(".FCS|.fcs","_beadGate.png",
                                              filename)),
                    plot = p, limitsize = FALSE)

    p <- dat_norm$lines
    ggplot2::ggsave(filename = file.path(plot_dir,
                                         gsub(".FCS|.fcs","_beadLines.png",
                                              filename)),
                    plot = p, limitsize = FALSE)

  }

  f@description$FILENAME <- basename(flow_frame@description$FILENAME)
  f@description$FIL <- basename(flow_frame@description$FILENAME)

  return(f)
}


