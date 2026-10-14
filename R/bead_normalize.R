#' Create a baseline file for bead normalization
#'
#' @description
#' Creates a reference flow frame from multiple FCS files. Mean bead
#' intensities from this reference are used for subsequent bead normalization.
#'
#' @param fcs_files Character vector containing paths to the FCS files used
#' to create the reference.
#' @param beads Numeric vector specifying the bead masses used for
#' normalization, as in \code{CATALYST::normCytof}. Default is
#' \code{c(89, 115, 140, 159, 175, 209)}.
#' @param to_plot Logical indicating whether diagnostic plots should be
#' generated. Default is \code{FALSE}.
#' @param out_dir Character specifying the directory where plots should be
#' saved when \code{to_plot = TRUE}. Default is the working directory.
#' @param k Integer specifying the width of the median window used for bead
#' smoothing, as in \code{CATALYST::normCytof}. Default is 80.
#' @param ncells Number of cells sampled from each FCS file to construct the
#' reference flow frame. Default is 25000 cells per file.
#' @param ... Additional arguments passed to \code{CATALYST::normCytof}.
#'
#' @return
#' A \code{flowFrame} containing the aggregated reference cells used for
#' bead normalization.
#'
#' @examples
#' fcs_files <- system.file(
#'   "extdata",
#'   paste0("Center", 1:4, "_batch1.fcs"),
#'   package = "CytoQP"
#' )
#'
#' set.seed(2)
#'
#' ref_file <- baseline_file(
#'   fcs_files = fcs_files,
#'   out_dir = file.path(tempdir(), "BeadNorm")
#' )
#'
#' @export
baseline_file <- function(fcs_files, beads = c(89, 115, 140, 159, 175, 209), to_plot = FALSE,
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


#' Perform bead-based normalization
#'
#' @description
#' Performs bead-based normalization of FCS files using normalization beads
#' spiked into the samples. Normalization is performed using
#' \code{CATALYST::normCytof}. A reference flow frame generated with
#' \code{baseline_file()} can be supplied to normalize multiple files to the
#' same bead reference.
#'
#' @param files Character vector or list containing paths to the FCS files
#' to be normalized.
#' @param cores Number of cores used for processing. Default is 1.
#' @param markers_to_keep Character vector specifying marker names or patterns
#' to retain after normalization, for example \code{"CD45"} or \code{"CD"}.
#' If \code{NULL}, all channels are retained.
#' @param non_mass_channel Character vector specifying additional non-mass
#' channels to retain. \code{Time} and event-length channels are retained
#' automatically.
#' @param beads Numeric vector specifying the bead masses used for
#' normalization, as in \code{CATALYST::normCytof}. Default is
#' \code{c(89, 115, 140, 159, 175, 209)}.
#' @param norm_to_ref Reference \code{flowFrame}, typically generated with
#' \code{baseline_file()}, to which the input files are normalized.
#' If \code{NULL}, each file is normalized to its own bead means.
#' @param remove_beads Logical indicating whether normalization beads should
#' be removed after normalization. Default is \code{TRUE}.
#' @param to_plot Logical indicating whether bead-gating and normalization
#' diagnostic plots should be generated. Default is \code{TRUE}.
#' @param out_dir Character specifying the directory where normalized FCS
#' files and diagnostic plots are saved. If \code{NULL}, a directory named
#' \code{"BeadNorm"} is created in the current working directory.
#' @param k Integer specifying the width of the median window used for bead
#' smoothing, as in \code{CATALYST::normCytof}. Default is 80.
#' @param ... Additional arguments passed to \code{CATALYST::normCytof}.
#'
#' @return
#' A list containing the results of the normalization for each input FCS file.
#' Normalized FCS files are written to \code{out_dir}.
#'
#' @examples
#' fcs_files <- system.file(
#'   "extdata",
#'   paste0("Center", 1:4, "_batch1.fcs"),
#'   package = "CytoQP"
#' )
#'
#' set.seed(2)
#'
#' ref_file <- baseline_file(
#'   fcs_files = fcs_files,
#'   out_dir = file.path(tempdir(), "BeadNorm")
#' )
#'
#' files_beadnorm <- bead_normalize(
#'   files = fcs_files,
#'   cores = 1,
#'   non_mass_channel = NULL,
#'   norm_to_ref = ref_file,
#'   to_plot = TRUE,
#'   remove_beads = TRUE,
#'   k = 80,
#'   out_dir = file.path(tempdir(), "BeadNorm"),
#'   markers_to_keep = c(
#'     "CD",
#'     "HLA",
#'     "IgD",
#'     "TCR",
#'     "Siglec",
#'     "Ir",
#'     "Viability",
#'     "IL",
#'     "IFNa",
#'     "TNF",
#'     "TGF",
#'     "MIP",
#'     "MCP",
#'     "Granz"
#'   )
#' )
#'
#' @export
bead_normalize <- function(
    files,
    cores = 1,
    markers_to_keep = NULL,
    non_mass_channel = NULL,
    beads = c(89, 115, 140, 159, 175, 209),
    norm_to_ref = NULL,
    remove_beads = TRUE,
    to_plot = TRUE,
    out_dir = NULL,
    k = 80,
    ...) {
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
                           beads,
                           non_mass_channel,
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
                                beads = c(89, 115, 140, 159, 175, 209),
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
                                  plot = to_plot,
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


