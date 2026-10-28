#' Clean Flow Rate and Signal
#'
#' @description
#' Cleans abnormal flow-rate regions and unstable signal regions in cytometry
#' data using flowAI-based flow-rate quality control and flowCut-based signal
#' cleaning.
#'
#' @param files Character vector or list containing paths to FCS files.
#' @param cores Number of cores to use. Default is 1.
#' @param to_plot Character specifying which diagnostic plots should be
#' generated. Options are \code{"All"}, \code{"Flagged Only"}, and
#' \code{"None"}.
#' @param clean_flow_rate Logical indicating whether flow-rate cleaning should
#' be performed. Default is \code{TRUE}.
#' @param clean_signal Logical indicating whether signal cleaning should be
#' performed. Default is \code{TRUE}.
#' @param out_dir Directory where cleaned FCS files and diagnostic plots are
#' written. If \code{NULL}, a \code{"Cleaned"} directory is created in the
#' current working directory.
#' @param alpha Statistical significance level used for flow-rate anomaly
#' detection. Default is 0.01.
#' @param data_type Character specifying the cytometry data type:
#' \code{"MC"} for mass cytometry or \code{"FC"} for flow cytometry.
#' @param channels_to_clean Optional character vector specifying channels used
#' for signal cleaning.
#' @param Segment Number of events per segment used by flowCut.
#' Default is 1000.
#' @param arcsine_transform Logical indicating whether the data should be
#' transformed before signal cleaning. Default is \code{TRUE}.
#' @param non_used_bead_ch Optional character vector specifying bead channels
#' that should be excluded from signal cleaning.
#' @param MaxPercCut Maximum proportion of events that flowCut is allowed to
#' remove. Default is 0.5.
#' @param UseOnlyWorstChannels Logical passed to flowCut indicating whether
#' cleaning should be based on automatically detected worst channels.
#' @param AllowFlaggedRerun Logical passed to flowCut indicating whether
#' flagged files can be analyzed a second time.
#' @param AlwaysClean Logical passed to flowCut indicating whether cleaning
#' should also be performed for relatively stable files.
#' @param ... Additional arguments passed to the signal-cleaning procedure.
#'
#' @return
#' A list of paths to the cleaned FCS files (V1-compatible).
#'
#' @examples
#' fcs_file <- system.file(
#'   "extdata",
#'   "Center1_batch1.fcs",
#'   package = "CytoQP"
#' )
#'
#' clean_files(
#'   files = fcs_file,
#'   cores = 1,
#'   to_plot = "None",
#'   data_type = "MC",
#'   out_dir = file.path(tempdir(), "Cleaned")
#' )
#'
#' @import ggplot2
#'
#' @export
clean_files <- function(files,
                        cores = 1,
                        to_plot = "All",
                        clean_flow_rate = TRUE,
                        clean_signal = TRUE,
                        out_dir = NULL,
                        alpha = 0.01,
                        data_type = "MC",
                        channels_to_clean = NULL,
                        Segment = 1000,
                        arcsine_transform = TRUE,
                        non_used_bead_ch = NULL,
                        MaxPercCut = 0.5,
                        UseOnlyWorstChannels = TRUE,
                        AllowFlaggedRerun = TRUE,
                        AlwaysClean = TRUE,
                        ...) {
  # Check parameters
  if(!is(files, "character") & !is(files, "list")) {
    stop("files must be a character vector or a list")
  }

  if(is(files, "list")){
    files <- unlist(files)
  }

  if (any(!is(cores, "numeric") | cores < 1)){
    stop("cores must be a positive number")
  }

  if(!is(clean_flow_rate, "logical")){
    stop("clean_flow_rate must be logical")
  }

  if(!is(clean_signal, "logical")){
    stop("clean_signal must be logical")
  }

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Cleaned")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

  if (!all(file.exists(files))){
    stop("incorrect file path, the fcs file does not exist")
  }


  # Analysis with a single core
  if (cores == 1) {

    output_files <- lapply(
      files,
      function(x) {
        .save_bead_clean(
          file = x,
          to_plot = to_plot,
          clean_flow_rate = clean_flow_rate,
          clean_signal = clean_signal,
          out_dir = out_dir,
          alpha = alpha,
          data_type = data_type,
          channels_to_clean = channels_to_clean,
          Segment = Segment,
          arcsine_transform = arcsine_transform,
          non_used_bead_ch = non_used_bead_ch,
          MaxPercCut = MaxPercCut,
          UseOnlyWorstChannels = UseOnlyWorstChannels,
          AllowFlaggedRerun = AllowFlaggedRerun,
          AlwaysClean = AlwaysClean,
          ...
        )
      }
    )

  } else {

    # Windows-compatible parallelization
    output_files <- BiocParallel::bplapply(
      files,
      function(x) {
        .save_bead_clean(
          file = x,
          to_plot = to_plot,
          clean_flow_rate = clean_flow_rate,
          clean_signal = clean_signal,
          out_dir = out_dir,
          alpha = alpha,
          data_type = data_type,
          channels_to_clean = channels_to_clean,
          Segment = Segment,
          arcsine_transform = arcsine_transform,
          non_used_bead_ch = non_used_bead_ch,
          MaxPercCut = MaxPercCut,
          UseOnlyWorstChannels = UseOnlyWorstChannels,
          AllowFlaggedRerun = AllowFlaggedRerun,
          AlwaysClean = AlwaysClean,
          ...
        )
      },
      BPPARAM = BiocParallel::MulticoreParam(
        workers = cores
      )
    )
  }

  # V1 returned the list produced by lapply()/bplapply().
  # Preserve the list contract for scripts that use files_clean[[i]].
  return(output_files)
}



#' Clean Signal with PeacoQC
#'
#' @description
#' Cleans unstable signal regions in cytometry data using PeacoQC.
#' The selected channels are transformed before quality control when
#' \code{arcsine_transform = TRUE}, while the original untransformed events
#' corresponding to PeacoQC good cells are written to the output FCS files.
#'
#' @param files Character vector or list containing paths to FCS files.
#' @param channels_to_clean Character vector specifying the channels used
#' for PeacoQC signal cleaning.
#' @param cores Number of cores to use. Default is 1.
#' @param plot Passed to \code{PeacoQC::PeacoQC()}. Controls generation of
#' diagnostic plots. Default is 20.
#' @param out_dir Directory where cleaned FCS files and PeacoQC output are
#' written. Default is the current working directory.
#' @param name_dir Name of the PeacoQC results directory.
#' Default is \code{"Cleaned"}.
#' @param MAD MAD threshold used by PeacoQC. Increasing this value makes
#' the algorithm less strict. Default is 6.
#' @param IT_limit IsolationTree threshold used by PeacoQC. Increasing this
#' value makes the algorithm less strict. Default is 0.6.
#' @param remove_zeros Logical indicating whether zero values should be removed
#' before peak detection. Default is \code{TRUE}.
#' @param arcsine_transform Logical indicating whether the data should be
#' transformed before PeacoQC analysis. Default is \code{TRUE}.
#' @param data_type Character specifying the cytometry data type:
#' \code{"MC"} for mass cytometry or \code{"FC"} for flow cytometry.
#' @param suffix_old Character string to replace in input filenames.
#' Default is \code{"beadNorm"} (as in V1).
#' @param suffix_new Character string used in cleaned output filenames.
#' Default is \code{"cleaned"}.
#'
#' @return
#' A list of paths to the cleaned FCS files (V1-compatible).
#'
#' @examples
#' fcs_file <- system.file(
#'   "extdata",
#'   "Center1_batch1.fcs",
#'   package = "CytoQP"
#' )
#'
#' ff <- flowCore::read.FCS(
#'   fcs_file,
#'   transformation = FALSE
#' )
#'
#' channels_to_clean <- find_mass_ch(
#'   ff,
#'   value = TRUE
#' )
#'
#' files_clean <- clean_files_PQC(
#'   files = fcs_file,
#'   channels_to_clean = channels_to_clean,
#'   cores = 1,
#'   plot = FALSE,
#'   data_type = "MC",
#'   out_dir = file.path(tempdir(), "PeacoQC")
#' )
#'
#' @export
clean_files_PQC <- function(
    files,
    channels_to_clean,
    cores = 1,
    plot = 20,
    out_dir = ".",
    name_dir = "Cleaned",
    MAD = 6,
    IT_limit = 0.6,
    remove_zeros = TRUE,
    arcsine_transform = TRUE,
    data_type = "MC",
    suffix_old = "beadNorm",
    suffix_new = "cleaned"
) {

  # Check input files
  if (
    !is(files, "character") &&
    !is(files, "list")
  ) {
    stop(
      "files must be a character vector or a list.",
      call. = FALSE
    )
  }

  if (is(files, "list")) {
    files <- unlist(
      files,
      use.names = FALSE
    )
  }

  if (
    length(cores) != 1 ||
    !is.numeric(cores) ||
    cores < 1
  ) {
    stop(
      "cores must be a positive number.",
      call. = FALSE
    )
  }

  if (!all(file.exists(files))) {
    stop(
      "One or more FCS files do not exist.",
      call. = FALSE
    )
  }

  if (!data_type %in% c("MC", "FC")) {
    stop(
      "data_type must be either 'MC' or 'FC'.",
      call. = FALSE
    )
  }

  # V1 used the current working directory when out_dir was NULL.
  if (is.null(out_dir)) out_dir <- getwd()
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Function applied to each FCS file
  clean_one_file <- function(file) {

    message(
      "Cleaning signal with PeacoQC: ",
      basename(file)
    )

    ff <- flowCore::read.FCS(
      file,
      transformation = FALSE
    )

    channels_to_transform <- find_mass_ch(
      ff,
      value = FALSE
    )

    # Transform data for PeacoQC
    if (arcsine_transform) {

      if (data_type == "MC") {

        ff_t <- flowCore::transform(
          ff,
          flowCore::transformList(
            flowCore::colnames(ff)[channels_to_transform],
            CytoNorm::cytofTransform
          )
        )

      } else {

        ff_t <- flowCore::transform(
          ff,
          flowCore::transformList(
            flowCore::colnames(ff)[channels_to_transform],
            flowCore::arcsinhTransform(
              a = 0,
              b = 1 / 150,
              c = 0
            )
          )
        )
      }

    } else {

      ff_t <- ff
    }

    # Run PeacoQC
    pqc <- PeacoQC::PeacoQC(
      ff_t,
      channels = channels_to_clean,
      plot = plot,
      output_directory = out_dir,
      name_directory = name_dir,
      MAD = MAD,
      IT_limit = IT_limit,
      remove_zeros = remove_zeros,
      save_fcs = FALSE,
      time_unit = 10000
    )

    # Define output path
    input_name <- basename(file)

    if (grepl(suffix_old, input_name, fixed = TRUE)) {

      output_name <- sub(
        suffix_old,
        suffix_new,
        input_name,
        fixed = TRUE
      )

    } else {
      # V1 did not append a suffix when suffix_old was absent.
      output_name <- input_name
    }

    output_file <- file.path(
      out_dir,
      output_name
    )

    # Never overwrite the input FCS file in place.
    if (normalizePath(file, mustWork = TRUE) ==
        normalizePath(output_file, mustWork = FALSE)) {
      stop("PeacoQC output would overwrite the input FCS file: ",
           file, call. = FALSE)
    }

    # Save original-scale good cells
    flowCore::write.FCS(
      ff[pqc$GoodCells, ],
      filename = output_file
    )

    return(output_file)
  }

  # Run cleaning
  if (cores == 1) {

    output_files <- lapply(
      files,
      clean_one_file
    )

  } else {

    output_files <- BiocParallel::bplapply(
      files,
      clean_one_file,
      BPPARAM = BiocParallel::MulticoreParam(
        workers = cores
      )
    )
  }

  # V1 returned the list produced by lapply()/bplapply().
  # Preserve the list contract for scripts that use files_clean[[i]].
  return(output_files)
}

.save_bead_clean <- function(
    file,
    to_plot = "All",
    clean_flow_rate = TRUE,
    clean_signal = TRUE,
    out_dir = getwd(),
    alpha = 0.01,
    data_type = "MC",
    channels_to_clean = NULL,
    Segment = 1000,
    arcsine_transform = TRUE,
    non_used_bead_ch = NULL,
    MaxPercCut = 0.5,
    UseOnlyWorstChannels = TRUE,
    AllowFlaggedRerun = TRUE,
    AlwaysClean = TRUE,
    ...
) {

  # Read FCS file
  ff <- flowCore::read.FCS(
    filename = file,
    transformation = FALSE
  )

  # Clean flow rate
  if (clean_flow_rate) {
    message(
      "Cleaning flow rate: ",
      basename(file)
    )

    ff <- .clean_flow_rate_ind(
      flow_frame = ff,
      to_plot = to_plot,
      out_dir = out_dir,
      alpha = alpha,
      data_type = data_type
    )
  }

  # Clean signal
  if (clean_signal) {
    message(
      "Cleaning signal: ",
      basename(file)
    )

    ff <- .clean_signal_ind(
      flow_frame = ff,
      channels_to_clean = channels_to_clean,
      to_plot = to_plot,
      Segment = Segment,
      out_dir = out_dir,
      arcsine_transform = arcsine_transform,
      non_used_bead_ch = non_used_bead_ch,
      MaxPercCut = MaxPercCut,
      UseOnlyWorstChannels = UseOnlyWorstChannels,
      AllowFlaggedRerun = AllowFlaggedRerun,
      AlwaysClean = AlwaysClean,
      data_type = data_type,
      ...
    )
  }

  # Define output file
  # V1 accepted *_beadNorm.fcs; newer V2 workflows often supply
  # *_compensated.fcs. Preserve both filename conventions.
  output_name <- basename(file)
  output_name <- gsub("_beadNorm", "_cleaned", output_name, fixed = TRUE)
  output_name <- gsub("_compensated", "_cleaned", output_name, fixed = TRUE)
  output_file <- file.path(out_dir, output_name)

  # Write cleaned FCS file
  flowCore::write.FCS(
    ff,
    filename = output_file
  )

  return(output_file)
}


.clean_flow_rate_ind <- function(flow_frame, to_plot = TRUE,
                                 out_dir = getwd(), alpha = 0.01,
                                 data_type = "MC") {

  if (data_type == "MC"){
    time_division <- 100
    timestep <- 0.01
  }
  else if (data_type == "FC") {
    time_division <- 1
    fcs_keywords <- flowCore::keyword(flow_frame)
    word <- grep("TIMESTEP", names(fcs_keywords), ignore.case = TRUE)
    if (length(word) == 0L) {
      stop("FC FCS file does not contain a TIMESTEP keyword.", call. = FALSE)
    }
    timestep <- suppressWarnings(as.numeric(fcs_keywords[[word[1]]]))
    if (length(timestep) != 1L || !is.finite(timestep) || timestep <= 0) {
      stop("Invalid TIMESTEP keyword in FC FCS file.", call. = FALSE)
    }
  }
  else{
    stop("type of data MC or FC needs to be specified")
  }

  if (to_plot == "None") {
    to_plot <- FALSE
  }
  else {
    to_plot <- TRUE

  }

  flow_frame@exprs[, "Time"] <- flow_frame@exprs[, "Time"]/time_division


  FlowRateData <- .flow_rate_bin_adapted(flow_frame,
                                         timeCh = "Time",
                                         timestep = timestep)

  FlowRateQC <- .flow_rate_check_adapted(x = flow_frame,
                                         FlowRateData = FlowRateData,
                                         alpha = alpha,
                                         use_decomp = TRUE)

  if (to_plot == TRUE){

    out_dir <- file.path(out_dir, "FlowRateCleaning")
    if(!dir.exists(out_dir)){
      dir.create(out_dir)
    }

    png(file.path(out_dir,
                  gsub(".fcs", "_flowAI.png",
                       basename(flow_frame@description$FILENAME),
                       ignore.case = TRUE)),
        width = 800,
        height = 600)
    if(data_type == "MC"){
      FlowRateQC$res_fr_QC[,1] <- timestep
    }
    p <- .plot_flowrate(FlowRateQC, data_type = data_type)
    print(p)
    dev.off()
  }

  flow_frame_cl <- flow_frame[FlowRateQC$goodCellIDs,]
  flow_frame_cl@exprs[,"Time"] <- flow_frame_cl@exprs[,"Time"]*time_division

  return(flow_frame_cl)
}

#' flow_rate_bin_adapted
#'
#' @description flow_rate_bin_adapted
#'
#' @param x flow frame
#' @param second_fraction the fraction of the seconds used in the data.
#' @param timeCh Time channel.
#' @param timestep Step in time measurement
#'
#' @return timeFlowData, the cell assignment to the bins.
#'
#' @references this code is strongly based on flowAI::flow_rate_bin.
.flow_rate_bin_adapted <- function (x, second_fraction = 0.1, timeCh = timeCh,
                                    timestep = timestep)
{
  xx <- flowCore::exprs(x)[, timeCh]
  idx <- c(1:nrow(x))
  endsec <- ceiling(timestep * max(xx))
  lenx <- length(xx)
  secbegin <- as.numeric(gsub("(.*)(\\.)(.{0}).*", "\\1\\2\\3", xx[1]))
  tbins <- seq(secbegin, endsec/timestep, by = as.numeric(second_fraction)/timestep)
  if (tail(tbins, n=1) < endsec/timestep){
    tbins <- c(tbins, tail(tbins, n=1) + 10)
  }
  if (secbegin == 0){
    secbegin2 <- 0
  } else {
    secbegin2 <- as.numeric(gsub("(.*)(\\.)(.{1}).*", "\\1\\2\\3", xx[1]/100))
  }

  secbin <- seq(secbegin2, endsec, by = as.numeric(second_fraction))
  minbin <- round(secbin/60, 3)
  nrBins <- length(tbins) - 1
  tbCounts <- c(0, hist(xx, tbins, plot = FALSE)$counts)
  expEv <- lenx/(nrBins)
  binID <- do.call(c, mapply(rep, x = 1:length(tbCounts),
                             times = tbCounts, SIMPLIFY = FALSE))
  if (length(idx) != length(binID))
    stop("length of cell ID not equal length of bin ID")
  timeFlowData <- list(frequencies = cbind(tbins, minbin,
                                           secbin, tbCounts),
                       cellBinID = data.frame(cellID = idx,
                                              binID = binID),
                       info = data.frame(second_fraction = second_fraction,
                                         expFrequency = expEv, bins = nrBins))
  return(timeFlowData)
}


#' Plots flow rate
#'
#' @description plots flow rate for .fcs files
#'
#' @param data_type Character, if "MC" (mass cytometry) or "FC" data are analyzed
#' @param FlowRateQC list obtained using flowAI:::flow_rate_check function
#'
#' @return xgraph plot
#'
#' @references this code is adapted from the flowAI:::flow_rate_plot()
#' Monaco, G., Chen, H., Poidinger, M., Chen, J., de Magalhães, J.P.,
#' and Larbi, A. (2016). flowAI: automatic and interactive anomaly discerning
#' tools for flow cytometry data. Bioinformatics 32, 2473–2480.
.plot_flowrate <- function (FlowRateQC, data_type = "MC")
{
  if (data_type == "MC"){
    lab <- "Time (10 * Seconds)"
  } else {
    lab <- "Time (Seconds)"
  }
  second_fraction <- FlowRateQC$res_fr_QC$second_fraction
  num_obs = FlowRateQC$res_fr_QC$num_obs
  frequencies = as.data.frame(FlowRateQC$frequencies)
  anoms = as.data.frame(FlowRateQC$anoms)
  anoms_points = as.data.frame(cbind(sec_anom = frequencies$secbin[anoms$index],
                                     count_anom = anoms$anoms))
  xgraph <- ggplot2::ggplot(frequencies, ggplot2::aes_string(x = "secbin", y = "tbCounts")) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), text = element_text(size = 30)) +
    ggplot2::geom_line(colour = "darkblue")
  xgraph <- xgraph + ggplot2::labs(x = lab, y = paste0("Number of events per 1/",
                                                       1/second_fraction, " of a second"))
  if (!is.null(anoms_points)) {
    xgraph <- xgraph + ggplot2::geom_point(data = anoms_points, aes_string(x = "sec_anom",
                                                                           y = "count_anom"),
                                           color = "green4", size = 5,
                                           shape = 1, stroke = 3)
  }
  return(xgraph)
}

.clean_signal_ind <- function(flow_frame,
                              channels_to_clean = NULL,
                              to_plot = "All",
                              Segment = 1000,
                              out_dir = getwd(),
                              arcsine_transform = TRUE,
                              non_used_bead_ch = NULL,
                              MaxPercCut = 0.5,
                              UseOnlyWorstChannels = TRUE,
                              AllowFlaggedRerun = TRUE,
                              AlwaysClean = TRUE,
                              data_type = "MC",
                              ...){

  channels_to_transform <- find_mass_ch(flow_frame, value = FALSE)

  if (arcsine_transform){

    if(data_type == "MC"){
      ff_t <- flowCore::transform(flow_frame,
                                  flowCore::transformList(flowCore::colnames(flow_frame)[channels_to_transform],
                                                          CytoNorm::cytofTransform))
    }
    else if (data_type == "FC"){
      ff_t <- flowCore::transform(flow_frame,
                                  flowCore::transformList(flowCore::colnames(flow_frame)[channels_to_transform],
                                                          flowCore::arcsinhTransform(a = 0, b = 1/150, c = 0)))

    }
    else {
      stop("specify data type MC or FC")
    }

  }
  else {
    ff_t <- flow_frame
  }

  if (!is.null(channels_to_clean)) {

    ch_to_clean <- which(
      flowCore::colnames(flow_frame) %in% channels_to_clean
    )

    if (!("TIME" %in% toupper(
      flowCore::colnames(flow_frame)[ch_to_clean]
    ))) {

      ind_time <- grep(
        "TIME",
        toupper(flowCore::colnames(flow_frame))
      )

      channels <- unique(
        sort(c(ch_to_clean, ind_time))
      )

    } else {

      channels <- ch_to_clean
    }
  }
  else {

    if (is.null(non_used_bead_ch)) {
      non_bead_ch <- "140"
    }
    else {
      non_bead_ch <- paste(non_used_bead_ch, collapse="|")
    }

    ind_Time <- grep("TIME", flowCore::colnames(flow_frame), value = T, ignore.case = T)
    ch_to_clean <- c(ind_Time, find_mass_ch(flow_frame, value = TRUE))
    ind_nonbeads <- grep(non_bead_ch, flowCore::colnames(flow_frame), value = TRUE)
    channels <- ch_to_clean[!(ch_to_clean %in% ind_nonbeads)]
    channels <- grep(paste(channels, collapse = "|"), flowCore::colnames(flow_frame))
  }

  out_dir <- file.path(out_dir, "SignalCleaning")
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }

  cleaned_data <- flowCut::flowCut(f = ff_t,
                                   Segment = Segment,
                                   MaxPercCut = MaxPercCut,
                                   Channels = channels,
                                   FileID = gsub("_beadNorm|_compensated", "_flowCutCleaned",
                                                 basename(flow_frame@description$FILENAME)),
                                   Plot = to_plot,
                                   Directory = out_dir,
                                   UseOnlyWorstChannels = UseOnlyWorstChannels,
                                   AllowFlaggedRerun = AllowFlaggedRerun,
                                   AlwaysClean = AlwaysClean, ...)

  ff_t_clean <- cleaned_data$frame

  if (arcsine_transform){

    if(data_type == "MC"){
      ff_clean <- flowCore::transform(ff_t_clean,
                                      flowCore::transformList(flowCore::colnames(ff_t_clean)[channels_to_transform],
                                                              CytoNorm::cytofTransform.reverse))
    }
    else if (data_type == "FC"){
      ff_clean <- flowCore::transform(ff_t_clean,
                                      flowCore::transformList(flowCore::colnames(ff_t_clean)[channels_to_transform],
                                                              CytoNorm::cytofTransform.reverse(x = 150)))

    }
    else {
      stop("specify data type MC or FC")
    }

  }
  else {
    ff_clean <- ff_t_clean
  }

  return(ff_clean)
}

#' @references this code uses internal functions from flowAI package
#' Monaco, G., Chen, H., Poidinger, M., Chen, J., de Magalhães, J.P.,
#' and Larbi, A. (2016). flowAI: automatic and interactive anomaly discerning
#' tools for flow cytometry data. Bioinformatics 32, 2473–2480.
.flow_rate_check_adapted <- function (x, FlowRateData, alpha = alpha,
                                      use_decomp = use_decomp)
{
  fr_frequences <- FlowRateData$frequencies
  fr_cellBinID <- FlowRateData$cellBinID
  second_fraction <- FlowRateData$info["second_fraction"]
  if (length(unique(fr_frequences[, 2])) == 1) {
    fr_autoqc <- NULL
  }
  else {
    fr_autoqc <- .anomaly_detection_addapted(fr_frequences[, "tbCounts"],
                                             alpha = alpha,
                                             use_decomp = use_decomp)
  }
  if (is.null(fr_autoqc) || is.null(fr_autoqc$anoms)) {
    badPerc <- 0
    newx <- x
    goodCellIDs <- fr_cellBinID$cellID
    badCellIDs <- NULL
  }
  else {
    goodCellIDs <- fr_cellBinID$cellID[!(fr_cellBinID$binID %in%
                                           fr_autoqc$anoms$index)]
    badCellIDs <- setdiff(fr_cellBinID$cellID, goodCellIDs)
    badPerc <- round(1 - (length(goodCellIDs)/nrow(fr_cellBinID)),
                     4)
    params <- flowCore::parameters(x)
    keyval <- flowCore::keyword(x)
    sub_exprs <- flowCore::exprs(x)
    sub_exprs <- sub_exprs[goodCellIDs, ]
    newx <- flowCore::flowFrame(exprs = sub_exprs, parameters = params,
                                description = keyval)
  }
  cat(paste0(100 * badPerc, "% of anomalous cells detected in the flow rate
             check. \n"))
  return(list(anoms = fr_autoqc$anoms, frequencies = fr_frequences,
              FRnewFCS = newx, goodCellIDs = goodCellIDs,
              badCellIDs = badCellIDs,
              res_fr_QC = data.frame(second_fraction = second_fraction,
                                     num_obs = fr_autoqc$num_obs,
                                     badPerc = badPerc)))
}


.anomaly_detection_addapted <- function (x, max_anoms = 0.49,
                                         direction = "both", alpha = 0.01,
                                         use_decomp = TRUE, period = 1,
                                         verbose = FALSE)
{
  if (is.vector(x) && is.numeric(x)) {
    x <- ts(x, frequency = period)
  }
  else if (is.ts(x)) {
  }
  else {
    stop("data must be a time series object or a vector that holds numeric
         values.")
  }
  if (length(rle(is.na(c(NA, x, NA)))$values) > 3) {
    stop("Data contains non-leading NAs. We suggest replacing NAs with
         interpolated values (see na.approx in Zoo package).")
  }
  else {
    x <- na.omit(x)
  }
  if (max_anoms > 0.49) {
    stop(paste("max_anoms must be less than 50% of the data points (max_anoms =",
               round(max_anoms * length(x), 0), " data_points =",
               length(x), ")."))
  }
  if (!direction %in% c("pos", "neg", "both")) {
    stop("direction options are: pos | neg | both.")
  }
  if (!(0.01 <= alpha || alpha <= 0.1)) {
    print("Warning: alpha is the statistical significance level, and is usually
          between 0.01 and 0.1")
  }
  if (is.null(period)) {
    stop("Period must be set to the number of data points in a single period")
  }
  if (use_decomp) {
    x_cf <- .cffilter_adapted(x)
    med_t <- trunc(median(x_cf$trend))
    sign_n <- sign(x_cf$trend - med_t)
    sign_n[which(sign_n == 0)] <- 1
    x_2 <- as.vector(trunc(abs(x - med_t) + abs(x_cf$cycle)) *
                       sign_n)
    trend <- x_cf$trend
  }
  else {
    x_2 <- as.vector(x - median(x))
    trend <- x
  }
  anomaly_direction = switch(direction, pos = data.frame(one_tail = TRUE,
                                                         upper_tail = TRUE),
                             neg = data.frame(one_tail = TRUE,
                                              upper_tail = FALSE),
                             both = data.frame(one_tail = FALSE,
                                               upper_tail = TRUE))
  n <- length(x_2)
  data_det <- data.frame(index = 1:length(x), values = x_2,
                         or_values = trend)
  max_outliers <- trunc(n * max_anoms)
  func_ma <- match.fun(median)
  func_sigma <- match.fun(IQR)
  R_idx <- 1L:max_outliers
  num_anoms <- 0L
  one_tail <- anomaly_direction$one_tail
  upper_tail <- anomaly_direction$upper_tail
  for (i in 1L:max_outliers) {
    if (verbose)
      message(paste(i, "/", max_outliers, "completed"))
    if (one_tail) {
      if (upper_tail) {
        ares <- data_det[[2L]] - func_ma(data_det[[2L]])
      }
      else {
        ares <- func_ma(data_det[[2L]]) - data_det[[2L]]
      }
    }
    else {
      ares = abs(data_det[[2L]] - func_ma(data_det[[2L]]))
    }
    data_sigma <- func_sigma(ares)
    if (data_sigma == 0)
      break
    ares <- ares/data_sigma
    R <- max(ares)
    temp_max_idx <- which(ares == R)[1L]
    R_idx[i] <- data_det[[1L]][temp_max_idx]
    data_det <- data_det[-which(data_det[[1L]] == R_idx[i]),
    ]
    if (one_tail) {
      p <- 1 - alpha/(n - i + 1)
    }
    else {
      p <- 1 - alpha/(2 * (n - i + 1))
    }
    t <- qt(p, (n - i - 1L))
    lam <- t * (n - i)/sqrt((n - i - 1 + t^2) * (n - i +
                                                   1))
    if (R > lam)
      num_anoms <- i
  }
  if (num_anoms > 0) {
    R_idx <- R_idx[1L:num_anoms]
    all_data <- data.frame(index = 1:length(x), anoms = x)
    anoms_data <- subset(all_data, (all_data[[1]] %in% R_idx))
  }
  else {
    anoms_data <- NULL
  }
  return(list(anoms = anoms_data, num_obs = n))
}

.cffilter_adapted <- function (x, pl = NULL, pu = NULL, root = FALSE,
                               drift = FALSE,
                               type = c("asymmetric", "symmetric",
                                        "fixed", "baxter-king",
                                        "trigonometric"), nfix = NULL,
                               theta = 1)
{
  type = match.arg(type)
  if (is.null(root))
    root <- FALSE
  if (is.null(drift))
    drift <- FALSE
  if (is.null(theta))
    theta <- 1
  if (is.null(type))
    type <- "asymmetric"
  if (is.ts(x))
    freq = frequency(x)
  else freq = 1
  if (is.null(pl)) {
    if (freq > 1)
      pl = trunc(freq * 1.5)
    else pl = 2
  }
  if (is.null(pu))
    pu = trunc(freq * 8)
  if (is.null(nfix))
    nfix = freq * 3
  nq = length(theta) - 1
  b = 2 * pi/pl
  a = 2 * pi/pu
  xname = deparse(substitute(x))
  xo = x
  x = as.matrix(x)
  n = nrow(x)
  nvars = ncol(x)
  if (n < 5)
    warning("# of observations < 5")
  if (n < (2 * nq + 1))
    stop("# of observations must be at least 2*q+1")
  if (pu <= pl)
    stop("pu must be larger than pl")
  if (pl < 2) {
    warning("pl less than 2 , reset to 2")
    pl = 2
  }
  if (root != 0 && root != 1)
    stop("root must be 0 or 1")
  if (drift < 0 || drift > 1)
    stop("drift must be 0 or 1")
  if ((type == "fixed" || type == "baxter-king") && nfix <
      1)
    stop("fixed lag length must be >= 1")
  if (type == "fixed" & nfix < nq)
    stop("fixed lag length must be >= q")
  if ((type == "fixed" || type == "baxter-king") && nfix >=
      n/2)
    stop("fixed lag length must be < n/2")
  if (type == "trigonometric" && (n - 2 * floor(n/2)) != 0)
    stop("trigonometric regressions only available for even n")
  theta = as.matrix(theta)
  m1 = nrow(theta)
  m2 = ncol(theta)
  if (m1 > m2)
    th = theta
  else th = t(theta)
  g = convolve(th, th, type = "open")
  cc = g[(nq + 1):(2 * nq + 1)]
  j = 1:(2 * n)
  B = as.matrix(c((b - a)/pi, (sin(j * b) - sin(j * a))/(j *
                                                           pi)))
  R = matrix(0, n, 1)
  if (nq > 0) {
    R0 = B[1] * cc[1] + 2 * t(B[2:(nq + 1)]) * cc[2:(nq +
                                                       1)]
    R[1] = pi * R0
    for (i in 2:n) {
      dj = Bge(i - 2, nq, B, cc)
      R[i] = R[i - 1] - dj
    }
  }
  else {
    R0 = B[1] * cc[1]
    R[1] = pi * R0
    for (j in 2:n) {
      dj = 2 * pi * B[j - 1] * cc[1]
      R[j] = R[j - 1] - dj
    }
  }
  AA = matrix(0, n, n)
  if (type == "asymmetric") {
    if (nq == 0) {
      for (i in 1:n) {
        AA[i, i:n] = t(B[1:(n - i + 1)])
        if (root)
          AA[i, n] = R[n + 1 - i]/(2 * pi)
      }
      AA[1, 1] = AA[n, n]
      AAu = AA
      AAu[!upper.tri(AAu)] <- 0
      AA = AA + .flipud_adapted(.fliplr_adapted(AAu))
    }
    else {
      A = Abuild(n, nq, g, root)
      Ainv = solve(A)
      for (np in 0:ceiling(n/2 - 1)) {
        d = matrix(0, n, 1)
        ii = 0
        for (jj in (np - root):(np + 1 + root - n)) {
          ii = ii + 1
          d[ii] = Bge(jj, nq, B, cc)
        }
        if (root == 1)
          d[n - 1] = R[n - np]
        Bhat = Ainv %*% d
        AA[np + 1, ] = t(Bhat)
      }
      AA[(ceiling(n/2) + 1):n, ] = .flipud_adapted(.fliplr_adapted(AA[1:floor(n/2),
      ]))
    }
  }
  if (type == "symmetric") {
    if (nq == 0) {
      for (i in 2:ceiling(n/2)) {
        np = i - 1
        AA[i, i:(i + np)] = t(B[1:(1 + np)])
        if (root)
          AA[i, i + np] = R[np + 1]/(2 * pi)
        AA[i, (i - 1):(i - np)] = AA[i, (i + 1):(i +
                                                   np)]
      }
      AA[(ceiling(n/2) + 1):n, ] = .flipud_adapted(.fliplr_adapted(AA[1:floor(n/2),
      ]))
    }
    else {
      for (np in nq:ceiling(n/2 - 1)) {
        nf = np
        nn = 2 * np + 1
        A = Abuild(nn, nq, g, root)
        Ainv = solve(A)
        d = matrix(0, nn, 1)
        ii = 0
        for (jj in (np - root):(-nf + root)) {
          ii = ii + 1
          d[ii] = Bge(jj, nq, B, cc)
        }
        if (root)
          d[nn - 1] = R[nf + 1]
        Bhat = Ainv %*% d
        AA[np + 1, 1:(2 * np + 1)] = t(Bhat)
      }
      AA[(ceiling(n/2) + 1):n, ] = .flipud_adapted(.fliplr_adapted(AA[1:floor(n/2),
      ]))
    }
  }
  if (type == "fixed") {
    if (nq == 0) {
      bb = matrix(0, 2 * nfix + 1, 1)
      bb[(nfix + 1):(2 * nfix + 1)] = B[1:(nfix + 1)]
      bb[nfix:1] = B[2:(nfix + 1)]
      if (root) {
        bb[2 * nfix + 1] = R[nfix + 1]/(2 * pi)
        bb[1] = R[nfix + 1]/(2 * pi)
      }
      for (i in (nfix + 1):(n - nfix)) AA[i, (i - nfix):(i +
                                                           nfix)] = t(bb)
    }
    else {
      nn = 2 * nfix + 1
      A = Abuild(nn, nq, g, root)
      Ainv = solve(A)
      d = matrix(0, nn, 1)
      ii = 0
      for (jj in (nfix - root):(-nfix + root)) {
        ii = ii + 1
        d[ii] = Bge(jj, nq, B, cc)
      }
      if (root)
        d[nn - 1] = R[nn - nfix]
      Bhat = Ainv %*% d
      for (ii in (nfix + 1):(n - nfix)) AA[ii, (ii - nfix):(ii +
                                                              nfix)] = t(Bhat)
    }
  }
  if (type == "baxter-king") {
    bb = matrix(0, 2 * nfix + 1, 1)
    bb[(nfix + 1):(2 * nfix + 1)] = B[1:(nfix + 1)]
    bb[nfix:1] = B[2:(nfix + 1)]
    bb = bb - sum(bb)/(2 * nfix + 1)
    for (i in (nfix + 1):(n - nfix)) AA[i, (i - nfix):(i +
                                                         nfix)] = t(bb)
  }
  if (type == "trigonometric") {
    jj = 1:(n/2)
    jj = jj[((n/pu) <= jj & jj <= (n/pl) & jj < (n/2))]
    if (!any(jj))
      stop("frequency band is empty in trigonometric regression")
    om = 2 * pi * jj/n
    if (pl > 2) {
      for (t in 1:n) {
        for (k in n:1) {
          l = t - k
          tmp = sum(cos(om * l))
          AA[t, k] = tmp
        }
      }
    }
    else {
      for (t in 1:n) {
        for (k in n:1) {
          l = t - k
          tmp = sum(cos(om * l))
          tmp2 = (cos(pi * (t - l)) * cos(pi * t))/2
          AA[t, k] = tmp + tmp2
        }
      }
    }
    AA = AA * 2/n
  }
  if (root) {
    tst = max(abs(c(apply(AA, 1, sum))))
    if ((tst > 1e-09) && root) {
      warning("Bhat does not sum to 0 ")
      cat("test =", tst, "\n")
    }
  }
  if (drift)
    x = undrift(x)
  x.cycle = AA %*% as.matrix(x)
  if (type == "fixed" || type == "symmetric" || type == "baxter-king") {
    if (nfix > 0)
      x.cycle[c(1:nfix, (n - nfix + 1):n)] = NA
  }
  x.trend = x - x.cycle
  if (is.ts(xo)) {
    tsp.x = tsp(xo)
    x.cycle = ts(x.cycle, start = tsp.x[1], frequency = tsp.x[3])
    x.trend = ts(x.trend, start = tsp.x[1], frequency = tsp.x[3])
    x = ts(x, start = tsp.x[1], frequency = tsp.x[3])
  }
  if (type == "asymmetric")
    title = "Chiristiano-Fitzgerald Asymmetric Filter"
  if (type == "symmetric")
    title = "Chiristiano-Fitzgerald Symmetric Filter"
  if (type == "fixed")
    title = "Chiristiano-Fitzgerald Fixed Length Filter"
  if (type == "baxter-king")
    title = "Baxter-King Fixed Length Filter"
  if (type == "trigonometric")
    title = "Trigonometric Regression Filter"
  res <- list(cycle = x.cycle, trend = x.trend, fmatrix = AA,
              title = title, xname = xname, call = as.call(match.call()),
              type = type, pl = pl, pu = pu, nfix = nfix, root = root,
              drift = drift, theta = theta, method = "cffilter", x = x)
  return(structure(res, class = "mFilter"))
}

.flipud_adapted <- function (x)
{
  apply(as.matrix(x), 2, rev)
}

.fliplr_adapted <- function (x)
{
  t(apply(as.matrix(x), 1, rev))
}
