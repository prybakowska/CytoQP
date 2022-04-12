#' Gate intact cells
#'
#' @description Performs gating of intact cells using flowDensity package and
#' deGate function.
#'
#' @param flow_frame A flowframe that contains cytometry data.
#' @param file_name Character, the file name used for saving the flow frame
#' (if save_gated_flow_frame = TRUE) and for plotting, if NULL (default)
#' the file name stored in keyword FIL will be used.
#' @param tinypeak_removal_head Numeric from 0-1, as in deGate to exclude/include
#' tiny peaks in the head of the density distribution curve for both Iridium
#' channels. Default set to 0.8.
#' @param tinypeak_removal_tail The same as tinypeak_removal1 but for the tail
#' in the density distribution curve.Default set to 0.8.
#' @param alpha_head Numeric, 0-1, as in deGate specify the significance of change
#' in the slope being detected at the head of the density distribution curve.
#' Default 0.05.
#' @param alpha_tail The same as in alpha_head but for the tail of the density
#' distribution curve.Default 0.1.
#' @param arcsine_transform Logical, if the data should be transformed
#' with arcsine transformation and cofactor 5. If FALSE the data won't be
#' transformed, thus transformed flow frame should be used if needed.
#' Default TRUE.
#' @param save_gated_flow_frame Logical, if gated flow frame should be saved.
#' Only cells falling into intact cell region will be saved. Default set to FALSE.
#' @param suffix Character, suffix placed in the name of saved fcs file, only
#' if save_gated_flow_frame = TRUE.Defult is "_intact_gated".
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), Gated).
#' @param ... Additional parameters to pass to flowDensity::deGate()
#'
#' @examples
#' # Set input directory
#' aggregate_dir <- file.path(dir, "Aggregated")
#'
#' # List files for gating
#' files <- list.files(path = aggregate_dir,
#'                     pattern = ".fcs$",
#'                     full.names = TRUE)
#'
#' # Create directory to store plot
#' gate_dir <- file.path(getwd(), "Gated")
#' if(!dir.exists(gate_dir)){dir.create(gate_dir)}
#'
#' # Gate the files and plot the gating strategy for each file
#' n_plots <- 1
#' png(file.path(gate_dir, paste0("gating.png")),
#'     width = n_plots * 300, height = length(files) * 300)
#' layout(matrix(1:(length(files) * n_plots), ncol = n_plots, byrow = TRUE))
#'
#' for (file in files){
#'
#'  ff <- flowCore::read.FCS(filename = file,
#'                           transformation = FALSE)
#'
#'  ff <- gate_intact_cells(flow_frame = ff,
#'                          file_name = basename(file), save_gated_flow_frame = TRUE)
#'}
#'
#'dev.off()
#'
#' @export
#'
#' @return An untransformed flow frame with intact cells only
gate_intact_cells <- function(flow_frame,
                              file_name = NULL,
                              tinypeak_removal_head = 0.8,
                              tinypeak_removal_tail = 0.8,
                              alpha_head = 0.05,
                              alpha_tail = 0.1,
                              arcsine_transform = TRUE,
                              save_gated_flow_frame = FALSE,
                              out_dir = NULL,
                              suffix = "_intact_gated",
                              ...){

  # Check parameters
  if(!is(flow_frame, "flowFrame")) {
    stop("flow_frame must be a flow frame" )
  }

  if (is.null(file_name)){
    file_name <- ff@description$FIL
  } else {
    file_name
  }

  if(arcsine_transform == TRUE){

    ff_t <- flowCore::transform(ff,
                                flowCore::transformList(
                                  flowCore::colnames(ff)[grep("Di", flowCore::colnames(ff))],
                                  CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }

  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("intact")))

  tr <- list()
  for(m in c("Ir193Di", "Ir191Di")){

    tr[[m]] <- c(flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal_head,
                                     upper = FALSE, use.upper = TRUE,
                                     alpha = alpha_head, verbose = F, count.lim = 3, ...),
                 flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal_tail,
                                     upper = TRUE, use.upper = TRUE,
                                     alpha = alpha_tail, verbose = F, count.lim = 3, ...))
  }

  for(m in c("Ir193Di", "Ir191Di")){
    selection[ff_t@exprs[,m] < tr[[m]][1], "intact"] <- FALSE
    selection[ff_t@exprs[,m] > tr[[m]][2], "intact"] <- FALSE
  }

  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c("Ir193Di", "Ir191Di"),
                        main = paste0(basename(file_name)," ( ", format(round(percentage, 2),
                                                                        nsmall = 2), "% )"))

  graphics::abline(h = c(tr[["Ir191Di"]]))
  graphics::abline(v = c(tr[["Ir193Di"]]))
  graphics::points(ff_t@exprs[!selection[,"intact"], c("Ir193Di", "Ir191Di")], pch = ".")

  ff <- ff[selection[,"intact"], ]

  if(save_gated_flow_frame){

    .save_flowframe(ff, out_dir, suffix, file_name)
  }

  return(ff)
}


#' Gate singlet cells
#'
#' @description Detects outliers in the selected channel(s) using MAD
#' (mean absolute deviation).
#'
#' @param flow_frame A flowframe that contains cytometry data.
#' @param file_name Character, the file name used for saving the flow frame
#' (if save_gated_flow_frame = TRUE) and for plotting, if NULL (default)
#' the file name stored in keyword FIL will be used.
#' @param channels character, channels name to be used for gating, default is
#' to Event_length.
#' @param n_mad Numeric, number of MADs to detect outliers.Default set to 2.
#' @param arcsine_transform Logical, if the data should be transformed
#' with arcsine transformation and cofactor 5. If FALSE the data won't be
#' transformed, thus transformed flow frame should be used if needed.
#' Default TRUE.
#' @param save_gated_flow_frame Logical, if gated flow frame should be saved.
#' Only cells falling into intact cell region will be saved. Default set to FALSE.
#' @param suffix Character, suffix placed in the name of saved fcs file, only
#' if save_gated_flow_frame = TRUE. Default is "_singlets_gated".
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), Gated).
#' @param ... Additional arguments to pass to flowDensity::plotDens().
#'
#' @examples
#' #' # Set input directory
#' aggregate_dir <- file.path(dir, "Aggregated")
#'
#' # List files for gating
#' files <- list.files(path = aggregate_dir,
#'                     pattern = ".fcs$",
#'                     full.names = TRUE)
#'
#' # Create directory to store plot
#' gate_dir <- file.path(getwd(), "Gated")
#' if(!dir.exists(gate_dir)){dir.create(gate_dir)}
#'
#' # Gate the files and plot the gating strategy for each file
#' n_plots <- 1
#' png(file.path(gate_dir, paste0("gating.png")),
#'     width = n_plots * 300, height = length(files) * 300)
#' layout(matrix(1:(length(files) * n_plots), ncol = n_plots, byrow = TRUE))
#'
#' for (file in files){
#'
#'  ff <- flowCore::read.FCS(filename = file,
#'                           transformation = FALSE)
#'
#'  ff <- gate_singlet_cells(flow_frame = ff,
#'                           channels = "Event_length",
#'                          file_name = basename(file), save_gated_flow_frame = TRUE)
#'}
#'
#'dev.off()
#'
#' @export
#'
#' @return An untransformed flow frame with singlets only
gate_singlet_cells <- function(flow_frame,
                               file_name = NULL,
                               channels = "Event_length",
                               arcsine_transform = TRUE,
                               save_gated_flow_frame = NULL,
                               suffix = "_singlets_gated",
                               n_mad = 2,
                               out_dir = NULL,
                               ...){

  # Check parameters
  if(!is(flow_frame, "flowFrame")) {
    stop("flow_frame must be a flow frame" )
  }

  if (is.null(file_name)){
    file_name <- flow_frame@description$FIL
  } else {
    file_name
  }

  if(arcsine_transform == TRUE){

    flow_frame_t <- flowCore::transform(flow_frame,
                                        flowCore::transformList(
                                          flowCore::colnames(flow_frame)[grep("Di", flowCore::colnames(flow_frame))],
                                          CytoNorm::cytofTransform))
  } else {
    flow_frame_t <- flow_frame
  }

  selection <- matrix(TRUE,
                      nrow = nrow(flow_frame),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("singlets")))

  selection[, "singlets"] <- .remove_mad_outliers(flow_frame = flow_frame_t,
                                                 channels = channels,
                                                 main = paste("Singlets", file_name),
                                                 n_mad = n_mad,
                                                 xlim = c(0, 100), ylim = c(0, 8), ...)

  flow_frame <- flow_frame[selection[,"singlets"], ]

  if(save_gated_flow_frame){

    .save_flowframe(flow_frame, out_dir, suffix, file_name)
  }

  return(flow_frame)

}


#' Gate live cells
#'
#' @description Performs gating of live cells using flowDensity::deGate
#'
#' @param flow_frame A flowframe that contains cytometry data.
#' @param file_name Character, the file name used for saving the flow frame
#' (if save_gated_flow_frame = TRUE) and for plotting, if NULL (default)
#' the file name stored in keyword FIL will be used,
#' @param viability_channel Character, the channel name used for viability staining
#' @param tinypeak_removal_viability Numeric from 0-1, as in deGate to exclude/include
#' tiny peaks in the tail of the density distribution curve for both viability channel
#' @param tinypeak_removal_Iridium The same as tinypeak_removal_viability but for
#' the head and tail of the density distribution curve in Iridium channel
#' @param alpha_viability Numeric, 0-1, as in deGate specify the significance of change
#' in the slope of viability channel
#' @param alpha_Iridium The same as in alpha_viability but for the Iridium
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5.
#' @param save_gated_flow_frame Logical, if gated flow frame should be saved.
#' Only cells falling into intact cell region will be saved. Default set to FALSE.
#' @param suffix Character, suffix placed in the name of saved fcs file, only
#' if save_gated_flow_frame = TRUE.Defult is "_intact_gated".
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), Gated).
#' @param ... Arguments to pass to flowDensity::plotDens().
#'
#' @examples
#'
#' #' Set input directory
#' aggregate_dir <- file.path(dir, "Aggregated")
#'
#' # List files for gating
#' files <- list.files(path = aggregate_dir,
#'                     pattern = ".fcs$",
#'                     full.names = TRUE)
#'
#' # Create directory to store plot
#' gate_dir <- file.path(getwd(), "Gated")
#' if(!dir.exists(gate_dir)){dir.create(gate_dir)}
#'
#' # Gate the files and plot the gating strategy for each file
#' n_plots <- 1
#' png(file.path(gate_dir, paste0("gating.png")),
#' width = n_plots * 300,
#' height = length(files) * 300)
#' layout(matrix(1:(length(files) * n_plots), ncol = n_plots, byrow = TRUE))
#'
#' for (file in files){
#'
#'  ff <- flowCore::read.FCS(filename = file,
#'                           transformation = FALSE)
#'
#'  ff <- gate_live_cells(flow_frame = ff,
#'                        viability_channel = "Pt195Di", save_gated_flow_frame = TRUE,
#'                        file_name = basename(file))
#'}
#'
#'dev.off()
#'
#' @export
#'
#' @return An untransformed flow frame with live cells only
gate_live_cells <- function(flow_frame,
                            file_name = NULL,
                            viability_channel,
                            tinypeak_removal_viability = 0.8,
                            alpha_viability = 0.1,
                            tinypeak_removal_Iridium = 0.8,
                            alpha_Iridium = 0.05,
                            arcsine_transform = TRUE,
                            save_gated_flow_frame = FALSE,
                            suffix = "_live_gated",
                            out_dir = NULL,... ){

  # Check parameters
  if(!is(flow_frame, "flowFrame")) {
    stop("flow_frame must be a flow frame" )
  }

  ff <- flow_frame

  if (is.null(file_name)){
    file_name <- ff@description$FIL
  } else {
    file_name
  }

  if(arcsine_transform == TRUE){

    ff_t <- flowCore::transform(ff,
                                flowCore::transformList( flowCore::colnames(ff)[grep("Di",  flowCore::colnames(ff))],
                                                         CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }

  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("live")))


  v_ch <- grep(viability_channel,  flowCore::colnames(ff), value = T)

  tr <- list()
  for(m in c("Ir191Di", v_ch)){
    if (m == v_ch) {
      upper = TRUE
      alpha = alpha_viability
      tr[[m]] <- flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal_viability,
                                     upper = upper, use.upper = TRUE,
                                     alpha = alpha, verbose = F, count.lim = 3)

    } else {
      alpha = alpha_Iridium
      tr[[m]] <- c(flowDensity::deGate(ff_t, m,
                                       tinypeak.removal = tinypeak_removal_Iridium,
                                       upper = FALSE, use.upper = TRUE,
                                       alpha = alpha,  verbose = F, count.lim = 3),
                   flowDensity::deGate(ff_t, m,
                                       tinypeak.removal = tinypeak_removal_Iridium,
                                       upper = TRUE, use.upper = TRUE,
                                       alpha = alpha, verbose = F, count.lim = 3))

    }
  }

  for(m in c(v_ch, "Ir191Di")){
    if (m == v_ch) {
      selection[ff_t@exprs[,m] > tr[[m]][1], "live"] <- FALSE
    } else {
      selection[ff_t@exprs[,m] < tr[[m]][1], "live"] <- FALSE
      selection[ff_t@exprs[,m] > tr[[m]][2], "live"] <- FALSE
    }
  }
  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c(v_ch, "Ir191Di"),
                        main = paste0(file_name," ( ", format(round(percentage, 2),
                                                              nsmall = 2), "% )"),
                        xlim = c(0, 8), ylim = c(0, 8), ...)

  graphics::abline(h = tr[["Ir191Di"]])
  graphics::abline(v = tr[[v_ch]])

  graphics::points(ff_t@exprs[!selection[,"live"], c(v_ch, "Ir191Di")], pch = ".")

  ff <- ff[selection[,"live"], ]

  if(save_gated_flow_frame){

    .save_flowframe(flow_frame, out_dir, suffix, file_name)
  }

  return(ff)

}


.save_flowframe <- function(flow_frame,
                            out_dir,
                            suffix,
                            file_name){
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Gated")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  fname <- gsub(pattern = ".fcs", replacement = paste0(suffix, ".fcs"),
                x = file_name, ignore.case = TRUE)

  flowCore::write.FCS(x = flow_frame, filename = file.path(out_dir, fname))
}


#' remove_mad_outliers
#'
#' @description detects outliers in the selected channel(s) using MAD
#' (mean absolute deviation).
#'
#' @param flow_frame A flowframe that contains cytometry data.
#' @param channels character, channel names used for gating, default is set to
#' "Event_length".
#' @param n_mad numeric, how many MAD should be use to detect outliers
#' @param mad_f function used to compute deviation, default set to "mad"
#' @param plot logicle, if to plot the data, default TRUE
#' @param center the centre: defaults to the median
#' @param main character, title of the plot, default set to ""
#' @param ... other arguments to pass plotDens
#'
#' @return matrix with the selected cells
.remove_mad_outliers <- function(flow_frame,
                                 channels = "Event_length",
                                 n_mad = 2,
                                 mad_f = mad,
                                 plot = TRUE,
                                 center = "center",
                                 main = "",
                                 ...){
  boundaries <- matrix(NA,
                       nrow = 5,
                       ncol = length(channels),
                       dimnames = list(c("median", "center", "mad", "l_lim", "u_lim"),
                                       channels))
  for (channel in channels) {
    x <- flow_frame@exprs[, channel]
    boundaries["median", channel] <- stats::median(x)
    boundaries["center", channel] <- stats::density(x)$x[which.max(stats::density(x)$y)]
    boundaries["mad", channel] <- mad_f(x,
                                        center = boundaries[center, channel] )
    boundaries["l_lim", channel] <- boundaries[center, channel] - n_mad * boundaries["mad", channel]
    boundaries["u_lim", channel] <- boundaries[center, channel] + n_mad * boundaries["mad", channel]
  }

  selection <- rep(TRUE, nrow(flow_frame))
  for (channel in channels) {
    selection <- selection & (flow_frame@exprs[, channel] > boundaries["l_lim", channel])
    selection <- selection & (flow_frame@exprs[, channel] < boundaries["u_lim", channel])
  }
  percentage <- (sum(selection)/length(selection))*100
  if (plot) {
    flowDensity::plotDens(flow_frame,
                          c(channels, "Ir191Di"),
                          main = paste0(main, " ( ", format(round(percentage, 2),
                                                            nsmall = 2), "% )"),
                          ...)
    if(length(channels) == 2) {
      graphics::points(flow_frame@exprs[!selection, channels], col = "red", pch = ".")
      graphics::abline(v = boundaries[c("l_lim", "u_lim"), channels[1]], col = "black")
      graphics::abline(h = boundaries[c("l_lim", "u_lim"), channels[2]], col = "black")
    } else if(length(channels) == 1) {
      graphics::points(flow_frame@exprs[!selection, c(channels, "Ir191Di")], pch = ".")
      graphics::abline(v = boundaries[c("l_lim", "u_lim"), channels[1]], col = "black")
    }
  }

  return(selection)
}
