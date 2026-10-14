#' Internal package imports
#'
#' @importFrom methods is
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom graphics hist
#' @importFrom stats IQR convolve frequency is.ts mad median na.omit qt sd ts tsp
#' @importFrom utils capture.output tail write.csv
#' @noRd
NULL
#' Finds all the mass/fluorochrome channels for the flow frame
#'
#' @description Finds all the mass channels
#'
#' @param flow_frame Untransformed flow frame
#' @param channels Pattern for non-mass channels, default is
#' "Time|Event_length|Center|Offset|Width|Residual|SSC|FSC|File_scattered"
#' @param ... Additional arguments to pass to grep
#'
#' @return Logical vector with TRUE values for mass channels
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
#' mass_channels <- find_mass_ch(
#'   flow_frame = ff
#' )
#'
#' @export
find_mass_ch <- function(flow_frame,
                         channels = "Time|Event_length|Center|Offset|Width|Residual|SSC|FSC|File_scattered",
                         ...){
  non_mass_ch <- grep(c(channels),
                      flowCore::colnames(flow_frame),
                      invert = TRUE, ...)
  return(non_mass_ch)
}


#' Tests the order of the files
#'
#' @param x the name of the file 1 (after removing uncommon prefix/sufix)
#' @param y the name of the file 2 (after removing uncommon prefix/sufix)
#'
#' @return warrning if he files are not in order
#' @export
#'
#' @examples
#' test_match_order(
#'   x = c("sample1.fcs", "sample2.fcs"),
#'   y = c("sample1.fcs", "sample2.fcs")
#' )
test_match_order <- function(x,y) {

  if (isTRUE(all.equal(x,y))) print('Files are ordered')

  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y))))
    warning('Perfect match but wrong order. Please order the files or define uncommon_prefix parameter')

  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y))))
    warning('No match, please make sure that files are in the same order or define uncommon_prefix parameter')
}
