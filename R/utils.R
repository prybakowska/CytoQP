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


test_match_order <- function(x,y) {

  if (isTRUE(all.equal(x,y))) print('Files are ordered')

  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y))))
    warning('Perfect match but wrong order. Please order the files or define uncommon_prefix parameter')

  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y))))
    warning('No match, please make sure that files are in the same order or define uncommon_prefix parameter')
}
