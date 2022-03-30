#' Split big flow frames into smaller fcs files
#'
#' @description Split the big flow frames into smaller ones.
#'
#' @param flow_frame Flow frame containing cytometry data.
#' @param event_per_flow_frame Numeric, the number of events to be split to
#' small flow frames, default is set to 500000.
#' @param min_cell_per_fcs Numeric, minimal number of cells in flow frame
#' to save fcs file, default 20000.
#' @param out_dir Character, pathway to where the files should be saved,
#' default is set to file.path(getwd(), Splitted).
#'
#' @return Save splitted fcs files in out_dir.
#'
#' @examples
#' ff <- flowCore::read.FCS(list.files(path = raw_data_dir,
#' pattern = ".FCS", full.names = TRUE)[1])
#'
#' split_big_flowFrames(flow_frame = ff, event_per_flow_frame = 100000)
#'
#' @export
split_big_flowFrames <- function(flow_frame,
                                 event_per_flow_frame = 500000,
                                 out_dir = NULL,
                                 min_cell_per_fcs = 20000){

  total_events <- nrow(flow_frame)
  res_breaks <- .make_breaks(event_per_flow_frame, total_events)

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Splitted")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  for (i in as.numeric((names(res_breaks$breaks)))) {
    id <- res_breaks$breaks[[i]]

    if (i < 10) {
      num <- paste0("0", i)
    }

    if(nrow(flow_frame[id, ]) > min_cell_per_fcs){
      flowCore::write.FCS(flow_frame[id, ],
                          file.path(out_dir,
                                    gsub(".fcs|.FCS",
                                         paste0("_", num, ".fcs"), flow_frame@description$ORIGINALGUID)))
    }

  }
}


#'Calculates breaks for flow frame splitting
.make_breaks <- function(event_per_flow_frame, total_events){
  breaks <- .split_flowFrames(seq_len(total_events),
                              event_per_flow_frame)

  names(breaks) <- seq_along(breaks)

  return(list("breaks"=breaks, "events_per_flowframe"=event_per_flow_frame))
}


#'Calculates beginning and end of each flow frame.
#'Code from Annelies Emmaneel (2021). PeacoQC:
#'Peak-based selection of high quality cytometry data. R package
#'version 1.4.0.
.split_flowFrames <- function(vec, seg.length) {
  starts=seq(1, length(vec), by=seg.length)
  ends  = starts + seg.length - 1
  ends[ends > length(vec)]=length(vec)

  lapply(seq_along(starts), function(i) vec[starts[i]:ends[i]])
}
