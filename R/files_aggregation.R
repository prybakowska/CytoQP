#' Deconvolute and aggregate debarcoded files
#'
#' @description Performs aggregation of debarcoded files, assigning user defined
#' name to each file.
#'
#' @param fcs_files Character, full path to the fcs_files.
#' @param md Metadata. Must contain the following columns:
#' batch_column: defines to which batch each file belongs;
#' barcode_name: defines to which barcode each file belongs
#' fcs_new_name: a name for the fcs file that will be given after deconvolution
#' and aggregation.
#' @param barcode_column Character string specifying the column in \code{md}
#' containing the barcode identifiers.
#' @param batch_column Character string specifying the column in \code{md}
#' containing the batch identifiers.
#' @param cores Number of cores to be used.
#' @param channels_to_keep Character vector with channel names to be kept.
#' Default NULL.
#' @param maxcells Numeric, maximum cells to randomly aggregate from each file,
#' default is set to NULL, which means that all the cells will be aggregated.
#' @param write_agg_file Logical, if the fcs files should be saved, if TRUE
#' files will be saved in out_dir. Default set to TRUE.
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), Aggregated).
#'
#' @return List of the pathways to aggregated files.
#'
#' @examples
#' fcs_files <- system.file(
#'   "extdata",
#'   paste0("Center", 1:4, "_batch1.fcs"),
#'   package = "CytoQP"
#' )
#'
#' md <- data.frame(
#'   barcode_name = "Center1",
#'   batch = "batch1",
#'   fcs_new_name = "sample_aggregated.fcs"
#' )
#'
#' files_aggregated <- aggregate_files(
#'   fcs_files = fcs_files,
#'   md = md,
#'   barcode_column = "barcode_name",
#'   batch_column = "batch",
#'   cores = 1,
#'   write_agg_file = FALSE,
#'   out_dir = file.path(tempdir(), "Aggregated")
#' )
#'
#' @export
aggregate_files <- function(fcs_files,
                            md,
                            barcode_column,
                            batch_column,
                            cores = 1,
                            channels_to_keep = NULL,
                            maxcells = NULL,
                            write_agg_file = TRUE,
                            out_dir = NULL){

  # Check parameters
  if(!is(fcs_files, "character") & !is(fcs_files, "list")) {
    stop("files must be a character vector or a list")
  }

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if (any(!is(cores, "numeric") | cores < 1)){
    stop("cores must be a positive number")
  }

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Aggregated")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  # Parallelized analysis
  aggregatedFiles <- BiocParallel::bplapply(seq_len(nrow(md)), function(i) {
    patterns <- as.character(md[i, c(barcode_column, batch_column)])

    files_to_agg <- grep(pattern = patterns[2],
                         grep(pattern = patterns[1],
                              fcs_files, value = TRUE),
                         value = TRUE)

    print(paste0("Creating ", md[[i, "fcs_new_name"]]))

    outputFile <- md[[i, "fcs_new_name"]]

    .aggregate_ind(fcs_files = files_to_agg,
                   channels_to_keep = channels_to_keep,
                   outputFile = outputFile,
                   maxcells = maxcells,
                   write_agg_file = write_agg_file,
                   out_dir = out_dir)
  },
  BPPARAM = BiocParallel::MulticoreParam(workers = cores))

  return(unlist(aggregatedFiles))

}


.aggregate_ind <- function(fcs_files,
                           channels_to_keep = NULL,
                           outputFile = "aggregate.fcs",
                           maxcells = NULL,
                           write_agg_file = FALSE,
                           out_dir = getwd()) {
  nFiles <- length(fcs_files)
  flowFrame <- NULL

  if(!dir.exists(out_dir))(dir.create(out_dir))

  for (i in seq_len(nFiles)) {
    f <- flowCore::read.FCS(fcs_files[i])

    if(!is.null(maxcells)){
      cell_idx <- sample(
        seq_len(nrow(f)),
        min(nrow(f), maxcells)
      )

      f <- f[cell_idx, ]
    }

    m <- matrix(rep(i, nrow(f)))
    m2 <- m + stats::rnorm(length(m), 0, 0.1)
    m <- cbind(m, m2)
    colnames(m) <- c("File", "File_scattered")
    prev_agg <- length(grep("File[0-9]*$", colnames(f)))
    if (prev_agg > 0) {
      colnames(m) <- paste0(colnames(m), prev_agg + 1)
    }
    if(is.null(channels_to_keep)){
      f <- flowCore::fr_append_cols(f, m)
    } else {
      f <- flowCore::fr_append_cols(f[ , channels_to_keep], m)
    }
    if (is.null(flowFrame)) {
      flowFrame <- f
      flowFrame@description$`$FIL` <- gsub(".*/", "",
                                           outputFile)
      flowFrame@description$FILENAME <- gsub(".*/", "",
                                             outputFile)
    }
    else {
      f@exprs[, "Time"] <- f@exprs[, "Time"] + max(flowFrame@exprs[,"Time"]) + 1000
      flowCore::exprs(flowFrame) <- rbind(flowCore::exprs(flowFrame),
                                          flowCore::exprs(f))
    }
  }
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame) -
                                 1, "Rmin", sep = "")]] <- 0
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame) -
                                 1, "Rmax", sep = "")]] <- nFiles + 1
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame),
                               "Rmin", sep = "")]] <- 0
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame),
                               "Rmax", sep = "")]] <- nFiles + 1
  flowFrame@description[[paste("$P", ncol(flowFrame) - 1,
                               "B", sep = "")]] <- 32
  flowFrame@description[[paste("$P", ncol(flowFrame), "B",
                               sep = "")]] <- 32
  flowFrame@description$FIL <- gsub(".*/", "", outputFile)

  if(write_agg_file == TRUE){

    flowCore::write.FCS(x = flowFrame, filename = file.path(out_dir, outputFile), endian = "big")
  }

  return(file.path(out_dir, outputFile))

  # return(flowFrame)
}
