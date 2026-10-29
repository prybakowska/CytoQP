#' Debarcodes files
#'
#' @description Performs sample debarcoding.
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param cores Number of cores to be used.
#' @param file_batch_id Character vector containing one batch label for each
#' input FCS file.
#' @param file_score Optional data frame containing file quality information.
#' Must contain \code{file_names} and \code{quality} columns. Only files with
#' \code{quality == "good"} are processed. Default is \code{NULL}.
#' @param out_dir Character, path where output files and plots are saved.
#' @param file_dir Optional directory where debarcoded FCS files are saved.
#' @param min_threshold Logical, whether a minimum debarcoding threshold
#'   should be applied.
#' @param threshold Numeric minimum threshold for debarcoding.
#' @param to_plot Logical, whether plots for yields and debarcoding quality
#'   should be generated.
#' @param barcodes_used Character vector containing barcode names to use, or a
#' named list of barcode names for each batch. List names must match the
#' corresponding values in \code{file_batch_id}. Default is \code{NULL}.
#' @param less_than_th Logical, whether files with cutoffs below threshold
#'   should be recorded.
#' @param barcode_key Barcode key accepted by CATALYST::assignPrelim, such as
#' a data.frame or matrix with sample names as row names, or a numeric
#' vector of barcode masses when no subset is specified.
#' @return Character vector containing paths to debarcoded FCS files.
#' @examples
#' fcs_file <- system.file(
#'   "extdata",
#'   "Center1_batch1.fcs",
#'   package = "CytoQP"
#' )
#'
#' sample_key <- CATALYST::sample_key
#'
#' barcodes_used <- rownames(sample_key)[
#'   c(20, 11, 2)
#' ]
#'
#' files_debarcoded <- debarcode_files(
#'   fcs_files = fcs_file,
#'   file_batch_id = "batch1",
#'   barcodes_used = list(
#'     batch1 = barcodes_used
#'   ),
#'   barcode_key = sample_key,
#'   to_plot = FALSE,
#'   out_dir = file.path(tempdir(), "Debarcoded")
#' )
#'
#' @export
debarcode_files <- function(fcs_files,
                            cores = 1,
                            file_batch_id,
                            file_score = NULL,
                            out_dir = NULL,
                            file_dir = NULL,
                            min_threshold = TRUE,
                            threshold = 0.18,
                            to_plot = TRUE,
                            barcodes_used = NULL,
                            less_than_th = FALSE,
                            barcode_key = NULL) {

  if (is(fcs_files, "list")) {
    fcs_files <- unlist(fcs_files)
  }

  if (anyDuplicated(fcs_files) != 0) {
    stop("names of fcs files are duplicated")
  }

  if (!all(file.exists(fcs_files))) {
    stop("incorrect file path, the fcs file does not exist")
  }

  if (length(file_batch_id) != length(fcs_files)) {
    stop("the length of file_batch_id is not equal to the length of fcs_files")
  }

  if (!is.null(file_score)) {
    if (!inherits(file_score, "data.frame")) {
      stop("file_score is not a data frame")
    } else {
      good_files <- file_score$file_names[file_score$quality == "good"]
      good_idx <- basename(fcs_files) %in% good_files
      file_batch_id <- file_batch_id[good_idx]
      fcs_files <- fcs_files[good_idx]
    }
  }

  if (length(fcs_files) == 0L) {
    stop("No FCS files remain for debarcoding")
  }

  if (is.null(out_dir)) {
    out_dir <- file.path(getwd(), "Debarcoded")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  outFiles <- BiocParallel::bplapply(
    fcs_files,
    function(file) {
      .debarcode_ind(
        file,
        fcs_files = fcs_files,
        file_batch_id = file_batch_id,
        file_score = file_score,
        out_dir = out_dir,
        min_threshold = min_threshold,
        threshold = threshold,
        to_plot = to_plot,
        barcodes_used = barcodes_used,
        less_than_th = less_than_th,
        barcode_key = barcode_key,
        file_dir = file_dir
      )
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = cores)
  )

  lessFiles <- lapply(outFiles, function(x) x[[1]])
  debarcodedFiles <- lapply(outFiles, function(x) x[[2]])

  if (less_than_th) {
    saveRDS(
      unlist(lessFiles),
      file.path(out_dir, "files_with_lower_debarcoding_threshold.RDS")
    )
  }

  return(unlist(debarcodedFiles))
}

.debarcode_ind <- function(file,
                           fcs_files,
                           file_batch_id,
                           file_score,
                           out_dir,
                           min_threshold,
                           threshold,
                           to_plot,
                           barcodes_used,
                           less_than_th,
                           barcode_key,
                           file_dir = NULL) {

  print(paste0("   ", Sys.time()))
  print(paste0("   Debarcoding ", file))

  ff <- flowCore::read.FCS(file, transformation = FALSE)
  file_id <- which(file == fcs_files)
  batch_id <- as.character(file_batch_id[file_id])

  if (!is.null(barcodes_used)) {
    if (is.list(barcodes_used)) {
      s_key <- barcode_key[
        rownames(barcode_key) %in% barcodes_used[[batch_id]],
        , drop = FALSE
      ]
    } else {
      s_key <- barcode_key[
        rownames(barcode_key) %in% barcodes_used,
        , drop = FALSE
      ]
    }
  } else {
    s_key <- barcode_key
  }

  if (is.null(s_key) || NROW(s_key) == 0L) {
    stop("No matching barcodes found for batch: ", batch_id)
  }

  dat <- CATALYST::prepData(ff)
  dat <- CATALYST::assignPrelim(dat, bc_key = s_key)
  dat <- CATALYST::estCutoffs(dat)

  fileOut <- NULL
  sep_cutoffs <- S4Vectors::metadata(dat)$sep_cutoffs

  if (min_threshold) {
    if (any(sep_cutoffs < threshold, na.rm = TRUE)) {
      warning(paste0(
        "cutoff lower than ", threshold, " detected for ",
        basename(file), "; cutoff will be set to ", threshold
      ))
      fileOut <- basename(file)
    }
    sep_cutoffs[which(sep_cutoffs < threshold)] <- threshold
  } else {
    if (any(sep_cutoffs < threshold, na.rm = TRUE)) {
      warning(paste0(
        "cutoff lower than ", threshold, " detected for ", basename(file)
      ))
      fileOut <- basename(file)
    }
  }

  sep_cutoffs[is.na(sep_cutoffs)] <- 1
  S4Vectors::metadata(dat)$sep_cutoffs <- sep_cutoffs

  if (to_plot) {
    p <- CATALYST::plotYields(dat, which = rownames(s_key))
    grDevices::pdf(file.path(
      out_dir,
      sub("(?i)\\.fcs$", "_yields.pdf", basename(file), perl = TRUE)
    ))
    for (name in names(p)) print(p[[name]])
    grDevices::dev.off()
  }

  dat <- CATALYST::applyCutoffs(dat)

  if (to_plot) {
    p <- CATALYST::plotEvents(dat, n = 500)
    grDevices::pdf(file.path(
      out_dir,
      sub("(?i)\\.fcs$", "_debarcode_quality.pdf", basename(file), perl = TRUE)
    ))
    for (name in names(p)) print(p[[name]])
    grDevices::dev.off()
  }

  dat <- dat[, dat$bc_id != 0]
  fs <- CATALYST::sce2fcs(dat, split_by = "bc_id")

  if (is.null(file_dir)) {
    file_dir <- file.path(out_dir, batch_id)
  }
  if (!dir.exists(file_dir)) dir.create(file_dir, recursive = TRUE)

  file_name <- sub("(?i)_cleaned\\.fcs$|\\.fcs$", "", basename(file), perl = TRUE)
  output_names <- paste0(
    rownames(fs@phenoData), "_", file_name, "_debarcoded.fcs"
  )

  flowCore::write.flowSet(
    fs,
    outdir = file_dir,
    filename = output_names
  )

  filePaths <- file.path(file_dir, output_names)
  return(list(fileOut, filePaths))
}
