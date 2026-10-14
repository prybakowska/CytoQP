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
#' @param barcode_key Matrix as in CATALYST::assignPrelim.
#' @param suffix_old Character string to replace in input filenames.
#' Default is \code{"cleaned"}.
#' @param suffix_new Character string used in debarcoded output filenames.
#' Default is \code{"debarcoded"}.
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
                            barcode_key = NULL,
                            suffix_old = "cleaned",
                            suffix_new = "debarcoded") {

  if (is(fcs_files, "list")) {
    fcs_files <- unlist(fcs_files)
  }

  if (anyDuplicated(fcs_files) != 0) {
    stop(
      "names of fcs files are duplicated",
      call. = FALSE
    )
  }

  if (!all(file.exists(fcs_files))) {
    stop(
      "incorrect file path, the fcs file does not exist",
      call. = FALSE
    )
  }

  if (length(file_batch_id) != length(fcs_files)) {
    stop(
      "the length of file_batch_id is not equal to the length of fcs_files",
      call. = FALSE
    )
  }

  if (is.null(barcode_key)) {
    stop(
      "barcode_key must be provided",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Optional file quality filtering
  # ---------------------------------------------------------------------------

  if (!is.null(file_score)) {

    if (!inherits(file_score, "data.frame")) {
      stop(
        "file_score is not a data frame",
        call. = FALSE
      )
    }

    good_files <- file_score$file_names[
      file_score$quality == "good"
    ]

    good_idx <- basename(fcs_files) %in% good_files

    fcs_files <- fcs_files[
      good_idx
    ]

    file_batch_id <- file_batch_id[
      good_idx
    ]
  }


  # ---------------------------------------------------------------------------
  # Output directory
  # ---------------------------------------------------------------------------

  if (is.null(out_dir)) {
    out_dir <- file.path(
      getwd(),
      "Debarcoded"
    )
  }

  if (!dir.exists(out_dir)) {
    dir.create(
      out_dir,
      recursive = TRUE
    )
  }


  # ---------------------------------------------------------------------------
  # Parallel backend
  # ---------------------------------------------------------------------------

  if (.Platform$OS.type == "windows") {

    BPPARAM <- BiocParallel::SnowParam(
      workers = cores
    )

  } else {

    BPPARAM <- BiocParallel::MulticoreParam(
      workers = cores
    )
  }


  # ---------------------------------------------------------------------------
  # Debarcode files
  # ---------------------------------------------------------------------------

  outFiles <- BiocParallel::bplapply(
    fcs_files,
    function(file) {

      .debarcode_ind(
        file = file,
        fcs_files = fcs_files,
        file_batch_id = file_batch_id,
        out_dir = out_dir,
        file_dir = file_dir,
        min_threshold = min_threshold,
        threshold = threshold,
        to_plot = to_plot,
        barcodes_used = barcodes_used,
        less_than_th = less_than_th,
        barcode_key = barcode_key,
        suffix_old = suffix_old,
        suffix_new = suffix_new
      )
    },
    BPPARAM = BPPARAM
  )


  # ---------------------------------------------------------------------------
  # Collect output
  # ---------------------------------------------------------------------------

  lessFiles <- lapply(
    outFiles,
    function(x) {
      x[[1]]
    }
  )

  debarcodedFiles <- lapply(
    outFiles,
    function(x) {
      x[[2]]
    }
  )

  if (less_than_th) {
    saveRDS(
      unlist(lessFiles),
      file.path(
        out_dir,
        "files_with_lower_debarcoding_threshold.RDS"
      )
    )
  }

  return(
    unlist(debarcodedFiles)
  )
}


# ==============================================================================
# Internal helper
# ==============================================================================

.debarcode_ind <- function(file,
                           fcs_files,
                           file_batch_id,
                           out_dir,
                           file_dir = NULL,
                           min_threshold,
                           threshold,
                           to_plot,
                           barcodes_used,
                           less_than_th,
                           barcode_key,
                           suffix_old,
                           suffix_new
                           ) {

  print(
    paste0(
      "   ",
      Sys.time()
    )
  )

  print(
    paste0(
      "   Debarcoding ",
      file
    )
  )


  # ---------------------------------------------------------------------------
  # Read FCS
  # ---------------------------------------------------------------------------

  ff <- flowCore::read.FCS(
    file,
    transformation = FALSE
  )

  file_id <- which(
    file == fcs_files
  )

  batch_id <- file_batch_id[
    file_id
  ]

  if (length(batch_id) != 1L) {
    stop(
      paste0(
        "Could not uniquely determine batch for file: ",
        basename(file)
      ),
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Select barcode key
  # ---------------------------------------------------------------------------

  if (!is.null(barcodes_used)) {

    if (is.list(barcodes_used)) {

      if (is.null(barcodes_used[[batch_id]])) {
        stop(
          paste0(
            "No barcode information found for batch: ",
            batch_id
          ),
          call. = FALSE
        )
      }

      s_key <- barcode_key[
        rownames(barcode_key) %in%
          barcodes_used[[batch_id]],
        ,
        drop = FALSE
      ]

    } else {

      s_key <- barcode_key[
        rownames(barcode_key) %in%
          barcodes_used,
        ,
        drop = FALSE
      ]
    }

  } else {

    s_key <- barcode_key
  }

  if (nrow(s_key) == 0L) {
    stop(
      paste0(
        "No matching barcodes were found for batch: ",
        batch_id
      ),
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Preliminary barcode assignment
  # ---------------------------------------------------------------------------

  dat <- CATALYST::prepData(
    ff
  )

  dat <- CATALYST::assignPrelim(
    dat,
    bc_key = s_key
  )

  dat <- CATALYST::estCutoffs(
    dat
  )


  # ---------------------------------------------------------------------------
  # Barcode separation thresholds
  # ---------------------------------------------------------------------------

  fileOut <- NULL

  sep_cutoffs <- S4Vectors::metadata(dat)$sep_cutoffs

  if (min_threshold) {

    if (
      any(
        sep_cutoffs < threshold,
        na.rm = TRUE
      )
    ) {

      warning(
        paste0(
          "cutoff lower than ",
          threshold,
          " has been detected for ",
          basename(file),
          ", cutoff will be set to ",
          threshold
        )
      )

      fileOut <- basename(file)
    }

    id <- sep_cutoffs < threshold

    id[
      is.na(id)
    ] <- FALSE

    sep_cutoffs[
      id
    ] <- threshold

  } else {

    if (
      any(
        sep_cutoffs < threshold,
        na.rm = TRUE
      )
    ) {

      warning(
        paste0(
          "cutoff lower than ",
          threshold,
          " detected for ",
          basename(file)
        )
      )

      fileOut <- basename(file)
    }
  }


  # ---------------------------------------------------------------------------
  # Replace missing cutoffs
  # ---------------------------------------------------------------------------

  sep_cutoffs[
    is.na(sep_cutoffs)
  ] <- 1

  S4Vectors::metadata(dat)$sep_cutoffs <-
    sep_cutoffs


  # ---------------------------------------------------------------------------
  # Plot barcode yields
  # ---------------------------------------------------------------------------

  if (to_plot) {

    p <- CATALYST::plotYields(
      dat,
      which = rownames(s_key)
    )

    grDevices::pdf(
      file.path(
        out_dir,
        gsub(
          "(?i)\\.fcs$",
          "_yields.pdf",
          basename(file),
          perl = TRUE
        )
      )
    )

    for (name in names(p)) {
      print(
        p[[name]]
      )
    }

    grDevices::dev.off()
  }


  # ---------------------------------------------------------------------------
  # Apply barcode cutoffs
  # ---------------------------------------------------------------------------

  dat <- CATALYST::applyCutoffs(
    dat
  )


  # ---------------------------------------------------------------------------
  # Plot debarcoding quality
  # ---------------------------------------------------------------------------

  if (to_plot) {

    p <- CATALYST::plotEvents(
      dat,
      n = 500
    )

    grDevices::pdf(
      file.path(
        out_dir,
        gsub(
          "(?i)\\.fcs$",
          "_debarcode_quality.pdf",
          basename(file),
          perl = TRUE
        )
      )
    )

    for (name in names(p)) {
      print(
        p[[name]]
      )
    }

    grDevices::dev.off()
  }


  # ---------------------------------------------------------------------------
  # Remove unassigned cells
  # ---------------------------------------------------------------------------

  dat <- dat[
    ,
    dat$bc_id != 0
  ]


  # ---------------------------------------------------------------------------
  # Convert barcode populations to FCS
  # ---------------------------------------------------------------------------

  fs <- CATALYST::sce2fcs(
    dat,
    split_by = "bc_id"
  )


  # ---------------------------------------------------------------------------
  # Output directory
  # ---------------------------------------------------------------------------

  if (is.null(file_dir)) {

    file_dir_current <- file.path(
      out_dir,
      batch_id
    )

  } else {

    file_dir_current <- file_dir
  }

  if (!dir.exists(file_dir_current)) {
    dir.create(
      file_dir_current,
      recursive = TRUE
    )
  }


  # ---------------------------------------------------------------------------
  # Output file names
  # ---------------------------------------------------------------------------

  file_name <- sub(
    "(?i)\\.fcs$",
    "",
    basename(file),
    perl = TRUE
  )

  if (grepl(suffix_old, file_name, fixed = TRUE)) {

    file_name <- sub(
      suffix_old,
      suffix_new,
      file_name,
      fixed = TRUE
    )

  } else {

    file_name <- paste0(
      file_name,
      "_",
      suffix_new
    )
  }

  output_names <- paste0(
    rownames(fs@phenoData),
    "_",
    file_name,
    ".fcs"
  )


  # ---------------------------------------------------------------------------
  # Write FCS files
  # ---------------------------------------------------------------------------

  flowCore::write.flowSet(
    fs,
    outdir = file_dir_current,
    filename = output_names
  )

  filePaths <- file.path(
    file_dir_current,
    output_names
  )


  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------

  return(
    list(
      fileOut,
      filePaths
    )
  )
}
