# CytoQP: V1-compatible spillover visualization and CyTOF compensation.
# Keep the original public functions and their argument names/defaults:
#   plotspillmat(), compensation().
# No additional public compensation interfaces are introduced.

#' Plot spillover values (internal)
#'
#' @param data_df Numeric spillover matrix expressed as percentages.
#' @param custom_colors Four colors for the fill gradient.
#' @param title Plot title.
#' @param ... Retained for compatibility with original plotting calls.
#' @return A ggplot object.
#' @keywords internal
.plot_spill <- function(data_df,
                        custom_colors = c("white", "lightcoral", "red2", "darkred"),
                        title = "Spillover Matrix", ...) {
  data_df <- as.matrix(data_df)
  if (!is.numeric(data_df) || length(dim(data_df)) != 2L ||
      any(dim(data_df) == 0L)) {
    stop("data_df must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (length(custom_colors) != 4L) {
    stop("custom_colors must contain exactly four colors.", call. = FALSE)
  }

  # V1 used row/column sums <= 100 to drop channels. That accidentally
  # removed valid channels with a 100% diagonal and no measured spillover.
  # Drop only entirely zero/NA channels, never valid diagonal-only channels.
  keep_rows <- rowSums(abs(data_df), na.rm = TRUE) > 0
  keep_cols <- colSums(abs(data_df), na.rm = TRUE) > 0
  if (!any(keep_rows) || !any(keep_cols)) {
    stop("Spillover matrix contains no non-zero channels.", call. = FALSE)
  }
  data_df <- data_df[keep_rows, keep_cols, drop = FALSE]

  off_diag_values <- data_df[is.finite(data_df) & data_df != 100]
  if (length(off_diag_values) > 0L) {
    scale_data <- range(off_diag_values)
  } else {
    scale_data <- c(0, 1)
  }
  if (diff(scale_data) == 0) {
    scale_data <- scale_data + c(-0.01, 0.01)
  }

  plot_df <- as.data.frame(as.table(data_df), stringsAsFactors = FALSE)
  ggplot2::ggplot(plot_df, ggplot2::aes(x = Var2, y = Var1, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(
      colors = c(custom_colors[1:3], "grey", custom_colors[4]),
      limits = scale_data, guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 22)
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::geom_hline(
      yintercept = seq(0.5, nrow(data_df) + 0.5),
      color = "black", linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = seq(0.5, ncol(data_df) + 0.5),
      color = "black", linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(Freq != 0 & Freq != 100,
                                  as.character(Freq), "")),
      color = "black", size = 2
    ) +
    ggplot2::coord_fixed(ratio = 1)
}


#' Extend isotope mappings (internal)
#'
#' @param ff A flowFrame, spillover matrix, or data.frame.
#' @return A named numeric isotope list compatible with CATALYST.
#' @keywords internal
.fixed_isotypes <- function(ff) {
  isotope_list <- CATALYST::isotope_list
  if (inherits(ff, "flowFrame")) {
    marker_names <- flowCore::colnames(ff)
  } else if (is.matrix(ff) || is.data.frame(ff)) {
    marker_names <- colnames(ff)
  } else {
    stop("ff must be a flowFrame, matrix, or data.frame.", call. = FALSE)
  }
  if (is.null(marker_names) || length(marker_names) == 0L) {
    return(isotope_list)
  }
  # Only bona fide mass-channel labels; V1's broad regex could include CD3,
  # Time, and other non-isotope metadata as isotope symbols.
  marker_names <- marker_names[grepl("^[A-Za-z]+[0-9]+Di$", marker_names)]
  if (!length(marker_names)) return(isotope_list)
  symbols <- sub("^([A-Za-z]+)[0-9]+Di$", "\\1", marker_names)
  masses <- as.numeric(sub("^[A-Za-z]+([0-9]+)Di$", "\\1", marker_names))
  for (symbol in setdiff(unique(symbols), names(isotope_list))) {
    possible_masses <- unique(masses[symbols == symbol])
    if (length(possible_masses) != 1L) {
      warning("Cannot infer a unique isotope mass for ", symbol,
              "; supply isotope_list explicitly.", call. = FALSE)
      next
    }
    isotope_list[symbol] <- possible_masses
  }
  isotope_list
}


#' Estimate the spillover matrix (internal)
#'
#' @param input A flowFrame containing single-stained CyTOF beads.
#' @param choosen_channels Optional character vector of channels (V1 spelling).
#' @param isotope_list Retained V1 parameter; custom mappings are used by
#'   compensation(), not the CATALYST spillover-estimation steps.
#' @param ... Additional V1 arguments.
#' @return A spillover matrix of proportions.
#' @keywords internal
.catalyst_compensation <- function(input, choosen_channels = NULL,
                                   isotope_list = NULL, ...) {
  if (!inherits(input, "flowFrame")) {
    stop("input must be a flowFrame.", call. = FALSE)
  }
  available_channels <- flowCore::colnames(input)
  if (is.null(choosen_channels)) {
    choosen_channels <- available_channels[grepl("Di", available_channels)]
  }
  if (!is.character(choosen_channels) || length(choosen_channels) == 0L ||
      anyNA(choosen_channels)) {
    stop("choosen_channels must contain at least one channel.", call. = FALSE)
  }
  missing_channels <- setdiff(choosen_channels, available_channels)
  if (length(missing_channels)) {
    stop("Channels not found in compensation beads: ",
         paste(missing_channels, collapse = ", "), call. = FALSE)
  }
  bc_key <- as.numeric(gsub("\\D", "", choosen_channels))
  if (anyNA(bc_key)) {
    stop("Cannot extract mass numbers from choosen_channels.", call. = FALSE)
  }
  bc_key <- sort(unique(bc_key))

  # Preserve the V1 CATALYST algorithm; no changes to the estimation method.
  sce <- CATALYST::prepData(input)
  sce <- CATALYST::assignPrelim(sce, bc_key)
  sce <- CATALYST::applyCutoffs(CATALYST::estCutoffs(sce))
  sce <- CATALYST::computeSpillmat(sce)
  spillover_matrix <- S4Vectors::metadata(sce)$spillover_matrix
  if (is.null(spillover_matrix)) {
    stop("CATALYST did not return a spillover matrix.", call. = FALSE)
  }
  spillover_matrix
}


#' Plot Spillover Matrix
#'
#' @description
#' Creates a spillover matrix from single-stained bead FCS data and plots it,
#' or plots an existing spillover matrix. This is the original CytoQP V1 entry
#' point, with all original parameters retained.
#'
#' @param input A numeric spillover matrix, data.frame, or compensation-bead
#'   FCS file path. Matrices must contain proportions (diagonal near 1).
#' @param choosen_channels Optional character vector of mass channels.
#' @param transformation V1 parameter: when FALSE, inverse CytoNorm
#'   transformation is applied to compensation-bead data before estimation.
#'   Default TRUE.
#' @param isotope_list Optional isotope list retained for V1 compatibility.
#' @param savespillmat Optional CSV path to save the estimated spillover matrix
#'   (as proportions) when the input is an FCS file.
#' @param ... Optional plot parameters such as custom_colors and title.
#' @return A ggplot object.
#' @examples
#' m <- matrix(c(1, 0.01, 0.02, 0.01, 1, 0.03, 0.02, 0.03, 1),
#'             nrow = 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' plotspillmat(m)
#' @export
plotspillmat <- function(input, choosen_channels = NULL,
                         transformation = TRUE, isotope_list = NULL,
                         savespillmat = NULL, ...) {
  if (is.matrix(input) || is.data.frame(input)) {
    sm <- as.matrix(input)
    if (!is.numeric(sm)) {
      stop("The spillover matrix must be numeric.", call. = FALSE)
    }
  } else if (is.character(input) && length(input) == 1L &&
             file.exists(input) && flowCore::isFCSfile(input)) {
    ff <- flowCore::read.FCS(input, transformation = FALSE,
                             truncate_max_range = FALSE)
    if (isFALSE(transformation)) {
      channels <- flowCore::colnames(ff)
      channels <- channels[grepl("Di", channels)]
      if (length(channels)) {
        ff <- flowCore::transform(
          ff, flowCore::transformList(channels,
                                      CytoNorm::cytofTransform.reverse)
        )
      }
    }
    # Do not pass plotting arguments (...) to the spillover estimator.
    sm <- .catalyst_compensation(
      input = ff, choosen_channels = choosen_channels,
      isotope_list = isotope_list
    )
    if (!is.null(savespillmat)) {
      utils::write.csv(sm, file = savespillmat, row.names = TRUE)
    }
  } else {
    stop("input must be a numeric matrix, data.frame, or valid FCS file.",
         call. = FALSE)
  }
  # CATALYST spillover values are proportions; V1 displays percentages.
  .plot_spill(round(sm * 100, digits = 1), ...)
}


#' Data Compensation
#'
#' @description
#' Compensates CyTOF data using CATALYST. This is the original CytoQP V1
#' function and argument interface. V2 fixes validation, file output handling,
#' and isotope detection without renaming public arguments or changing defaults.
#'
#' @param samples_input A flowFrame or character vector/list of FCS paths.
#' @param compensation_input A spillover matrix/data.frame of proportions, or
#'   compensation-bead FCS file path.
#' @param choosen_channels Optional channels used for bead-based estimation.
#' @param transformation Original V1 behavior: when FALSE for flowFrame input,
#'   applies inverse CytoNorm transformation before compensation. Default TRUE.
#' @param isotope_list Optional named isotope-to-mass mapping.
#' @param out_dir Directory for saved FCS files. NULL creates Compensated_files.
#' @param suffix_old Substring replaced in source filenames. Default "cleaned".
#' @param suffix_new Replacement substring. Default "compensated".
#' @param ... Original V1 argument retained; unused additional arguments
#'   produce a warning instead of being silently discarded.
#' @return For flowFrame input returns a compensated flowFrame. For FCS paths,
#'   writes compensated files and returns invisibly NULL, as in V1.
#' @examples
#' \dontrun{
#' compensation(
#'   samples_input = c("sample_cleaned.fcs"),
#'   compensation_input = "Compensation_beads.fcs",
#'   out_dir = "Compensated_files"
#' )
#' }
#' @export
compensation <- function(samples_input, compensation_input,
                         choosen_channels = NULL, transformation = TRUE,
                         isotope_list = NULL, out_dir = NULL,
                         suffix_old = "cleaned", suffix_new = "compensated", ...) {
  if (is.list(samples_input) && !inherits(samples_input, "flowFrame")) {
    samples_input <- unlist(samples_input, use.names = FALSE)
  }
  sample_is_frame <- inherits(samples_input, "flowFrame")
  if (!sample_is_frame &&
      (!is.character(samples_input) || length(samples_input) == 0L ||
       anyNA(samples_input))) {
    stop("samples_input must be a flowFrame or FCS file paths.", call. = FALSE)
  }
  if (!is.logical(transformation) || length(transformation) != 1L ||
      is.na(transformation)) {
    stop("transformation must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(list(...))) {
    warning("Additional arguments in '...' are not used by the CATALYST ",
            "compensation workflow.", call. = FALSE)
  }

  # Estimate or accept the existing spillover matrix.
  if (is.matrix(compensation_input) || is.data.frame(compensation_input)) {
    sm <- as.matrix(compensation_input)
  } else if (is.character(compensation_input) &&
             length(compensation_input) == 1L &&
             file.exists(compensation_input) &&
             flowCore::isFCSfile(compensation_input)) {
    bead_file <- flowCore::read.FCS(
      compensation_input, transformation = FALSE, truncate_max_range = FALSE
    )
    sm <- .catalyst_compensation(
      input = bead_file, choosen_channels = choosen_channels,
      isotope_list = isotope_list
    )
  } else {
    stop("compensation_input must be a spillover matrix, data.frame, ",
         "or valid compensation-bead FCS path.", call. = FALSE)
  }
  if (!is.numeric(sm) || length(dim(sm)) != 2L || any(dim(sm) == 0L) ||
      any(!is.finite(sm))) {
    stop("The spillover matrix must contain finite numeric values.",
         call. = FALSE)
  }

  if (is.null(isotope_list)) {
    if (sample_is_frame) {
      isotope_list <- .fixed_isotypes(samples_input)
    } else {
      missing_files <- samples_input[!file.exists(samples_input)]
      if (length(missing_files)) {
        stop("FCS files not found: ", paste(missing_files, collapse = ", "),
             call. = FALSE)
      }
      first_sample <- flowCore::read.FCS(
        samples_input[1], transformation = FALSE, truncate_max_range = FALSE
      )
      isotope_list <- .fixed_isotypes(first_sample)
    }
  }
  # Explicit user-provided isotope_list must not be silently overwritten.
  if (is.null(names(isotope_list)) || anyNA(names(isotope_list))) {
    stop("isotope_list must be named by isotope symbol.", call. = FALSE)
  }
  if (!"BCKG" %in% names(isotope_list)) isotope_list[["BCKG"]] <- 190

  if (!sample_is_frame) {
    # V1 did not initialize out_dir for matrix input, causing path/NULL bugs.
    if (is.null(out_dir)) out_dir <- "Compensated_files"
    if (!is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir)) {
      stop("out_dir must be a non-empty directory path.", call. = FALSE)
    }
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(out_dir)) stop("Could not create out_dir.", call. = FALSE)
    missing_files <- samples_input[!file.exists(samples_input)]
    if (length(missing_files)) {
      stop("FCS files not found: ", paste(missing_files, collapse = ", "),
           call. = FALSE)
    }
    # Preserve V1 name replacement; do not silently add V2-only suffixes.
    output_names <- gsub(suffix_old, suffix_new, basename(samples_input),
                         fixed = TRUE)
    output_paths <- file.path(out_dir, output_names)
    if (anyDuplicated(output_paths)) {
      stop("Output names are not unique; check sample basenames.",
           call. = FALSE)
    }
    input_paths <- normalizePath(samples_input, mustWork = TRUE)
    for (i in seq_along(output_paths)) {
      output_path <- normalizePath(output_paths[i], mustWork = FALSE)
      if (output_path %in% input_paths) {
        stop("Refusing to overwrite an input FCS file: ", output_path,
             call. = FALSE)
      }
    }
  }

  if (sample_is_frame) {
    sample_files <- list(samples_input)
  } else {
    sample_files <- as.list(samples_input)
  }
  for (i in seq_along(sample_files)) {
    sample_ff <- if (sample_is_frame) {
      sample_files[[i]]
    } else {
      flowCore::read.FCS(sample_files[[i]], transformation = FALSE,
                         truncate_max_range = FALSE)
    }
    original_desc <- sample_ff@parameters$desc
    if (!is.null(original_desc)) {
      sample_ff@parameters$desc <- stringr::str_split_fixed(
        original_desc, "_", 2
      )[, 1]
    }

    # Retain V1 flowFrame transformation behavior; see real-data test plan.
    if (sample_is_frame && isFALSE(transformation)) {
      channels <- flowCore::colnames(sample_ff)
      channels <- channels[grepl("Di", channels)]
      if (length(channels)) {
        sample_ff <- flowCore::transform(
          sample_ff, flowCore::transformList(
            channels, CytoNorm::cytofTransform.reverse
          )
        )
      }
    }

    sce_sample <- CATALYST::prepData(sample_ff)
    sample_comp <- CATALYST::compCytof(
      sce_sample, sm, overwrite = FALSE, isotope_list = isotope_list
    )
    out_fcs <- CATALYST::sce2fcs(sample_comp, assay = "compexprs")
    channels <- flowCore::colnames(out_fcs)
    channels <- channels[grepl("Di", channels)]
    if (length(channels)) {
      out_fcs <- flowCore::transform(
        out_fcs, flowCore::transformList(
          channels, CytoNorm::cytofTransform.reverse
        )
      )
    }
    if (!is.null(original_desc)) out_fcs@parameters$desc <- original_desc

    if (sample_is_frame) return(out_fcs)
    flowCore::write.FCS(out_fcs, filename = output_paths[i])
  }
  invisible(NULL)
}
