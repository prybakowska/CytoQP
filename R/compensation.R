# ==============================================================================
# CytoQP V2
# Spillover estimation, visualization, and CyTOF compensation
# ==============================================================================


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

#' Plot a Spillover Matrix
#'
#' Internal helper used by plot_spillover().
#'
#' @param spillover_matrix Numeric matrix or data.frame containing spillover
#'   values expressed as percentages.
#' @param colors Character vector of four colors used for the gradient.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @keywords internal
.plot_spillover_matrix <- function(
    spillover_matrix,
    colors = c("white", "lightcoral", "red2", "darkred"),
    title = "Spillover Matrix"
) {

  spillover_matrix <- as.matrix(spillover_matrix)

  if (!is.numeric(spillover_matrix)) {
    stop(
      "spillover_matrix must contain numeric values.",
      call. = FALSE
    )
  }

  if (length(colors) != 4L) {
    stop(
      "colors must contain exactly four colors.",
      call. = FALSE
    )
  }

  row_totals <- rowSums(
    spillover_matrix,
    na.rm = TRUE
  )

  col_totals <- colSums(
    spillover_matrix,
    na.rm = TRUE
  )

  rows_to_remove <- which(row_totals <= 100)
  cols_to_remove <- which(col_totals <= 100)

  if (
    length(rows_to_remove) > 0L &&
    length(rows_to_remove) < nrow(spillover_matrix)
  ) {
    spillover_matrix <- spillover_matrix[
      -rows_to_remove,
      ,
      drop = FALSE
    ]
  }

  if (
    length(cols_to_remove) > 0L &&
    length(cols_to_remove) < ncol(spillover_matrix)
  ) {
    spillover_matrix <- spillover_matrix[
      ,
      -cols_to_remove,
      drop = FALSE
    ]
  }

  n_rows <- nrow(spillover_matrix)
  n_cols <- ncol(spillover_matrix)

  off_diagonal_values <- spillover_matrix[
    spillover_matrix != 100
  ]

  off_diagonal_values <- off_diagonal_values[
    is.finite(off_diagonal_values)
  ]

  if (length(off_diagonal_values) == 0L) {

    fill_limits <- range(
      spillover_matrix,
      na.rm = TRUE
    )

  } else {

    fill_limits <- range(
      off_diagonal_values,
      na.rm = TRUE
    )
  }

  if (
    length(fill_limits) != 2L ||
    any(!is.finite(fill_limits))
  ) {
    fill_limits <- c(0, 1)
  }

  if (diff(fill_limits) == 0) {
    fill_limits <- fill_limits + c(-0.01, 0.01)
  }

  plot_data <- as.data.frame(
    as.table(spillover_matrix),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Var2,
      y = Var1,
      fill = Freq
    )
  ) +
    ggplot2::geom_tile(
      color = "white"
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c(
        colors[1:3],
        "grey",
        colors[4]
      ),
      limits = fill_limits,
      guide = "none"
    ) +
    ggplot2::geom_hline(
      yintercept = seq(
        0.5,
        n_rows + 0.5
      ),
      color = "black",
      linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = seq(
        0.5,
        n_cols + 0.5
      ),
      color = "black",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(
          Freq != 0 & Freq != 100,
          as.character(Freq),
          ""
        )
      ),
      color = "black",
      size = 2
    ) +
    ggplot2::coord_fixed(
      ratio = 1
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 22
      )
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    )
}


#' Extend the CATALYST Isotope List
#'
#' Detects isotopes present in the input but absent from the default
#' CATALYST isotope list and adds them automatically.
#'
#' @param x A flowFrame, matrix, or data.frame.
#'
#' @return A named numeric isotope list.
#' @keywords internal
.get_isotope_list <- function(x) {

  isotope_list <- CATALYST::isotope_list

  if (inherits(x, "flowFrame")) {

    marker_names <- flowCore::colnames(x)

    marker_names <- marker_names[
      grepl("Di", marker_names)
    ]

  } else if (
    is.matrix(x) ||
    is.data.frame(x)
  ) {

    marker_names <- colnames(x)

  } else {

    stop(
      "x must be a flowFrame, matrix, or data.frame.",
      call. = FALSE
    )
  }

  if (
    is.null(marker_names) ||
    length(marker_names) == 0L
  ) {
    return(isotope_list)
  }

  input_isotopes <- unique(
    stringr::str_extract(
      marker_names,
      "[A-Za-z]+"
    )
  )

  input_isotopes <- input_isotopes[
    !is.na(input_isotopes)
  ]

  missing_isotopes <- setdiff(
    input_isotopes,
    names(isotope_list)
  )

  for (isotope in missing_isotopes) {

    matched_channel <- marker_names[
      grepl(
        isotope,
        marker_names,
        fixed = TRUE
      )
    ][1]

    isotope_mass <- suppressWarnings(
      as.numeric(
        gsub(
          "\\D",
          "",
          matched_channel
        )
      )
    )

    if (
      length(isotope_mass) == 1L &&
      !is.na(isotope_mass)
    ) {
      isotope_list[isotope] <- isotope_mass
    }
  }

  isotope_list
}


#' Estimate a CyTOF Spillover Matrix
#'
#' Internal helper using the CATALYST compensation workflow.
#'
#' @param input A flowFrame containing single-stained compensation beads.
#' @param channels Optional character vector containing channels to use.
#'
#' @return A numeric spillover matrix.
#' @keywords internal
.compute_spillover <- function(
    input,
    channels = NULL
) {

  if (!inherits(input, "flowFrame")) {
    stop(
      "input must be a flowFrame.",
      call. = FALSE
    )
  }

  available_channels <- flowCore::colnames(input)

  if (is.null(channels)) {

    channels <- available_channels[
      grepl(
        "Di",
        available_channels
      )
    ]
  }

  if (length(channels) == 0L) {
    stop(
      "No CyTOF channels were detected.",
      call. = FALSE
    )
  }

  missing_channels <- setdiff(
    channels,
    available_channels
  )

  if (length(missing_channels) > 0L) {
    stop(
      paste0(
        "The following channels were not found: ",
        paste(
          missing_channels,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  barcode_key <- as.numeric(
    gsub(
      "\\D",
      "",
      channels
    )
  )

  if (anyNA(barcode_key)) {
    stop(
      "Mass numbers could not be extracted from all selected channels.",
      call. = FALSE
    )
  }

  barcode_key <- sort(
    unique(barcode_key)
  )

  sce <- CATALYST::prepData(
    input
  )

  sce <- CATALYST::assignPrelim(
    sce,
    barcode_key
  )

  sce <- CATALYST::applyCutoffs(
    CATALYST::estCutoffs(
      sce
    )
  )

  sce <- CATALYST::computeSpillmat(
    sce
  )

  spillover_matrix <-
    S4Vectors::metadata(sce)$spillover_matrix

  if (is.null(spillover_matrix)) {
    stop(
      "CATALYST did not return a spillover matrix.",
      call. = FALSE
    )
  }

  spillover_matrix
}


# ------------------------------------------------------------------------------
# Exported functions
# ------------------------------------------------------------------------------

#' Calculate a CyTOF Spillover Matrix
#'
#' @description
#' Calculates a spillover matrix from single-stained CyTOF compensation beads
#' using the CATALYST compensation workflow.
#'
#' @param input Path to a compensation-bead FCS file or a \code{flowFrame}.
#' @param channels Optional character vector specifying the CyTOF channels
#' used for spillover estimation. If \code{NULL}, channels containing
#' \code{"Di"} are used.
#' @param save_matrix Optional path where the calculated spillover matrix
#' should be saved as a CSV file. If \code{NULL}, the matrix is not written
#' to disk.
#'
#' @return
#' A numeric matrix containing the estimated spillover values.
#'
#' @examples
#' compensation_file <- system.file(
#'   "extdata",
#'   "Compensation_beads.fcs",
#'   package = "CytoQP"
#' )
#'
#' spillover_matrix <- calculate_spillover(
#'   input = compensation_file
#' )
#'
#' @export
calculate_spillover <- function(
    input,
    channels = NULL,
    save_matrix = NULL
) {

  if (inherits(input, "flowFrame")) {

    compensation_beads <- input

  } else if (
    is.character(input) &&
    length(input) == 1L &&
    file.exists(input) &&
    flowCore::isFCSfile(input)
  ) {

    compensation_beads <- flowCore::read.FCS(
      input,
      transformation = FALSE,
      truncate_max_range = FALSE
    )

  } else {

    stop(
      "input must be a flowFrame or a valid FCS file.",
      call. = FALSE
    )
  }

  spillover_matrix <- .compute_spillover(
    input = compensation_beads,
    channels = channels
  )

  if (!is.null(save_matrix)) {

    utils::write.csv(
      spillover_matrix,
      file = save_matrix,
      row.names = TRUE
    )
  }

  spillover_matrix
}


#' Plot a CyTOF Spillover Matrix
#'
#' @description
#' Visualizes an existing CyTOF spillover matrix.
#'
#' @param input A numeric matrix or data.frame containing a spillover matrix.
#' @param colors Character vector of four colors used for the heatmap.
#' @param title Plot title.
#'
#' @return
#' A \code{ggplot} object containing the spillover heatmap.
#'
#' @examples
#' compensation_file <- system.file(
#'   "extdata",
#'   "Compensation_beads.fcs",
#'   package = "CytoQP"
#' )
#'
#' spillover_matrix <- calculate_spillover(
#'   input = compensation_file
#' )
#'
#' plot_spillover(
#'   input = spillover_matrix
#' )
#'
#' @export
plot_spillover <- function(
    input,
    colors = c(
      "white",
      "lightcoral",
      "red2",
      "darkred"
    ),
    title = "Spillover Matrix"
) {

  if (
    !is.matrix(input) &&
    !is.data.frame(input)
  ) {
    stop(
      "input must be a spillover matrix or data.frame.",
      call. = FALSE
    )
  }

  spillover_matrix <- as.matrix(
    input
  )

  if (!is.numeric(spillover_matrix)) {
    stop(
      "The spillover matrix must contain numeric values.",
      call. = FALSE
    )
  }

  # CATALYST spillover values are proportions.
  # Convert them to percentages for visualization.
  plot_matrix <- round(
    spillover_matrix * 100,
    digits = 1
  )

  .plot_spillover_matrix(
    spillover_matrix = plot_matrix,
    colors = colors,
    title = title
  )
}


#' Compensate CyTOF Data
#'
#' @description
#' Applies a previously calculated spillover matrix to CyTOF FCS files using
#' the CATALYST compensation workflow.
#'
#' @param samples A \code{flowFrame} or character vector containing paths to
#' the FCS files to be compensated.
#' @param compensation Numeric matrix or data.frame containing the spillover
#' matrix, typically generated with \code{calculate_spillover()}.
#' @param transformed Logical indicating whether the supplied sample data are
#' already transformed. Default is \code{TRUE}.
#' @param isotope_list Optional named vector containing isotope masses. If
#' \code{NULL}, it is generated automatically from the input data.
#' @param output_dir Character specifying the directory where compensated FCS
#' files are written. Default is \code{"Compensated_files"}.
#' @param suffix_old Character string to replace in input filenames.
#' @param suffix_new Character string used in compensated output filenames.
#' Default is \code{"compensated"}.
#'
#' @return
#' If \code{samples} is a \code{flowFrame}, a compensated \code{flowFrame} is
#' returned. If \code{samples} contains FCS file paths, compensated files are
#' written to \code{output_dir} and their paths are returned invisibly.
#'
#' @examples
#' fcs_files <- system.file(
#'   "extdata",
#'   paste0("Center", 1:4, "_batch1.fcs"),
#'   package = "CytoQP"
#' )
#'
#' compensation_file <- system.file(
#'   "extdata",
#'   "Compensation_beads.fcs",
#'   package = "CytoQP"
#' )
#'
#' spillover_matrix <- calculate_spillover(
#'   input = compensation_file
#' )
#'
#' files_compensated <- compensate_cytof(
#'   samples = fcs_files,
#'   compensation = spillover_matrix,
#'   output_dir = file.path(tempdir(), "Compensated_files")
#' )
#'
#' @export
compensate_cytof <- function(
    samples,
    compensation,
    transformed = TRUE,
    isotope_list = NULL,
    output_dir = "Compensated_files",
    suffix_old = "beadNorm",
    suffix_new = "compensated"
) {

  valid_samples <-
    inherits(samples, "flowFrame") ||
    (
      is.character(samples) &&
        length(samples) > 0L
    )

  if (!valid_samples) {
    stop(
      "samples must be a flowFrame or character vector of FCS paths.",
      call. = FALSE
    )
  }

  if (
    !is.logical(transformed) ||
    length(transformed) != 1L ||
    is.na(transformed)
  ) {
    stop(
      "transformed must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Prepare spillover matrix
  # ---------------------------------------------------------------------------

  if (
    !is.matrix(compensation) &&
    !is.data.frame(compensation)
  ) {
    stop(
      "compensation must be a spillover matrix or data.frame.",
      call. = FALSE
    )
  }

  spillover_matrix <- as.matrix(
    compensation
  )

  if (!is.numeric(spillover_matrix)) {
    stop(
      "The spillover matrix must contain numeric values.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Prepare isotope list
  # ---------------------------------------------------------------------------

  if (is.null(isotope_list)) {

    if (inherits(samples, "flowFrame")) {

      isotope_list <- .get_isotope_list(
        samples
      )

    } else {

      missing_files <- samples[
        !file.exists(samples)
      ]

      if (length(missing_files) > 0L) {
        stop(
          paste0(
            "The following FCS files do not exist:\n",
            paste(
              missing_files,
              collapse = "\n"
            )
          ),
          call. = FALSE
        )
      }

      first_sample <- flowCore::read.FCS(
        samples[1],
        transformation = FALSE,
        truncate_max_range = FALSE
      )

      isotope_list <- .get_isotope_list(
        first_sample
      )
    }
  }

  # Required for CATALYST compensation compatibility
  if (!"BCKG" %in% names(isotope_list)) {
    isotope_list[["BCKG"]] <- 190
  }


  # ---------------------------------------------------------------------------
  # Helper: compensate one flowFrame
  # ---------------------------------------------------------------------------

  compensate_one <- function(sample_ff) {

    if (!inherits(sample_ff, "flowFrame")) {
      stop(
        "sample_ff must be a flowFrame.",
        call. = FALSE
      )
    }

    original_desc <- sample_ff@parameters$desc

    # CATALYST requires clean marker descriptions
    if (!is.null(original_desc)) {

      cleaned_desc <- stringr::str_split_fixed(
        original_desc,
        "_",
        2
      )[, 1]

      sample_ff@parameters$desc <- cleaned_desc
    }

    # Preserve original CytoQP transformation behavior
    if (!transformed) {

      cytof_channels <- flowCore::colnames(sample_ff)[
        grepl(
          "Di",
          flowCore::colnames(sample_ff)
        )
      ]

      if (length(cytof_channels) > 0L) {

        sample_ff <- flowCore::transform(
          sample_ff,
          flowCore::transformList(
            cytof_channels,
            CytoNorm::cytofTransform.reverse
          )
        )
      }
    }

    sce <- CATALYST::prepData(
      sample_ff
    )

    compensated <- CATALYST::compCytof(
      sce,
      spillover_matrix,
      overwrite = FALSE,
      isotope_list = isotope_list
    )

    output_ff <- CATALYST::sce2fcs(
      compensated,
      assay = "compexprs"
    )

    cytof_channels <- flowCore::colnames(output_ff)[
      grepl(
        "Di",
        flowCore::colnames(output_ff)
      )
    ]

    if (length(cytof_channels) > 0L) {

      output_ff <- flowCore::transform(
        output_ff,
        flowCore::transformList(
          cytof_channels,
          CytoNorm::cytofTransform.reverse
        )
      )
    }

    if (!is.null(original_desc)) {
      output_ff@parameters$desc <- original_desc
    }

    output_ff
  }


  # ---------------------------------------------------------------------------
  # Single flowFrame
  # ---------------------------------------------------------------------------

  if (inherits(samples, "flowFrame")) {

    return(
      compensate_one(
        samples
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Multiple FCS files
  # ---------------------------------------------------------------------------

  missing_files <- samples[
    !file.exists(samples)
  ]

  if (length(missing_files) > 0L) {
    stop(
      paste0(
        "The following FCS files do not exist:\n",
        paste(
          missing_files,
          collapse = "\n"
        )
      ),
      call. = FALSE
    )
  }

  if (!dir.exists(output_dir)) {

    dir.create(
      output_dir,
      recursive = TRUE
    )
  }

  output_files <- character(
    length(samples)
  )

  for (i in seq_along(samples)) {

    file <- samples[i]

    message(
      sprintf(
        "Compensating [%d/%d]: %s",
        i,
        length(samples),
        basename(file)
      )
    )

    sample_ff <- flowCore::read.FCS(
      file,
      transformation = FALSE,
      truncate_max_range = FALSE
    )

    output_ff <- compensate_one(
      sample_ff
    )

    output_name <- gsub(
      suffix_old,
      suffix_new,
      basename(file),
      fixed = TRUE
    )

    if (identical(
      output_name,
      basename(file)
    )) {

      output_name <- sub(
        "(?i)\\.fcs$",
        paste0(
          "_",
          suffix_new,
          ".fcs"
        ),
        basename(file),
        perl = TRUE
      )
    }

    output_path <- file.path(
      output_dir,
      output_name
    )

    flowCore::write.FCS(
      output_ff,
      filename = output_path
    )

    output_files[i] <- output_path
  }

  invisible(
    output_files
  )
}

