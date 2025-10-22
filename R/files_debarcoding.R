#' Debarcodes files
#'
#' @description Performs sample debarcoding.
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param cores Number of cores to be used
#' @param file_batch_id Character vector with batch label for each fcs_file,
#' the order and the length needs to be the same as in fcs_files. If only batch
#' is processed can be prepared as e.g. file_batch_id <- rep("batch", length(files))
#' @param file_score Data frame with quality scores obtained from
#' file_quality_check.Default set to NULL.
#' @param out_dir Character, pathway to where the plots should be saved,
#' only if argument to_plot = TRUE, default is set to working directory
#' @param min_threshold Logical, if the minimal threshold for barcoding
#' should be applied.Default set to TRUE.
#' @param threshold Numeric, value for the minimum threshold for debarcoding,
#' default is set to 0.18, only if min_threshold set to TRUE.
#' @param to_plot Logical, if plots for yields and debarcoding quality should
#' be plotted.
#' @param barcodes_used Character vector with the names of the barcodes that
#' were used, eg. barcode 1 is the same as A1. Or a list with the barcodes
#' name per batch. If NULL (default) all the barcodes contained in
#' sample_key will be used, regarding the batch.
#' @param less_than_th Logical, if the name of the files for which lower threshold
#' than set in parameter threshold was detected. Default is set to FALSE.
#' @param barcode_key matrix as in CATALYST::assignPrelim, the debarcoding scheme.
#' A binary matrix with sample names as row names and numeric masses as column names
#' OR a vector of numeric masses corresponding to barcode channels.
#' When the latter is supplied, 'assignPrelim' will create a scheme of the
#' appropriate format internally.
#'
#' @return Save debarcoded fcs files in out_dir. If parameter to_plot set
#' to TRUE, save plots for yields and debarcodig quality in out_dir. If less_than_th
#' set to TRUE, save file names for which threshold lower than in parameter threshold
#' was detected "files_with_lower_debarcoding_threshold.RDS" in out_dir.
#'
#' @examples
#' # Set input directory
#'clean_dir <- file.path(dir, "Cleaned")
#'
#'# Define files for debarcoding
#' files <- list.files(clean_dir,
#'                    pattern = "_cleaned.fcs$",
#'                    full.names = TRUE)
#'
#'# Read in file scores if calculated
#'file_scores <- readRDS(list.files(path = dir,
#'                                  recursive = TRUE,
#'                                  full.names = TRUE,
#'                                  pattern = "Quality_AOF_score.RDS"))
#'
#'# Define file batch ID for each file
#'file_batch_id <- stringr::str_match(basename(files),
#'                                    "(day[0-9]*).*.fcs")[,2]
#'
#'# Read in metadata
#'md <- utils::read.csv(file.path(dir, "RawFiles", "meta_data.csv"))
#'
#'# read in barcode key
#'sample_key <- CATALYST::sample_key

#'# Extract information about barcodes used in each batch
#'barcodes_list <- list()
#'for (batch in unique(file_batch_id)){
#'  idx <- md[md[,"BATCH"] == batch, "BARCODE"]
#'  barcodes_list[[batch]] <- rownames(sample_key)[idx]
#'}
#'
#'# Debarcode files
#'debarcode_files(fcs_files = files,
#'                out_dir = NULL,
#'                file_score = file_scores,
#'                min_threshold = TRUE,
#'                barcodes_used = barcodes_list,
#'                file_batch_id = file_batch_id,
#'                less_than_th = TRUE,
#'                barcode_key = sample_key)
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
                            barcode_key = NULL){


  if(is(fcs_files, "list")){
    fcs_files <- unlist(fcs_files)
  }

  if(anyDuplicated(fcs_files) != 0){
    stop("names of fcs files are duplicated")
  }

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if(length(file_batch_id) != length(fcs_files)){
    stop("the lenght of the file_batch_id is not equal to the lenght of fcs_files")
  }

  if(!is.null(file_score)){
    if(!inherits(file_score, "data.frame")) {
      stop("file_scores is not a data frame")
    } else {
      # Select good quality files
      good_files <- file_scores$file_names[file_scores$quality == "good"]
      fcs_files_clean <- fcs_files[basename(fcs_files) %in% good_files]
      file_batch_id <- file_batch_id[basename(fcs_files) %in% good_files]
      fcs_files <- fcs_files_clean
    }
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Debarcoded")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  outFiles <-  BiocParallel::bplapply(fcs_files, function(file) {
    .debarcode_ind(file,
                   fcs_files = fcs_files,
                   file_batch_id = file_batch_id,
                   file_score = file_score,
                   out_dir = out_dir,
                   min_threshold = min_threshold,
                   threshold = threshold,
                   to_plot = to_plot,
                   barcodes_used = barcodes_used,
                   less_than_th = less_than_th,
                   barcode_key = barcode_key)
  },
  BPPARAM = BiocParallel::MulticoreParam(workers = cores))

  lessFiles <- lapply(outFiles, function(x) {return(x[[1]])})
  debarcodedFiles <- lapply(outFiles, function(x) {return(x[[2]])})

  if(less_than_th){
    saveRDS(unlist(lessFiles), file.path(out_dir, "files_with_lower_debarcoding_threshold.RDS"))
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
                           barcode_key) {

  print(paste0("   ", Sys.time()))
  print(paste0("   Debarcoding ", file))
  ff <- flowCore::read.FCS(file, transformation = FALSE)

  file_id <- which(file == fcs_files)
  batch_id <- file_batch_id[file_id]

  if(!is.null(barcodes_used)){
    if(is.list(barcodes_used)){
      s_key <- barcode_key[rownames(barcode_key) %in% barcodes_list[[batch_id]],]
    } else {
      s_key <- barcode_key[rownames(barcode_key) %in% barcodes_used,]
    }

  } else {
    s_key <- barcode_key
  }

  dat <- CATALYST::prepData(ff)
  dat <- CATALYST::assignPrelim(dat, bc_key = s_key)
  rownames(dat)[SummarizedExperiment::rowData(dat)$is_bc]
  # table(colData(dat)$bc_id)
  dat <- CATALYST::estCutoffs(dat)

  if (min_threshold){
    if(any(S4Vectors::metadata(dat)$sep_cutoffs < threshold)){
      warning(paste0("cutoff lower than 0.18 has been detected for ", basename(file),
                     ", cutoff will be set to 0.18"))
      fileOut <- basename(file)
    }

    id <- S4Vectors::metadata(dat)$sep_cutoffs < threshold
    S4Vectors::metadata(dat)$sep_cutoffs[id] <- threshold

  } else {
    if(any(S4Vectors::metadata(dat)$sep_cutoffs < threshold)){
      warning(paste0("cutoff lower than ", threshold, " detected for ", basename(file)))
      fileOut <- basename(file)
    }
  }

  id <- is.na(S4Vectors::metadata(dat)$sep_cutoffs)
  S4Vectors::metadata(dat)$sep_cutoffs[id] <- 1

  if (to_plot){
    p <- CATALYST::plotYields(dat, which = rownames(s_key))

    pdf(file.path(out_dir, paste(gsub(".fcs", "_yields.pdf", basename(file)))))
    for (name in names(p)){
      print(p[[name]])
    }
    dev.off()
  }

  dat <- CATALYST::applyCutoffs(dat)

  if (to_plot){
    p <- CATALYST::plotEvents(dat, n = 500)

    pdf(file.path(out_dir, paste(gsub(".fcs", "_debarcode_quality.pdf",
                                      basename(file)))))
    for (name in names(p)){
      print(p[[name]])
    }
    dev.off()
  }

  dat <- dat[, dat$bc_id !=0]
  fs <- CATALYST::sce2fcs(dat, split_by = "bc_id")

  if(is.null(file_dir)){
    file_dir <- file.path(out_dir, batch_id)
  }
  if(!dir.exists(file_dir)){dir.create(file_dir)}

  file_name <- gsub("_cleaned.fcs|.fcs", "", basename(file))

  flowCore::write.flowSet(fs, outdir = file_dir,
                          filename = paste0(rownames(fs@phenoData), "_", file_name,
                                            "_debarcoded.fcs"))

  filePaths <- file.path(file_dir, paste0(rownames(fs@phenoData), "_", file_name,
                                         "_debarcoded.fcs"))

  return(list(fileOut, filePaths))
}
