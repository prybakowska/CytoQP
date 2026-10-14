#' Visualize batch using umap dimensional reduction
#'
#' @description Plots batch effect using UMAP and clustering markers
#'
#' @param files_before_norm Character, full path to the unnormalized fcs_files.
#' @param files_after_norm Character, full path to the normalized fcs_files.
#' @param cores Number of cores to be used
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), CytoNormed).
#' @param clustering_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted.
#' These markers are used for building and plotting UMAP.
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5, default is set to TRUE.
#' @param batch_labels Character vector containing the batch label for each FCS file.
#' @param batch_pattern Character, batch pattern to be match in the fcs file name.
#' @param manual_colors Character, vector of the colors to be used,
#' the number of colors needs to be equal to the length of batch_pattern.
#' @param cells_total Number of cells to plot per each file.
#' @param transform_list Transformation list to pass to the flowCore
#' transform function, see flowCore::transformList(), if different transformation
#' than arcsine is needed. Only if arcsine_transform is FALSE. If NULL and
#' arcsine_transform = FALSE no transformation will be applied.
#' @param n_neighbors The size of local neighborhood in UMAP analysis, default
#' set to 15, as in uwot::umap().
#' It is recommended to set it to the number of files in each batch.
#'
#' @import ggplot2
#'
#' @examples
#' files_before_norm <- system.file(
#'   "extdata",
#'   paste0("Center", 1:4, "_batch1.fcs"),
#'   package = "CytoQP"
#' )
#'
#' files_after_norm <- files_before_norm
#'
#' batch_labels <- paste0(
#'   "Center",
#'   1:4
#' )
#'
#' set.seed(789)
#'
#' plot_batch(
#'   files_before_norm = files_before_norm,
#'   files_after_norm = files_after_norm,
#'   batch_labels = batch_labels,
#'   cores = 1,
#'   clustering_markers = c(
#'     "CD",
#'     "HLA",
#'     "IgD",
#'     "TCR",
#'     "Siglec"
#'   ),
#'   out_dir = file.path(tempdir(), "CytoNormed")
#' )
#'
#' @export
#'
#' @return save plots for batch effect in the out_dir

plot_batch <- function(files_before_norm,
                       files_after_norm,
                       batch_labels = NULL,
                       batch_pattern = NULL,
                       cores = 1,
                       out_dir = NULL,
                       clustering_markers = "CD|HLA|IgD|PD|BAFF|TCR",
                       arcsine_transform = TRUE,
                       manual_colors = NULL,
                       cells_total = 1000,
                       transform_list = NULL,
                       n_neighbors = length(files_before_norm)){

  if(!(length(files_after_norm) == length(files_before_norm))){
    stop("files_before and after does not have the same length")
  }

  fcs_files <- c(files_after_norm, files_before_norm)

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
                   basename(files_before_norm))

  files_list <- list("files_before_norm" = files_before_norm,
                     "files_after_norm" = files_after_norm)


  if(is.null(batch_labels) & is.null(batch_pattern)){
    stop("define batch_labels or batch_pattern")
  } else if (!(is.null(batch_labels)) & !(is.null(batch_pattern))){
    stop("both batch_labels and batch_pattern are defined, desellect one option by  setting to NULL")
  }

  if(!is.null(batch_labels)){
    if(length(files_after_norm) != length(batch_labels)){
      stop("The lenght of batch labels is not equal to the lenght of files")
    }
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "CytoNormed")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}


  # Parallelized analysis
  plots <- BiocParallel::bplapply(names(files_list), function(x) {
    .plot_batch_ind(name = x,
                    files = files_list[[x]],
                    batch_labels = batch_labels,
                    batch_pattern = batch_pattern,
                    out_dir = out_dir,
                    clustering_markers = clustering_markers,
                    arcsine_transform = arcsine_transform,
                    manual_colors = manual_colors,
                    cells_total = cells_total,
                    transform_list = transform_list,
                    n_neighbors = n_neighbors)},
    BPPARAM = BiocParallel::MulticoreParam(workers = cores))



  png(file.path(out_dir, "batch_effect.png"),
      width = length(plots)*1500,
      height = 1500, res = 300)
  gridExtra::grid.arrange(grobs = plots, ncol = 2)
  dev.off()
}


.plot_batch_ind <- function(name,
                            files,
                            batch_labels,
                            batch_pattern,
                            out_dir,
                            clustering_markers,
                            arcsine_transform,
                            manual_colors,
                            cells_total,
                            transform_list,
                            n_neighbors) {

  ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = files,
                                         cTotal = length(files) * cells_total,
                                         verbose = TRUE,
                                         writeMeta = FALSE,
                                         writeOutput = FALSE,
                                         outputFile = file.path(out_dir,
                                                                paste0("aggregated_for_batch_plotting.fcs")))

  if(arcsine_transform){
    ff_agg <- flowCore::transform(ff_agg,
                                  flowCore::transformList(grep("Di", flowCore::colnames(ff_agg),
                                                               value = TRUE),
                                                          CytoNorm::cytofTransform))
  } else if (!is.null(transform_list)){
    ff_agg <- flowCore::transform(ff_agg, transform_list)
  } else {
    ff_agg <- ff_agg
  }

  markers <- FlowSOM::GetMarkers(ff_agg, flowCore::colnames(ff_agg))

  cl_markers <- paste(clustering_markers, collapse="|")
  cl_markers <- grep(cl_markers, markers, value = T)

  ff_agg@exprs[, names(cl_markers)] <- apply(ff_agg@exprs[, names(cl_markers)],
                                             2, function(x){
                                               q <- stats::quantile(x, 0.9999)
                                               x[x > q] <- q
                                               x
                                             })

  samp <- length(files)
  ff_samp <- ff_agg@exprs[sample(nrow(ff_agg@exprs), samp*cells_total), ]

  dimred_res <- uwot::umap(X = ff_samp[, names(cl_markers)],
                           n_neighbors = n_neighbors, scale = TRUE)

  dimred_df <- data.frame(dim1 = dimred_res[,1], dim2= dimred_res[,2],
                          ff_samp[, names(cl_markers)])

  dimred_df$file_id <- ff_samp[,"File"]

  if(!(is.null(batch_pattern))){
    dimred_df$batch <- sapply(files[dimred_df$file_id], function(file) {
      stringr::str_match(file, batch_pattern)[,1]})
  }

  if (!is.null(batch_labels)) {
    dimred_df$batch <- batch_labels[
      as.integer(dimred_df$file_id)
    ]
  }


  p <- ggplot2::ggplot(dimred_df,  aes_string(x = "dim1", y = "dim2", color = "batch")) +
    geom_point(aes(color = batch), size = 3, position="jitter") +
    ggtitle(name)+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    theme(panel.background = element_rect(fill = "white", colour = "black",
                                          size = 2, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(color="black", size=26,
                                       hjust = 0.95, face = "bold"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text(size = 23, color = "black"),
          strip.background = element_rect(fill = "white"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22),
          legend.position = "bottom",
          # legend.key.size = unit(3,"point"),
          legend.key = element_blank())

  if (!is.null(manual_colors)){
    p <- p+scale_color_manual(values = manual_colors)
  }

  return(p)
}


#' Extracts percentages and MSI for cell populations
#'
#' @description Performs FlowSOM clustering and extracts cluster and metacluster
#' frequency and MSI. It is imputing 0 values when NAs are detected in MSI for
#' clusters and metaclusters.
#'
#' @param files Pathway to the files
#' @param nCells Numeric, number of cells to be cluster per each file,
#' default is set to 50 000.
#' @param phenotyping_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted.
#' @param functional_markers Character vector, marker names to be used for
#' functional markers, can be full marker name
#' e.g. "IL-6" or "IL" if all IL-markers needs to be plotted.
#' @param xdim Numeric, parameter to pass to FlowSOM, width of the SOM grid,
#' default is set to 10.
#' @param ydim Numeric, parameter to pass to FlowSOM, geight of the SOM grid,
#' default is set to 10.
#' @param n_metaclusters Numeric, exact number of clusters for metaclustering
#' in FlowSOM, default is set to 35.
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#' be saved, default is set to file.path(getwd(), "CytoNormed").
#' @param seed Numeric, set to obtain reproducible results, default is set to NULL.
#' @param arcsine_transform arcsine_transform Logical, if the data should
#' be transformed with arcsine transformation and cofactor 5, default is set to TRUE
#' @param save_matrix Logical, if the results should be saved, if TRUE (default)
#' matrice will be saved in out_dir.
#' @param file_name The prefix/name to use for saving the either flowsom result
#' or frequency matrix, only if save_matrix or save_flowsom_result set to TRUE.
#' deafult NULL.
#' @param transform_list Transformation list to pass to the flowCore
#' transform function, see flowCore::transformList, if different transformation
#' than arcsine is needed. Only if arcsine_transform is FALSE. If NULL and
#' arcsine_transform = FALSE no transformation will be applied.Default set to NULL.
#' @param save_flowsom_result Logical, if FlowSOM and FlowSOM plots should be
#' saved. If TRUE (default) files will be saved in out_dir.
#' @param impute_0_values Logical, if 0 values should be imputed for MSI values
#' if NAs are present.
#'
#' @import ggplot2
#'
#' @examples
#' fcs_files <- system.file(
#'   "extdata",
#'   paste0("Center", 1:4, "_batch1.fcs"),
#'   package = "CytoQP"
#' )
#'
#' mx <- extract_pctgs_msi_per_flowsom(
#'   files = fcs_files,
#'   nCells = 1000,
#'   phenotyping_markers = c(
#'     "CD",
#'     "HLA",
#'     "IgD"
#'   ),
#'   functional_markers = NULL,
#'   xdim = 5,
#'   ydim = 5,
#'   n_metaclusters = 5,
#'   save_matrix = FALSE,
#'   save_flowsom_result = FALSE,
#'   seed = 343,
#'   out_dir = file.path(tempdir(), "CytoNormed")
#' )
#'
#' @return A list containing cluster frequencies, metacluster frequencies,
#' cluster MSIs, and metacluster MSIs.
#'
#' @export

extract_pctgs_msi_per_flowsom <- function(files,
                                          nCells = 50000,
                                          phenotyping_markers = c("CD", "HLA", "IgD"),
                                          functional_markers = NULL,
                                          xdim = 10,
                                          ydim = 10,
                                          n_metaclusters = 35,
                                          out_dir = NULL,
                                          seed = NULL,
                                          arcsine_transform = TRUE,
                                          save_matrix = TRUE,
                                          file_name = NULL,
                                          transform_list = NULL,
                                          save_flowsom_result = TRUE,
                                          impute_0_values = TRUE) {

  if (!all(file.exists(files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "CytoNormed")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if(is.null(file_name)){
    f <- file_name
  } else {
    f <- paste0(file_name, "_")
  }

  nCells <- length(files) * nCells
  print(paste("aggregating files for normalization"))

  ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = files,
                                         cTotal = nCells,
                                         writeOutput = FALSE,
                                         outputFile = file.path(out_dir, paste0(f, "flowsom_aggregated_file.fcs")))


  if(arcsine_transform){
    ff_agg <- flowCore::transform(ff_agg,
                                  flowCore::transformList(grep("Di", flowCore::colnames(ff_agg),
                                                               value = TRUE),
                                                          CytoNorm::cytofTransform))
  } else if (!is.null(transform_list)){
    ff_agg <- flowCore::transform(ff_agg, transform_list)
  } else {
    ff_agg <- ff_agg
  }


  markers <- FlowSOM::GetMarkers(ff_agg, flowCore::colnames(ff_agg))
  phenotyping_channels <- grep(paste(phenotyping_markers,
                                     collapse = ("|")), markers, value = TRUE)
  functional_channels <- grep(paste(functional_markers,
                                    collapse = ("|")), markers, value = TRUE)

  # Define parameters for FlowSOM analysis
  xdim <- xdim
  ydim <- ydim
  nClus <- n_metaclusters
  s <- seed

  print(paste("building FlowSOM"))
  fsom <- FlowSOM::FlowSOM(ff_agg,
                           colsToUse = names(phenotyping_channels),
                           scale = FALSE,
                           nClus = nClus,
                           seed = s,
                           xdim = xdim,
                           ydim = ydim)

  if(save_flowsom_result){
    fsomPlot <- FlowSOM::PlotStars(fsom, backgroundValues = fsom$metaclustering)
    fsomTsne <- FlowSOM::PlotDimRed(fsom = fsom, cTotal = 5000, seed = s)

    figure <- suppressWarnings(ggpubr::ggarrange(fsomPlot, fsomTsne,
                                                 # labels = c("FlowSOM clustering", "tsne"),
                                                 ncol = 2, nrow = 1))

    ggplot2::ggsave(filename = paste0(f, "FlowSOM_Figures.pdf"), plot = figure, device = "pdf", path = out_dir,
                    width =24, height = 10)

    saveRDS(object = fsom, file = file.path(out_dir, paste0(f, "FlowSOM_model.RDS")))
    rm(figure)
  }


  # Define matrices for frequency (pctgs) calculation and MSI (msi). These calculation is performed
  # for clusters (cl) and metaclusters (mcl)
  cl_pctgs <- matrix(data = NA, nrow = length(files),
                     ncol = xdim * ydim,
                     dimnames = list(basename(files), 1:(xdim*ydim)))

  mcl_pctgs <- matrix(data = NA, nrow = length(files),
                      ncol = nClus,
                      dimnames = list(basename(files), 1:nClus))
  mfi_cl_names <- apply(expand.grid(paste0("Cl", seq_len(fsom$map$nNodes)),
                                    FlowSOM::GetMarkers(ff_agg,
                                                        unique(c(phenotyping_channels,functional_channels)))),
                        1, paste, collapse = "_")
  mfi_mc_names <- apply(expand.grid(paste0("MC", 1:nClus),
                                    FlowSOM::GetMarkers(ff_agg,
                                                        unique(c(phenotyping_channels,functional_channels)))),
                        1, paste, collapse = "_")
  cl_msi <- matrix(NA,
                   nrow = length(files),
                   ncol = fsom$map$nNodes * length(unique(names(c(phenotyping_channels,
                                                                  functional_channels)))),
                   dimnames = list(basename(files), mfi_cl_names))
  mcl_msi <- matrix(NA,
                    nrow = length(files),
                    ncol =  length(mfi_mc_names),
                    dimnames = list(basename(files), mfi_mc_names))

  rm(ff_agg)

  print(paste("calculating frequency and MSI"))

  for (i in unique(fsom$data[,"File"])){

    file <- basename(files[i])

    id <- which(fsom$data[,"File"] == i)
    fsom_subset <- FlowSOM::FlowSOMSubset(fsom = fsom, ids = id)

    cl_counts <- rep(0, xdim * ydim)
    counts_tmp <- table(FlowSOM::GetClusters(fsom_subset))
    cl_counts[as.numeric(names(counts_tmp))] <- counts_tmp

    cl_pctgs[file,] <- (cl_counts/sum(cl_counts, na.rm = T))*100

    mcl_counts <- tapply(cl_counts, fsom$metaclustering, sum)
    mcl_pctgs[file,] <- tapply(cl_pctgs[file,], fsom$metaclustering, sum)

    cluster_mfis <- FlowSOM::GetClusterMFIs(fsom_subset)
    cl_msi[file,] <- as.numeric(cluster_mfis[,unique(names(c(phenotyping_channels,functional_channels)))])
    mcluster_mfis <- as.matrix(FlowSOM::GetMetaclusterMFIs(list(FlowSOM = fsom_subset,
                                                                metaclustering = fsom$metaclustering)))
    mcl_msi[file,] <- as.numeric(mcluster_mfis[,unique(names(c(phenotyping_channels,functional_channels)))])

    rm(fsom_subset)
  }

  rm(fsom)

  if(impute_0_values){
    cl_msi <- apply(cl_msi, 2,
                    function(x){
                      missing <- which(is.na(x))
                      x[missing] <- 0
                      x
                    })

    mcl_msi <- apply(mcl_msi, 2,
                     function(x){
                       missing <- which(is.na(x))
                       x[missing] <- 0
                       x
                     })
  }

  # store the matrices in the list for convenient plotting
  all_mx <- list("Cluster_frequencies" = cl_pctgs,
                 "Metacluster_frequencies" = mcl_pctgs,
                 "Cluster_MSIs" = cl_msi,
                 "Metacluster_MSIs" = mcl_msi)

  if(save_matrix){
    saveRDS(object = all_mx,
            file = file.path(out_dir,
                             paste0(f, "cell_frequency_and_msi_using_FlowSOM.RDS")))
  }

  res <- all_mx

  # clean memory
  gc()

  return(res)

}


#' Prepares data for plotting cell frequency and MSI
#'
#' @description Performs dimensional reduction and constructs
#' data frame for plotting cell frequencies and MSI
#' per clusters and metaclusters obtained from extract_pctgs_msi_per_flowsom function.
#' For the MSI it is taking the MSI > 1.
#'
#' @param frequency_msi_list List containing matrices with cell frequency and msi
#' obtained in step extract_pctgs_msi_per_flowsom.
#' @param matrix_type The name of the matrix to be plotted.
#' @param seed Numeric set to obtain reproducible results, default NULL.
#' @param n_neighbours The size of local neighborhood in UMAP analysis, default
#' set to 15, as in uwot::umap().
#' It is recommended to set it to the number of files in each batch.
#'
#' @return data frame for plotting
#'
#' @examples
#' frequency_msi_list <- list(
#'   before = list(
#'     Metacluster_frequencies = matrix(
#'       runif(40),
#'       nrow = 10
#'     )
#'   ),
#'   after = list(
#'     Metacluster_frequencies = matrix(
#'       runif(40),
#'       nrow = 10
#'     )
#'   )
#' )
#'
#' df_plot <- prepare_data_for_plotting(
#'   frequency_msi_list = frequency_msi_list,
#'   matrix_type = "Metacluster_frequencies",
#'   n_neighbours = 3,
#'   seed = 1
#' )
#'
#' @export
prepare_data_for_plotting <- function(frequency_msi_list,
                                      matrix_type,
                                      n_neighbours = 15,
                                      seed = NULL){
  print(matrix_type)

  # set the number of closest neighborhoods
  n <- n_neighbours

  # process files before normalization
  df_b <- frequency_msi_list[["before"]][[matrix_type]]

  # select the columns for which MSI is higher than 0.2 SD
  if(grepl("msi", matrix_type, ignore.case = TRUE)){
    id_cols <-  which(apply(df_b, 2, sd) > 0.2)
    df_b <- df_b[,id_cols]
  }

  # build UMAP for files before the normalization
  if(!is.null(seed)){
    set.seed(seed)
  }
  df_b_umap <- data.frame(suppressMessages(uwot::umap(df_b, n_neighbors = n, scale = TRUE)))

  # process files after normalization
  df_a <- frequency_msi_list[["after"]][[matrix_type]]

  # select the columns for which MSI is higher than 0.2 SD
  if(grepl("msi", matrix_type, matrix_type, ignore.case = TRUE)){
    id_cols <-  which(apply(df_a, 2, sd) > 0.2)
    df_a <- df_a[,id_cols]
  }

  # build UMAP for files after the normalization
  if(!is.null(seed)){
    set.seed(seed)
  }
  df_a_umap <- data.frame(suppressMessages(uwot::umap(df_a, n_neighbors = n, scale = TRUE)))

  # extract rownames to use the for ggplot annotation
  rnmes <- c(rownames(df_b), rownames(df_a))

  #join two UMAP data frames
  dr <- data.frame(rbind(df_b_umap, df_a_umap), check.names = F)
  colnames(dr) <- c("dim1", "dim2")

  dr$normalization <- c(rep("Raw", length(rownames(df_b_umap))),
                        rep("Normalized", length(rownames(df_a_umap))))

  return(dr)
}

#' Plot the distribution of the features across batches in two dimensional space
#'
#' @description Uses ggplot draw the distribution of the samples across batches
#'
#' @param df_plot Data frame that contains dim1 and dim2 obtained upon dimensional
#' reduction. This data frame is generated by prepare_data_for_plotting.
#' @param fill As in ggplot2, defines by which variable the points are colored.
#' @param shape As in ggplot2, defines by which variable the shapes are plotted.
#' @param color As in ggplot2, defines by which variable the borderd of the
#' points are colored.
#' @param split_by_normalization Logical, if TRUE (default), plots are split by normalization status.
#' @param fill_legend_name Character, specifies the name of the legend for fill.
#' @param color_legend_name Character, specifies the name of the legend for fill.
#' @param shape_legend_name Character, specifies the name of the legend for fill.
#' @param title Character, the title of the plot.
#' @param manual_colors Character, vector of the colors to be used,
#' the number of colors needs to be equal to the length of batch_pattern.
#' @examples
#' df_plot <- data.frame(
#'   dim1 = rnorm(20),
#'   dim2 = rnorm(20),
#'   normalization = rep(
#'     c("Raw", "Normalized"),
#'     each = 10
#'   ),
#'   batch = rep(
#'     c("batch1", "batch2"),
#'     10
#'   )
#' )
#'
#' p <- plot_batch_using_freq_msi(
#'   df_plot = df_plot,
#'   color = df_plot$batch,
#'   split_by_normalization = TRUE,
#'   title = "Batch effect"
#' )
#'
#' @return ggplot
#' @export
plot_batch_using_freq_msi <- function(df_plot,
                                      fill = NULL,
                                      shape = NULL,
                                      color = NULL,
                                      split_by_normalization = TRUE,
                                      fill_legend_name = NULL,
                                      color_legend_name = NULL,
                                      shape_legend_name = NULL,
                                      title = NULL,
                                      manual_colors = NULL){

  if (split_by_normalization) {
    df_plot$normalization <- factor(df_plot$normalization,
                                    levels = c("Raw", "Normalized"))
    p <- ggplot(df_plot, aes(x = dim1, y = dim2)) + geom_point(data = df_plot,
                                                               aes(x = dim1, y = dim2, fill = fill, shape = shape,
                                                                   color = color), size = 3) + ggtitle(title) +
      facet_wrap(~normalization)
    if(!is.null(manual_colors)){
      p <- p + scale_color_manual(values = manual_colors)
    }
  }
  else {
    p <- ggplot(df_plot, aes(x = dim1, y = dim2)) +
      geom_point(data = df_plot,
                 aes(x = dim1, y = dim2, fill = fill, shape = shape,
                     color = color), size = 3) + ggtitle(title)
    if(!is.null(manual_colors)){
      p <- p + scale_color_manual(values = manual_colors)
    }
  }
  p <- p + theme(panel.background = element_rect(fill = "white",
                                                 colour = "black", size = 1, linetype = "solid"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), legend.position = "right",
                 legend.key = element_blank(), strip.background = element_rect(fill = "white",
                                                                               colour = "black"))
  p <- p + labs(fill = fill_legend_name, color = color_legend_name,
                shape = shape_legend_name)
  return(p)
}


#' Plot listed plots with one common legend
#'
#' @param plot_lists The list of ggplots
#' @param nrow Numeric, number of rows to arrange the plots
#' @param ncol Numeric, number of columns to arrange the plots.
#' @param position Character, Position of the legend, either "bottom" or "right".
#'
#' @return combained ggplot
#' @export
#'
#' @import gridExtra
#' @import ggplot2
#'
#' @examples
#' p1 <- ggplot2::ggplot(
#'   mtcars,
#'   ggplot2::aes(x = wt, y = mpg, color = factor(cyl))
#' ) +
#'   ggplot2::geom_point()
#'
#' p2 <- ggplot2::ggplot(
#'   mtcars,
#'   ggplot2::aes(x = hp, y = mpg, color = factor(cyl))
#' ) +
#'   ggplot2::geom_point()
#'
#' combined <- grid_arrange_common_legend(
#'   plot_lists = list(p1, p2),
#'   nrow = 1,
#'   ncol = 2,
#'   position = "right"
#' )
#'
grid_arrange_common_legend <- function(plot_lists,
                                       nrow = 1,
                                       ncol = length(plot_lists),
                                       position = c("bottom", "right")) {

  position <- match.arg(position)
  g <- ggplotGrob(plot_lists[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plot_lists, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = grid::unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = grid::unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid::grid.newpage()
  grid::grid.draw(combined)
  return(combined)

}

