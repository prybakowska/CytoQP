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
#' # Define files before normalization
#' gate_dir <- file.path(dir, "Gated")
#' files_before_norm <- list.files(gate_dir,
#'                                 pattern = ".fcs",
#'                                 full.names = T)
#'
#' # Define files after normalization
#' norm_dir <- file.path(dir, "CytoNormed")
#' files_after_norm <- list.files(norm_dir,
#'                                pattern = ".fcs",
#'                                full.names = T)
#'
#' # files needs to be in the same order, check and order if needed
#' test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
#'                  basename(files_before_norm))
#'
#' batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]
#'
#' # Plot batch effect
#' set.seed(789)
#' plot_batch(files_before_norm = files_before_norm,
#'            files_after_norm = files_after_norm,
#'            batch_labels = batch_labels,
#'            cores = 1,
#'            out_dir = norm_dir,
#'            clustering_markers = c("CD", "IgD", "HLA"),
#'            manual_colors = c("darkorchid4", "darkorange", "chartreuse4"))
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



  png(file.path(norm_dir, "batch_effect.png"),
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

  dimred_df$file_id <- ff_samp[,"File2"]

  if(!(is.null(batch_pattern))){
    dimred_df$batch <- sapply(files[dimred_df$file_id], function(file) {
      stringr::str_match(file, batch_pattern)[,1]})
  }

  if(!(is.null(batch_labels))){
    dimred_df$batch <- batch_labels
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
#' @param file_list List, pathway to the files before and after normalization
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
#' list of matrices will be saved in out_dir.
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
#'
#' # Define files before normalization
#' gate_dir <- file.path(dir, "Gated")
#' files_before_norm <- list.files(gate_dir,
#'                                 pattern = ".fcs",
#'                                 full.names = T)
#'
#' # Define files after normalization
#' norm_dir <- file.path(dir, "CytoNormed")
#' files_after_norm <- list.files(norm_dir,
#'                                pattern = ".fcs",
#'                                full.names = T)
#'
#' # files needs to be in the same order, check and order if needed
#' test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
#'                  basename(files_before_norm))
#'
#' batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]
#'
#' mx <- extract_pctgs_msi_per_flowsom(files_after_norm = files_after_norm,
#'                                     files_before_norm = files_before_norm,
#'                                     nCells = 50000,
#'                                     phenotyping_markers = c("CD", "HLA", "IgD"),
#'                                     functional_markers = c("MIP", "MCP", "IL",
#'                                                            "IFNa", "TNF", "TGF",
#'                                                            "Gr"),
#'                                     xdim = 10,
#'                                     ydim = 10,
#'                                     n_metaclusters = 35,
#'                                     out_dir = norm_dir,
#'                                     arcsine_transform = TRUE,
#'                                     save_matrix = TRUE,
#'                                     seed = 343)
#'
#' @return list of matrices that contain calculation for
#' cl_pctgs (cluster percentages), mcl_pctgs (metaclusters percentages),
#' cl_msi (cluster MSIs for selected markers), mcl_msi (metaclusters MSI
#' for selected markers). If save_matrix = TRUE, saves this matrices in out_dir.
#' FlowSOM objects for normalized and unnormalized data,
#' if save_flowsom_result set to TRUE.
#'
#' @export
extract_pctgs_msi_per_flowsom <- function(files_before_norm,
                                          files_after_norm,
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
                                          transform_list = NULL,
                                          save_flowsom_result = TRUE,
                                          impute_0_values = TRUE) {

  if (!all(file.exists(c(files_after_norm, files_before_norm)))){
    stop("incorrect file path, the fcs file does not exist")
  }

  file_list <- list("before" = files_before_norm,
                    "after" = files_after_norm)

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "CytoNormed")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  res <- list()
  for (f in names(file_list)){

    nCells <- length(file_list[[f]]) * 50000
    print(paste("aggregating files for", f, "normalization"))
    ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = file_list[[f]],
                                           cTotal = nCells,
                                           writeOutput = F,
                                           outputFile = file.path(out_dir, paste0(f, "_flowsom_agg.fcs")))


    if(arcsine_transform){
      ff_aggt <- flowCore::transform(ff_agg,
                                     flowCore::transformList(grep("Di", flowCore::colnames(ff_agg),
                                                                  value = TRUE),
                                                             CytoNorm::cytofTransform))
    } else if (!is.null(transform_list)){
      ff_aggt <- flowCore::transform(ff_agg, transform_list)
    } else {
      ff_aggt <- ff_agg
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

    print(paste("building FlowSOM for", f, "normalization"))
    fsom <- FlowSOM::FlowSOM(ff_aggt,
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

      ggplot2::ggsave(filename = paste0(f, "_FlowSOM.pdf"), plot = figure, device = "pdf", path = out_dir,
                      width =24, height = 10)

      saveRDS(object = fsom, file = file.path(out_dir, paste0(f, "_flowsom.RDS")))
    }


    # Define matrices for frequency (pctgs) calculation and MSI (msi). These calculation is performed
    # for clusters (cl) and metaclusters (mcl)
    cl_pctgs <- matrix(data = NA, nrow = length(file_list[[f]]),
                       ncol = xdim * ydim,
                       dimnames = list(basename(file_list[[f]]), 1:(xdim*ydim)))

    mcl_pctgs <- matrix(data = NA, nrow = length(file_list[[f]]),
                        ncol = nClus,
                        dimnames = list(basename(file_list[[f]]), 1:nClus))
    mfi_cl_names <- apply(expand.grid(paste0("Cl", seq_len(fsom$map$nNodes)),
                                      FlowSOM::GetMarkers(ff_agg,
                                                          unique(c(phenotyping_channels,functional_channels)))),
                          1, paste, collapse = "_")
    mfi_mc_names <- apply(expand.grid(paste0("MC", 1:nClus),
                                      FlowSOM::GetMarkers(ff_agg,
                                                          unique(c(phenotyping_channels,functional_channels)))),
                          1, paste, collapse = "_")
    cl_msi <- matrix(NA,
                     nrow = length(file_list[[f]]),
                     ncol = fsom$map$nNodes * length(unique(names(c(phenotyping_channels,
                                                                    functional_channels)))),
                     dimnames = list(basename(file_list[[f]]), mfi_cl_names))
    mcl_msi <- matrix(NA,
                      nrow = length(file_list[[f]]),
                      ncol =  length(mfi_mc_names),
                      dimnames = list(basename(file_list[[f]]), mfi_mc_names))

    print(paste("calculating frequency and msi for:", f, "normalization"))

    for (i in unique(fsom$data[,"File2"])){

      file <- basename(file_list[[f]][i])

      id <- which(fsom$data[,"File2"] == i)
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

    }

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

    res[[f]] <- all_mx
  }

  if(save_matrix){
    saveRDS(object = res, file = file.path(out_dir, "cell_frequency_and_msi_list_using_FlowSOM.RDS"))
  }

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
#' # Define files before normalization
#' gate_dir <- file.path(dir, "Gated")
#' files_before_norm <- list.files(gate_dir,
#'                                 pattern = ".fcs",
#'                                 full.names = T)
#'
#' # Define files after normalization
#' norm_dir <- file.path(dir, "CytoNormed")
#' files_after_norm <- list.files(norm_dir,
#'                                pattern = ".fcs",
#'                                full.names = T)
#'
#' # files needs to be in the same order, check and order if needed
#' test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
#'                  basename(files_before_norm))
#'
#' batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]
#'
#' mx <- extract_pctgs_msi_per_flowsom(files_after_norm = files_after_norm,
#'                                     files_before_norm = files_before_norm,
#'                                     nCells = 50000,
#'                                     phenotyping_markers = c("CD", "HLA", "IgD"),
#'                                     functional_markers = c("MIP", "MCP", "IL",
#'                                                            "IFNa", "TNF", "TGF",
#'                                                            "Gr"),
#'                                     xdim = 10,
#'                                     ydim = 10,
#'                                     n_metaclusters = 35,
#'                                     out_dir = norm_dir,
#'                                     arcsine_transform = TRUE,
#'                                     save_matrix = TRUE,
#'                                     seed = 343)
#' # create the list to store the plots
#'  plots <- list()
#'  for (name in names(mx[[1]])){
#'  df_plot <- prepare_data_for_plotting(frequency_msi_list = mx,
#'                                        matrix_type = name,
#'                                       n_neighbours = 11, seed = 1)
#'
#'
#'  batch <- stringr::str_match(rownames(df_plot), "day[0-9]*")[,1]
#'  samples_id <- ifelse(grepl("p1", rownames(df_plot)),"p1",
#'                       ifelse(grepl("p2", rownames(df_plot)), "p2", "ref"))
#'  stimulation <- stringr::str_match(rownames(df_plot), "UNS|RSQ|IMQ|LPS|CPG")[,1]
#'
#'  plots[[name]] <- plot_batch_using_freq_msi(df_plot = df_plot, fill = batch,
#'                                             shape = samples_id, color = batch,
#'                                             split_by_normalization = TRUE, title = name)
#'
#'}
#'
#' gg_a <- grid_arrange_common_legend(plot_lists = plots, nrow = 2, ncol = 2,
#'                                   position = "right")
#'
#'ggplot2::ggsave(filename ="batch_effect_frequency_MSI.png",
#'                device = "png",
#'               path = norm_dir,
#'                plot = gg_a,
#'                units = "cm",
#'                width = 22,
#'                height = 14, dpi = 300)
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
#' @param split_by_batch Logical, if TRUE (defult) the plots are split
#' (facet_wrap) by normalization.
#' @param fill_legend_name Character, specifies the name of the legend for fill.
#' @param color_legend_name Character, specifies the name of the legend for fill.
#' @param shape_legend_name Character, specifies the name of the legend for fill.
#' @param title Character, the title of the plot.
#' @param manual_colors Character, vector of the colors to be used,
#' the number of colors needs to be equal to the length of batch_pattern.
#'
#' @return ggplot
#' @export
#'
#' @examples
#' # Define files before normalization
#' gate_dir <- file.path(dir, "Gated")
#' files_before_norm <- list.files(gate_dir,
#'                                 pattern = ".fcs",
#'                                 full.names = T)
#'
#' # Define files after normalization
#' norm_dir <- file.path(dir, "CytoNormed")
#' files_after_norm <- list.files(norm_dir,
#'                                pattern = ".fcs",
#'                                full.names = T)
#'
#' # files needs to be in the same order, check and order if needed
#' test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
#'                  basename(files_before_norm))
#'
#' batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]
#'
#' mx <- extract_pctgs_msi_per_flowsom(files_after_norm = files_after_norm,
#'                                     files_before_norm = files_before_norm,
#'                                     nCells = 50000,
#'                                     phenotyping_markers = c("CD", "HLA", "IgD"),
#'                                     functional_markers = c("MIP", "MCP", "IL",
#'                                                            "IFNa", "TNF", "TGF",
#'                                                            "Gr"),
#'                                     xdim = 10,
#'                                     ydim = 10,
#'                                     n_metaclusters = 35,
#'                                     out_dir = norm_dir,
#'                                     arcsine_transform = TRUE,
#'                                     save_matrix = TRUE,
#'                                     seed = 343)
#' # create the list to store the plots
#'  plots <- list()
#'  for (name in names(mx[[1]])){
#'  df_plot <- prepare_data_for_plotting(frequency_msi_list = mx,
#'                                        matrix_type = name,
#'                                       n_neighbours = 11, seed = 1)
#'
#'
#'  batch <- stringr::str_match(rownames(df_plot), "day[0-9]*")[,1]
#'  samples_id <- ifelse(grepl("p1", rownames(df_plot)),"p1",
#'                       ifelse(grepl("p2", rownames(df_plot)), "p2", "ref"))
#'  stimulation <- stringr::str_match(rownames(df_plot), "UNS|RSQ|IMQ|LPS|CPG")[,1]
#'
#'  plots[[name]] <- plot_batch_using_freq_msi(df_plot = df_plot, fill = batch,
#'                                             shape = samples_id, color = batch,
#'                                             split_by_normalization = TRUE, title = name)
#'
#'}
#'
#' gg_a <- grid_arrange_common_legend(plot_lists = plots, nrow = 2, ncol = 2,
#'                                   position = "right")
#'
#'ggplot2::ggsave(filename ="batch_effect_frequency_MSI.png",
#'                device = "png",
#'               path = norm_dir,
#'                plot = gg_a,
#'                units = "cm",
#'                width = 22,
#'                height = 14, dpi = 300)
#'
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
      scale_color_manual(values = manual_colors)+
      facet_wrap(~normalization)
  }
  else {
    p <- ggplot(df_plot, aes(x = dim1, y = dim2)) +
      geom_point(data = df_plot,
                 aes(x = dim1, y = dim2, fill = fill, shape = shape,
                     color = color), size = 3) + ggtitle(title)+
      scale_color_manual(values = manual_colors)
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
#' gg_a <- grid_arrange_common_legend(plot_lists = plots, nrow = 2, ncol = 2,
#' position = "right")
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

