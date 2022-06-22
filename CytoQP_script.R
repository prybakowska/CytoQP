
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("prybakowska/CytoQP")
library(CytoQP)

# ------------------------------------------------------------------------------
# Bead normalization -----------------------------------------------------------
#-------------------------------------------------------------------------------

# define full pathway to the files that you want to normalize
files_raw <- list.files(path = getwd(),
                        pattern = ".FCS",
                        full.names = TRUE)

# create baseline file to which all the files will be normalized
set.seed(2)
ref_sample <- baseline_file(fcs_files = files_raw,
                            beads = "dvs")

# Normalize files
files_beadnorm <- bead_normalize(files = files_raw,
                                 cores = 1,
                                 non_mass_channel = NULL,
                                 norm_to_ref = ref_sample,
                                 to_plot = TRUE,
                                 remove_beads = TRUE,
                                 k = 80,
                                 markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir",
                                                     "Viability","IL", "IFNa",
                                                     "TNF", "TGF", "MIP", "MCP", "Granz"))

# ------------------------------------------------------------------------------
# Visualized files after bead normalization  -----------------------------------
#-------------------------------------------------------------------------------

# Define batch id and sample id for each file
batch_pattern <- stringr::str_match(basename(files_raw),
                                    "(?i).*(day[0-9]*).*.FCS")[,2]

plot_marker_quantiles(files_after_norm = files_beadnorm,
                      files_before_norm = files_raw,
                      batch_pattern = batch_pattern,
                      plot_name = "Marker_distribution_across_aliquots.pdf",
                      arcsine_transform = TRUE,
                      remove_beads = TRUE,
                      bead_channel = "140",
                      uncommon_prefix = "_beadNorm.fcs|.FCS",
                      markers_to_plot = c("CD", "HLA", "IgD", "IL",
                                          "TNF","TGF", "GR", "IFNa"),
                      manual_colors = c("darkorchid4", "darkorange",
                                        "darkgreen"))

# ------------------------------------------------------------------------------
# Signal Cleaning --------------------------------------------------------------
#-------------------------------------------------------------------------------

files_clean <- clean_files(files = files_beadnorm,
                           cores = 1,
                           to_plot = "All",
                           data_type = "MC",
                           Segment = 1000,
                           arcsine_transform = TRUE,
                           non_used_bead_ch = "140")

# ------------------------------------------------------------------------------
# File outliers detection -----------------------------------------------------
#-------------------------------------------------------------------------------

# Define batch_id for each file
file_batch_id <- stringr::str_match(basename(unlist(files_clean)),
                                    "(day[0-9]*).*.fcs")[,2]

file_scores <- file_quality_check(fcs_files = files_clean,
                                  file_batch_id = file_batch_id,
                                  phenotyping_markers = c("Ir","CD",
                                                          "HLA", "IgD"),
                                  arcsine_transform = TRUE,
                                  nClus = 10,
                                  sd = 3, out_dir = outdir)

# ------------------------------------------------------------------------------
# Files debarcoding ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Define file batch ID for each file
file_batch_id <- stringr::str_match(basename(unlist(files_clean)),
                                    "(day[0-9]*).*.fcs")[,2]

# Load metadata or read in
md <- meta_data

# read in or create barcode key
sample_key <- CATALYST::sample_key

# Extract information about barcodes used in each batch
barcodes_list <- list()
for (batch in unique(file_batch_id)){
  idx <- md[md[,"BATCH"] == batch, "BARCODE"]
  barcodes_list[[batch]] <- rownames(sample_key)[idx]
}

# Debarcode files
files_debarcoded <- debarcode_files(fcs_files = files_clean,
                                    file_score = file_scores,
                                    min_threshold = TRUE,
                                    barcodes_used = barcodes_list,
                                    file_batch_id = file_batch_id,
                                    less_than_th = TRUE,
                                    barcode_key = sample_key)

# ------------------------------------------------------------------------------
# Files aggregation and file name deconvolution --------------------------------
# ------------------------------------------------------------------------------

# Assign barcodes names
md$barcode_name <- paste0(rownames(CATALYST::sample_key)[md$BARCODE])

# Assign new sample names specifying patient id and its batch name
md$fcs_new_name <- paste0(md$ID, "_", md$STIM, "_", md$BATCH, ".fcs")

# Aggregate and deconvolute file names
files_agg <- aggregate_files(fcs_files = files_debarcoded,
                             md,
                             barcode_column = "barcode_name",
                             batch_column = "BATCH",
                             cores = 1,
                             write_agg_file = TRUE)

# ------------------------------------------------------------------------------
# Files gating -----------------------------------------------------------------
#-------------------------------------------------------------------------------

# Create directory to store plot
gate_dir <- file.path(getwd(), "Gated")
if(!dir.exists(gate_dir)){dir.create(gate_dir)}

# Gate the files and plot the gating strategy for each file
n_plots <- 3
png(file.path(gate_dir, paste0("gating.png")),
    width = n_plots * 300, height = length(files_agg) * 300)
layout(matrix(1:(length(files_agg) * n_plots),
              ncol = n_plots, byrow = TRUE))

for (file in files_agg){

  print(file)

  ff <- flowCore::read.FCS(filename = file,
                           transformation = FALSE)

  ff <- gate_intact_cells(flow_frame = ff,
                          file_name = basename(file),
                          save_gated_flow_frame = FALSE)

  ff <- gate_singlet_cells(flow_frame = ff,
                           channels = "Event_length",
                           file_name = basename(file),
                           save_gated_flow_frame = FALSE)

  ff <- gate_live_cells(flow_frame = ff,
                        viability_channel = "Pt195Di",
                        save_gated_flow_frame = TRUE,
                        file_name = basename(file), suffix = "_gated")
}

dev.off()

# ------------------------------------------------------------------------------
# Normalization using reference sample -----------------------------------------
#-------------------------------------------------------------------------------

# Set input directory
gate_dir <- file.path(getwd(), "Gated")

# Define files for normalization including reference sample
files_ref <- list.files(gate_dir, pattern = "*_gated.fcs$",
                        full.names = TRUE, recursive = TRUE)

df <- data.frame("file_paths" = files_ref,
                 "batch_labels" = stringr::str_match(
                   files_ref, "day[0-9]*")[,1],
                 "ref_ids" = grepl("REF", files_ref))

# Build normalization model
model <- train_REF_model(df = df,
                         markers_to_normalize =
                           c("CD", "HLA", "IgD","IL", "TN",
                             "MCP", "MIP","Gran", "IFNa", "TG"),
                         arcsine_transform = TRUE,
                         nQ = 2, limit = c(0,8),
                         quantileValues = c(0.05, 0.95),
                         goal = "mean",norm_with_clustering = FALSE,
                         save_model = TRUE)

# Normalize files
normalize_REF(model = model, df = df, arcsine_transform = TRUE,
              norm_with_clustering = FALSE)

# ------------------------------------------------------------------------------
# Plot batch effect ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Define files before normalization
gate_dir <- file.path(getwd(), "Gated")
files_before_norm <- list.files(gate_dir,pattern = ".fcs", full.names = TRUE)

# Define files after normalization
norm_dir <- file.path(getwd(), "CytoNormed")
files_after_norm <- list.files(norm_dir, pattern = ".fcs",full.names = TRUE)

# files needs to be in the same order, check and order if needed
test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
                 basename(files_before_norm))

batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]

# Plot batch effect
set.seed(789)
plot_batch(files_before_norm = files_before_norm,
           files_after_norm = files_after_norm,
           batch_labels = batch_labels,
           cores = 1, out_dir = norm_dir,
           clustering_markers = c("CD", "IgD", "HLA"),
           manual_colors =
             c("darkorchid4", "darkorange", "chartreuse4"))

# Plot quantiles
batch_pattern <- "day[0-9]*"
plot_marker_quantiles(files_after_norm = files_after_norm,
                      files_before_norm = files_before_norm,
                      batch_labels = batch_labels,
                      arcsine_transform = TRUE,
                      markers_to_plot =
                        c("CD", "HLA", "IgD", "IL", "TNF",
                          "TGF", "GR", "IFNa", "MCP", "MIP"),
                      manual_colors =
                        c("darkorchid4", "darkorange", "darkgreen"),
                      out_dir = norm_dir,
                      plot_name = "Marker_distribution_across_batches.pdf")

# Extract cell frequency and MSI
files_list <- list("after" = files_after_norm, "before" = files_before_norm)

results <- list()
for (name in names(files_list)){

  files <- files_list[[name]]
  mx <- extract_pctgs_msi_per_flowsom(files = files,
                                      nCells = 50000,
                                      phenotyping_markers =
                                        c("CD", "HLA", "IgD"),
                                      functional_markers =
                                        c("MIP", "MCP", "IL",
                                          "IFNa", "TNF", "TGF",
                                          "Gr"),
                                      xdim = 10, ydim = 10, n_metaclusters = 35,
                                      out_dir = norm_dir,
                                      arcsine_transform = TRUE,
                                      save_matrix = TRUE, file_name = name,
                                      seed = 343, impute_0_values = TRUE)
  results[[name]] <- mx
}

# create the list to store the plots
plots <- list()
for (name in names(results[[1]])){
  df_plot <- prepare_data_for_plotting(frequency_msi_list = results,
                                       matrix_type = name,
                                       n_neighbours = 11, seed = 35)

  batch <- stringr::str_match(rownames(df_plot), "day[0-9]*")[,1]
  samples_id <- ifelse(grepl("REF", rownames(df_plot)),"REF",
                       ifelse(grepl("p2", rownames(df_plot)), "p2", "p1"))

  plots[[name]] <- plot_batch_using_freq_msi(df_plot = df_plot, fill = batch,
                                             shape = samples_id, color = batch,
                                             split_by_normalization  = TRUE,
                                             title = name,
                                             manual_colors =
                                               c("darkorchid4", "darkorange", "darkgreen"))

}

gg_a <- grid_arrange_common_legend(plot_lists = plots, nrow = 2, ncol = 2,
                                   position = "right")

ggplot2::ggsave(filename ="batch_effect_frequency_MSI.png",
                device = "png",path = norm_dir,plot = gg_a,
                units = "cm",width = 22, height = 14, dpi = 300)


