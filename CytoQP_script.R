# # ==============================================================================
# # CytoQP V2
# # CyTOF preprocessing workflow
# # ==============================================================================
#
# library(CytoQP)


# ==============================================================================
# 0. INPUT FILES
# ==============================================================================

# Example barcoded FCS files
files_raw <- list.files(
  path = "inst/extdata",
  pattern = "^Center[1-4]_batch1\\.fcs$",
  full.names = TRUE
)

# Single-stained compensation beads
compensation_file <- file.path(
  "inst/extdata",
  "Compensation_beads.fcs"
)


# ==============================================================================
# 1. CHANNEL SELECTION
# ==============================================================================

# Read one representative FCS file
ff <- flowCore::read.FCS(
  files_raw[1],
  transformation = FALSE,
  truncate_max_range = FALSE
)

# Marker groups used for signal cleaning
markers_to_clean <- c(
  "CD",
  "Siglec",
  "TCR",
  "IgD",
  "BAFF",
  "LOX",
  "OLF",
  "HLA",
  "ICAM",
  "PD",
  "DNA",
  "Viability"
)

# Extract marker descriptions
marker_desc <- flowCore::pData(
  flowCore::parameters(ff)
)$desc

# Identify corresponding FCS channels
channels_to_clean <- FlowSOM::GetChannels(
  ff,
  grep(
    paste(markers_to_clean, collapse = "|"),
    marker_desc,
    value = TRUE
  )
)

# Compensation uses the same marker channels, excluding channels
# that should not be compensated
channels_to_compensate <- channels_to_clean[
  !grepl(
    "Ir|Pt195Di|Y89Di|Pt198Di",
    channels_to_clean
  )
]


# ==============================================================================
# 2. BASELINE REFERENCE
# ==============================================================================

# Select/create the baseline reference used for bead normalization
set.seed(2)

ref_file <- baseline_file(
  fcs_files = files_raw,
  beads = c(89, 115, 140, 159, 175, 209)
)


# ==============================================================================
# 3. BEAD NORMALIZATION
# ==============================================================================

files_beadnorm <- bead_normalize(
  files = files_raw,
  cores = 1,
  non_mass_channel = NULL,
  norm_to_ref = ref_file,
  to_plot = TRUE,
  remove_beads = TRUE,
  k = 80,
  beads = c(89, 115, 140, 159, 175, 209)
)


# ==============================================================================
# 4. SPILLOVER MATRIX
# ==============================================================================

# Calculate the spillover matrix once from the compensation beads
spillover_matrix <- calculate_spillover(
  input = compensation_file,
  channels = channels_to_compensate
)

# Visualize the calculated spillover matrix
plot_spillover(
  input = spillover_matrix
)


# ==============================================================================
# 5. COMPENSATION
# ==============================================================================

# Apply the same spillover matrix to all bead-normalized FCS files
files_compensated <- compensate_cytof(
  samples = unlist(files_beadnorm),
  compensation = spillover_matrix,
  output_dir = "Compensated_files"
)


# ==============================================================================
# 6. SIGNAL CLEANING
# ==============================================================================

# Remove acquisition anomalies using PeacoQC
files_clean <- clean_files(
  files = files_compensated,
  cores = 1,
  to_plot = "All",
  data_type = "MC",
  Segment = 1000,
  arcsine_transform = TRUE,
  non_used_bead_ch = "140",
  channels_to_clean = channels_to_clean
)


# ==============================================================================
# 7. DEBARCODING
# ==============================================================================

# CATALYST palladium barcode key
sample_key <- CATALYST::sample_key

# Example dataset contains three barcoded samples
barcode_idx <- c(20, 11, 2)

barcodes_used <- rownames(sample_key)[
  barcode_idx
]

# All example files belong to batch1
file_batch_id <- rep(
  "batch1",
  length(files_clean)
)

# Same three barcodes are used in all example pooled FCS files
barcodes_list <- list(
  batch1 = barcodes_used
)

# Debarcode pooled FCS files
files_debarcoded <- debarcode_files(
  fcs_files = files_clean,
  min_threshold = TRUE,
  barcodes_used = barcodes_list,
  file_batch_id = file_batch_id,
  less_than_th = TRUE,
  barcode_key = sample_key
)


# ==============================================================================
# 8. GATING
# ==============================================================================

gate_dir <- "Gated"

if (!dir.exists(gate_dir)) {
  dir.create(
    gate_dir,
    recursive = TRUE
  )
}

# Three gating steps:
# intact cells -> singlets -> live cells
n_plots <- 3

grDevices::png(
  filename = file.path(
    gate_dir,
    "gating.png"
  ),
  width = n_plots * 300,
  height = length(files_debarcoded) * 300
)

graphics::layout(
  matrix(
    seq_len(
      length(files_debarcoded) * n_plots
    ),
    ncol = n_plots,
    byrow = TRUE
  )
)

for (file in files_debarcoded) {

  message(
    "Gating: ",
    basename(file)
  )

  ff <- flowCore::read.FCS(
    filename = file,
    transformation = FALSE
  )

  # Intact cells
  ff <- gate_intact_cells(
    flow_frame = ff,
    file_name = basename(file),
    save_gated_flow_frame = FALSE
  )

  # Singlets
  ff <- gate_singlet_cells(
    flow_frame = ff,
    channels = "Event_length",
    file_name = basename(file),
    save_gated_flow_frame = FALSE
  )

  # Live cells
  ff <- gate_live_cells_update(
    flow_frame = ff,
    viability_channel = "Pt195Di",
    eosinophil_channel = "In113Di",
    save_gated_flow_frame = TRUE,
    file_name = basename(file),
    suffix = "_gated"
  )
}

grDevices::dev.off()
