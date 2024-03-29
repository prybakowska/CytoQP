% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bead_normalize.R
\name{bead_normalize}
\alias{bead_normalize}
\title{Bead-based normalization}
\usage{
bead_normalize(
  files,
  cores = 1,
  markers_to_keep = NULL,
  non_mass_channel = NULL,
  beads = "dvs",
  norm_to_ref = NULL,
  remove_beads = TRUE,
  to_plot = TRUE,
  out_dir = NULL,
  k = 80,
  ...
)
}
\arguments{
\item{files}{Character vector or list with the paths of the raw files.}

\item{cores}{Number of cores to be used.Works only for not-Widows users.}

\item{markers_to_keep}{Character vector, marker names to be kept after
the normalization, can be full marker name e.g. "CD45" or "CD".
If NULL (default) all markers will be normalized and kept in flowframe.
Selection of the markers will reduce file volume and speedup the analysis.
Non-mass channels like Time, Event_length, Gaussian parameters and in addition
palladium barcoding channels are kept if non_mass_ch set to NULL.}

\item{non_mass_channel}{Character vector, non-mass channels to keep for
further analysis. Can be full channel name like Eu151Di or 151.
By default "Time" and "event_length" will be always kept in the flow frame.}

\item{beads}{Character, as in CATALYST::normCytof, "dvs"
(for bead masses 140, 151, 153 ,165, 175)
or "beta" (for bead masses 139, 141, 159, 169, 175)
or a numeric vector of masses. Default is set to "dvs".}

\item{norm_to_ref}{flow frame, created by baseline_file function to which
input data will be normalized, default is set to NULL.}

\item{remove_beads}{Logical, as in CATALYST::normCytof if beads should be
removed from fcs files. Default set to TRUE. Note, should be set to FALSE if
none of the channels is beads-specific.}

\item{to_plot}{Logical if to plot bead gate and bead normalization lines
for each file.Defaults is set to TRUE.}

\item{out_dir}{Character, pathway to where the bead normalized fcs files
and plots should be saved, for plots only if argument to_plot = TRUE,
default is set to file.path(getwd(), BeadNorm).}

\item{k}{The same as in CATALYST::normCytof, integer width of the
median window used for bead smoothing (affects visualizations only).}

\item{...}{Additional arguments to pass to normCytof.}
}
\value{
Save bead-normalized fcs files and plots to out_dir.
}
\description{
Performs bead-based normalization using beads spiked in
the sample. It is based on functions from CATALYST package.
}
\examples{
# set input directory (pathway to the files that are going to be normalized)
raw_data_dir <- file.path(dir, "RawFiles")

# set a directory where bead-normalized fcs files and plots will be saved
bead_norm_dir <- file.path(dir, "BeadNorm")

# define full pathway to the files that you want to normalize
files <- list.files(raw_data_dir,
                    pattern = ".FCS$",
                    full.names = TRUE)

# create baseline file to which all the files will be normalized
set.seed(2)
ref_sample <- baseline_file(fcs_files = files,
                            beads = "dvs",
                            out_dir = bead_norm_dir)

# Normalize files
bead_normalize(files, cores = 1,
               out_dir = bead_norm_dir,
               non_mass_channel = NULL,
               norm_to_ref = ref_sample,
               to_plot = TRUE,
               remove_beads = TRUE,
               k = 80,
               markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir",
                                 "Viability","IL", "IFNa",
                                   "TNF", "TGF", "MIP", "MCP", "Granz"))

}
