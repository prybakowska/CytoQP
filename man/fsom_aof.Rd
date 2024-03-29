% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/files_quality_check.R
\name{fsom_aof}
\alias{fsom_aof}
\title{Prepares FlowSOM}
\usage{
fsom_aof(
  fcs_files,
  phenotyping_markers,
  markers_to_score,
  nCells = length(fcs_files) * 10000,
  xdim = 10,
  ydim = 10,
  nClus = 10,
  out_dir = NULL,
  batch = NULL,
  arcsine_transform = TRUE,
  transform_list = NULL,
  my_colors = NULL,
  seed = NULL,
  to_plot = TRUE
)
}
\arguments{
\item{fcs_files}{Character, full path to fcs_files.}

\item{phenotyping_markers}{Character vector, marker names to be used for
clustering, can be full marker name e.g. "CD45" or "CD" if all CD-markers
needs to be plotted.}

\item{markers_to_score}{Character vector, marker names to be used for
flowsom clustering including DNA marker Iridium and viability staining
if available. Can be full marker name e.g. "CD45" or pattern "CD" if
all CD-markers needs to be plotted. Default is set to NULL, thus the aof scores
will be calculated for the markers included in phenotyping_markers.}

\item{nCells}{Numeric, the total number of cells, to use for FlowSOM clustering.
This number is determined by total number of fcs files, as a default 10000 cells
is used per file}

\item{xdim}{Numeric, as in FlowSOM, parameter to pass to FlowSOM,
width of the SOM grid}

\item{ydim}{Numeric, as in FlowSOM, parameter to pass to FlowSOM,
height of the SOM grid}

\item{nClus}{Numeric, exact number of clusters for metaclustering}

\item{out_dir}{Character, pathway to where the FlowSOM clustering plot should
be saved, default is set to working directory.}

\item{batch}{Character, the name of the acquisition batch for each fcs file.
This name is passed to FlowSOM plot name, default is set to NULL}

\item{arcsine_transform}{Logical, if the data should be transformed with
arcsine transformation and cofactor 5. Default set to TRUE.}

\item{transform_list}{Transformation list to pass to the flowCore
transform function see flowCore::transformList.}

\item{my_colors}{An array specifying colors to be used for the background
coloring of metaclusters in FlowSOM and t-SNE plot. Must have a length equal
to the nClus.If NULL (default) colors defined in FlowSOM package will be used.}

\item{seed}{numeric, set to obtain reproducible results, default set to NULL.}

\item{to_plot}{Logical, if FlowSOM tree and t-SNE map should be plotted,
default set to TRUE.}
}
\value{
fsom object
}
\description{
Builds FlowSOM object
}
