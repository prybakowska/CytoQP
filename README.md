# CytoQP
Cytof Quality Pipeline

## Introduction
The CytoQP package provide a workflow for data preprocessing. The preprocessing tools include: bead-based normalization; 
flow rate and signal cleaning; files quality check; debracoding; files deconvolution and aggregation; gating of singlets, 
intact and viable cells; normalization with the reference sample, and feature extraction using FlowSOM algorithm. Additionally 
multiple visualization functions are provided. 
This pipeline is designed for blood samples with multibatch and large scale studies. It is best suited for samples acquired in aliquots 
(to avoid long water or CAS exposure) in each batch.

## Installation
You can install this package using the devtools library.

devtools::install_github(prybakowska/CytoQP)
