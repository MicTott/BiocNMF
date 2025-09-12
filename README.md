# NMFscape

Fast Non-negative Matrix Factorization for Single Cell and Spatial Data using RcppML

## Overview

NMFscape provides high-performance non-negative matrix factorization (NMF) methods for SingleCellExperiment and SpatialExperiment objects. The package features:

- **Fast NMF**: Uses the RcppML backend for high-performance matrix factorization
- **Consensus NMF**: Implements true consensus clustering methodology from Kotliar et al. (eLife 2019)
- **Bioconductor Integration**: Seamless workflow with SingleCellExperiment and SpatialExperiment objects
- **Rich Diagnostics**: Comprehensive stability metrics and visualization functions

## Features

### Standard NMF
- `runNMFscape()` - Fast NMF using RcppML backend
- Results stored in `reducedDims()` with basis in `metadata()`

### Consensus NMF (cNMF)
- `runConsensusNMF()` - True consensus clustering implementation
- Density filtering to remove outlier gene programs
- K-means clustering of gene expression programs across runs  
- Clustering-based stability metrics for optimal K selection
- Non-negative least squares refitting for optimal reconstruction

### Accessor Functions
- `getConsensusGEPs()` - Extract consensus gene expression programs
- `getGEPUsage()` - Get program usage per cell
- `getStabilityMetrics()` - Clustering stability metrics
- `getTopGEPFeatures()` - Top contributing genes per program

### Visualization
- `plotStability()` - Stability metrics for K selection
- `plotGEPs()` - Gene expression program heatmaps
- `plotGEPUsage()` - Program usage on reduced dimensions

## Installation

```r
# Install from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("MicTott/NMFscape")

# Load package
library(NMFscape)
```

## Quick Start

```r
library(NMFscape)
library(scuttle)

# Create example data
sce <- mockSCE(ngenes = 1000, ncells = 500)
sce <- logNormCounts(sce)

# Standard NMF
sce <- runNMFscape(sce, k = 10)
nmf_coords <- reducedDim(sce, "NMF")

# Consensus NMF for robust gene expression programs
sce <- runConsensusNMF(sce, k_range = 5:15, n_runs = 100)

# Access results
geps <- getConsensusGEPs(sce)           # Gene expression programs
usage <- getGEPUsage(sce)               # Cell usage
optimal_k <- getOptimalK(sce)           # Selected K value
top_genes <- getTopGEPFeatures(sce)     # Top genes per program

# Visualize results
plotStability(sce)                      # K selection plot
plotGEPs(sce, programs = 1:5)           # Program heatmap
```

## Methodology

### Consensus NMF Algorithm

NMFscape implements the consensus clustering approach from [Kotliar et al. (2019)](https://elifesciences.org/articles/43803):

1. **Multiple NMF runs**: Generate many gene expression programs
2. **Density filtering**: Remove outlier programs using local density
3. **Consensus clustering**: Group similar programs with k-means
4. **Stability assessment**: Clustering-based metrics for K selection
5. **Matrix refitting**: Optimal usage via non-negative least squares

This approach provides much more robust and stable gene expression programs compared to simple averaging methods.

## Requirements

- R (>= 4.3.0)
- Bioconductor packages: SingleCellExperiment, SummarizedExperiment
- RcppML for fast NMF computation
- ggplot2 for visualization (suggested)

## Citation

If you use NMFscape in your research, please cite:

- The original consensus NMF paper: Kotliar et al. (2019) "Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq" eLife
- RcppML: DeBruine et al. (2021) "High-performance non-negative matrix factorization for large single cell data"

## License

This package is licensed under the Artistic-2.0 license.

## Issues

Please report issues on [GitHub Issues](https://github.com/MicTott/NMFscape/issues).