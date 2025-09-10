#' Visualize NMF programs on UMAP coordinates
#'
#' Convenience function to plot NMF program usage on UMAP coordinates.
#' This is a wrapper around vizDimRed() specifically for UMAP plotting.
#'
#' @param x A SingleCellExperiment object with NMF results and UMAP coordinates
#' @param nmf_name Character, name of NMF result to use (default "NMF")
#' @param program Character or integer, which NMF program to plot (e.g., "NMF_1", "NMF_4", or 4)
#' @param point_size Numeric, size of points (default 0.8)
#' @param alpha Numeric, transparency of points (default 1.0)
#' @param color_scale Character, color scale to use: "viridis", "plasma", "inferno", "magma" (default "viridis")
#' @param title Character, plot title (if NULL, auto-generates from program name)
#'
#' @return A ggplot object
#' @export
#' @examples
#' # Basic usage
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' sce <- runBiocNMF(sce, k = 5)
#' sce <- runUMAP(sce)
#' 
#' # Plot NMF program 1 on UMAP
#' vizUMAP(sce, program = 1)
#' 
#' # Plot specific program by name with custom styling
#' vizUMAP(sce, program = "NMF_4", color_scale = "plasma", point_size = 1.2)
vizUMAP <- function(x, nmf_name = "NMF", program = 1, point_size = 0.8, 
                    alpha = 1.0, color_scale = "viridis", title = NULL) {
    
    # Call vizDimRed with UMAP-specific settings
    vizDimRed(x = x, dimred = "UMAP", nmf_name = nmf_name, program = program,
              point_size = point_size, alpha = alpha, color_scale = color_scale,
              title = title, dims = c(1, 2))
}