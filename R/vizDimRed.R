#' Visualize NMF programs on dimension reduction plots
#'
#' Plot NMF program usage on dimension reduction coordinates (UMAP, PCA, etc.)
#' with flexible color mapping and customization options.
#'
#' @param x A SingleCellExperiment object with NMF results and dimension reduction
#' @param dimred Character, name of dimension reduction to use (default "UMAP")
#' @param nmf_name Character, name of NMF result to use (default "NMF")
#' @param program Character or integer, which NMF program to plot (e.g., "NMF_1", "NMF_4", or 4)
#' @param point_size Numeric, size of points (default 0.8)
#' @param alpha Numeric, transparency of points (default 1.0)
#' @param color_scale Character, color scale to use: "viridis", "plasma", "inferno", "magma", or custom (default "viridis")
#' @param title Character, plot title (if NULL, auto-generates from program name)
#' @param dims Integer vector of length 2, which dimensions to plot (default c(1,2))
#'
#' @return A ggplot object
#' @export
#' @examples
#' # Basic usage
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 5)
#' sce <- runUMAP(sce)
#' 
#' # Plot NMF program 1 on UMAP
#' vizDimRed(sce, program = 1)
#' 
#' # Plot specific program by name
#' vizDimRed(sce, program = "NMF_4")
#' 
#' # Use PCA instead of UMAP
#' sce <- runPCA(sce)
#' vizDimRed(sce, dimred = "PCA", program = 2)
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c scale_color_gradientn
#' @importFrom ggplot2 ggtitle theme_minimal labs xlab ylab
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
vizDimRed <- function(x, dimred = "UMAP", nmf_name = "NMF", program = 1,
                      point_size = 0.8, alpha = 1.0, color_scale = "viridis",
                      title = NULL, dims = c(1, 2)) {
    
    # Input validation
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (!dimred %in% reducedDimNames(x)) {
        stop("Dimension reduction '", dimred, "' not found in reducedDims(x). ",
             "Available: ", paste(reducedDimNames(x), collapse = ", "))
    }
    
    if (!nmf_name %in% reducedDimNames(x)) {
        stop("NMF result '", nmf_name, "' not found in reducedDims(x). ",
             "Run runNMFscape() first.")
    }
    
    # Get dimension reduction coordinates
    dimred_coords <- reducedDim(x, dimred)
    if (ncol(dimred_coords) < max(dims)) {
        stop("Requested dimensions ", paste(dims, collapse = ", "), 
             " not available in ", dimred, " (has ", ncol(dimred_coords), " dimensions)")
    }
    
    # Get NMF coordinates
    nmf_coords <- reducedDim(x, nmf_name)
    
    # Handle program specification (by index or name)
    if (is.numeric(program)) {
        if (program < 1 || program > ncol(nmf_coords)) {
            stop("Program index ", program, " out of range. ",
                 "Available programs: 1 to ", ncol(nmf_coords))
        }
        program_idx <- program
        program_name <- colnames(nmf_coords)[program_idx]
        if (is.null(program_name)) {
            program_name <- paste0("Program_", program_idx)
        }
    } else {
        # Program specified by name
        if (!program %in% colnames(nmf_coords)) {
            stop("Program '", program, "' not found. ",
                 "Available programs: ", paste(colnames(nmf_coords), collapse = ", "))
        }
        program_idx <- which(colnames(nmf_coords) == program)
        program_name <- program
    }
    
    # Create data frame for plotting
    df <- data.frame(
        x = dimred_coords[, dims[1]],
        y = dimred_coords[, dims[2]],
        program_usage = nmf_coords[, program_idx]
    )
    
    # Set column names for axes
    dimred_labels <- colnames(dimred_coords)
    if (is.null(dimred_labels)) {
        dimred_labels <- paste0(dimred, c("1", "2"))
    }
    x_label <- dimred_labels[dims[1]]
    y_label <- dimred_labels[dims[2]]
    
    # Create base plot
    p <- ggplot(df, aes(x = x, y = y, color = program_usage)) +
        geom_point(size = point_size, alpha = alpha) +
        theme_minimal() +
        labs(x = x_label, y = y_label, color = "Usage")
    
    # Apply color scale
    if (color_scale == "viridis") {
        p <- p + scale_color_viridis_c(option = "viridis")
    } else if (color_scale == "plasma") {
        p <- p + scale_color_viridis_c(option = "plasma")
    } else if (color_scale == "inferno") {
        p <- p + scale_color_viridis_c(option = "inferno")
    } else if (color_scale == "magma") {
        p <- p + scale_color_viridis_c(option = "magma")
    } else {
        # Assume custom color scale - user can modify the returned plot
        p <- p + scale_color_viridis_c()
    }
    
    # Add title
    if (is.null(title)) {
        title <- paste("NMF Program:", program_name)
    }
    p <- p + ggtitle(title)
    
    return(p)
}