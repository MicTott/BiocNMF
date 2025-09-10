#' Plot stability metrics from consensus NMF
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result (default "cNMF")
#' @param metrics Character vector, which metrics to plot (default all)
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # plotStability(sce)
plotStability <- function(x, name = "cNMF", 
                         metrics = c("stability", "silhouette", "reproducibility")) {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package required for plotting")
    }
    
    stability_data <- getStabilityMetrics(x, name = name)
    optimal_k <- getOptimalK(x, name = name)
    
    # Reshape data for plotting
    plot_data <- data.frame()
    for (metric in metrics) {
        if (metric %in% colnames(stability_data)) {
            metric_data <- data.frame(
                k = stability_data$k,
                value = stability_data[[metric]],
                metric = metric
            )
            plot_data <- rbind(plot_data, metric_data)
        }
    }
    
    if (nrow(plot_data) == 0) {
        stop("No valid metrics found to plot")
    }
    
    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = value, color = metric)) +
        ggplot2::geom_line(size = 1) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_x_continuous(breaks = unique(plot_data$k)) +
        ggplot2::labs(
            title = "Consensus NMF Stability Metrics",
            x = "Number of factors (k)",
            y = "Metric value",
            color = "Metric"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "bottom"
        )
    
    # Add optimal k line if available
    if (!is.null(optimal_k)) {
        p <- p + ggplot2::geom_vline(xintercept = optimal_k, 
                                   linetype = "dashed", 
                                   color = "red", 
                                   alpha = 0.7) +
            ggplot2::annotate("text", x = optimal_k, y = max(plot_data$value) * 0.9,
                            label = paste("Optimal k =", optimal_k),
                            color = "red", hjust = -0.1)
    }
    
    return(p)
}

#' Plot gene expression programs (GEPs) heatmap
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result (default "cNMF")
#' @param k Integer, which k value to plot. If NULL, uses optimal k
#' @param programs Integer vector, which programs to plot (default all)
#' @param n_genes Integer, number of top genes per program to show (default 20)
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # plotGEPs(sce, programs = 1:3, n_genes = 15)
plotGEPs <- function(x, name = "cNMF", k = NULL, programs = NULL, n_genes = 20) {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package required for plotting")
    }
    
    geps <- getConsensusGEPs(x, name = name, k = k)
    
    if (is.null(programs)) {
        programs <- seq_len(ncol(geps))
    }
    
    programs <- programs[programs <= ncol(geps)]
    if (length(programs) == 0) {
        stop("No valid programs specified")
    }
    
    # Get top genes for each program
    top_genes_list <- getTopGEPFeatures(x, name = name, k = k, n = n_genes)
    
    # Create data for heatmap
    plot_data <- data.frame()
    for (i in programs) {
        program_name <- names(top_genes_list)[i]
        if (is.null(program_name)) {
            program_name <- paste0("GEP_", i)
        }
        
        top_genes <- top_genes_list[[i]]
        gene_data <- data.frame(
            gene = top_genes,
            program = program_name,
            weight = geps[top_genes, i],
            rank = seq_along(top_genes)
        )
        plot_data <- rbind(plot_data, gene_data)
    }
    
    # Create heatmap
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = program, y = reorder(gene, weight), 
                                                fill = weight)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "white", high = "red", 
                                    name = "Weight") +
        ggplot2::labs(
            title = "Top Genes per Gene Expression Program",
            x = "Gene Expression Program",
            y = "Genes"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            axis.text.y = ggplot2::element_text(size = 8),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    
    return(p)
}

#' Plot GEP usage in reduced dimension space
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result (default "cNMF")
#' @param k Integer, which k value to use. If NULL, uses optimal k
#' @param programs Integer vector, which programs to plot (default 1:4)
#' @param reduction Character, which reducedDim to use for coordinates (default "UMAP")
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # sce <- runUMAP(sce)
#' # plotGEPUsage(sce, programs = 1:4, reduction = "UMAP")
plotGEPUsage <- function(x, name = "cNMF", k = NULL, programs = 1:4, reduction = "UMAP") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package required for plotting")
    }
    
    if (!reduction %in% reducedDimNames(x)) {
        stop("Reduction '", reduction, "' not found. Available: ", 
             paste(reducedDimNames(x), collapse = ", "))
    }
    
    # Get coordinates and usage
    coords <- reducedDim(x, reduction)
    usage <- getGEPUsage(x, name = name, k = k)
    
    programs <- programs[programs <= ncol(usage)]
    if (length(programs) == 0) {
        stop("No valid programs specified")
    }
    
    # Prepare plotting data
    plot_data <- data.frame(
        x = coords[, 1],
        y = coords[, 2]
    )
    
    # Add usage for each program
    for (i in programs) {
        program_name <- colnames(usage)[i]
        if (is.null(program_name)) {
            program_name <- paste0("GEP_", i)
        }
        plot_data[[program_name]] <- usage[, i]
    }
    
    # Reshape for faceting
    coord_data <- plot_data[, c("x", "y")]
    usage_data <- plot_data[, -(1:2), drop = FALSE]
    
    plot_long <- data.frame()
    for (program_name in colnames(usage_data)) {
        program_data <- data.frame(
            x = coord_data$x,
            y = coord_data$y,
            usage = usage_data[[program_name]],
            program = program_name
        )
        plot_long <- rbind(plot_long, program_data)
    }
    
    # Create plot
    p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = x, y = y, color = usage)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.7) +
        ggplot2::scale_color_gradient(low = "lightgrey", high = "red", 
                                    name = "Usage") +
        ggplot2::facet_wrap(~ program, scales = "free") +
        ggplot2::labs(
            title = "Gene Expression Program Usage",
            x = paste0(reduction, "_1"),
            y = paste0(reduction, "_2")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            axis.ticks = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank()
        )
    
    return(p)
}

#' Plot comparison of stability metrics across different k values
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result (default "cNMF")
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:15)
#' # plotStabilityComparison(sce)
plotStabilityComparison <- function(x, name = "cNMF") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package required for plotting")
    }
    
    stability_data <- getStabilityMetrics(x, name = name)
    optimal_k <- getOptimalK(x, name = name)
    
    # Normalize metrics for comparison
    normalized_data <- stability_data
    for (col in c("stability", "silhouette", "reproducibility")) {
        if (col %in% colnames(normalized_data)) {
            col_data <- normalized_data[[col]]
            normalized_data[[col]] <- (col_data - min(col_data)) / 
                (max(col_data) - min(col_data) + 1e-10)
        }
    }
    
    # Reshape for plotting
    plot_data <- data.frame()
    for (metric in c("stability", "silhouette", "reproducibility")) {
        if (metric %in% colnames(normalized_data)) {
            metric_data <- data.frame(
                k = normalized_data$k,
                value = normalized_data[[metric]],
                metric = metric
            )
            plot_data <- rbind(plot_data, metric_data)
        }
    }
    
    # Create radar-like plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = value, fill = metric)) +
        ggplot2::geom_area(alpha = 0.6, position = "identity") +
        ggplot2::scale_x_continuous(breaks = unique(plot_data$k)) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(
            title = "Normalized Stability Metrics Comparison",
            x = "Number of factors (k)",
            y = "Normalized metric value",
            fill = "Metric"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "bottom"
        )
    
    # Add optimal k line
    if (!is.null(optimal_k)) {
        p <- p + ggplot2::geom_vline(xintercept = optimal_k, 
                                   linetype = "dashed", 
                                   color = "black", 
                                   size = 1) +
            ggplot2::annotate("text", x = optimal_k, y = 0.9,
                            label = paste("Optimal k =", optimal_k),
                            hjust = -0.1)
    }
    
    return(p)
}