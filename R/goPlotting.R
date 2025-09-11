# GO Enrichment Plotting Functions
# Functions to visualize GO enrichment results from runProgramGOEnrichment

#' Plot GO Enrichment Bar Chart
#' 
#' Creates a horizontal bar chart of the top enriched GO terms
#' 
#' @param go_result A data.frame from clusterProfiler enrichGO result
#' @param n_terms Number of top terms to display (default: 10)
#' @param title Plot title
#' @param color_by Column to color bars by ("pvalue", "p.adjust", "Count")
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_viridis_c labs theme_minimal theme element_text
#' @importFrom dplyr arrange head mutate
plotGOBarChart <- function(go_result, n_terms = 10, title = NULL, color_by = "p.adjust") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required for data manipulation")
  }
  
  if (nrow(go_result) == 0) {
    stop("No GO terms found in results")
  }
  
  # Prepare data
  plot_data <- go_result %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n_terms) %>%
    dplyr::mutate(
      Description = factor(Description, levels = rev(Description)),
      neg_log_pval = -log10(p.adjust)
    )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = neg_log_pval, y = Description)) +
    ggplot2::geom_col(ggplot2::aes_string(fill = color_by), alpha = 0.8) +
    ggplot2::scale_fill_viridis_c(name = gsub("\\.", " ", color_by)) +
    ggplot2::labs(
      x = "-log10(adjusted p-value)",
      y = "GO Terms",
      title = title %||% "Top Enriched GO Terms"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  return(p)
}

#' Plot GO Enrichment Dot Plot
#' 
#' Creates a dot plot showing enrichment significance and gene count
#' 
#' @param go_result A data.frame from clusterProfiler enrichGO result
#' @param n_terms Number of top terms to display (default: 15)
#' @param title Plot title
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c scale_size_continuous labs theme_minimal theme element_text
#' @importFrom dplyr arrange head mutate
#' @importFrom scales scientific
plotGODotPlot <- function(go_result, n_terms = 15, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required for data manipulation")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("scales package required for formatting")
  }
  
  if (nrow(go_result) == 0) {
    stop("No GO terms found in results")
  }
  
  # Prepare data
  plot_data <- go_result %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n_terms) %>%
    dplyr::mutate(
      Description = factor(Description, levels = rev(Description)),
      neg_log_pval = -log10(p.adjust)
    )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = neg_log_pval, y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust), alpha = 0.8) +
    ggplot2::scale_color_viridis_c(name = "Adj. p-value", trans = "log10", 
                         labels = scales::scientific) +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
    ggplot2::labs(
      x = "-log10(adjusted p-value)",
      y = "GO Terms",
      title = title %||% "GO Enrichment Overview"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  return(p)
}

#' Plot Multiple Programs GO Comparison
#' 
#' Creates a comparison plot of GO enrichment across multiple NMF programs
#' 
#' @param go_results_list Named list of GO enrichment results
#' @param n_terms Number of top terms per program (default: 5)
#' @param ncol Number of columns for facet layout
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap scale_fill_viridis_d labs theme_minimal theme element_text
#' @importFrom dplyr arrange head mutate bind_rows
#' @importFrom purrr map_dfr
plotMultipleProgramsGO <- function(go_results_list, n_terms = 5, ncol = 2) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required for data manipulation")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("purrr package required for list operations")
  }
  
  # Combine data from all programs
  combined_data <- purrr::map_dfr(names(go_results_list), function(program_name) {
    go_result <- go_results_list[[program_name]]
    
    if (nrow(go_result) == 0) {
      return(data.frame())
    }
    
    go_result %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = n_terms) %>%
      dplyr::mutate(
        Program = program_name,
        neg_log_pval = -log10(p.adjust)
      )
  })
  
  if (nrow(combined_data) == 0) {
    stop("No GO terms found in any program")
  }
  
  # Create plot
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = neg_log_pval, y = reorder(Description, neg_log_pval))) +
    ggplot2::geom_col(ggplot2::aes(fill = Program), alpha = 0.8) +
    ggplot2::facet_wrap(~Program, scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::labs(
      x = "-log10(adjusted p-value)",
      y = "GO Terms",
      title = "GO Enrichment Comparison Across NMF Programs"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  return(p)
}

#' Plot GO Network
#' 
#' Creates a network plot showing relationships between GO terms
#' 
#' @param go_result A data.frame from clusterProfiler enrichGO result
#' @param n_terms Number of top terms to include (default: 20)
#' @param title Plot title
#' @return ggplot2 object (requires enrichplot package)
#' @export
plotGONetwork <- function(go_result, n_terms = 20, title = NULL) {
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("enrichplot package required for network plots. Install with: BiocManager::install('enrichplot')")
  }
  
  if (nrow(go_result) == 0) {
    stop("No GO terms found in results")
  }
  
  # Convert to enrichResult object if needed
  if (!inherits(go_result, "enrichResult")) {
    # This is a simplified conversion - in practice you might need the full enrichResult object
    message("Note: Network plot works best with full enrichResult objects from clusterProfiler")
  }
  
  # Create network plot
  p <- enrichplot::emapplot(go_result, showCategory = n_terms) +
    ggplot2::labs(title = title %||% "GO Terms Network")
  
  return(p)
}

#' Save GO Plots to PDF
#' 
#' Convenience function to save multiple GO plots to a single PDF
#' 
#' @param go_results_list Named list of GO enrichment results
#' @param filename Output PDF filename
#' @param width PDF width in inches
#' @param height PDF height in inches
#' @export
#' @importFrom grDevices pdf dev.off
saveGOPlotsToPDF <- function(go_results_list, filename = "go_enrichment_plots.pdf", 
                            width = 12, height = 8) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package required for combining plots")
  }
  
  grDevices::pdf(filename, width = width, height = height)
  
  # Plot for each program
  for (program_name in names(go_results_list)) {
    go_result <- go_results_list[[program_name]]
    
    if (nrow(go_result) == 0) {
      cat("No GO terms for", program_name, "- skipping\n")
      next
    }
    
    cat("Creating plots for", program_name, "\n")
    
    # Create bar chart and dot plot
    p1 <- tryCatch({
      plotGOBarChart(go_result, title = paste(program_name, "- Bar Chart"))
    }, error = function(e) {
      ggplot2::ggplot() + ggplot2::ggtitle(paste("Error plotting", program_name, ":", e$message))
    })
    
    p2 <- tryCatch({
      plotGODotPlot(go_result, title = paste(program_name, "- Dot Plot"))
    }, error = function(e) {
      ggplot2::ggplot() + ggplot2::ggtitle(paste("Error plotting", program_name, ":", e$message))
    })
    
    # Combine and print
    combined_plot <- p1 / p2
    print(combined_plot)
  }
  
  # Add comparison plot if multiple programs
  if (length(go_results_list) > 1) {
    cat("Creating comparison plot\n")
    tryCatch({
      comparison_plot <- plotMultipleProgramsGO(go_results_list)
      print(comparison_plot)
    }, error = function(e) {
      cat("Error creating comparison plot:", e$message, "\n")
    })
  }
  
  grDevices::dev.off()
  cat("GO plots saved to", filename, "\n")
}