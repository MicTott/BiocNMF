#' Find All Differentially Expressed Programs (1-vs-all analysis)
#'
#' Uses scran::findMarkers to identify NMF programs that are differentially
#' expressed in each cell type compared to all other cell types. Treats each
#' NMF program as a "gene" and performs 1-vs-all comparisons.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param cell_type_col Character, column name in colData containing cell type labels
#' @param nmf_name Character, name of NMF result to use (default "NMF")
#' @param test Character, statistical test to use: "wilcox" (default), "t", "binom"
#' @param pval.type Character, p-value type: "some", "any", "all" (default "some")
#' @param log.p Logical, whether to log-transform p-values (default FALSE)
#' @param min.prop Numeric, minimum proportion of cells that must express program (default 0.1)
#'
#' @return List of DataFrames (one per cell type) with differential expression statistics
#' @export
#' @examples
#' # Run NMF and find differentially expressed programs
#' sce <- runBiocNMF(sce, k = 8)
#' deps <- FindAllDEPs(sce, "cell_type")
#' 
#' # View results for first cell type
#' head(deps[[1]])
#'
#' @importFrom scran findMarkers
#' @importFrom SingleCellExperiment reducedDim colData
FindAllDEPs <- function(x, cell_type_col, nmf_name = "NMF",
                       test = "wilcox", pval.type = "some", 
                       log.p = FALSE, min.prop = 0.1) {
    
    # Input validation
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (!nmf_name %in% reducedDimNames(x)) {
        stop("NMF result '", nmf_name, "' not found. Run runBiocNMF() first.")
    }
    
    if (!cell_type_col %in% colnames(colData(x))) {
        stop("Cell type column '", cell_type_col, "' not found in colData(x)")
    }
    
    # Get NMF program usage (cells x programs)
    program_usage <- reducedDim(x, nmf_name)
    
    # Get cell type labels
    cell_types <- colData(x)[[cell_type_col]]
    
    # Check for missing cell type labels
    if (any(is.na(cell_types))) {
        warning("Removing ", sum(is.na(cell_types)), " cells with missing cell type labels")
        keep_cells <- !is.na(cell_types)
        program_usage <- program_usage[keep_cells, , drop = FALSE]
        cell_types <- cell_types[keep_cells]
    }
    
    # Transpose for findMarkers (programs x cells, like genes x cells)
    program_matrix <- t(program_usage)
    
    # Set row names to program names
    if (is.null(rownames(program_matrix))) {
        rownames(program_matrix) <- paste0("Program_", seq_len(nrow(program_matrix)))
    }
    
    # Run findMarkers
    markers <- scran::findMarkers(program_matrix, groups = cell_types,
                                  test.type = test, pval.type = pval.type,
                                  log.p = log.p, min.prop = min.prop)
    
    # Add fold change calculations to each result
    cell_type_names <- names(markers)
    for (ct in cell_type_names) {
        marker_df <- markers[[ct]]
        
        # Calculate fold enrichment for each program
        in_celltype <- cell_types == ct
        fold_enrichment <- numeric(nrow(marker_df))
        mean_usage_in <- numeric(nrow(marker_df))
        mean_usage_out <- numeric(nrow(marker_df))
        
        for (i in seq_len(nrow(marker_df))) {
            program_idx <- rownames(marker_df)[i]
            usage_in <- program_matrix[program_idx, in_celltype]
            usage_out <- program_matrix[program_idx, !in_celltype]
            
            mean_in <- mean(usage_in, na.rm = TRUE)
            mean_out <- mean(usage_out, na.rm = TRUE)
            
            mean_usage_in[i] <- mean_in
            mean_usage_out[i] <- mean_out
            
            # Avoid division by zero
            fold_enrichment[i] <- if (mean_out == 0) {
                if (mean_in == 0) 1 else Inf
            } else {
                mean_in / mean_out
            }
        }
        
        # Add fold enrichment and mean usage columns
        markers[[ct]]$fold_enrichment <- fold_enrichment
        markers[[ct]]$mean_usage_in <- mean_usage_in
        markers[[ct]]$mean_usage_out <- mean_usage_out
        
        # Reorder columns for better readability
        col_order <- c("fold_enrichment", "mean_usage_in", "mean_usage_out", 
                      setdiff(names(markers[[ct]]), c("fold_enrichment", "mean_usage_in", "mean_usage_out")))
        markers[[ct]] <- markers[[ct]][, col_order]
    }
    
    return(markers)
}


#' Visualize Differentially Expressed Programs with volcano plots
#'
#' Creates volcano plots showing fold enrichment vs. statistical significance
#' for differentially expressed NMF programs across cell types.
#'
#' @param deps_results List of DataFrames from FindAllDEPs()
#' @param cell_types Character vector, which cell types to plot (default: all)
#' @param fold_threshold Numeric, fold enrichment threshold for significance (default 1.5)
#' @param pval_threshold Numeric, p-value threshold for significance (default 0.05)
#' @param top_n Integer, number of top programs to label (default 5)
#' @param ncol Integer, number of columns for faceting (default 3)
#'
#' @return A ggplot object
#' @export
#' @examples
#' deps <- FindAllDEPs(sce, "cell_type") 
#' plotDEPsVolcano(deps)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline
#' @importFrom ggplot2 scale_color_manual facet_wrap theme_minimal labs
#' @importFrom ggplot2 theme element_text
plotDEPsVolcano <- function(deps_results, cell_types = NULL,
                           fold_threshold = 1.5, pval_threshold = 0.05,
                           top_n = 5, ncol = 3) {
    
    # Select cell types to plot
    if (is.null(cell_types)) {
        cell_types <- names(deps_results)
    }
    
    # Combine all results into one data frame - handle DataFrame structure properly
    plot_data <- do.call(rbind, lapply(cell_types, function(ct) {
        deps_df <- deps_results[[ct]]
        
        # Create clean data frame with only needed columns
        df <- data.frame(
            cell_type = ct,
            program = rownames(deps_df),
            fold_enrichment = deps_df$fold_enrichment,
            p.value = deps_df$p.value,
            FDR = deps_df$FDR,
            stringsAsFactors = FALSE
        )
        
        df$neg_log10_pval <- -log10(df$p.value)
        df$log2_fold <- log2(df$fold_enrichment)
        
        # Handle infinite values properly
        if (any(is.infinite(df$log2_fold))) {
            max_finite <- max(df$log2_fold[is.finite(df$log2_fold)], na.rm = TRUE)
            min_finite <- min(df$log2_fold[is.finite(df$log2_fold)], na.rm = TRUE)
            
            if (is.finite(max_finite)) {
                df$log2_fold[is.infinite(df$log2_fold) & df$log2_fold > 0] <- max_finite + 0.5
            } else {
                df$log2_fold[is.infinite(df$log2_fold) & df$log2_fold > 0] <- 5
            }
            
            if (is.finite(min_finite)) {
                df$log2_fold[is.infinite(df$log2_fold) & df$log2_fold < 0] <- min_finite - 0.5
            } else {
                df$log2_fold[is.infinite(df$log2_fold) & df$log2_fold < 0] <- -5
            }
        }
        
        return(df)
    }))
    
    # Create significance categories
    plot_data$significance <- "Not Significant"
    plot_data$significance[plot_data$fold_enrichment >= fold_threshold & 
                          plot_data$p.value <= pval_threshold] <- "Enriched"
    plot_data$significance[plot_data$fold_enrichment <= (1/fold_threshold) & 
                          plot_data$p.value <= pval_threshold] <- "Depleted"
    
    # Identify top programs to label for each cell type
    plot_data$label <- ""
    for (ct in cell_types) {
        ct_data <- plot_data[plot_data$cell_type == ct, ]
        ct_data <- ct_data[order(-ct_data$neg_log10_pval, -abs(ct_data$log2_fold)), ]
        top_programs <- head(ct_data$program, top_n)
        plot_data$label[plot_data$cell_type == ct & plot_data$program %in% top_programs] <- 
            plot_data$program[plot_data$cell_type == ct & plot_data$program %in% top_programs]
    }
    
    # Create volcano plot
    p <- ggplot(plot_data, aes(x = log2_fold, y = neg_log10_pval)) +
        geom_point(aes(color = significance), alpha = 0.7, size = 2) +
        geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray50", alpha = 0.8) +
        geom_vline(xintercept = log2(fold_threshold), linetype = "dashed", color = "gray50", alpha = 0.8) +
        geom_vline(xintercept = log2(1/fold_threshold), linetype = "dashed", color = "gray50", alpha = 0.8) +
        scale_color_manual(values = c("Not Significant" = "gray70", 
                                     "Enriched" = "red3", 
                                     "Depleted" = "blue3")) +
        facet_wrap(~cell_type, ncol = ncol, scales = "free") +
        theme_minimal() +
        labs(x = "Log2 Fold Enrichment", 
             y = "-Log10 P-value",
             title = "Differentially Expressed Programs Across Cell Types",
             subtitle = "Each panel shows 1-vs-all comparison for each cell type",
             color = "Significance") +
        theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(size = 11, hjust = 0.5),
              strip.text = element_text(size = 11, face = "bold"),
              axis.text = element_text(size = 9),
              legend.position = "bottom")
    
    # Add labels for top programs if ggrepel is available
    if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(aes(label = ifelse(label != "", label, "")),
                                         size = 3, max.overlaps = 10)
    }
    
    return(p)
}


#' Visualize Differentially Expressed Programs with dot plot
#'
#' Creates a dot plot showing fold enrichment and significance of differentially
#' expressed NMF programs across cell types, similar to gene expression dot plots.
#'
#' @param deps_results List of DataFrames from FindAllDEPs()
#' @param fold_threshold Numeric, minimum fold enrichment to display (default 1.2)
#' @param pval_threshold Numeric, maximum p-value to display (default 0.05)
#' @param top_programs Integer, maximum number of programs to show (default 20)
#'
#' @return A ggplot object  
#' @export
#' @examples
#' deps <- FindAllDEPs(sce, "cell_type")
#' plotDEPsDots(deps)
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous
#' @importFrom ggplot2 scale_color_viridis_c theme_minimal labs theme
#' @importFrom ggplot2 element_text
plotDEPsDots <- function(deps_results, fold_threshold = 1.2, 
                        pval_threshold = 0.05, top_programs = 20) {
    
    # Combine all results and filter for significant enrichments
    plot_data <- do.call(rbind, lapply(names(deps_results), function(ct) {
        deps_df <- deps_results[[ct]]
        
        # Create clean data frame
        df <- data.frame(
            cell_type = ct,
            program = rownames(deps_df),
            fold_enrichment = deps_df$fold_enrichment,
            p.value = deps_df$p.value,
            FDR = deps_df$FDR,
            stringsAsFactors = FALSE
        )
        
        # Filter for significant enrichments
        return(df[df$fold_enrichment >= fold_threshold & df$p.value <= pval_threshold, ])
    }))
    
    if (nrow(plot_data) == 0) {
        stop("No significant enrichments found with current thresholds. ",
             "Try lowering fold_threshold or increasing pval_threshold.")
    }
    
    # Select top programs by significance
    program_importance <- aggregate((-log10(plot_data$p.value)), 
                                   by = list(plot_data$program), 
                                   FUN = max)
    names(program_importance) <- c("program", "max_neg_log10_pval")
    program_importance <- program_importance[order(-program_importance$max_neg_log10_pval), ]
    
    if (nrow(program_importance) > top_programs) {
        keep_programs <- head(program_importance$program, top_programs)
        plot_data <- plot_data[plot_data$program %in% keep_programs, ]
    }
    
    # Add negative log10 p-value for color mapping
    plot_data$neg_log10_pval <- -log10(plot_data$p.value)
    
    # Create dot plot
    p <- ggplot(plot_data, aes(x = cell_type, y = program)) +
        geom_point(aes(size = fold_enrichment, color = neg_log10_pval)) +
        scale_size_continuous(name = "Fold\nEnrichment", range = c(1, 8),
                             breaks = c(1.5, 2, 3, 5), labels = c("1.5", "2", "3", "5+")) +
        scale_color_viridis_c(name = "-Log10\nP-value", option = "plasma") +
        theme_minimal() +
        labs(x = "Cell Type", y = "NMF Program",
             title = "Differentially Expressed Programs Across Cell Types") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 14, hjust = 0.5),
              legend.position = "right")
    
    return(p)
}