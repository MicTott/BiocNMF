#' Visualize top genes per NMF program with heatmap
#'
#' Creates a heatmap showing expression of top unique genes per NMF program
#' across different cell types. Helps interpret the biological meaning of programs.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param cell_type_col Character, column name in colData containing cell type labels
#' @param nmf_name Character, name of NMF result to use (default "NMF")
#' @param n_genes Integer, number of top genes per program to include (default 5)
#' @param max_genes Integer, maximum total genes to show (default 50)
#' @param assay_name Character, which assay to use for expression values (default "logcounts")
#' @param scale_rows Logical, whether to scale rows (genes) for visualization (default TRUE)
#' @param show_programs Logical, whether to show program annotations on heatmap (default TRUE)
#'
#' @return A pheatmap object
#' @export
#' @examples
#' sce <- runBiocNMF(sce, k = 8)
#' plotProgramGenes(sce, "cell_type", n_genes = 8, max_genes = 40)
#'
#' @importFrom pheatmap pheatmap
plotProgramGenes <- function(x, cell_type_col, nmf_name = "NMF", n_genes = 5, 
                            max_genes = 50, assay_name = "logcounts",
                            scale_rows = TRUE, show_programs = TRUE) {
    
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
    
    if (!assay_name %in% assayNames(x)) {
        stop("Assay '", assay_name, "' not found in x")
    }
    
    # Check if pheatmap is available
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required but not available. Install with: install.packages('pheatmap')")
    }
    
    # Get top genes per program
    top_genes_list <- getTopFeatures(x, name = nmf_name, n = n_genes * 3)  # Get extra to ensure uniqueness
    
    # Get unique genes across all programs (prioritize by max loading)
    basis_matrix <- getBasis(x, name = nmf_name)
    all_genes <- unique(unlist(top_genes_list))
    
    # Calculate max loading for each gene across all programs
    max_loadings <- apply(basis_matrix[all_genes, , drop = FALSE], 1, max)
    
    # Select top genes by max loading, ensuring representation from each program
    selected_genes <- c()
    for (i in seq_along(top_genes_list)) {
        program_genes <- head(top_genes_list[[i]], n_genes)
        program_genes <- program_genes[!program_genes %in% selected_genes]  # Remove already selected
        selected_genes <- c(selected_genes, program_genes)
        
        if (length(selected_genes) >= max_genes) break
    }
    
    # Trim to max_genes if needed
    if (length(selected_genes) > max_genes) {
        gene_loadings <- max_loadings[selected_genes[1:max_genes]]
        selected_genes <- names(sort(gene_loadings, decreasing = TRUE))
    }
    
    # Calculate mean expression per cell type
    cell_types <- colData(x)[[cell_type_col]]
    expr_matrix <- assay(x, assay_name)[selected_genes, , drop = FALSE]
    
    # Aggregate by cell type
    cell_type_means <- t(sapply(selected_genes, function(gene) {
        tapply(expr_matrix[gene, ], cell_types, mean, na.rm = TRUE)
    }))
    
    # Handle any missing values
    cell_type_means[is.na(cell_type_means)] <- 0
    
    # Create program annotations for genes
    gene_annotations <- data.frame(
        Program = factor(rep(NA, length(selected_genes))),
        row.names = selected_genes,
        stringsAsFactors = FALSE
    )
    
    # Assign each gene to its top program
    for (gene in selected_genes) {
        gene_loadings <- basis_matrix[gene, ]
        top_program <- names(gene_loadings)[which.max(gene_loadings)]
        gene_annotations[gene, "Program"] <- top_program
    }
    
    # Create program colors
    n_programs <- ncol(basis_matrix)
    program_colors <- rainbow(n_programs, alpha = 0.7)
    names(program_colors) <- colnames(basis_matrix)
    
    annotation_colors <- list(Program = program_colors)
    
    # Create the heatmap
    if (show_programs) {
        pheat <- pheatmap::pheatmap(
            cell_type_means,
            scale = if (scale_rows) "row" else "none",
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            show_colnames = TRUE,
            annotation_row = gene_annotations,
            annotation_colors = annotation_colors,
            fontsize_row = max(8, min(12, 100/length(selected_genes))),  # Adjust font size based on number of genes
            fontsize_col = 10,
            main = paste("Top Genes per NMF Program Across Cell Types"),
            color = colorRampPalette(c("blue", "white", "red"))(100)
        )
    } else {
        pheat <- pheatmap::pheatmap(
            cell_type_means,
            scale = if (scale_rows) "row" else "none",
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            show_colnames = TRUE,
            fontsize_row = max(8, min(12, 100/length(selected_genes))),
            fontsize_col = 10,
            main = paste("Top Genes per NMF Program Across Cell Types"),
            color = colorRampPalette(c("blue", "white", "red"))(100)
        )
    }
    
    return(pheat)
}


#' Prepare gene sets for enrichment analysis from NMF programs
#'
#' Extracts top genes from each NMF program and prepares them for gene set
#' enrichment analysis. Returns both gene lists and a summary table.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param nmf_name Character, name of NMF result to use (default "NMF")
#' @param n_genes Integer, number of top genes per program (default 50)
#' @param min_loading Numeric, minimum loading value to include gene (default 0.01)
#'
#' @return List with gene_sets (named list of gene vectors) and summary table
#' @export
#' @examples
#' sce <- runBiocNMF(sce, k = 8)
#' gene_sets <- prepareProgramGeneSets(sce, n_genes = 100)
#' # Use with clusterProfiler, fgsea, etc.
#'
prepareProgramGeneSets <- function(x, nmf_name = "NMF", n_genes = 50, min_loading = 0.01) {
    
    # Input validation
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (!nmf_name %in% reducedDimNames(x)) {
        stop("NMF result '", nmf_name, "' not found. Run runBiocNMF() first.")
    }
    
    # Get basis matrix and top genes
    basis_matrix <- getBasis(x, name = nmf_name)
    top_genes_list <- getTopFeatures(x, name = nmf_name, n = n_genes)
    
    # Filter genes by minimum loading
    filtered_gene_sets <- list()
    summary_info <- data.frame(
        Program = character(),
        N_Genes = integer(),
        Top_Gene = character(),
        Top_Loading = numeric(),
        stringsAsFactors = FALSE
    )
    
    for (i in seq_along(top_genes_list)) {
        program_name <- names(top_genes_list)[i]
        genes <- top_genes_list[[i]]
        
        # Filter by minimum loading
        loadings <- basis_matrix[genes, i]
        keep_genes <- genes[loadings >= min_loading]
        
        if (length(keep_genes) > 0) {
            filtered_gene_sets[[program_name]] <- keep_genes
            
            # Add summary info
            top_gene <- keep_genes[1]
            top_loading <- loadings[top_gene]
            
            summary_info <- rbind(summary_info, data.frame(
                Program = program_name,
                N_Genes = length(keep_genes),
                Top_Gene = top_gene,
                Top_Loading = round(top_loading, 4),
                stringsAsFactors = FALSE
            ))
        }
    }
    
    return(list(
        gene_sets = filtered_gene_sets,
        summary = summary_info
    ))
}


#' Run Gene Ontology enrichment on NMF program gene sets
#'
#' Wrapper function for running GO enrichment analysis on top genes from each
#' NMF program using clusterProfiler (if available).
#'
#' @param x A SingleCellExperiment object with NMF results  
#' @param nmf_name Character, name of NMF result to use (default "NMF")
#' @param n_genes Integer, number of top genes per program (default 50)
#' @param organism Character, organism for GO analysis: "human", "mouse" (default "human")
#' @param ont Character, GO ontology: "BP", "MF", "CC", "ALL" (default "BP")
#' @param pvalue_cutoff Numeric, p-value cutoff for significance (default 0.05)
#' @param qvalue_cutoff Numeric, q-value cutoff for significance (default 0.1)
#'
#' @return List of clusterProfiler enrichResult objects (one per program)
#' @export
#' @examples  
#' \dontrun{
#' sce <- runBiocNMF(sce, k = 8)
#' go_results <- runProgramGOEnrichment(sce, organism = "mouse")
#' # View results for first program
#' head(go_results[[1]])
#' }
#'
runProgramGOEnrichment <- function(x, nmf_name = "NMF", n_genes = 50,
                                  organism = "human", ont = "BP", 
                                  pvalue_cutoff = 0.05, qvalue_cutoff = 0.1) {
    
    # Check if clusterProfiler is available
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required but not available.\n",
             "Install with: BiocManager::install('clusterProfiler')")
    }
    
    # Set organism database
    if (organism == "human") {
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            stop("Package 'org.Hs.eg.db' is required for human analysis.\n",
                 "Install with: BiocManager::install('org.Hs.eg.db')")
        }
        orgdb <- "org.Hs.eg.db"
    } else if (organism == "mouse") {
        if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
            stop("Package 'org.Mm.eg.db' is required for mouse analysis.\n",
                 "Install with: BiocManager::install('org.Mm.eg.db')")
        }
        orgdb <- "org.Mm.eg.db"
    } else {
        stop("Organism must be 'human' or 'mouse'")
    }
    
    # Prepare gene sets
    gene_sets_data <- prepareProgramGeneSets(x, nmf_name = nmf_name, n_genes = n_genes)
    gene_sets <- gene_sets_data$gene_sets
    
    # Run GO enrichment for each program
    go_results <- list()
    
    for (program in names(gene_sets)) {
        cat("Running GO enrichment for", program, "...\n")
        
        tryCatch({
            enrichment <- clusterProfiler::enrichGO(
                gene = gene_sets[[program]],
                OrgDb = orgdb,
                ont = ont,
                pAdjustMethod = "BH",
                pvalueCutoff = pvalue_cutoff,
                qvalueCutoff = qvalue_cutoff,
                readable = TRUE
            )
            
            go_results[[program]] <- enrichment
            
        }, error = function(e) {
            warning("GO enrichment failed for ", program, ": ", e$message)
            go_results[[program]] <- NULL
        })
    }
    
    # Remove failed results
    go_results <- go_results[!sapply(go_results, is.null)]
    
    if (length(go_results) == 0) {
        warning("No successful GO enrichment results")
        return(NULL)
    }
    
    cat("GO enrichment completed for", length(go_results), "programs\n")
    return(go_results)
}