#' Run Consensus Non-negative Matrix Factorization on SingleCellExperiment objects
#'
#' Performs consensus NMF (cNMF) by running multiple NMF iterations and combining
#' results to identify stable gene expression programs. This approach increases
#' robustness and accuracy compared to single NMF runs.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k_range Integer vector, range of K values to test (default 5:15)
#' @param n_runs Integer, number of NMF runs per K value (default 100)
#' @param assay Character or integer, which assay to use (default "logcounts")
#' @param name Character, name for storing results (default "cNMF")
#' @param subset_row Vector specifying which features to use
#' @param tol Numeric, tolerance for convergence (default 1e-5)
#' @param maxit Integer, maximum iterations per run (default 100)
#' @param L1 Numeric vector of length 2, L1 regularization for [w, h] (default c(0,0))
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param n_cores Integer, number of cores for parallel processing (default 1)
#' @param ... Additional arguments passed to RcppML::nmf
#'
#' @return The input object with cNMF results stored in metadata(x)$cNMF
#'   and optimal usage matrix in reducedDims
#'
#' @examples
#' # Single-cell example
#' library(scuttle)
#' sce <- mockSCE(ngenes = 500, ncells = 200)
#' sce <- logNormCounts(sce)
#' sce <- runConsensusNMF(sce, k_range = 5:10, n_runs = 20)
#' 
#' # Access consensus results
#' geps <- getConsensusGEPs(sce)
#' usage <- getGEPUsage(sce)
#'
#' @export
#' @importFrom RcppML nmf
#' @importFrom S4Vectors metadata<-
#' @importFrom parallel mclapply
#' @importFrom Matrix Matrix
runConsensusNMF <- function(x, k_range = 5:15, n_runs = 100, 
                           assay = "logcounts", name = "cNMF",
                           subset_row = NULL, tol = 1e-5, maxit = 100,
                           L1 = c(0, 0), seed = NULL, verbose = TRUE,
                           n_cores = 1, ...) {
    
    # Ensure Matrix package is available for RcppML
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required but not available")
    }
    
    # Attach Matrix to search path if not already attached (RcppML requirement)
    if (!"package:Matrix" %in% search()) {
        attachNamespace("Matrix")
    }
    
    # Input validation
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }
    
    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }
    
    if (min(k_range) < 2) {
        stop("k_range must contain values >= 2")
    }
    
    if (n_runs < 2) {
        stop("n_runs must be >= 2 for consensus analysis")
    }
    
    # Extract and prepare data
    mat <- assay(x, assay)
    
    # Subset features if requested
    if (!is.null(subset_row)) {
        mat <- mat[subset_row, , drop = FALSE]
        feature_names <- rownames(x)[subset_row]
    } else {
        feature_names <- rownames(x)
    }
    
    # Ensure non-negative values
    if (any(mat < 0)) {
        warning("Negative values detected. Setting to zero.")
        mat[mat < 0] <- 0
    }
    
    if (verbose) {
        message("Running consensus NMF with k_range=", min(k_range), ":", max(k_range), 
                ", n_runs=", n_runs)
    }
    
    # Set up parallel processing
    if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
        use_parallel <- TRUE
    } else {
        use_parallel <- FALSE
        if (n_cores > 1) {
            warning("parallel package not available, using single core")
        }
    }
    
    # Set seed for reproducibility
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # Initialize results storage
    all_results <- list()
    stability_metrics <- data.frame(
        k = integer(),
        stability = numeric(),
        silhouette = numeric(),
        reproducibility = numeric(),
        cophenetic_correlation = numeric()
    )
    
    # Run consensus NMF for each K
    for (k in k_range) {
        if (verbose) {
            message("Processing k=", k, " (", which(k_range == k), "/", length(k_range), ")")
        }
        
        # Generate seeds for reproducible parallel runs
        run_seeds <- sample.int(n = 10000, size = n_runs)
        
        # Function to run single NMF
        run_single_nmf <- function(run_seed) {
            set.seed(run_seed)
            tryCatch({
                RcppML::nmf(A = mat, k = k, tol = tol, maxit = maxit,
                           L1 = L1, verbose = FALSE, ...)
            }, error = function(e) {
                warning("NMF run failed for k=", k, ", seed=", run_seed, ": ", e$message)
                return(NULL)
            })
        }
        
        # Run multiple NMF iterations
        if (use_parallel) {
            nmf_runs <- parallel::mclapply(run_seeds, run_single_nmf, 
                                         mc.cores = n_cores)
        } else {
            nmf_runs <- lapply(run_seeds, run_single_nmf)
        }
        
        # Remove failed runs
        nmf_runs <- nmf_runs[!sapply(nmf_runs, is.null)]
        
        if (length(nmf_runs) < n_runs * 0.5) {
            warning("More than 50% of runs failed for k=", k, ". Consider adjusting parameters.")
        }
        
        if (length(nmf_runs) == 0) {
            warning("All runs failed for k=", k, ". Skipping this k value.")
            next
        }
        
        # Calculate consensus using clustering approach and stability metrics
        consensus_result <- calculateConsensusMatrices(nmf_runs, k, feature_names, 
                                                      colnames(x), 
                                                      density_threshold = 0.5)
        stability_result <- calculateStabilityMetrics(nmf_runs, consensus_result)
        
        # Refit usage matrix using consensus spectra (non-negative least squares)
        refitted_usage <- refitUsageMatrix(mat, consensus_result$consensus_w)
        
        # Store results
        all_results[[as.character(k)]] <- list(
            k = k,
            n_successful_runs = length(nmf_runs),
            consensus_w = consensus_result$consensus_w,
            consensus_h = refitted_usage,  # Use refitted usage matrix
            consensus_matrix = consensus_result$consensus_matrix,
            cluster_assignments = consensus_result$cluster_assignments,
            cluster_sizes = consensus_result$cluster_sizes,
            n_filtered = consensus_result$n_filtered,
            n_total_spectra = consensus_result$n_total_spectra,
            stability = stability_result$stability,
            silhouette = stability_result$silhouette,
            reproducibility = stability_result$reproducibility,
            cophenetic_correlation = stability_result$cophenetic_correlation,
            raw_runs = if (n_runs <= 20) nmf_runs else NULL  # Store raw runs only for small n_runs
        )
        
        # Add to stability metrics
        stability_metrics <- rbind(stability_metrics, data.frame(
            k = k,
            stability = stability_result$stability,
            silhouette = stability_result$silhouette,
            reproducibility = stability_result$reproducibility,
            cophenetic_correlation = stability_result$cophenetic_correlation
        ))
        
        if (verbose) {
            message("  Stability: ", round(stability_result$stability, 3),
                   ", Silhouette: ", round(stability_result$silhouette, 3))
        }
    }
    
    if (length(all_results) == 0) {
        stop("No successful consensus NMF results. Check input data and parameters.")
    }
    
    # Select optimal K
    optimal_k <- selectOptimalK(stability_metrics, verbose = verbose)
    
    # Store results in metadata
    cnmf_metadata <- list(
        k_range = k_range,
        n_runs = n_runs,
        optimal_k = optimal_k,
        stability_metrics = stability_metrics,
        all_results = all_results,
        assay_used = assay,
        subset_row = subset_row,
        parameters = list(
            tol = tol,
            maxit = maxit,
            L1 = L1,
            n_cores = n_cores
        ),
        call = match.call()
    )
    
    # Store in metadata
    if (is.null(metadata(x)$cNMF)) {
        metadata(x)$cNMF <- list()
    }
    metadata(x)$cNMF[[name]] <- cnmf_metadata
    
    # Store optimal usage in reducedDims
    if (!is.null(optimal_k) && as.character(optimal_k) %in% names(all_results)) {
        optimal_usage <- t(all_results[[as.character(optimal_k)]]$consensus_h)
        rownames(optimal_usage) <- colnames(x)
        colnames(optimal_usage) <- paste0("cNMF_", seq_len(optimal_k))
        reducedDim(x, name) <- optimal_usage
    }
    
    if (verbose) {
        message("Consensus NMF completed. Optimal k=", optimal_k)
        message("Results stored in metadata(x)$cNMF$", name)
        message("Optimal usage stored in reducedDim(x, '", name, "')")
    }
    
    return(x)
}