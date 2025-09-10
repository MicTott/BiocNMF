#' Extract consensus gene expression programs (GEPs) from cNMF results
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result to extract (default "cNMF")
#' @param k Integer, which k value to extract. If NULL, uses optimal k
#'
#' @return Matrix with features x programs (consensus basis matrix)
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # geps <- getConsensusGEPs(sce)
getConsensusGEPs <- function(x, name = "cNMF", k = NULL) {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (is.null(metadata(x)$cNMF) || is.null(metadata(x)$cNMF[[name]])) {
        stop("cNMF results '", name, "' not found. Run runConsensusNMF() first.")
    }
    
    cnmf_data <- metadata(x)$cNMF[[name]]
    
    # Determine k value to use
    if (is.null(k)) {
        k <- cnmf_data$optimal_k
        if (is.null(k)) {
            stop("No optimal k found and k not specified")
        }
    }
    
    k_str <- as.character(k)
    if (!k_str %in% names(cnmf_data$all_results)) {
        stop("Results for k=", k, " not found. Available k values: ", 
             paste(names(cnmf_data$all_results), collapse = ", "))
    }
    
    return(cnmf_data$all_results[[k_str]]$consensus_w)
}

#' Extract GEP usage matrix from cNMF results  
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result to extract (default "cNMF")
#' @param k Integer, which k value to extract. If NULL, uses optimal k
#'
#' @return Matrix with cells x programs (consensus usage matrix)
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # usage <- getGEPUsage(sce)
getGEPUsage <- function(x, name = "cNMF", k = NULL) {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (is.null(metadata(x)$cNMF) || is.null(metadata(x)$cNMF[[name]])) {
        stop("cNMF results '", name, "' not found. Run runConsensusNMF() first.")
    }
    
    cnmf_data <- metadata(x)$cNMF[[name]]
    
    # Determine k value to use
    if (is.null(k)) {
        k <- cnmf_data$optimal_k
        if (is.null(k)) {
            stop("No optimal k found and k not specified")
        }
        # For optimal k, can also get from reducedDims
        if (name %in% reducedDimNames(x)) {
            return(reducedDim(x, name))
        }
    }
    
    k_str <- as.character(k)
    if (!k_str %in% names(cnmf_data$all_results)) {
        stop("Results for k=", k, " not found. Available k values: ", 
             paste(names(cnmf_data$all_results), collapse = ", "))
    }
    
    # Return transposed H matrix (cells x programs)
    return(t(cnmf_data$all_results[[k_str]]$consensus_h))
}

#' Get stability metrics from cNMF results
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result to extract (default "cNMF")
#'
#' @return Data frame with stability metrics for each k value
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # metrics <- getStabilityMetrics(sce)
getStabilityMetrics <- function(x, name = "cNMF") {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (is.null(metadata(x)$cNMF) || is.null(metadata(x)$cNMF[[name]])) {
        stop("cNMF results '", name, "' not found. Run runConsensusNMF() first.")
    }
    
    return(metadata(x)$cNMF[[name]]$stability_metrics)
}

#' Get optimal k value from cNMF results
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result (default "cNMF")
#'
#' @return Integer, optimal k value
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # k_opt <- getOptimalK(sce)
getOptimalK <- function(x, name = "cNMF") {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (is.null(metadata(x)$cNMF) || is.null(metadata(x)$cNMF[[name]])) {
        stop("cNMF results '", name, "' not found. Run runConsensusNMF() first.")
    }
    
    return(metadata(x)$cNMF[[name]]$optimal_k)
}

#' Get top features for each consensus gene expression program
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result (default "cNMF")
#' @param k Integer, which k value to use. If NULL, uses optimal k
#' @param n Integer, number of top features per program (default 20)
#'
#' @return List of character vectors, top features for each GEP
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # top_genes <- getTopGEPFeatures(sce, n = 20)
getTopGEPFeatures <- function(x, name = "cNMF", k = NULL, n = 20) {
    geps <- getConsensusGEPs(x, name = name, k = k)
    
    top_features <- apply(geps, 2, function(col) {
        idx <- order(col, decreasing = TRUE)[seq_len(min(n, length(col)))]
        rownames(geps)[idx]
    })
    
    # Convert to named list
    factor_names <- colnames(geps)
    if (is.null(factor_names)) {
        factor_names <- paste0("GEP_", seq_len(ncol(geps)))
    }
    
    if (is.matrix(top_features)) {
        result <- as.list(as.data.frame(top_features, stringsAsFactors = FALSE))
    } else {
        result <- list(top_features)
    }
    names(result) <- factor_names
    
    return(result)
}

#' Get consensus NMF information and parameters
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result (default "cNMF")
#'
#' @return List with cNMF run information (excluding large matrices)
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # info <- getConsensusNMFInfo(sce)
getConsensusNMFInfo <- function(x, name = "cNMF") {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (is.null(metadata(x)$cNMF) || is.null(metadata(x)$cNMF[[name]])) {
        stop("cNMF results '", name, "' not found. Run runConsensusNMF() first.")
    }
    
    cnmf_data <- metadata(x)$cNMF[[name]]
    
    # Create clean info without large matrices
    info <- list(
        k_range = cnmf_data$k_range,
        n_runs = cnmf_data$n_runs,
        optimal_k = cnmf_data$optimal_k,
        assay_used = cnmf_data$assay_used,
        subset_row = cnmf_data$subset_row,
        parameters = cnmf_data$parameters,
        available_k_values = names(cnmf_data$all_results),
        stability_summary = summary(cnmf_data$stability_metrics),
        call = cnmf_data$call
    )
    
    return(info)
}

#' Get available k values from cNMF results
#'
#' @param x A SingleCellExperiment object with cNMF results
#' @param name Character, name of the cNMF result (default "cNMF")
#'
#' @return Integer vector of available k values
#' @export
#' @examples
#' # sce <- runConsensusNMF(sce, k_range = 5:10)
#' # k_values <- getAvailableK(sce)
getAvailableK <- function(x, name = "cNMF") {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (is.null(metadata(x)$cNMF) || is.null(metadata(x)$cNMF[[name]])) {
        stop("cNMF results '", name, "' not found. Run runConsensusNMF() first.")
    }
    
    cnmf_data <- metadata(x)$cNMF[[name]]
    return(as.integer(names(cnmf_data$all_results)))
}