#' Calculate consensus matrices using clustering approach (true cNMF)
#'
#' This implements the proper consensus clustering approach from Kotliar et al.
#' where gene expression programs from multiple runs are clustered rather than averaged.
#'
#' @param nmf_runs List of NMF results from RcppML::nmf
#' @param k Integer, number of factors
#' @param feature_names Character vector of feature names
#' @param sample_names Character vector of sample names
#' @param density_threshold Numeric, threshold for density filtering (default 0.5)
#' @param max_iter Integer, maximum iterations for k-means (default 1000)
#'
#' @return List containing consensus matrices and clustering information
#' @keywords internal
calculateConsensusMatrices <- function(nmf_runs, k, feature_names, sample_names,
                                      density_threshold = 0.5, max_iter = 1000) {
    n_runs <- length(nmf_runs)
    n_features <- length(feature_names)
    n_samples <- length(sample_names)
    
    if (n_runs < 2) {
        stop("Need at least 2 successful runs for consensus clustering")
    }
    
    # Step 1: Collect all gene expression programs (spectra) from all runs
    all_spectra <- list()
    all_usage <- list()
    run_ids <- c()
    
    for (i in seq_along(nmf_runs)) {
        run <- nmf_runs[[i]]
        if (!is.null(run) && !is.null(run$w) && !is.null(run$h)) {
            # Normalize spectra to probability distributions (sum to 1)
            w_normalized <- apply(run$w, 2, function(x) x / (sum(x) + 1e-10))
            h_normalized <- run$h
            
            for (j in seq_len(k)) {
                all_spectra[[length(all_spectra) + 1]] <- w_normalized[, j]
                all_usage[[length(all_usage) + 1]] <- h_normalized[j, ]
                run_ids <- c(run_ids, i)
            }
        }
    }
    
    if (length(all_spectra) < k) {
        stop("Insufficient spectra for consensus clustering")
    }
    
    # Convert to matrix for clustering (spectra x features)
    spectra_matrix <- do.call(rbind, all_spectra)
    usage_matrix <- do.call(rbind, all_usage)
    
    # Step 2: Apply density filtering to remove outlier spectra
    filtered_results <- applyDensityFiltering(spectra_matrix, density_threshold)
    valid_indices <- filtered_results$valid_indices
    
    if (length(valid_indices) < k) {
        warning("Density filtering removed too many spectra, using all spectra")
        valid_indices <- seq_len(nrow(spectra_matrix))
    }
    
    filtered_spectra <- spectra_matrix[valid_indices, , drop = FALSE]
    filtered_usage <- usage_matrix[valid_indices, , drop = FALSE]
    filtered_run_ids <- run_ids[valid_indices]
    
    # Step 3: Cluster filtered spectra using k-means
    if (nrow(filtered_spectra) >= k) {
        # Use k-means clustering
        set.seed(42)  # For reproducibility
        km_result <- stats::kmeans(filtered_spectra, centers = k, 
                                 iter.max = max_iter, nstart = 20)
        cluster_assignments <- km_result$cluster
        cluster_centers <- km_result$centers
    } else {
        # Fallback: use all available spectra
        cluster_assignments <- seq_len(nrow(filtered_spectra))
        cluster_centers <- filtered_spectra
    }
    
    # Step 4: Build consensus matrix showing co-clustering frequency
    consensus_matrix <- buildConsensusMatrix(valid_indices, cluster_assignments, 
                                           length(all_spectra))
    
    # Step 5: Calculate consensus gene expression programs
    consensus_spectra <- matrix(0, nrow = n_features, ncol = k)
    consensus_usage <- matrix(0, nrow = k, ncol = n_samples)
    cluster_sizes <- numeric(k)
    
    for (cluster_id in seq_len(k)) {
        cluster_members <- which(cluster_assignments == cluster_id)
        
        if (length(cluster_members) > 0) {
            # Use median instead of mean for robustness
            if (length(cluster_members) == 1) {
                consensus_spectra[, cluster_id] <- filtered_spectra[cluster_members, ]
                consensus_usage[cluster_id, ] <- filtered_usage[cluster_members, ]
            } else {
                consensus_spectra[, cluster_id] <- apply(filtered_spectra[cluster_members, , drop = FALSE], 
                                                       2, median)
                consensus_usage[cluster_id, ] <- apply(filtered_usage[cluster_members, , drop = FALSE], 
                                                     2, median)
            }
            cluster_sizes[cluster_id] <- length(cluster_members)
        } else {
            # Empty cluster - use cluster center
            consensus_spectra[, cluster_id] <- cluster_centers[cluster_id, ]
            consensus_usage[cluster_id, ] <- 0  # Will be refitted later
            cluster_sizes[cluster_id] <- 0
        }
    }
    
    # Step 6: Reorder by total contribution (largest first)
    total_contributions <- colSums(consensus_spectra)
    reorder_idx <- order(total_contributions, decreasing = TRUE)
    
    consensus_spectra <- consensus_spectra[, reorder_idx, drop = FALSE]
    consensus_usage <- consensus_usage[reorder_idx, , drop = FALSE]
    cluster_sizes <- cluster_sizes[reorder_idx]
    
    # Set names
    rownames(consensus_spectra) <- feature_names
    colnames(consensus_spectra) <- paste0("cNMF_", seq_len(k))
    rownames(consensus_usage) <- paste0("cNMF_", seq_len(k))
    colnames(consensus_usage) <- sample_names
    
    return(list(
        consensus_w = consensus_spectra,
        consensus_h = consensus_usage,
        consensus_matrix = consensus_matrix,
        cluster_assignments = cluster_assignments,
        cluster_sizes = cluster_sizes,
        n_filtered = length(valid_indices),
        n_total_spectra = length(all_spectra),
        density_threshold = density_threshold
    ))
}

#' Calculate stability metrics for consensus NMF using clustering approach
#'
#' @param nmf_runs List of NMF results
#' @param consensus_result List with consensus matrices and clustering info
#'
#' @return List with stability metrics
#' @keywords internal
calculateStabilityMetrics <- function(nmf_runs, consensus_result) {
    n_runs <- length(nmf_runs)
    k <- ncol(consensus_result$consensus_w)
    
    if (n_runs < 2) {
        return(list(stability = 0, silhouette = 0, reproducibility = 0))
    }
    
    # Extract clustering information
    consensus_matrix <- consensus_result$consensus_matrix
    n_total_spectra <- consensus_result$n_total_spectra
    n_filtered <- consensus_result$n_filtered
    cluster_sizes <- consensus_result$cluster_sizes
    
    # Calculate clustering stability using consensus matrix
    # Stability = mean co-clustering probability
    if (!is.null(consensus_matrix) && nrow(consensus_matrix) > 1) {
        # Average of upper triangle (excluding diagonal)
        upper_tri <- upper.tri(consensus_matrix)
        if (sum(upper_tri) > 0) {
            stability <- mean(consensus_matrix[upper_tri])
        } else {
            stability <- 0
        }
    } else {
        stability <- 0
    }
    
    # Calculate silhouette-like score based on cluster quality
    if (length(cluster_sizes) > 1) {
        # Balance of cluster sizes (more balanced = better)
        size_balance <- 1 - stats::sd(cluster_sizes) / (mean(cluster_sizes) + 1e-10)
        size_balance <- max(0, min(1, size_balance))
        
        # Filter retention rate
        filter_retention <- n_filtered / n_total_spectra
        
        # Combined silhouette score
        silhouette <- 0.7 * size_balance + 0.3 * filter_retention
    } else {
        silhouette <- 0
    }
    
    # Reproducibility = fraction of spectra that survived filtering and clustering
    if (n_total_spectra > 0) {
        # How many spectra contributed to consensus vs total generated
        reproducibility <- n_filtered / n_total_spectra
        
        # Adjust based on cluster occupancy
        occupied_clusters <- sum(cluster_sizes > 0)
        cluster_occupancy <- occupied_clusters / k
        reproducibility <- reproducibility * cluster_occupancy
    } else {
        reproducibility <- 0
    }
    
    # Additional metric: Cophenetic correlation (measuring clustering quality)
    cophenetic_corr <- calculateCopheneticCorrelation(consensus_matrix)
    
    return(list(
        stability = max(0, min(1, stability)),
        silhouette = max(0, min(1, silhouette)), 
        reproducibility = max(0, min(1, reproducibility)),
        cophenetic_correlation = cophenetic_corr,
        n_filtered = n_filtered,
        n_total_spectra = n_total_spectra,
        cluster_sizes = cluster_sizes
    ))
}

#' Calculate cophenetic correlation for clustering quality assessment
#'
#' @param consensus_matrix Consensus matrix showing co-clustering frequencies
#'
#' @return Cophenetic correlation coefficient
#' @keywords internal
calculateCopheneticCorrelation <- function(consensus_matrix) {
    if (is.null(consensus_matrix) || nrow(consensus_matrix) < 3) {
        return(0)
    }
    
    # Convert consensus matrix to distance matrix
    distance_matrix <- 1 - consensus_matrix
    
    # Perform hierarchical clustering
    if (any(is.na(distance_matrix)) || any(is.infinite(distance_matrix))) {
        return(0)
    }
    
    tryCatch({
        hc <- stats::hclust(stats::as.dist(distance_matrix), method = "average")
        cophenetic_dist <- stats::cophenetic(hc)
        
        # Calculate correlation between original distances and cophenetic distances
        original_dist <- stats::as.dist(distance_matrix)
        corr <- stats::cor(as.vector(original_dist), as.vector(cophenetic_dist))
        
        return(ifelse(is.na(corr), 0, corr))
    }, error = function(e) {
        return(0)
    })
}

#' Find best matching between two factor matrices
#'
#' @param mat1 First matrix
#' @param mat2 Second matrix  
#'
#' @return List with permutation and correlation matrix
#' @keywords internal
findBestMatching <- function(mat1, mat2) {
    if (ncol(mat1) != ncol(mat2)) {
        stop("Matrices must have same number of columns")
    }
    
    k <- ncol(mat1)
    
    # Calculate correlation matrix
    cor_mat <- cor(mat1, mat2)
    
    # Find best permutation using Hungarian algorithm approximation
    # For small k, try all permutations; for large k, use greedy approach
    if (k <= 8) {
        # Try all permutations for small k
        perms <- permutations(k)
        best_score <- -Inf
        best_perm <- seq_len(k)
        
        for (i in seq_len(nrow(perms))) {
            perm <- perms[i, ]
            score <- sum(diag(cor_mat[, perm]))
            if (score > best_score) {
                best_score <- score
                best_perm <- perm
            }
        }
    } else {
        # Greedy approach for larger k
        best_perm <- greedyMatching(cor_mat)
    }
    
    return(list(
        perm = best_perm,
        correlation_matrix = cor_mat,
        score = sum(diag(cor_mat[, best_perm]))
    ))
}

#' Generate all permutations for small k
#'
#' @param k Integer, number of elements
#' @return Matrix of permutations
#' @keywords internal
permutations <- function(k) {
    if (k == 1) return(matrix(1, nrow = 1))
    if (k == 2) return(matrix(c(1,2,2,1), nrow = 2, byrow = TRUE))
    
    # For k > 2, use recursive approach (only for small k)
    if (k > 8) {
        stop("Too many permutations for k > 8")
    }
    
    # Simple implementation for small k
    factorial_k <- factorial(k)
    result <- matrix(0, nrow = factorial_k, ncol = k)
    
    # Generate using lexicographic order
    elements <- seq_len(k)
    idx <- 1
    
    # Helper function for generating permutations
    generate_perms <- function(arr, l, r) {
        if (l == r) {
            result[idx, ] <<- arr
            idx <<- idx + 1
        } else {
            for (i in l:r) {
                # Swap
                temp <- arr[l]
                arr[l] <- arr[i]
                arr[i] <- temp
                
                # Recurse
                generate_perms(arr, l + 1, r)
                
                # Backtrack
                temp <- arr[l]
                arr[l] <- arr[i]
                arr[i] <- temp
            }
        }
    }
    
    generate_perms(elements, 1, k)
    return(result)
}

#' Greedy matching for larger k
#'
#' @param cor_mat Correlation matrix
#' @return Best permutation vector
#' @keywords internal
greedyMatching <- function(cor_mat) {
    k <- ncol(cor_mat)
    used <- logical(k)
    perm <- integer(k)
    
    for (i in seq_len(k)) {
        # Find best available match for column i
        available <- which(!used)
        if (length(available) == 0) break
        
        best_idx <- available[which.max(cor_mat[i, available])]
        perm[i] <- best_idx
        used[best_idx] <- TRUE
    }
    
    return(perm)
}

#' Select optimal K from stability metrics
#'
#' @param stability_metrics Data frame with k, stability, silhouette, reproducibility
#' @param method Character, method for selection (default "combined")
#' @param verbose Logical, whether to print selection details
#'
#' @return Integer, optimal k value
#' @keywords internal
selectOptimalK <- function(stability_metrics, method = "combined", verbose = TRUE) {
    if (nrow(stability_metrics) == 0) {
        warning("No stability metrics available")
        return(NULL)
    }
    
    if (nrow(stability_metrics) == 1) {
        return(stability_metrics$k[1])
    }
    
    # Normalize metrics to 0-1 scale
    norm_stability <- (stability_metrics$stability - min(stability_metrics$stability)) / 
        (max(stability_metrics$stability) - min(stability_metrics$stability) + 1e-10)
    
    norm_silhouette <- (stability_metrics$silhouette - min(stability_metrics$silhouette)) / 
        (max(stability_metrics$silhouette) - min(stability_metrics$silhouette) + 1e-10)
    
    norm_reproducibility <- (stability_metrics$reproducibility - min(stability_metrics$reproducibility)) / 
        (max(stability_metrics$reproducibility) - min(stability_metrics$reproducibility) + 1e-10)
    
    # Combined score (weighted average)
    if (method == "combined") {
        combined_score <- 0.4 * norm_stability + 0.3 * norm_silhouette + 0.3 * norm_reproducibility
        optimal_idx <- which.max(combined_score)
    } else if (method == "stability") {
        optimal_idx <- which.max(stability_metrics$stability)
    } else if (method == "silhouette") {
        optimal_idx <- which.max(stability_metrics$silhouette)
    } else if (method == "reproducibility") {
        optimal_idx <- which.max(stability_metrics$reproducibility)
    } else {
        stop("Unknown method: ", method)
    }
    
    optimal_k <- stability_metrics$k[optimal_idx]
    
    if (verbose) {
        message("Selected optimal k=", optimal_k, " based on ", method, " method")
        message("  Stability: ", round(stability_metrics$stability[optimal_idx], 3))
        message("  Silhouette: ", round(stability_metrics$silhouette[optimal_idx], 3))
        message("  Reproducibility: ", round(stability_metrics$reproducibility[optimal_idx], 3))
    }
    
    return(optimal_k)
}

#' Apply density filtering to remove outlier spectra
#'
#' This implements local density filtering similar to the original cNMF
#'
#' @param spectra_matrix Matrix of spectra (rows = spectra, cols = features)
#' @param density_threshold Numeric, threshold for density filtering
#'
#' @return List with valid indices and density scores
#' @keywords internal
applyDensityFiltering <- function(spectra_matrix, density_threshold = 0.5) {
    n_spectra <- nrow(spectra_matrix)
    
    if (n_spectra <= 3) {
        # Too few spectra for meaningful density filtering
        return(list(
            valid_indices = seq_len(n_spectra),
            density_scores = rep(1, n_spectra)
        ))
    }
    
    # Calculate pairwise distances between spectra using cosine distance
    dist_matrix <- as.matrix(cosine_dist_matrix(spectra_matrix))
    
    # Calculate local density for each spectrum
    # Use k-nearest neighbors approach
    k_neighbors <- min(5, n_spectra - 1)
    density_scores <- numeric(n_spectra)
    
    for (i in seq_len(n_spectra)) {
        # Get distances to all other spectra
        distances <- dist_matrix[i, -i]
        
        # Find k nearest neighbors
        nearest_dists <- sort(distances)[seq_len(k_neighbors)]
        
        # Local density is inverse of mean distance to k nearest neighbors
        density_scores[i] <- 1 / (mean(nearest_dists) + 1e-10)
    }
    
    # Normalize density scores
    density_scores <- (density_scores - min(density_scores)) / 
        (max(density_scores) - min(density_scores) + 1e-10)
    
    # Filter based on threshold
    valid_indices <- which(density_scores >= density_threshold)
    
    return(list(
        valid_indices = valid_indices,
        density_scores = density_scores
    ))
}

#' Build consensus matrix showing co-clustering frequency
#'
#' @param valid_indices Indices of spectra that passed filtering
#' @param cluster_assignments Cluster assignments for filtered spectra
#' @param n_total_spectra Total number of spectra before filtering
#'
#' @return Matrix showing co-clustering frequency
#' @keywords internal
buildConsensusMatrix <- function(valid_indices, cluster_assignments, n_total_spectra) {
    n_filtered <- length(valid_indices)
    
    # Initialize consensus matrix
    consensus_matrix <- matrix(0, nrow = n_total_spectra, ncol = n_total_spectra)
    
    # Fill in co-clustering information for filtered spectra
    for (i in seq_len(n_filtered)) {
        for (j in seq_len(n_filtered)) {
            idx_i <- valid_indices[i]
            idx_j <- valid_indices[j]
            
            # Spectra are co-clustered if they have the same cluster assignment
            if (cluster_assignments[i] == cluster_assignments[j]) {
                consensus_matrix[idx_i, idx_j] <- 1
            }
        }
    }
    
    return(consensus_matrix)
}

#' Calculate cosine distance between vectors
#'
#' @param x First vector
#' @param y Second vector
#'
#' @return Cosine distance (1 - cosine similarity)
#' @keywords internal
cosine_distance <- function(x, y) {
    dot_product <- sum(x * y)
    norm_x <- sqrt(sum(x^2))
    norm_y <- sqrt(sum(y^2))
    
    if (norm_x == 0 || norm_y == 0) {
        return(1)  # Maximum distance
    }
    
    cosine_sim <- dot_product / (norm_x * norm_y)
    return(1 - cosine_sim)
}

# Custom distance function for stats::dist
cosine_dist_matrix <- function(x) {
    n <- nrow(x)
    dist_mat <- matrix(0, n, n)
    
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            if (i != j) {
                dist_mat[i, j] <- cosine_distance(x[i, ], x[j, ])
            }
        }
    }
    
    return(stats::as.dist(dist_mat))
}

#' Refit usage matrix using consensus spectra via non-negative least squares
#'
#' This step ensures the final usage matrix is optimally fitted to the original
#' data using the consensus gene expression programs
#'
#' @param data_matrix Original data matrix (features x samples)
#' @param consensus_spectra Consensus basis matrix (features x factors)
#'
#' @return Refitted usage matrix (factors x samples)
#' @keywords internal
refitUsageMatrix <- function(data_matrix, consensus_spectra) {
    n_factors <- ncol(consensus_spectra)
    n_samples <- ncol(data_matrix)
    
    # Initialize refitted usage matrix
    refitted_usage <- matrix(0, nrow = n_factors, ncol = n_samples)
    
    # Fit each sample separately using non-negative least squares
    for (i in seq_len(n_samples)) {
        sample_data <- data_matrix[, i]
        
        # Solve: sample_data ≈ consensus_spectra %*% usage[, i]
        # Subject to: usage[, i] >= 0
        usage_coeffs <- solveNNLS(consensus_spectra, sample_data)
        refitted_usage[, i] <- usage_coeffs
    }
    
    return(refitted_usage)
}

#' Solve non-negative least squares problem
#'
#' Solves min ||Ax - b||^2 subject to x >= 0
#'
#' @param A Matrix A (features x factors)
#' @param b Vector b (features)
#'
#' @return Non-negative solution vector x
#' @keywords internal
solveNNLS <- function(A, b) {
    k <- ncol(A)
    
    # Simple iterative algorithm for NNLS
    # Initialize with regular least squares solution, then project to non-negative
    tryCatch({
        # Regular least squares solution
        x <- solve(t(A) %*% A + diag(1e-10, k), t(A) %*% b)
        
        # Project to non-negative orthant
        x[x < 0] <- 0
        
        # Iterative refinement (simplified active set method)
        for (iter in seq_len(10)) {
            # Identify active set (non-zero variables)
            active <- which(x > 1e-10)
            
            if (length(active) == 0) {
                break
            }
            
            # Solve unconstrained problem on active set
            if (length(active) == 1) {
                A_active <- A[, active, drop = FALSE]
                x_new <- drop(solve(t(A_active) %*% A_active + 1e-10, t(A_active) %*% b))
            } else {
                A_active <- A[, active, drop = FALSE]
                x_new <- solve(t(A_active) %*% A_active + diag(1e-10, length(active)), 
                              t(A_active) %*% b)
            }
            
            # Update solution
            x_old <- x
            x[] <- 0
            x[active] <- pmax(0, x_new)
            
            # Check convergence
            if (sum(abs(x - x_old)) < 1e-8) {
                break
            }
        }
        
        return(as.vector(x))
        
    }, error = function(e) {
        # Fallback: return uniform positive values
        return(rep(1/k, k))
    })
}