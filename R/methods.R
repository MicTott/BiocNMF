#' Extract NMF basis matrix from SingleCellExperiment object
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param name Character, name of the NMF result to extract (default "NMF")
#'
#' @return Matrix with features x factors (basis matrix)
#' @export
#' @examples
#' # sce <- runNMF(sce, k = 10)
#' # basis <- getBasis(sce)
getBasis <- function(x, name = "NMF") {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    basis_name <- paste0(name, "_basis")
    if (!basis_name %in% names(metadata(x))) {
        stop("Basis matrix '", basis_name, "' not found in metadata. Run runNMF() first.")
    }
    
    return(metadata(x)[[basis_name]])
}

#' Extract NMF coefficient matrix from SingleCellExperiment object
#'
#' @param x A SingleCellExperiment object with NMF results  
#' @param name Character, name of the NMF result to extract (default "NMF")
#'
#' @return Matrix with cells x factors (coefficient matrix, same as reducedDim)
#' @export
#' @examples
#' # sce <- runNMF(sce, k = 10)
#' # coeffs <- getCoefficients(sce)
getCoefficients <- function(x, name = "NMF") {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (!name %in% reducedDimNames(x)) {
        stop("reducedDim '", name, "' not found. Run runNMF() first.")
    }
    
    return(reducedDim(x, name))
}

#' Calculate feature loadings (which genes contribute most to each factor)
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param name Character, name of the NMF result to use (default "NMF")
#' @param n Integer, number of top features to return per factor (default 10)
#'
#' @return List of character vectors, each containing top features for a factor
#' @export
#' @examples
#' # sce <- runNMF(sce, k = 10)
#' # top_genes <- getTopFeatures(sce, n = 20)
getTopFeatures <- function(x, name = "NMF", n = 10) {
    basis <- getBasis(x, name)
    
    top_features <- apply(basis, 2, function(col) {
        idx <- order(col, decreasing = TRUE)[seq_len(min(n, length(col)))]
        rownames(basis)[idx]
    })
    
    # Convert to named list
    factor_names <- colnames(basis)
    if (is.null(factor_names)) {
        factor_names <- paste0("Factor_", seq_len(ncol(basis)))
    }
    
    if (is.matrix(top_features)) {
        result <- as.list(as.data.frame(top_features, stringsAsFactors = FALSE))
    } else {
        result <- list(top_features)
    }
    names(result) <- factor_names
    
    return(result)
}

#' Reconstruct original matrix from NMF factors
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param name Character, name of the NMF result to use (default "NMF")
#'
#' @return Reconstructed matrix (basis %*% t(coefficients))
#' @export
#' @examples
#' # sce <- runNMF(sce, k = 10)
#' # reconstructed <- reconstructNMF(sce)
reconstructNMF <- function(x, name = "NMF") {
    basis <- getBasis(x, name)
    coeffs <- getCoefficients(x, name)
    
    return(basis %*% t(coeffs))
}