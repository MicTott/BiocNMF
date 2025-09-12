#' Run NMFscape Non-negative Matrix Factorization on SingleCellExperiment objects
#'
#' Performs NMF using RcppML on SingleCellExperiment or SpatialExperiment objects
#' and stores results in reducedDims slot with metadata
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer, number of factors for NMF (rank)
#' @param assay Character or integer, which assay to use (default "logcounts")
#' @param name Character, name for the reducedDim slot (default "NMF")
#' @param subset_row Vector specifying which features to use
#' @param tol Numeric, tolerance for convergence (default 1e-5)
#' @param maxit Integer, maximum iterations (default 100)
#' @param L1 Numeric vector of length 2, L1 regularization for [w, h] (default c(0,0))
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to RcppML::nmf
#'
#' @return The input object with NMF results stored in reducedDims(x, name)
#'   and metadata stored in metadata(x)$NMF
#'
#' @examples
#' # Single-cell example
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 10)
#' 
#' # Access NMF coefficients (cell loadings)
#' nmf_coords <- reducedDim(sce, "NMF")
#' 
#' # Access basis matrix from metadata
#' basis <- metadata(sce)$NMF_basis
#'
#' @export
#' @importFrom RcppML nmf
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom S4Vectors metadata<-
#' @importFrom Matrix Matrix
runNMFscape <- function(x, k, assay = "logcounts", name = "NMF", 
                   subset_row = NULL, tol = 1e-5, maxit = 100,
                   L1 = c(0, 0), seed = NULL, verbose = TRUE, ...) {
    
    # Ensure Matrix package is available and attached for RcppML
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required but not available")
    }
    
    # Load Matrix library if not already loaded (RcppML needs it on search path)
    if (!"package:Matrix" %in% search()) {
        library(Matrix, quietly = TRUE)
    }
    
    # Input validation
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }
    
    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }
    
    # Extract assay data
    mat <- assay(x, assay)
    
    # Subset features if requested
    if (!is.null(subset_row)) {
        mat <- mat[subset_row, , drop = FALSE]
    }
    
    # Ensure non-negative values (NMF requirement)
    if (any(mat < 0)) {
        warning("Negative values detected. Consider using log-transformed data.")
        mat[mat < 0] <- 0
    }
    
    # Set seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # Run NMF using RcppML
    if (verbose) {
        message("Running NMF with k=", k, " factors...")
    }
    
    nmf_result <- RcppML::nmf(A = mat, k = k, tol = tol, maxit = maxit,
                              L1 = L1, verbose = verbose, ...)
    
    # Extract results
    basis_matrix <- nmf_result$w  # Features x factors
    coeff_matrix <- t(nmf_result$h)  # Factors x cells -> Cells x factors
    
    # Set names
    if (!is.null(subset_row)) {
        rownames(basis_matrix) <- rownames(x)[subset_row]
    } else {
        rownames(basis_matrix) <- rownames(x)
    }
    colnames(basis_matrix) <- paste0("NMF_", seq_len(k))
    rownames(coeff_matrix) <- colnames(x)
    colnames(coeff_matrix) <- paste0("NMF_", seq_len(k))
    
    # Store coefficient matrix in reducedDims
    reducedDim(x, name) <- coeff_matrix
    
    # Store basis matrix in metadata
    basis_name <- paste0(name, "_basis")
    metadata(x)[[basis_name]] <- basis_matrix
    
    if (verbose) {
        message("NMF completed. Results stored in reducedDim(x, '", name, "')")
        message("Basis matrix stored in metadata(x)$", basis_name)
    }
    
    return(x)
}