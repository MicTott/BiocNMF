test_that("runNMF works with SingleCellExperiment", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    
    # Test basic functionality
    sce <- runNMF(sce, k = 5, verbose = FALSE)
    
    # Check that results are stored correctly
    expect_true("NMF" %in% reducedDimNames(sce))
    expect_true("NMF_basis" %in% names(metadata(sce)))
    
    # Check dimensions
    expect_equal(nrow(reducedDim(sce, "NMF")), ncol(sce))
    expect_equal(ncol(reducedDim(sce, "NMF")), 5)
    expect_equal(nrow(getBasis(sce)), nrow(sce))
    expect_equal(ncol(getBasis(sce)), 5)
    
    # Test accessor functions
    coeffs <- getCoefficients(sce)
    basis <- getBasis(sce)
    top_features <- getTopFeatures(sce, n = 10)
    
    expect_true(is.matrix(coeffs))
    expect_true(is.matrix(basis))
    expect_type(top_features, "list")
    expect_length(top_features, 5)
    
    # Test reconstruction
    reconstructed <- reconstructNMF(sce)
    expect_true(is.matrix(reconstructed))
    expect_equal(dim(reconstructed), c(nrow(sce), ncol(sce)))
})