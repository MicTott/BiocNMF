test_that("runConsensusNMF works with SingleCellExperiment", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    
    # Test basic functionality with small parameters for speed
    sce <- runConsensusNMF(sce, k_range = 3:5, n_runs = 5, verbose = FALSE)
    
    # Check that results are stored correctly
    expect_true("cNMF" %in% names(metadata(sce)$cNMF))
    expect_true("cNMF" %in% reducedDimNames(sce))
    
    # Check metadata structure
    cnmf_data <- metadata(sce)$cNMF$cNMF
    expect_true(is.list(cnmf_data))
    expect_true(all(c("k_range", "n_runs", "optimal_k", "stability_metrics", 
                      "all_results") %in% names(cnmf_data)))
    
    # Check stability metrics
    stability_metrics <- cnmf_data$stability_metrics
    expect_s3_class(stability_metrics, "data.frame")
    expect_true(all(c("k", "stability", "silhouette", "reproducibility", "cophenetic_correlation") %in% 
                   colnames(stability_metrics)))
    expect_equal(nrow(stability_metrics), length(3:5))
    
    # Check optimal k selection
    optimal_k <- cnmf_data$optimal_k
    expect_true(is.numeric(optimal_k))
    expect_true(optimal_k %in% 3:5)
    
    # Check reducedDim dimensions
    cnmf_coords <- reducedDim(sce, "cNMF")
    expect_equal(nrow(cnmf_coords), ncol(sce))
    expect_equal(ncol(cnmf_coords), optimal_k)
})

test_that("cNMF accessor functions work correctly", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce <- runConsensusNMF(sce, k_range = 3:5, n_runs = 5, verbose = FALSE)
    
    # Test getConsensusGEPs
    geps <- getConsensusGEPs(sce)
    expect_true(is.matrix(geps))
    expect_equal(nrow(geps), nrow(sce))
    expect_equal(ncol(geps), getOptimalK(sce))
    
    # Test getGEPUsage
    usage <- getGEPUsage(sce)
    expect_true(is.matrix(usage))
    expect_equal(nrow(usage), ncol(sce))
    expect_equal(ncol(usage), getOptimalK(sce))
    
    # Test getStabilityMetrics
    metrics <- getStabilityMetrics(sce)
    expect_s3_class(metrics, "data.frame")
    expect_equal(nrow(metrics), 3)
    
    # Test getOptimalK
    k_opt <- getOptimalK(sce)
    expect_true(is.numeric(k_opt))
    expect_true(k_opt %in% 3:5)
    
    # Test getTopGEPFeatures
    top_features <- getTopGEPFeatures(sce, n = 5)
    expect_type(top_features, "list")
    expect_length(top_features, k_opt)
    expect_true(all(sapply(top_features, length) == 5))
    
    # Test getConsensusNMFInfo
    info <- getConsensusNMFInfo(sce)
    expect_type(info, "list")
    expect_true(all(c("k_range", "optimal_k", "available_k_values") %in% names(info)))
    
    # Test getAvailableK
    available_k <- getAvailableK(sce)
    expect_true(is.integer(available_k))
    expect_true(all(available_k %in% 3:5))
})

test_that("cNMF works with different k values", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    sce <- mockSCE(ngenes = 50, ncells = 30)
    sce <- logNormCounts(sce)
    sce <- runConsensusNMF(sce, k_range = 2:4, n_runs = 3, verbose = FALSE)
    
    # Test accessing different k values
    geps_k2 <- getConsensusGEPs(sce, k = 2)
    geps_k3 <- getConsensusGEPs(sce, k = 3)
    geps_k4 <- getConsensusGEPs(sce, k = 4)
    
    expect_equal(ncol(geps_k2), 2)
    expect_equal(ncol(geps_k3), 3)
    expect_equal(ncol(geps_k4), 4)
    
    usage_k2 <- getGEPUsage(sce, k = 2)
    usage_k3 <- getGEPUsage(sce, k = 3)
    usage_k4 <- getGEPUsage(sce, k = 4)
    
    expect_equal(ncol(usage_k2), 2)
    expect_equal(ncol(usage_k3), 3)
    expect_equal(ncol(usage_k4), 4)
})

test_that("cNMF handles edge cases properly", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    sce <- mockSCE(ngenes = 50, ncells = 30)
    sce <- logNormCounts(sce)
    
    # Test error for invalid k range
    expect_error(runConsensusNMF(sce, k_range = 1:2, n_runs = 3),
                 "k_range must contain values >= 2")
    
    # Test error for insufficient runs
    expect_error(runConsensusNMF(sce, k_range = 3:4, n_runs = 1),
                 "n_runs must be >= 2")
    
    # Test accessor functions without cNMF results
    sce_empty <- mockSCE(ngenes = 50, ncells = 30)
    expect_error(getConsensusGEPs(sce_empty), "cNMF results .* not found")
    expect_error(getGEPUsage(sce_empty), "cNMF results .* not found")
    expect_error(getStabilityMetrics(sce_empty), "cNMF results .* not found")
})

test_that("cNMF parallel processing works", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("parallel")
    
    library(scuttle)
    sce <- mockSCE(ngenes = 50, ncells = 30)
    sce <- logNormCounts(sce)
    
    # Test with n_cores = 2 (if available)
    n_cores <- min(2, parallel::detectCores() - 1)
    if (n_cores >= 2) {
        sce <- runConsensusNMF(sce, k_range = 2:3, n_runs = 4, 
                              n_cores = n_cores, verbose = FALSE)
        
        # Should work the same as single core
        expect_true("cNMF" %in% reducedDimNames(sce))
        expect_true(is.matrix(getConsensusGEPs(sce)))
        expect_true(is.matrix(getGEPUsage(sce)))
    }
})