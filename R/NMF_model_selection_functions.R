# @title
# Get the reconstruction accuracy, \eqn{Q^{2}}, for a given choice of parameter
# values
#
# @param x The parameter values. Expected columns "k_vals", "alpha" and "fold".
# @param seed_val The seed to be used.
# @param verbose Default to \code{0} which will not print any messages, or can
# be set to \code{1} which will print messages.
#
# @return \eqn{Q^{2}}, a real value
# @importFrom MASS ginv
.get_q2_using_py <- function(x, seed_val, verbose = 0) {
    ##
    if (verbose == 1) {
        cat(paste("INFO-START", x["k_vals"], x["alpha"], "\n", sep = ","))
    } else {
        ## cat(".")
    }
    this_k <- as.numeric(x["k_vals"])
    this_alpha <- as.numeric(x["alpha"])
    #
    test_fold <- as.numeric(x["fold"])
    train_rows <-
        cvfolds$cvf_rows$subsets[cvfolds$cvf_rows$which != test_fold]
    train_cols <-
        cvfolds$cvf_cols$subsets[cvfolds$cvf_cols$which != test_fold]
    test_rows <-
        cvfolds$cvf_rows$subsets[cvfolds$cvf_rows$which == test_fold]
    test_cols <-
        cvfolds$cvf_cols$subsets[cvfolds$cvf_cols$which == test_fold]

    ## Split data matrix X -- separate the training and test parts
    ##        _      _
    ##   X = |  A  B  |
    ##       |_ C  D _|
    ##
    ## Reconstruct A, by performing NMF on D
    ## More details in Owen and Perry, Annals of Statistics, 2009
    ## for (itr in seq(100)) {
    new_ord <- sample(ncol(X), ncol(X), replace = FALSE)
    X <- X[, new_ord]
    ##
    submatrixD <- X[train_rows, train_cols]
    submatrixA <- X[test_rows, test_cols]
    submatrixB <- X[test_rows, train_cols]
    submatrixC <- X[train_rows, test_cols]

    ## NMF on submatrixD
    ## Setup params
    ## NMF call
    ## Using python/scikit-learn NMF

    ## new_ord <- sample(ncol(submatrixD), ncol(submatrixD))
    ## submatrixD <- submatrixD[ , new_ord]
    ##
    ## Sourcing perform_nmf_func, see:
    ## https://rstudio.github.io/reticulate/reference/py_is_null_xptr.html
    ##
    ## print(file.path(paths_list_src_path, 'perform_nmf.py'))
    ## source_python(file.path(paths_list_src_path, 'perform_nmf.py'))
    ## if(py_validate_xptr(perform_nmf_func)) {
    nmf_submatrixD <- perform_nmf_func(submatrixD,
                                        nPatterns = as.integer(this_k),
                                        nIter = as.integer(200),
                                        givenAlpha = this_alpha,
                                        givenL1_ratio = 1,
                                        seed_val = as.integer(seed_val)
                        )
    ## }else{
    ##   #Throw error
    ## }
    ##
    D_W <- nmf_submatrixD[[1]]
    D_H <- nmf_submatrixD[[2]]
    #
    reconstructed_submatrixA <-
        submatrixB %*% MASS::ginv(D_H) %*% MASS::ginv(D_W) %*% submatrixC
    ##
    if (verbose == 1) {
        cat(paste("INFO-END", x["k_vals"], x["alpha"], "\n", sep = ","))
    } else {
        ## if (itr %% 5 == 0) cat(",", itr)
        ## if (itr %% 50 == 0) cat(paste(itr, ", "))
        ## cat(",")
    }
    ##
    q2 <- .compute_q2(submatrixA, reconstructed_submatrixA)
    ## }
    ## return(mean(unlist(q2)))
    return (q2)
}
## =============================================================================




# @title Compute \eqn{Q^2} Value
#
# @description Computes the reconstruction accuracy, \eqn{Q^2}, given the
# original matrix and the recontructed matrix to check against.
#
# @param A The original matrix
# @param recA The reconstructed matrix
#
# @return Q^2, the reconstruction accuracy
.compute_q2 <- function(A, recA) {
    ## input: original matrix, reconstructed portion
    ## return: computed q2 val
    if (!is.matrix(A)) {
        stop("Original matrix is not of type matrix")
    }
    if (!is.matrix(recA)) {
        stop("Reconstructed matrix is not of type matrix")
    }
    if (sum(dim(A)) == 2 && is.na(A)) {
        stop("Empty: Original matrix")
    }
    if (sum(dim(recA)) == 2 && is.na(recA)) {
        stop("Empty: Reconstructed matrix")
    }
    second_term_num <- sum((A - recA) ^ 2)
    second_term_den <- sum(A ^ 2) # TO-DO: make sure this is not zero
    if (second_term_den != 0.0) {
        q2 <- 1 - (second_term_num / second_term_den)
    } else {
        q2 <- 0.0
    }
    return(q2)
}
## =============================================================================



# @title Perform Model Selection via Cross-Validation
#
# @description The function performs model selection via cross-validation for
# choosing the optimal values of parameters for NMF. These are K, the number of
#  factors in NMF, and \eqn{\alpha}, the sparsity coefficient.
#
# @param X The given data matrix
# @param param_ranges An object holding the range of values for parameters
# \code{k}, \code{alphaBase}, and \code{alphaPow}.
# @param kFolds Numeric The number of cross-validation folds.
# @param parallelDo Set to \code{1} if you want to parallelize, else \code{0}.
# @param nIterations Number of iterations for NMF.
# @param nCores If \code{parallel} is set to \code{1}
# @param seed_val The seed to be set.
# @param logfile The log file name.
# @param set_verbose Default to \code{0} which will not print any messages, or
# can be set to \code{1} which will print messages. Value passed to the
# \code{get_q2_val} function
#
# @return A tibble of grid_search_results
#
# @importFrom magrittr %>%
# @importFrom purrr cross_df
# @importFrom parallel makeCluster stopCluster detectCores clusterEvalQ
# @importFrom parallel clusterExport
.cv_model_select_pyNMF <- function(X,
                                    param_ranges,
                                    kFolds = 5,
                                    parallelDo = FALSE,
                                    nCores = NA,
                                    nIterations = 20,
                                    seed_val = 10208090,
                                    logfile = "outfile.txt",
                                    set_verbose = 1) {
    if (!is.matrix(X)) {
        stop("X not of type matrix")
    }
    ##
    if (kFolds < 3) {
        if (kFolds < 0) {
            stop("Number of cross-validation folds cannot be negative")
        } else {
            stop("Set at least 3 cross-validation folds")
        }
    } else if (kFolds > ncol(X)) {
        stop("CV folds should be less than or equal to #sequences.
                Standard values: 3, 5, 10.")
    }
    ## Check names in param_ranges list, the function relies on it below
    if (length(setdiff(names(param_ranges), c("alphaPow", "alphaBase",
                                                "k_vals"))) > 0) {
        stop(paste0(
            "Check param_ranges list, expecting three element names: ",
            c("alphaBase", "alphaPow", "k_vals")
        ))
    }
    ## Get cross-validation folds
    cvfolds <-
        .generate_folds(dim(X), kFolds, seed_val = seed_val)
    ## Convert to data frame grid, tidyverse approach
    ## Params to tune: alphaP, alphaA, #factors
    grid_search_params <- purrr::cross_df(list(
        k_vals = param_ranges$k_vals,
        alpha = param_ranges$alphaBase ^ param_ranges$alphaPow,
        fold = seq_len(kFolds),
        iteration = seq_len(nIterations),
        seed_val = seed_val,
        verbose = set_verbose
    ))
    ##
    cat(paste0("Grid search: ", nrow(grid_search_params), " combinations\n"))
    ##
    if (parallelDo) {
        cat(paste0("Opted: Parallel for grid search\n"))
        if (is.na(nCores)) {
            ## raise error or handle
            print("'nCores' not specified, performing serial grid search")
            stop("Number of cores to use not specified")
        } else {
            if (nCores <= parallel::detectCores()) {
                cat(paste0("No. of cores: ", nCores, "\n"))
            } else {
                stop("Specified more than available cores. Stopping ")
            }
            if (nCores > nrow(grid_search_params)) {
                stop("nCores > number of individual computations. Stopping ")
            }
        }
        cl <-
            parallel::makeCluster(nCores, type = "FORK", outfile = logfile)
        parallel::clusterEvalQ(cl, suppressWarnings(require(MASS)))
        ## ^for pseudo-inverse using function `ginv`
        parallel::clusterExport(
            cl = cl,
            varlist = c(".get_q2_using_py", ".compute_q2", "X", "cvfolds"),
            envir = environment()
        )
        ##
        q2_vals <-
            unlist(parallel::clusterApplyLB(cl = cl,
                                            seq_len(nrow(grid_search_params)),
                                            function(i) {
                                                .get_q2_using_py(
                                                    grid_search_params[i,],
                                                    seed_val,
                                                    verbose = set_verbose)
                                            }))
        cat("Stopping cluster...")
        parallel::stopCluster(cl)
        cat("done!\n")
    }
    if (!parallelDo) {
        # TO-DO: check params passed to rowLapply
        cat(paste0("Opted: Serial\n"))
        X <<- X
        cvfolds <<- cvfolds
        q2_vals <- unlist(BBmisc::rowLapply(
                            grid_search_params,
                            .get_q2_using_py,
                            seed_val,
                            verbose = set_verbose
                            )
                    )
    }
    ## The dummy copy for satisfying the constraint in rslurm.
    ## The filtered true copy maintained for downstream steps in the procedure
    grid_search_params <-
        dplyr::select(grid_search_params, k_vals, alpha, fold,
                        iteration)
    grid_search_results <-
        tibble::add_column(grid_search_params, q2_vals)
    return(grid_search_results)
}
## =============================================================================

# @title Generate Cross-Validation Data Splits
#
# @description This function generates the row and column indices for the
# cross-validation splits.
#
# @param Xdims Dimensions of matrix data matrix X.
# @param kFolds Number of cross-validation folds.
# @param seed_val The seed value.
#
# @return A list of two elements: rowIDs and columnIDs for different cross-
# validation folds.
# @importFrom cvTools cvFolds
#
.generate_folds <- function(Xdims, kFolds, seed_val = 10208090) {
    ## Argument seed_val is redundant given that the call set.seed which was
    ## here earlier is now removed.
    ## Xdims gives the dimensions of the matrix X
    cvf_rows <- cvTools::cvFolds(Xdims[1], K = kFolds, type = "random")
    cvf_cols <- cvTools::cvFolds(Xdims[2], K = kFolds, type = "random")
    return(list(cvf_rows = cvf_rows, cvf_cols = cvf_cols))
}
## =============================================================================

# @title Get the best performing value of K (number of factors in NMF)
#
# @param x A tibble, grid_search_results, as returned by
# \code{.cv_model_select_pyNMF}
#
# @return A number The best performing value of K.
.get_best_K <- function(x) {
    # Assumes, max q2_val is best
    # Returns simply the best performing K value
    # Check names in param_ranges list, the function relies on it below
    if (length(setdiff(
        names(x),
        c("k_vals", "alpha", "fold", "iteration",
            "q2_vals")
    )) > 0) {
        message(names(x))
        stop(paste0(
            "Check colnames in tibble, expecting five element names: ",
            c("k_vals", "alpha", "fold", "iteration", "q2_vals")
        ))
    }
    ##
    averages <-
        .get_q2_aggregates_chosen_var(x, chosen_var = x$k_vals, mean)
    idx_best <- as.numeric(which.max(unlist(averages["q2_vals"])))
    ##
    ## TO-DO: Check this
    ## q2_std <- unlist(aggregate(x, by=list(k = x$k_vals), sd)["q2_vals"])
    ## q2_threshold <- as.numeric( x[idx_best,"q2_vals"] - q2_std )
    ##
    best_K <- as.numeric(averages[idx_best, "rel_var"])
    return(best_K)
}
## =============================================================================


# @title Aggregate \eqn{Q^2} Values
#
# @description Aggregate the \eqn{Q^2} values from the grid search results.
#
# @param x The return object from \code{\link{.cv_model_select_pyNMF}}.
# @param chosen_var The variable to aggregate over.
# @param chosen_func The aggregate function to use (should be a function
# already existing wihtin R). Possible values are: \code{mean} and \code{sd}.
#
# @return The mean of $Q^2$ values per the chosen variable
# @importFrom stats aggregate
#
.get_q2_aggregates_chosen_var <- function(x, chosen_var, chosen_func) {
    ## Returns the mean of q2 values per the chosen variable
    averages <- stats::aggregate(x,
                                by = list(rel_var = chosen_var),
                                chosen_func,
                                simplify = TRUE)
    return(averages)
}
## =============================================================================


# @title Get Threshold Value for Selecting \eqn{\alpha}
#
# @description Get the threshold value for selection of \eqn{\alpha} by
# looking at cross-validation performance for K.
#
# @param model_selectK Cross-validation performance over K values.
#
# @return The \eqn{Q^2} threshold value.
# @importFrom stats sd
.get_q2_threshold_by_K <- function(model_selectK) {
    ##
    mean_by_K <-
        .get_q2_aggregates_chosen_var(model_selectK, model_selectK$k_vals,
                                        base::mean)
    sd_by_K <-
        .get_q2_aggregates_chosen_var(model_selectK, model_selectK$k_vals,
                                        stats::sd)
    se_by_K <- sd_by_K / sqrt(nrow(sd_by_K))
    ##
    best_K <- .get_best_K(model_selectK)
    idx_best_K <- which(sd_by_K$rel_var == best_K)
    ##
    q2_threshold <- list(mean = as.numeric(mean_by_K[idx_best_K, "q2_vals"]),
                        se = as.numeric(se_by_K[idx_best_K, "q2_vals"]))
    ##
    return(q2_threshold)
}
## =============================================================================


# @title Get Best \eqn{\alpha}.
#
# @description Get the best performing value of \eqn{\alpha}.
#
# @param for_alpha The return value from \code{\link{.cv_model_select_pyNMF}}.
# This is used for \eqn{\alpha}.
# @param for_k The return value from \code{\link{.cv_model_select_pyNMF}}.
# This is used for K.
# @param min_or_max Specify whether min or max is to be used as the condition
# to choose one when multiple values satisfy the threshold.
#
# @return The best performing value of \eqn{\alpha}.
.get_best_alpha <- function(for_alpha, for_k, min_or_max = min) {
    #
    # Use the one standard error rule
    # Value of alpha that has a reconstruction err not
    # more than 1 std.err of the mean q2 for best k at
    # alpha = 0
    #
    q2_threshold <- .get_q2_threshold_by_K(for_k)
    #
    cat(paste0("Q2 threshold, mean: ", q2_threshold$mean, "\n"))
    averagesA <-
        .get_q2_aggregates_chosen_var(for_alpha, for_alpha$alpha, mean)
    #
    idx_best_alpha <-
        which(unlist(averagesA["q2_vals"]) > q2_threshold$mean)
    #
    if (length(idx_best_alpha) > 1) {
        cat("IF")
        cat("Choosing the highest alpha: ")
        # Choose one (the highest -- max value) if many satisfy threshold
        cat(averagesA[idx_best_alpha, "alpha"])
        best_alpha <-
            min_or_max(as.numeric(averagesA[idx_best_alpha, "alpha"]))
    } else {
        cat("ELSE")
        print(idx_best_alpha)
        print(averagesA)
        best_alpha <- averagesA[idx_best_alpha, "alpha"]
    }
    return(best_alpha)
}
## =============================================================================

# @title Plot Cross-Validation Performance of K
#
# @description Plot showing performance of different values of K tested in
# cross-validation.
#
# @param averages The average of performance values for different combinations
# in grid search
#
# @return A ggplot object so you can simply call \code{print} or \code{save}
# on it later.
#
# @import ggplot2
.plot_cv_K <- function(averages) {
    # Using ggplot to plot
    if ("rel_var" %in% averages) {
        cat("Plotting Q2 vs. K: check colnames in the object
        (need 'rel_var' as 'k')")
        return(NULL)
    }
    if ("q2_vals" %in% averages) {
        cat("Plotting Q2 vs. K: check colnames in the object (need 'q2_vals')")
        return(NULL)
    }
    cat("Plotting Q2 as a function of K")
    p1 <- ggplot2::ggplot(averages, aes_string(x = "rel_var", y = "q2_vals")) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = averages$rel_var,
                            labels = averages$rel_var) +
        ggplot2::labs(
            title = "Reconstruction accuracy, Q\U00B2 = f(#Factors)",
            x = "#Factors (K)",
            y = paste0("Reconstruction accuracy (Q\U00B2)")
        )
    return(p1)
}
## =============================================================================


# @title Plot Cross-Validation Performance of \eqn{\alpha}
#
# @description Plot showing performance of different values of \eqn{\alpha}
# tested in cross-validation.
#
# @param model_selectAlpha The return value from
# \code{\link{.cv_model_select_pyNMF}}. This is used for \eqn{\alpha}.
# @param threshold The threshold to be applied when choosing the best
# performing value
#
# @return A ggplot object so you can simply call \code{print} or \code{save}
# on it later.
#
# @import ggplot2
# @importFrom rlang .data
# @importFrom stats sd
.plot_cv_Alpha <- function(model_selectAlpha, threshold = 0.0) {
    ## Using ggplot to plot
    mean_by_Alpha <- .get_q2_aggregates_chosen_var(model_selectAlpha,
                                                    model_selectAlpha$alpha,
                                                    base::mean)
    sd_by_Alpha <- .get_q2_aggregates_chosen_var(model_selectAlpha,
                                                    model_selectAlpha$alpha,
                                                    stats::sd)
    se_by_Alpha <- sd_by_Alpha / sqrt(nrow(sd_by_Alpha))
    cat("Plotting Q2 as a function of alpha")
    p1 <-
        ggplot2::ggplot(mean_by_Alpha, aes(x = log2(rlang::.data$rel_var),
                                            y = rlang::.data$q2_vals)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = threshold$mean - threshold$se,
                    linetype = "dashed") +
        ggplot2::geom_hline(yintercept = threshold$mean + threshold$se,
                    linetype = "dashed") +
        ggplot2::scale_x_continuous(
            breaks = log2(mean_by_Alpha$rel_var),
            labels = log2(mean_by_Alpha$rel_var)
        ) +
        ggplot2::labs(
            title = "Reconstruction accuracy, Q\U00B2 = f(\u03B1)",
            x = paste0("\u03B1", "=", expression(2 ^ x)),
            y = paste0("Reconstruction accuracy (Q\U00B2)")
        )
    return(p1)
}
## =============================================================================
