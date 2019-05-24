#' Get the reconstruction accuracy, \eqn{Q^2}, for a given choice of parameter
#' values
#'
#' @param x The parameter values. Expected columns "k_vals", "alpha" and "fold".
#' @param seed_val The seed to be used.
#' @param verbose Default to \code{0} which will not print any messages, or can
#' be set to \code{1} which will print messages.
#'
#' @return \eqn{Q^2} value
#' @export
#'
#' @examples
#'
#'
get_q2_using_py <- function(x, seed_val, verbose = 0){
      ##
      if(verbose == 1){
          print( paste("INFO-START", x["k_vals"], x["alpha"], sep=",") )
      }else{
          cat(".")
      }
      this_k <- as.numeric(x["k_vals"])
      this_alpha <- as.numeric(x["alpha"])
      #
      test_fold <- as.numeric(x["fold"])
      train_rows <- cvfolds$cvf_rows$subsets[cvfolds$cvf_rows$which != test_fold]
      train_cols <- cvfolds$cvf_cols$subsets[cvfolds$cvf_cols$which != test_fold]
      test_rows <- cvfolds$cvf_rows$subsets[cvfolds$cvf_rows$which == test_fold]
      test_cols <- cvfolds$cvf_cols$subsets[cvfolds$cvf_cols$which == test_fold]

      # Split data matrix X -- separate the training and test parts
      #        _      _
      #   X = |  A  B  |
      #       |_ C  D _|
      #
      # Reconstruct A, by performing NMF on D
      # More details in Owen and Perry, Annals of Statistics, 2009

      submatrixD <- X[train_rows, train_cols]
      submatrixA <- X[test_rows, test_cols]
      submatrixB <- X[test_rows, train_cols]
      submatrixC <- X[train_rows, test_cols]

      # NMF on submatrixD
      # Setup params
      # NMF call
      # Using python/scikit-learn NMF
      nmf_submatrixD <- perform_nmf(submatrixD,
                                    nPatterns = as.integer(this_k),
                                    nIter = as.integer(1000),
                                    givenAlpha = this_alpha,
                                    givenL1_ratio = 1,
                                    seed_val = as.integer(seed_val)
                                    )
      #
      D_W <- nmf_submatrixD[[1]]
      D_H <- nmf_submatrixD[[2]]
      #
      reconstructed_submatrixA <- submatrixB %*% ginv(D_H) %*% ginv(D_W) %*% submatrixC
      #
      if(verbose == 1){
          print( paste("INFO-END", x["k_vals"], x["alpha"], sep=",") )
      }else{
          cat(",")
      }
      #
      q2 <- compute_q2(submatrixA, reconstructed_submatrixA)
      return (q2)
}



#' Compute \eqn{Q^2} value
#'
#' @param A The original matrix
#' @param recA The reconstructed matrix
#'
#' @return Q^2, the reconstruction accuracy
#' @export
#'
#' @examples
#'
#'
compute_q2 <- function(A, recA){
  # input: original matrix, reconstructed portion
  # return: computed q2 val

  second_term_num <- sum((A - recA)^2)
  second_term_den <- sum(A^2) #TO-DO: make sure this is not zero
  if(second_term_den != 0.0){
    q2 <- 1 - (second_term_num/second_term_den)
  }else{
    q2 <- 0.0
  }
  return(q2)
}




#' Perform cross-validation and model selection
#'
#' @param X The given data matrix
#' @param param_ranges An object holding the range of values for parameters
#' \code{k}, \code{alphaBase}, and \code{alphaPow}.
#' @param kFolds The number of cross-validation folds.
#' @param parallel Set to \code{1} if you want to parallelize, else \code{0}.
#' @param nCores If \code{parallel} is set to \code{1}
#' @param seed_val The seed to be set.
#' @param logfile The log file name.
#' @param set_verbose Default to \code{0} which will not print any messages, or
#' can be set to \code{1} which will print messages. Value passed to the
#' \code{get_q2_val} function
#'
#' @return A tibble of grid_search_results
#' @export
#'
#' @examples
#'
#'
cv_model_select_pyNMF <- function(X,
                                  param_ranges,
                                  kFolds = 5,
                                  parallel = FALSE,
                                  nCores = NA,
                                  seed_val = 10208090,
                                  logfile='outfile.txt',
                                  set_verbose = 0){
      # For Moore-Penrose pseudoinverse, load them on the clusters directly
      #suppressPackageStartupMessages(library(MASS, quietly = TRUE))

      # Get cross-validation folds
      cvfolds <- generate_folds(dim(X), kFolds, seed_val = seed_val)

      # Params to tune: alphaP, alphaA, #factors
      # Tidyverse approach
      grid_search_params <- list(k_vals = param_ranges$k_vals,
                                 alpha = param_ranges$alphaBase^param_ranges$alphaPow,
                                 fold = 1:kFolds
                                  ) %>% cross_df() # Convert to data frame grid
      cat(paste0("Grid search: ", nrow(grid_search_params), " combinations\n"))
      #
      if(parallel){
        cat(paste0("Opted: Parallel for grid search\n"))
        if(is.na(nCores)){
          #raise error or handle
          print("'nCores' not specified, turning to serially perform grid search")
          stop("Number of cores to use not specified")
        }else{
          if(nCores <= detectCores()){
            cat(paste0("No. of cores: ", nCores, "\n"))
          }else{
            stop("Specified more than available cores. Stopping ")
          }
          if(nCores > nrow(grid_search_params)){
            stop("nCores more than number of individual computations. Stopping ")
          }
        }
        cl <- makeCluster(nCores, type="FORK", outfile=logfile)
        clusterEvalQ(cl, suppressWarnings(require(CoGAPS))) #for NMF using CoGAPS
        clusterEvalQ(cl, suppressWarnings(require(MASS))) #for pseudo-inverse using function `ginv`
        clusterExport(cl=cl, varlist=c("get_q2_using_py", "compute_q2", "X", "cvfolds"),
                      envir=environment())
        # q2_vals <- parallel::parRapply(cl=cl, grid_search_params, get_q2, seed_val)
        #
        # parRapply does not balance load dynamically.
        # We observed many of the nodes lying idle (process mode S in htop)
        # when a few of them were working. This happens because the indices
        # have been split statically beforehand and assigned to each node.
        # Using a load-balancing approach distributes jobs to each node after
        # one assigned to it earlier has been completed. This way one would
        # expect that nodes do not lie idle while there jobs running slowly
        # on few other nodes.
        #
        q2_vals <- unlist(parallel::clusterApplyLB(cl=cl, 1:nrow(grid_search_params),
                                                   function(i){get_q2_using_py(grid_search_params[i,], seed_val, verbose = set_verbose)}
        ))
        cat("Stopping cluster...")
        stopCluster(cl)
        cat("done!\n")
      }
      if(!parallel){
        #TO-DO: check params passed to rowLapply
        cat(paste0("Opted: Serial\n"))
        X <<- X
        cvfolds <<- cvfolds
        #nCoresUse <<- nCores
        q2_vals <- unlist(rowLapply(grid_search_params, get_q2_using_py, seed_val))
      }
      grid_search_results <- add_column(grid_search_params, q2_vals)
      return(grid_search_results)
}


#' Generate cross-validation data splits
#'
#' @param Xdims Dimensions of matrix data matrix X.
#' @param kFolds Number of cross-validation folds.
#' @param seed_val The seed value.
#'
#' @return A list of two elements: rowIDs and columnIDs for different cross-
#' validation folds.
#' @export
#'
#' @examples
#'
generate_folds <- function(Xdims, kFolds, seed_val){
      # Xdims gives the dimensions of the matrix X
      suppressPackageStartupMessages(require(cvTools, quietly = TRUE))
      set.seed(seed_val)
      cvf_rows <- cvFolds(Xdims[1], K = kFolds, type="random")
      cvf_cols <- cvFolds(Xdims[2], K = kFolds, type="random")
      return(list(cvf_rows = cvf_rows, cvf_cols = cvf_cols))
}


#' Get the best performing value of K (number of factors in NMF)
#'
#' @param x A tibble, grid_search_results, as returned by
#' \code{cv_model_select_pyNMF}
#'
#' @return A number The best performing value of K.
#' @export
#'
#' @examples
#'
get_best_k <- function(x){
      # Assumes, max q2_val is best
      # Returns simply the best performing K value
      averages <- get_q2_aggregates_chosen_var(x, chosen_var=x$k_vals, mean)
      idx_best <- as.numeric( which.max( unlist(averages["q2_vals"])) )
      #
      # TO-DO: Check this
      # q2_std <- unlist(aggregate(x, by=list(k = x$k_vals), sd)["q2_vals"])
      # q2_threshold <- as.numeric( x[idx_best,"q2_vals"] - q2_std )
      #
      best_k <- as.numeric( averages[idx_best, "rel_var"] )
      # best_vals <- list(averages = averages,
      #                   best_k = as.numeric( averages[idx_best, "k"] ),
      #                   q2_threshold = q2_threshold )
      return (best_k)
}


# Do something ---------------------------
#' Aggregate \eqn{Q^2} values fromthe grid search results
#'
#' @param x The grid_search_result return value.
#' @param chosen_var The variable to aggregate over.
#' @param chosen_func The aggregate function to use (should be something already
#'  existing wihtin R). Values such as \code{mean}, or \code{sd} are allowed.
#'
#' @return The mean of $Q^2$ values per the chosen variable
#' @export
#'
#' @examples
#'
get_q2_aggregates_chosen_var <- function(x, chosen_var, chosen_func){
      # Returns the mean of q2 values per the chosen variable
      #
      averages <- aggregate(x, by=list(rel_var = chosen_var), chosen_func, simplify = TRUE)
      return( averages )
}


# Do something ---------------------------
#' Get the threshold value for selection of \eqn{\alpha} by looking at
#' cross-validation performance for K.
#'
#' @param model_selectK Cross-validation performance over K values.
#'
#' @return The \eqn{Q^2} threshold value.
#' @export
#'
#' @examples
#'
get_q2_threshold_by_k <- function(model_selectK){
      #
      mean_by_k <- get_q2_aggregates_chosen_var(model_selectK, model_selectK$k_vals, mean)
      sd_by_k <- get_q2_aggregates_chosen_var(model_selectK, model_selectK$k_vals, sd)
      se_by_k <- sd_by_k/sqrt(nrow(sd_by_k))
      #
      best_k  <- get_best_k(model_selectK)
      idx_best_k <- which(sd_by_k$rel_var == best_k)
      #
      q2_threshold <- as.numeric( mean_by_k[idx_best_k, "q2_vals"] - se_by_k[idx_best_k,"q2_vals"])
      #
      return(q2_threshold)
}


# Do something ---------------------------
#' Get the best perfing value of \eqn{\alpha}.
#'
#' @param for_alpha grid_search_result copy for \eqn{\alpha}
#' @param for_k grid_search_search_result copy for K.
#' @param min_or_max Specify whether min or max is to be used as the condition to choose one when multiple values satisfy the threshold.
#'
#' @return The best performing value of \eqn{\alpha}.
#' @export
#'
#' @examples
#'
get_best_alpha <- function(for_alpha, for_k, min_or_max = min){
      #
      # Use the one standard error rule
      # Value of alpha that has a reconstruction err not
      # more than 1 std.err of the mean q2 for best k at
      # alpha = 0
      #
      q2_threshold <- get_q2_threshold_by_k(for_k)
      #
      cat(paste0("Q2 threshold: ", q2_threshold, "\n"))
      averagesA <- get_q2_aggregates_chosen_var(for_alpha, for_alpha$alpha, mean)
      #
      idx_best_alpha <- which(unlist(averagesA["q2_vals"]) > q2_threshold)
      #
      if(length(idx_best_alpha) > 1){
            cat("Choosing the highest alpha: ")
            # Choose one (the highest -- max value) if many satisfy threshold
            cat(averagesA[idx_best_alpha, "alpha"])
            best_alpha <- min_or_max(as.numeric(averagesA[idx_best_alpha, "alpha"]))
      }else{
            best_alpha <- averagesA[idx_best_alpha, "alpha"]
      }
      return( best_alpha )
}

# Do something ---------------------------
#' Title
#'
#' @param averages The average of performance values for different combinations
#' in grid_search
#'
#' @return A ggplot object so you can simply call \code{print} or \code{save}
#' on it later.
#' @export
#'
#' @examples
#'
plot_cv_K <- function(averages){
  # Using ggplot to plot
  if("rel_var" %in% averages){
      cat("Plotting Q2 vs. K: check colnames in the object (need 'rel_var' as 'k')")
      return(NULL)
  }
  if("q2_vals" %in% averages){
      cat("Plotting Q2 vs. K: check colnames in the object (need 'q2_vals')")
      return(NULL)
  }
  cat("Plotting Q2 as a function of K")
  p1 <- ggplot(averages, aes(x=rel_var, y=q2_vals)) +
            geom_point() +
            geom_line() +
            scale_x_continuous(breaks=averages$rel_var, labels=averages$rel_var) +
            labs(title = "Reconstruction accuracy, Q\U00B2 = f(#Factors)",
                 x = "#Factors (K)",
                 y = paste0("Reconstruction accuracy (Q\U00B2)")
            )
  return(p1)
}


# Do something ---------------------------
#' Plot cross-validation performance of alpha values
#'
#' @param averages The cross-validation averages
#' @param threshold The threshold to be applied when choosing the best
#' performing value
#'
#' @return A ggplot object so you can simply call \code{print} or \code{save}
#' on it later.
#' @export
#'
#' @examples
#'
plot_cv_Alpha <- function(averages, threshold = 0.0){
      # Using ggplot to plot
      if("rel_var" %in% averages){
        cat("Plotting Q2 vs. alpha: check colnames in the object (need 'a')")
        return(NULL)
      }
      cat("Plotting Q2 as a function of alpha")
      p1 <- ggplot(averages, aes(x=log2(rel_var), y=q2_vals) ) +
        geom_point() +
        geom_line() +
        geom_hline(yintercept = threshold, linetype="dashed") +
        scale_x_continuous(breaks=log2(averages$rel_var) , labels=log2(averages$rel_var)) +
        labs(title = "Reconstruction accuracy, Q\U00B2 = f(\u03B1)",
             x = paste0("\u03B1", "=", expression(2^x)),
             y = paste0("Reconstruction accuracy (Q\U00B2)")
        )
      return(p1)
}
