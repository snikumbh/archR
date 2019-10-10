#' @title Sparsify a matrix
#'
#' @description This function sparsifies a given matrix below a given a
#' threshold.
#'
#' @param given_mat The matrix to be sparsified.
#' @param threshold A value between 0.1 and 0.9 inclusive.
#'
#' @return Sparsified matrix of same dimension
#' @export
sparsify_mat <- function(given_mat, threshold = 0.5) {
    # Are the factors on the columns of the matrix? or on rows?
    # Currently assumed, on
    # the columns.  print(norm(given_mat, type='2'))
    if (!is.matrix(given_mat)) {
        stop("Expecting matrix")
    }
    if (sum(dim(given_mat)) == 2 && is.na(given_mat)) {
        stop("Matrix empty")
    } else {
        if (threshold < 0) {
            stop("Threshold should be non-negative")
        } else {
            sparsified_mat <- given_mat
            for (i in 1:ncol(given_mat)) {
                # print(sparsified_mat)
                threshold_i <- threshold * max(given_mat[, i])
                i.set0 <- which(sparsified_mat[, i] < threshold_i)
                sparsified_mat[i.set0, i] <- 0
                # print(sparsified_mat)
            }
            # re-scale
            sparsified_mat <- sparsified_mat *
                (norm(given_mat, type = "2")/norm(sparsified_mat,
                type = "2"))
            # print(norm(sparsified_mat, type='2'))
            return(sparsified_mat)
        }
    }
}
