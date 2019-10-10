#' @title Get a specified number of clusters from the given data.
#'
#' @description Get the specified number of clusters from the given data using
#' the chosen clusering method.
#'
#' @param givenMat A data matrix representing sequences along columns and their
#' features aling the rows.
#' @param clustMethod The method used for clustering the data. Default is
#' \code{kmeans}. The other value supported currently in \code{hclust}.
#' @param nCluster The numebr of clusters to be extracted.
#' @param distMethod.hclust If hclust is the chosen clustering method, provide
#' the distance function to be used.
#'
#' @return A list containing the clustering result \code{clust_sol}, and the
#' reordered indices \code{reordering_idx}.
#' @export
#' @importFrom stats dist kmeans hclust
#'
get_clusters <- function(givenMat, clustMethod, nCluster,
                         distMethod.hclust = "euclidean") {
    #
    # For consistency, the input argument matrix has sequences along columns,
    # but
    # we need sequences along rows here.
    givenMat2 <- t(givenMat)
    nSeqs <- ncol(givenMat)
    #
    #
    if (is.na(givenMat) && sum(dim(givenMat)) == 2) {
        stop("Empty matrix")
    }
    #
    if (length(nCluster) > 1) {
        stop("Expecting only one value for nCluster")
    }
    # else if (nCluster < 2) {
    #     stop("Ask for at least 2 clusters")
    # }
    else if (nCluster > ncol(givenMat)) {
        stop("nClusters more than #sequences")
    } else if (nCluster == 1) {
          ######
          # print("kmeans WITH ONLY 1 CLUSTER")
          kmeans_result <- NULL
          ######
          # print("=== kmeans in getCluster ===")
          # print(kmeans_result)
          # reorder sequences in the matrix by clusters
          reordering_idx <- vector("list", 1)
          reordering_idx[[1]] <- seq_len(nSeqs)
          #
          return(list(clustType = "kmeans", clust_sol = kmeans_result,
                      reordering_idx = reordering_idx))
          #
    } else {
          if (clustMethod == "kmeans") {
              # start <- Sys.time()
              kmeans_result <- suppressWarnings(
                                  stats::kmeans(givenMat2,
                                      centers = nCluster,
                                      iter.max = 1000,
                                      nstart = 50,
                                      algorithm = "Lloyd")
                                  )
              # print(Sys.time() - start)
              # print("=== kmeans in getCluster ===")
              # print(kmeans_result)
              # reorder sequences in the matrix by clusters
              reordering_idx <- lapply(seq_len(nCluster), function(x)
                which(kmeans_result$cluster == x))
              #
              return(list(clustType = "kmeans", clust_sol = kmeans_result,
                          reordering_idx = reordering_idx))
          }
          else {
              stop("Wrong clustMethod passed. Takes
                   'kmeans'/'hclust'/'kmedoids'")
          }
    }
    #
}
