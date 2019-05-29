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
get_clusters <- function(givenMat, clustMethod, nCluster, distMethod.hclust = "euclidean") {
  #
  # For consistency, the input argument matrix has sequences along columns, but
  # we need sequences along rows here.
  givenMat2 <- t(givenMat)
  #
  #
  if (is.na(givenMat) && sum(dim(givenMat)) == 2) {
    stop("Empty matrix")
  }
  #
  if (length(nCluster) > 1) {
    stop("Expecting only one value for nCluster")
  } else if (nCluster < 2) {
    stop("Ask for at least 2 clusters")
  } else if (nCluster > ncol(givenMat)) {
    stop("nClusters more than #sequences")
  } else {
    if (clustMethod == "kmeans") {
      print("kmeans")
      start <- Sys.time()
      kmeans_result <- suppressWarnings(stats::kmeans(givenMat2, centers = nCluster, iter.max = 1000, nstart = 50, algorithm = "Lloyd"))
      print(Sys.time() - start)
      # reorder sequences in the matrix by clusters
      reordering_idx <- sapply(1:nCluster, function(x) which(kmeans_result$cluster == x))
      #
      return(list(clust_sol = kmeans_result, reordering_idx = reordering_idx))
    } else if (clustMethod == "hclust") {
      # TO-DO check this part of code for option hclust
      print("hclust")
      start <- Sys.time()
      if (distMethod.hclust == "euclidean") {
        distMat <- stats::dist(givenMat2, method = distMethod.hclust)
      } else if (distMethod.hclust == "manhattan") {
        distMat <- stats::dist(givenMat2, method = distMethod.hclust)
      } else {
        stop("Wrong distMethod.hclust")
      }
      hclust_result <- stats::hclust(distMat)
      hclust_cl <- stats::cutree(hclust_result, nCluster)
      print(Sys.time() - start)
      # reorder sequences in the matrix by clusters
      reordering_idx <- hclust_result$order # this is still not cluster-wise ordering, but an ordering that avoids crossings of the branches in the dendrogram
      # sapply(1:nCluster, function(x) which(hclust_result$cluster == x))
      #
      return(list(clust_sol = hclust_cl, reordering_idx = reordering_idx))
    } else {
      stop("Wrong clustMethod passed. Takes 'kmeans'/'hclust'")
    }
  }
  #
}
