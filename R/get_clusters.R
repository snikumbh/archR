#' Get a specified number of clusters from the given data.
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
#'
#' @examples
#'
#'
get_clusters <- function(givenMat, clustMethod,  nCluster, distMethod.hclust = "euclidean"){

  # dist
  #

  #
  if(clustMethod == "kmeans"){
        print("kmeans")
        start <- Sys.time()
        kmeans_result <- suppressWarnings(kmeans(givenMat, centers = nCluster, iter.max = 1000, nstart = 50, algorithm = "Lloyd"))
        print(Sys.time()-start)
        # reorder sequences in the matrix by clusters
        reordering_idx <- sapply(1:nCluster, function(x) which(kmeans_result$cluster == x))
        #
        return (list(clust_sol = kmeans_result, reordering_idx = reordering_idx))
  }else if(clustMethod == "hclust"){
        # TO-DO check this part of code for option hclust
        print("hclust")
        start <- Sys.time()
        distMat <- dist(givenMat, method = distMethod.hclust)
        hclust_result <- hclust(distMat)
        hclust_cl <- cutree(hclust_res, nCluster)
        print(Sys.time()-start)
        # reorder sequences in the matrix by clusters
        reordering_idx <-  hclust_result$order # this is still not cluster-wise ordering, but an ordering that avoids crossings of the branches in the dendrogram
        # sapply(1:nCluster, function(x) which(hclust_result$cluster == x))
        #
        return (list(clust_sol = hclust_result, reordering_idx = reordering_idx))
  }
  #

}
