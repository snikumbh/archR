#' @title Find the appropriate numebr of clusters for the given data.
#'
#' @description Tests for different clusterings (different number of clusters)
#' of the given data. The average silhoutte value is used for choosing the best
#' suitable number of clusters.
#'
#' @param givenMat A data matrix representing sequences along columns and their
#' features aling the rows.
#' @param distMethod The distance method to be used. Default is \code{euclidean}.
#' Takes all values that can be passed to the \code{dist} method from the
#' \code{stats} package.
#' @param clustMethod The method used for clustering the data. Default is
#' \code{kmeans}. The other value supported currently in \code{hclust}.
#' @param nCluster_vals_test The different number of clusters to be tested.
#' Default is \code{\(3,\ldots,5\)}.
#'
#' @return The number of clusters most suitable for the given data.
#' @export
#'
#'
#'
#' @examples
#'
#'
choose_clusters <- function(givenMat,
                            distMethod = "euclidean",
                            clustMethod = "kmeans",
                            nCluster_vals_test = seq(3,5)
                            ){
  require(cluster)
  #
  # Default clustMethod: kmeans, other options: hclust?
  # Default distMethod: euclidean, other options: manhattan?
  #
  chosen_nClust_val <- 0
  # For consistency, the input argument matrix has sequences along columns, but
  # we need sequences along rows here.
  givenMat2 <- t(givenMat)
  #
  if(is.na(givenMat) && sum(dim(givenMat)) == 2){
        stop("Empty matrix")
  }
  #
  if(distMethod == "euclidean"){
        distMat <- dist(givenMat2, method = "euclidean")
  }
  else if(distMethod == "manhattan"){
        distMat <- dist(givenMat2, method = "manhattan")
  }
  else{
        stop("Wrong distMethod passed [Takes 'euclidean'/'manhattan'].")
  }
  if(max(nCluster_vals_test) < 2 ){
        stop("Ask for at least 2 clusters")
  }else if(any(nCluster_vals_test > ncol(givenMat))){
        stop("nClusters more than #sequences")
  }else{
        sil_vals <- matrix(rep(c(-1,-1), length(nCluster_vals_test) ), ncol = 2,
                           byrow = T)
        colnames(sil_vals) <- c("nClustVals", "Silhouette values")
        sil_vals[,1] <- nCluster_vals_test
        #
        if(clustMethod == "kmeans"){
              print("kmeans")
              start <- Sys.time()
              for (i in 1:length(nCluster_vals_test)){
                kmeans_res <- suppressWarnings(
                                    kmeans(givenMat2,
                                            centers = nCluster_vals_test[i],
                                            iter.max = 1000,
                                            nstart = 50,
                                            algorithm = "Lloyd"))
                sils <- silhouette(kmeans_res$cluster, dist=distMat)
                sil_vals[i,2] <- mean(sils[,"sil_width"])
              }
              print(Sys.time()-start)
              # Report best value
              print(sil_vals)
              chosen_nClust_val <- sil_vals[which.max(sil_vals[,2]),1]
              #
        }else if(clustMethod == "hclust"){
              print("hclust")
              start <- Sys.time()
              for (i in 1:length(nCluster_vals_test)){
                hclust_res <- hclust(distMat)
                sils <- silhouette(cutree(hclust_res, nCluster_vals_test[i]), dist=distMat)
                sil_vals[i,2] <- mean(sils[,"sil_width"])
              }
              print(Sys.time()-start)
              # Report best value
              print(sil_vals)
              chosen_nClust_val <- sil_vals[which.max(sil_vals[,2]),1]
              # cat("Number of clusters: ", sil_vals[which.max(sil_vals[,2]),1], "\n")
        }
        else{
              stop("Wrong clustMethod passed. Takes 'kmeans'/'hclust'")
        }
  }
  return (chosen_nClust_val)
}
