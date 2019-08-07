#' @title Find the appropriate numebr of clusters for the given data.
#'
#' @description Tests for different clusterings (different number of clusters)
#' of the given data. The average silhoutte value is used for choosing the best
#' suitable number of clusters.
#'
#' @param givenMatrix A data matrix representing sequences along columns and their
#' features along the rows.
#' @param distMethod The distance method to be used. Default is \code{euclidean}.
#' Takes all values that can be passed to the \code{dist} method from the
#' \code{stats} package.
#' @param clustMethod The method used for clustering the data. Default is
#' \code{kmeans}. The other value supported currently in \code{hclust}.
#' @param nClusterValsTest The different number of clusters to be tested.
#' Default is \code{\(3,\ldots,5\)}.
#'
#' @return The number of clusters most suitable for the given data.
#' @export
#' @importFrom stats dist kmeans hclust cutree
#' @importFrom cluster silhouette
#'
choose_clusters <- function(givenMatrix,
                            distMethod = "euclidean",
                            clustMethod = "kmeans",
                            nClusterValsTest = seq(3, 5)) {
  # require(cluster)
  #
  # Default clustMethod: kmeans, other options: hclust?
  # Default distMethod: euclidean, other options: manhattan?
  #
  chosenNClusterVal <- 0
  # For consistency, the input argument matrix has sequences along columns, but
  # we need sequences along rows here.
  givenMat2 <- t(givenMatrix)
  #
  if (is.na(givenMatrix) && sum(dim(givenMatrix)) == 2) {
      stop("Empty matrix")
  }
  #
  if (distMethod == "euclidean") {
      distMat <- stats::dist(givenMat2, method = "euclidean")
  }
  else if (distMethod == "manhattan") {
      distMat <- stats::dist(givenMat2, method = "manhattan")
  }
  else {
      stop("Wrong distMethod passed [Takes 'euclidean'/'manhattan'].")
  }
  if (max(nClusterValsTest) < 2) {
      stop("Ask for at least 2 clusters")
  } else if (any(nClusterValsTest > ncol(givenMatrix))) {
      stop("nClusters more than #sequences")
  } else {
      silhouetteVals <- matrix(rep(c(-1, -1), length(nClusterValsTest)),
                               ncol = 2, 
                               byrow = T
                          )
      colnames(silhouetteVals) <- c("nClustVals", "Silhouette values")
      silhouetteVals[, 1] <- nClusterValsTest
      #
      if (clustMethod == "kmeans") {
          # Perform k-means clustering  
          print("kmeans")
          start <- Sys.time()
          for (i in 1:length(nClusterValsTest)) {
            kmeansRes <- suppressWarnings(
              stats::kmeans(givenMat2,
                centers = nClusterValsTest[i],
                iter.max = 1000,
                nstart = 50,
                algorithm = "Lloyd"
              )
            )
            sils <- cluster::silhouette(kmeansRes$cluster, dist = distMat)
            silhouetteVals[i, 2] <- mean(sils[, "sil_width"])
          }
          print(Sys.time() - start)
          # Report best value
          print(silhouetteVals)
          chosenNClusterVal <- silhouetteVals[which.max(silhouetteVals[, 2]), 1]
          #
      } else if (clustMethod == "hclust") {
          # Perform hierarchical clustering  
          print("hclust")
          start <- Sys.time()
          for (i in 1:length(nClusterValsTest)) {
            hclustRes <- hclust(distMat)
            sils <- cluster::silhouette(cutree(hclustRes, nClusterValsTest[i]), dist = distMat)
            silhouetteVals[i, 2] <- mean(sils[, "sil_width"])
          }
          print(Sys.time() - start)
          # Report best value
          print(silhouetteVals)
          chosenNClusterVal <- silhouetteVals[
            which.max(silhouetteVals[, 2]), 1
          ]
      } else if (clustMethod == "kmedoids") {
          # Perform k-medoid clustering  
          print("kmedoids")
          start <- Sys.time()
          #
          for (i in 1:length(nClusterValsTest)) {
            set.seed(6666)
            kmedoidsRes <- suppressWarnings(
                              cluster::pam(x = givenMat2,
                                 k = nClusterValsTest[i],
                                 diss = FALSE,
                                 metric = distMethod,
                                 cluster.only = FALSE
                                 )
              
              # cluster::clara(givenMat2,
              #               k = nClusterValsTest[i],
              #               metric = distMethod,
              #               stand = FALSE,
              #               samples = 50,
              #               sampsize = 100,
              #               medoids.x = TRUE,
              #               rngR = TRUE,
              #               correct.d = TRUE
              #               )
                            
            )
            sils <- cluster::silhouette(kmedoidsRes$clustering, dist = distMat)
            silhouetteVals[i, 2] <- mean(sils[, "sil_width"])
            
          }
          #
          print(Sys.time() - start)
          # Report best value
          print(silhouetteVals)
          chosenNClusterVal <- silhouetteVals[which.max(silhouetteVals[, 2]), 1]
        
      }
      else {
        stop("Wrong clustMethod passed. Takes 'kmeans'/'hclust'")
      }
  }
  return(chosenNClusterVal)
}
