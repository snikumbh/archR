# Date created: 27/03/2019
# Author: snikumbh 
#

sparsify_mat <- function(given_mat){
  # Are the factors on the columns of the matrix? or on rows?
  # Currently assumed, on the columns.
  # print(norm(given_mat, type="2"))
  sparsified_mat <- given_mat
  for(i in 1:ncol(given_mat)){
    #print(sparsified_mat)
    threshold_i <- 0.5 * max(given_mat[,i])
    i.set0 <- which(sparsified_mat[,i] < threshold_i)
    sparsified_mat[i.set0, i] <- 0
    #print(sparsified_mat)
  }
  # re-scale
  sparsified_mat <- sparsified_mat *  (norm(given_mat, type="2")/norm(sparsified_mat, type="2"))
  # print(norm(sparsified_mat, type="2"))
  #
  return(sparsified_mat)
}