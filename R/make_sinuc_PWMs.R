#' @title Make a PWM-resembling matrix out of a given n-vector
#'
#' @description The given matrix (or simply a vector) is reshaped to have four rows for four nucleotides and a relevant number of columns.
#'
#' @param givenMatrix Actually a vector that will be reshaped into a (PWM)
#' matrix of DNA sequences.
#' @param add_pseudo_counts Boolean, taking values TRUE/T or FALSE/F, specifying
#' whether or not pseudocounts are added to the matrix.
#' @param scale Boolean, taking values TRUE/T or FALSE/F, specifying whether or
#' not the matrix is scaled column-wise, i.e., all columns summed to 1.
#'
#' @return A (PWM) matrix with 4 rows corresponding to the 4 nucleotides (A, C,
#' G, T) and the relevant number of columns (i.e., number of elements in given
#' vector/4)
#' @export
#'
#'
#' @examples
#'
#'
make_sinuc_PWMs <- function(givenMatrix, add_pseudo_counts = T, scale = T){
      # return PWM matrix
      #
      # TO-DO: Make more clear whether a matrix or a colum vector is expected
      #
      # if(!is.matrix(givenMatrix)){
      #   stop("givenMatrix not of type matrix")
      # }
      # if(sum(dim(givenMatrix)) == 2 && is.na(givenMatrix)){
      #   stop("Empty givenMatrix")
      # }
      sinuc <- c('A', 'C', 'G', 'T')
      if(add_pseudo_counts){
          givenMatrix <- givenMatrix + 10^-5
      }
      this_givenMatrix <- t(matrix(givenMatrix, ncol = length(sinuc), byrow=FALSE))
      rownames(this_givenMatrix) <- sinuc
      if(scale){
          scaled <- this_givenMatrix
          scaled <- sweep(this_givenMatrix, 2, colSums(this_givenMatrix), '/')
          return (scaled)
      }else{
          return(this_givenMatrix)
      }
}
