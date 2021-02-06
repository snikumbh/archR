#' @title Make a PWM-resembling matrix out of a given n-vector
#'
#' @description The given matrix (or simply a vector) is reshaped to have four
#' rows for four nucleotides and a relevant number of columns.
#'
#' @param mat Actually a vector that will be reshaped into a (PWM)
#' matrix of DNA sequences.
#' @param add_pseudo_counts Boolean, taking values TRUE/T or FALSE/F, specifying
#' whether or not pseudocounts are added to the matrix.
#' @param scale Boolean, taking values TRUE/T or FALSE/F, specifying whether or
#' not the matrix is scaled column-wise, i.e., all columns summed to 1.
#'
#' @return A (PWM) matrix with 4 rows corresponding to the 4 nucleotides (A, C,
#' G, T) and the relevant number of columns (i.e., number of elements in given
#' vector/4)
#'
#' @export
make_sinuc_PWMs <- function(mat, add_pseudo_counts = TRUE,
                            scale = TRUE) {
    ## return PWM matrix TO-DO: Make more clear whether a matrix or a column
    ## vector is expected
    ## if(!is.matrix(mat)){
    ## stop('mat not of type matrix') }
    ## if(sum(dim(mat)) == 2 && is.na(mat)){
    ## stop('Empty mat')
    ## }
    sinuc <- c("A", "C", "G", "T")
    if (add_pseudo_counts) {
        mat <- mat + 10^-5
    }
    this_mat <- t(matrix(mat, ncol = length(sinuc),
                                    byrow = FALSE))
    rownames(this_mat) <- sinuc
    if (scale) {
        scaled <- this_mat
        scaled <- base::sweep(this_mat, 2, colSums(this_mat),
                                "/")
        return(scaled)
    } else {
        return(this_mat)
    }
}


#' @title Similarly to the PWM-like matrix for mononucleotides, make one for 
#'  dinucleotides 
#' 
#' @description This function converts the basis matrix with basis vectors 
#' of dinucleotide information into matrix of dimension 
#' 16 x (sequence_length) for visualization.
#' 
#' @param mat Input basis matrix
#' @param add_pseudo_counts Whether pesudocounts are to be added. TRUE or FALSE.
#' @param scale Whether to perform per position scaling of the matrix. TRUE or 
#' FALSE
#' 
#' @return A (PWM) matrix with 4 rows corresponding to the 4 nucleotides (A, C,
#' G, T) and the relevant number of columns (i.e., number of elements in given
#' vector/4)
#' 
#' @export
make_dinuc_PWMs <- function(mat, add_pseudo_counts = TRUE,
                            scale = TRUE) {
    # return PWM matrix TO-DO: Make more clear whether a matrix or a colum
    # vector is expected
    # if(!is.matrix(mat)){
    # stop('mat not of type matrix')
    # }
    # if(sum(dim(mat)) == 2 && is.na(mat)){
    # stop('Empty mat')
    # }
    dna_alphabet <- c("A", "C", "G", "T")
    dinuc <- do.call(paste0, expand.grid(dna_alphabet, dna_alphabet))
    if (add_pseudo_counts) {
        mat <- mat + 10^-5
    }
    this_mat <- t(matrix(mat, ncol = length(dinuc),
                                    byrow = FALSE))
    rownames(this_mat) <- dinuc
    if (scale) {
        scaled <- this_mat
        scaled <- base::sweep(this_mat, 2, colSums(this_mat),
                                "/")
        return(scaled)
    } else {
        return(this_mat)
    }
}


## Collapse the dinucleotide matrix to single nucleotide
collapse_into_sinuc_matrix <- function(given_feature_mat, 
                                        dinuc_mat, 
                                        feature_names){
    # Collapse the matrix of dinuc features into sinuc features
    # When V is (#seqs x features)
    # sinuc_feature_mat <-  matrix(rep(0, 4*ncol(given_feature_mat)),
    #                              nrow=4, ncol=ncol(given_feature_mat))
    # When V is (features x #seqs)

    sinuc_feature_mat <-  matrix(rep(1, nrow(given_feature_mat)), nrow=4, 
                            ncol=nrow(given_feature_mat)/length(feature_names))

    ##useMat <- make_dinuc_PWMs(given_feature_mat, scale = FALSE, 
    ##add_pseudo_counts = FALSE)
    useMat <- dinuc_mat

    rownames(sinuc_feature_mat) <- c('A', 'C', 'G', 'T')
    splitted_names <- strsplit(feature_names, "")


    for(i in seq_along(splitted_names)){
        for (j in seq_len(ncol(dinuc_mat)-1)){
            sinuc_feature_mat[splitted_names[[i]][1], j] <- 
                sinuc_feature_mat[splitted_names[[i]][1], j] + 5*useMat[i,j]
            sinuc_feature_mat[splitted_names[[i]][2], j+1] <- 
                sinuc_feature_mat[splitted_names[[i]][2], j+1] + 5*useMat[i,j]
        }
    }
    return(sinuc_feature_mat)
}
