# @title Make a PWM-resembling matrix out of a given n-vector
#
# @description The given matrix (or simply a vector) is reshaped to have four
# rows for four nucleotides and a relevant number of columns.
#
# @param givenMatrix Actually a vector that will be reshaped into a (PWM)
# matrix of DNA sequences.
# @param add_pseudo_counts Boolean, taking values TRUE/T or FALSE/F, specifying
# whether or not pseudocounts are added to the matrix.
# @param scale Boolean, taking values TRUE/T or FALSE/F, specifying whether or
# not the matrix is scaled column-wise, i.e., all columns summed to 1.
#
# @return A (PWM) matrix with 4 rows corresponding to the 4 nucleotides (A, C,
# G, T) and the relevant number of columns (i.e., number of elements in given
# vector/4)
#
make_sinuc_PWMs <- function(givenMatrix, add_pseudo_counts = TRUE,
                            scale = TRUE) {
    ## return PWM matrix TO-DO: Make more clear whether a matrix or a column
    ## vector is expected
    ## if(!is.matrix(givenMatrix)){
    ## stop('givenMatrix not of type matrix') }
    ## if(sum(dim(givenMatrix)) == 2 && is.na(givenMatrix)){
    ## stop('Empty givenMatrix')
    ## }
    sinuc <- c("A", "C", "G", "T")
    if (add_pseudo_counts) {
        givenMatrix <- givenMatrix + 10^-5
    }
    this_givenMatrix <- t(matrix(givenMatrix, ncol = length(sinuc),
                                    byrow = FALSE))
    rownames(this_givenMatrix) <- sinuc
    if (scale) {
        scaled <- this_givenMatrix
        scaled <- base::sweep(this_givenMatrix, 2, colSums(this_givenMatrix),
                                "/")
        return(scaled)
    } else {
        return(this_givenMatrix)
    }
}




make_dinuc_PWMs <- function(givenMatrix, add_pseudo_counts = TRUE,
                            scale = TRUE) {
    # return PWM matrix TO-DO: Make more clear whether a matrix or a colum
    # vector is expected
    # if(!is.matrix(givenMatrix)){
    # stop('givenMatrix not of type matrix')
    # }
    # if(sum(dim(givenMatrix)) == 2 && is.na(givenMatrix)){
    # stop('Empty givenMatrix')
    # }
    dna_alphabet <- c("A", "C", "G", "T")
    dinuc <- do.call(paste0, expand.grid(dna_alphabet, dna_alphabet))
    if (add_pseudo_counts) {
        givenMatrix <- givenMatrix + 10^-5
    }
    this_givenMatrix <- t(matrix(givenMatrix, ncol = length(dinuc),
                                    byrow = FALSE))
    rownames(this_givenMatrix) <- dinuc
    if (scale) {
        scaled <- this_givenMatrix
        scaled <- base::sweep(this_givenMatrix, 2, colSums(this_givenMatrix),
                                "/")
        return(scaled)
    } else {
        return(this_givenMatrix)
    }
}


## Collapse the dinucleotide matrix to single nucleotide
collapse_into_sinuc_matrix <- function(given_feature_mat, dinuc_mat, feature_names){
    # Collapse the matrix of dinuc features into sinuc features
    # When V is (#seqs x features)
    # sinuc_feature_mat <-  matrix(rep(0, 4*ncol(given_feature_mat)),
    #                              nrow=4, ncol=ncol(given_feature_mat))
    # When V is (features x #seqs)

    sinuc_feature_mat <-  matrix(rep(1, nrow(given_feature_mat)),
                                 nrow=4, ncol=nrow(given_feature_mat)/length(feature_names) )

    ##useMat <- make_dinuc_PWMs(given_feature_mat, scale = FALSE, add_pseudo_counts = FALSE)
    useMat <- dinuc_mat

    rownames(sinuc_feature_mat) <- c('A', 'C', 'G', 'T')
    splitted_names <- strsplit(feature_names, "")


    for( i in 1:length(splitted_names)){
        for (j in 1:(ncol(dinuc_mat)-1)){
            sinuc_feature_mat[splitted_names[[i]][1], j] <- sinuc_feature_mat[splitted_names[[i]][1], j] + 5*useMat[i,j]
            sinuc_feature_mat[splitted_names[[i]][2], j+1] <- sinuc_feature_mat[splitted_names[[i]][2], j+1] + 5*useMat[i,j]
        }
    }
    return(sinuc_feature_mat)
}
