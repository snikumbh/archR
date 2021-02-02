# @title One-hot encode
#
# @description One-hot encode a given DNA sequence.
#
# @param givenSeq A single DNA sequence as character/string.
#
# @return The one-hot encoded sequence.
#
.one_hot_encode_sinuc <- function(givenSeq) {
    # Input: A DNA seq as a vector of caharacters (A/C/G/T)
    # Returns: A row matrix of
    # size 4*seqlen
    dna_alphabet <- c("A", "C", "G", "T")
    seqlen <- length(givenSeq)
    if (seqlen > 0) {
        one_hot_encoded <- matrix(rep(0, length(dna_alphabet) * seqlen),
                                    nrow = 1,
                                    byrow = TRUE)
        # characters to match
        one_hot_encoded[, 0 * seqlen + which(givenSeq == dna_alphabet[1])] <-
            1
        one_hot_encoded[, 1 * seqlen + which(givenSeq == dna_alphabet[2])] <-
            1
        one_hot_encoded[, 2 * seqlen + which(givenSeq == dna_alphabet[3])] <-
            1
        one_hot_encoded[, 3 * seqlen + which(givenSeq == dna_alphabet[4])] <-
            1
        colnames(one_hot_encoded) <- paste(rep(dna_alphabet, 
                            each = seqlen), seq_len(seqlen), sep=".")
        return(one_hot_encoded)
    } else {
        stop("Empty or NULL found")
    }
}
## =============================================================================


# @title One-hot encode TRInucleotide profiles
#
# @description One-hot encode the dinucleotide profile of a given DNA sequence.
#
# @param givenSeq A single sequence.
#
# @return The one-hot encoded sequence.
#
.one_hot_encode_trinuc <- function(givenSeq) {
    # Input: A DNA seq as a vector of caharacters (A/C/G/T)
    # Returns: A row matrix of
    # size 4*seqlen
    dna_alphabet <- c("A", "C", "G", "T")
    dna_alphabet_trinuc <- do.call(paste0, expand.grid(dna_alphabet,
                                                        dna_alphabet,
                                                        dna_alphabet))
    given_seqlen <- length(givenSeq)
    givenSeq_trinuc <- unlist(lapply(seq_len(given_seqlen - 2), function(x) {
        paste0(givenSeq[x], givenSeq[x + 1], givenSeq[x + 2])
    }))
    if (given_seqlen > 0) {
        one_hot_encoded_trinuc_profile <- matrix(
            rep(0, length(dna_alphabet_trinuc) *
                    given_seqlen), nrow = 1, byrow = TRUE)
        #
        for (i in seq_along(dna_alphabet_trinuc)) {
            one_hot_encoded_trinuc_profile[, (i - 1) * given_seqlen +
                    which(givenSeq_trinuc == dna_alphabet_trinuc[i])] <- 1
        }
        colnames(one_hot_encoded_trinuc_profile) <- 
            paste(rep(dna_alphabet_trinuc, each = given_seqlen), 
                seq_len(given_seqlen), sep=".")
        return(one_hot_encoded_trinuc_profile)
    } else {
        stop("Empty or NULL found")
    }
}
## =============================================================================



# @title One-hot encode dinucleotide profiles
#
# @description One-hot encode the dinucleotide profile of a given DNA sequence.
#
# @param givenSeq A single sequence.
#
# @return The one-hot encoded sequence.
#
.one_hot_encode_dinuc <- function(givenSeq) {
    # Input: A DNA seq as a vector of caharacters (A/C/G/T)
    # Returns: A row matrix of
    # size 4*seqlen
    dna_alphabet <- c("A", "C", "G", "T")
    dna_alphabet_dinuc <- do.call(paste0, expand.grid(dna_alphabet,
                                                        dna_alphabet))
    given_seqlen <- length(givenSeq)
    givenSeq_dinuc <- unlist(lapply(seq_len(given_seqlen - 1), function(x) {
        paste0(givenSeq[x], givenSeq[x + 1])
    }))
    if (given_seqlen > 0) {
        one_hot_encoded_dinuc_profile <- matrix(
            rep(0, length(dna_alphabet_dinuc) *
            given_seqlen), nrow = 1, byrow = TRUE)
        #
        for (i in seq_along(dna_alphabet_dinuc)) {
            one_hot_encoded_dinuc_profile[, (i - 1) * given_seqlen +
                    which(givenSeq_dinuc == dna_alphabet_dinuc[i])] <- 1
        }
        colnames(one_hot_encoded_dinuc_profile) <- 
            paste(rep(dna_alphabet_dinuc, each = given_seqlen), 
                seq_len(given_seqlen), sep=".")
        return(one_hot_encoded_dinuc_profile)
    } else {
        stop("Empty or NULL found")
    }
}
## =============================================================================

# @title One-hot decode
#
# @description One-hot decode a given one-hot encoded DNA sequence.
#
# @param oneHotEncodedSeqV A single one-hot encoded sequence vector.
#
# @return The one-hot decoded sequence of ACGTs.
#
.one_hot_decode <- function(oneHotEncodedSeqV) {

    dna_alphabet <- c("A", "C", "G", "T")
    seqlen <- length(oneHotEncodedSeqV)/length(dna_alphabet)
    decodedSeq <- rep("Q", seqlen)
    for (alpha_char in seq_along(dna_alphabet)) {
        cutp <- seqlen
        startp <- (alpha_char * cutp) - cutp + 1
        endp <- (alpha_char * cutp)
        decodedSeq[which(oneHotEncodedSeqV[startp:endp] == 1)] <-
            dna_alphabet[alpha_char]
    }
    return(paste0(decodedSeq, collapse = ""))
}
## =============================================================================

#' @title Get one-hot encoded sequences
#'
#' @description Get the one-hot encoding representation of the given sequences.
#'
#' @param seqs A DNAStringSet object holding the given DNA sequences
#' @param sinuc_or_dinuc 'sinuc' or 'dinuc'
#'
#' @return A sparse matrix of sequences represented with one-hot-encoding
#' @family input functions
#' @seealso \code{\link{prepare_data_from_FASTA}} for generating one-hot 
#' encoding of sequences from a FASTA file
#' @importFrom Matrix Matrix
#' 
#' @examples 
#' 
#' fname <- system.file("extdata", "example_data.fa", 
#'                         package = "archR", mustWork = TRUE)
#' 
#'                         
#' rawFasta <- prepare_data_from_FASTA(inputFastaFilename = fname,
#'                         rawSeq = TRUE)
#' 
#' get_one_hot_encoded_seqs(seqs = rawFasta, sinuc_or_dinuc = "dinuc")
#'                         
#' @export
get_one_hot_encoded_seqs <- function(seqs, sinuc_or_dinuc = "sinuc") {
    #
    if (!is.null(seqs) && length(seqs) > 0) {
        seqs_split_as_list <-
            base::strsplit(as.character(seqs), split = NULL)
        if (length(seqs_split_as_list) > 0) {
            if (sinuc_or_dinuc == "sinuc") {
                encoded_seqs <- lapply(seqs_split_as_list,
                                        .one_hot_encode_sinuc)
            } else if (sinuc_or_dinuc == "dinuc") {
                message("Generating dinucleotide profiles")
                encoded_seqs <-
                    lapply(seqs_split_as_list, .one_hot_encode_dinuc)
            }  else if (sinuc_or_dinuc == "trinuc") {
                message("Generating trinucleotide profiles")
                encoded_seqs <-
                    lapply(seqs_split_as_list, .one_hot_encode_trinuc)
            }

            encoded_seqs <- do.call(rbind, encoded_seqs)
            encoded_seqs <- t(encoded_seqs)
            ## Use Matrix package, store as sparse matrix
            ## Helps reduce object size and computation time
            encoded_seqs <- Matrix::Matrix(encoded_seqs, sparse = TRUE)
            ##
            return(encoded_seqs)
        }
    } else {
        stop("Empty or NULL found")
    }
}
## =============================================================================

# @title Assert attributes of sequences
#
# @description Assert the attributes of the sequences provided. This includes
# checking for (1) the length of the sequences, (2) characters in the
# sequences.
#
# @param givenSeqs DNA sequences as a list.
#
#
# @return nothing. Only prints a warning to the screen.
# @importFrom Biostrings width
#
.assert_seq_attributes <- function(givenSeqs) {
    # Check that all sequences are of same length
    seqs_split_as_list <-
        base::strsplit(as.character(givenSeqs), split = NULL)
    # length_vals <- unlist(lapply(seqs_split_as_list, length))
    length_vals <- levels(as.factor(Biostrings::width(givenSeqs)))
    char_levels <- levels(as.factor(unlist(seqs_split_as_list)))
    dna_alphabet <- c("A", "C", "G", "T")
    if (0 %in% length_vals) {
        # Checking sequences of length 0
        stop(paste0("Found ", which(0 == length_vals), " sequence(s) of length
                    zero\n"))
        #
    } else if (length(levels(as.factor(length_vals))) > 1) {
        # Checking all sequences are of same length
        stop("Sequences are of different length\n")
        #
    } else if (any(!(char_levels %in% dna_alphabet))) {
        # Check for non-alphabet characters
        # Raise either an error or just warn!
        warning(c("Non DNA-alphabet character in the sequences: ",
                char_levels), immediate. = TRUE)

    } else {
        # All OK!
        message("Sequences OK, ", levels(length_vals)[1])
    }
}
## =============================================================================

#' @title
#' Generate one-hot encoding of sequences given as FASTA file
#'
#' @description
#' Given a set of sequences in a FASTA file this function returns a sparse 
#' matrix with one-hot encoded sequences.
#' In this matrix, the sequence features are along rows, and sequences along 
#' columns. Currently, mono- and dinucleotide features for DNA sequences are 
#' supported. Therefore, the length of the feature vector is 4 and 16 times 
#' the length of the sequences (since the DNA alphabet is four characters) 
#' for mono- and dinucleotide features respectively.
#'
#' @param inputFastaFilename Provide the name (with complete path) of the input
#' FASTA file.
#'
#' @param rawSeq TRUE or FALSE, set this to TRUE if you want the raw sequences.
#'
#' @param sinuc_or_dinuc character string, 'sinuc' or 'dinuc' to select for
#' mono- or dinucleotide profiles.
#'
#' @return A sparse matrix of sequences represented with one-hot-encoding.
#' @family input functions
#' @seealso \code{\link{get_one_hot_encoded_seqs}} for directly using a 
#' DNAStringSet object
#' @importFrom Biostrings DNAStringSet
#' 
#' @examples 
#' 
#' fname <- system.file("extdata", "example_data.fa", 
#'                         package = "archR", mustWork = TRUE)
#' 
#' # mononucleotides feature matrix
#' prepare_data_from_FASTA(inputFastaFilename = fname,
#'                         sinuc_or_dinuc = "sinuc")
#' 
#' # dinucleotides feature matrix
#' prepare_data_from_FASTA(inputFastaFilename = fname,
#'                         sinuc_or_dinuc = "dinuc")
#'                        
#' # FASTA sequences as a Biostrings::DNAStringSet
#' prepare_data_from_FASTA(inputFastaFilename = fname,
#'                         rawSeq = TRUE)
#' 
#' @export
prepare_data_from_FASTA <- function(inputFastaFilename, rawSeq = FALSE,
                                    sinuc_or_dinuc = "sinuc") {
    if (!file.exists(inputFastaFilename)) {
        stop("File not found, please check if it exists")
    }
    
    if (rawSeq) {
        givenSeqs <- Biostrings::readDNAStringSet(
            filepath = inputFastaFilename, format = "fasta", use.names = TRUE)
        givenSeqs <- Biostrings::DNAStringSet(toupper(givenSeqs))
        return(givenSeqs)
    } else {
        #
        givenSeqs <- Biostrings::readDNAStringSet(
            filepath = inputFastaFilename, format = "fasta", use.names = FALSE)
        givenSeqs <- Biostrings::DNAStringSet(toupper(givenSeqs))
        .assert_seq_attributes(givenSeqs)
        message("Read ", length(givenSeqs), " sequences")
        #
        oheSeqs <- get_one_hot_encoded_seqs(givenSeqs,
                                            sinuc_or_dinuc = sinuc_or_dinuc)
        return(oheSeqs)
    }
}
## =============================================================================

