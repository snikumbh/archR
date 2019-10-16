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
    # size 4*seq_len
    dna_alphabet <- c("A", "C", "G", "T")
    seq_len <- length(givenSeq)
    if (seq_len > 0) {
        one_hot_encoded <- matrix(rep(0, length(dna_alphabet) * seq_len),
                                  nrow = 1,
                                  byrow = TRUE)
        # characters to match
        one_hot_encoded[, 0 * seq_len + which(givenSeq == dna_alphabet[1])] <-
            1
        one_hot_encoded[, 1 * seq_len + which(givenSeq == dna_alphabet[2])] <-
            1
        one_hot_encoded[, 2 * seq_len + which(givenSeq == dna_alphabet[3])] <-
            1
        one_hot_encoded[, 3 * seq_len + which(givenSeq == dna_alphabet[4])] <-
            1
        return(one_hot_encoded)
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
    # size 4*seq_len
    dna_alphabet <- c("A", "C", "G", "T")
    dna_alphabet_dinuc <- do.call(paste0, expand.grid(dna_alphabet,
                                                      dna_alphabet))
    given_seq_len <- length(givenSeq)
    givenSeq_dinuc <- unlist(lapply(seq_len(given_seq_len - 1), function(x) {
        paste0(givenSeq[x], givenSeq[x + 1])
    }))
    if (given_seq_len > 0) {
        one_hot_encoded_dinuc_profile <- matrix(
            rep(0, length(dna_alphabet_dinuc) *
            given_seq_len), nrow = 1, byrow = TRUE)
        #
        for (i in seq_along(dna_alphabet_dinuc)) {
            one_hot_encoded_dinuc_profile[, (i - 1) * given_seq_len +
                    which(givenSeq_dinuc == dna_alphabet_dinuc[i])] <- 1
        }
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
    seq_len <- length(oneHotEncodedSeqV)/length(dna_alphabet)
    decodedSeq <- rep("Q", seq_len)
    # print(decodedSeq)
    for (alpha_char in seq_along(dna_alphabet)) {
        cutp <- seq_len
        startp <- (alpha_char * cutp) - cutp + 1
        endp <- (alpha_char * cutp)
        decodedSeq[which(oneHotEncodedSeqV[startp:endp] == 1)] <-
            dna_alphabet[alpha_char]
        # print(decodedSeq)
    }
    return(paste0(decodedSeq, collapse = ""))
}
## =============================================================================

#' @title Get One-Hot Encoded Sequences
#'
#' @description Get the one-hot encoding representation of the given sequences.
#'
#' @param givenFastaSeqs A DNAStringSet object holding the given DNA sequences.
#' @param sinuc_or_dinuc 'sinuc' or 'dinuc'
#'
#' @return One-hot encoded sequences.
#' @export
get_one_hot_encoded_seqs <- function(givenFastaSeqs, sinuc_or_dinuc = "sinuc") {
    #
    if (!is.null(givenFastaSeqs) && length(givenFastaSeqs) > 0) {
        seqs_split_as_list <-
            base::strsplit(as.character(givenFastaSeqs), split = NULL)
        if (length(seqs_split_as_list) > 0) {
            if (sinuc_or_dinuc == "sinuc") {
                encoded_seqs <- lapply(seqs_split_as_list,
                                       .one_hot_encode_sinuc)
            } else if (sinuc_or_dinuc == "dinuc") {
                message("Generating dinucleotide profiles")
                encoded_seqs <-
                    lapply(seqs_split_as_list, .one_hot_encode_dinuc)
            }
            encoded_seqs <- do.call(rbind, encoded_seqs)
            return(t(encoded_seqs))
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
        warning("Non DNA-alphabet character in the sequences: ",
                char_levels, "\n")

    } else {
        # All OK!
        message("Sequences OK, ", levels(length_vals)[1])
    }
}
## =============================================================================

#' @title
#' Given a FASTA file, this function return a matrix with one-hot encoded
#' sequences along the columns or the simply the raw sequences.
#'
#' @description
#' In the data matrix, the features are along rows, and sequences along columns.
#' Currently, only mono-character features for DNA sequences are supported.
#' Therefore, the length of the feature vector is four times the
#' length of the sequences (since the DNA alphabet is four characters).
#'
#' @param inputFastaFilename Provide the name (with complete path) of the input
#' FASTA file.
#'
#' @param rawSeq TRUE or FALSE, set this to TRUE if you want the raw sequences.
#'
#' @param sinuc_or_dinuc character string, 'sinuc' or 'dinuc' to select for
#' single or dinucleotide profiles.
#'
#' @return A matrix of sequences represented with one-hot-encoding. Dimensions
#' of the matrix: 4*(sequence length) x number of sequences.
#' @importFrom Biostrings DNAStringSet
#' @export
prepare_data_from_FASTA <- function(inputFastaFilename, rawSeq = FALSE,
                                    sinuc_or_dinuc = "sinuc") {
    if (file.exists(inputFastaFilename)) {
        givenSeqs <-
            Biostrings::readDNAStringSet(filepath = inputFastaFilename,
                                         format = "fasta",
                                         use.names = FALSE)
        # givenSeqs <- seqinr::read.fasta(inputFastaFilename, seqtype = "DNA",
        #                                 as.string = TRUE)
    } else {
        stop("File not found, please check if it exists")
    }
    #
    givenSeqs <- Biostrings::DNAStringSet(toupper(givenSeqs))
    if (rawSeq) {
        return(givenSeqs)
    } else {
        #
        .assert_seq_attributes(givenSeqs)
        message("Read ", length(givenSeqs), " sequences")
        #
        oheSeqs <- get_one_hot_encoded_seqs(givenSeqs,
                                            sinuc_or_dinuc = sinuc_or_dinuc)
        return(oheSeqs)
    }
}
## =============================================================================
