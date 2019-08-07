#' @title One-hot encode
#'
#' @description One-hot encode a given DNA sequence.
#'
#' @param givenSeq A single sequence.
#'
#' @return The one-hot encoded sequence.
#'
one_hot_encode <- function(givenSeq) {
  # Input: A DNA seq as a vector of caharacters (A/C/G/T)
  # Returns: A row matrix of size 4*seq_len
  dna_alphabet <- c("A", "C", "G", "T")
  seq_len <- length(givenSeq)
  if (seq_len > 0) {
    one_hot_encoded <- matrix(rep(0, length(dna_alphabet) * seq_len), nrow = 1, byrow = T)
    # characters to match
    one_hot_encoded[, 0 * seq_len + which(givenSeq == dna_alphabet[1])] <- 1
    one_hot_encoded[, 1 * seq_len + which(givenSeq == dna_alphabet[2])] <- 1
    one_hot_encoded[, 2 * seq_len + which(givenSeq == dna_alphabet[3])] <- 1
    one_hot_encoded[, 3 * seq_len + which(givenSeq == dna_alphabet[4])] <- 1
    return(one_hot_encoded)
  } else {
    stop("One-hot encoding: empty sequence found")
  }
}

#' @title One-hot decode
#'
#' @description One-hot decode a given one-hot encoded DNA sequence.
#'
#' @param oneHotEncodedSeqV A single one-hot encoded sequence vector.
#'
#' @return The one-hot decoded sequence of ACGTs.
#'
one_hot_decode <- function(oneHotEncodedSeqV) {
  
  dna_alphabet <- c("A", "C", "G", "T")
  seq_len <- length(oneHotEncodedSeqV)/length(dna_alphabet)
  decodedSeq <- rep("Q", seq_len)
  # print(decodedSeq)
  for (alpha_char in 1:length(dna_alphabet)){
    cutp <- seq_len
    startp <- (alpha_char*cutp) - cutp + 1
    endp <- (alpha_char*cutp)
    decodedSeq[which(oneHotEncodedSeqV[startp:endp] == 1)] <- dna_alphabet[alpha_char]
    # print(decodedSeq)
  }
  return(paste0(decodedSeq, collapse=""))
}

#' @title Get One-Hot Encoded Sequences
#'
#' @description Get the one-hot encoding representation of the given sequences.
#'
#' @param givenFastaSeqs List of sequences.
#'
#' @return One-hot encoded sequences.
#'
get_one_hot_encoded_seqs <- function(givenFastaSeqs) {
  #
  if (length(givenFastaSeqs) > 0) {
    seqs_split_as_list <- strsplit(givenFastaSeqs, split = NULL)
    if (length(seqs_split_as_list) > 0) {
      encoded_seqs <- lapply(seqs_split_as_list, one_hot_encode)
      encoded_seqs <- do.call(rbind, encoded_seqs)
      return(encoded_seqs)
    }
  } else {
    stop("List of sequences empty")
  }
}


#' @title Assert attributes of sequences
#'
#' @description Assert the attributes of the sequences provided. This includes checking for
#' (1) the length of the sequences, (2) characters in the sequences.
#'
#' @param givenSeqs DNA sequences as a list.
#'
#'
#' @return nothing. Only prints a warning to the screen.
#'
#'
assert_attributes <- function(givenSeqs) {
  # Check that all sequences are of same length
  seqs_split_as_list <- strsplit(givenSeqs, split = NULL)
  length_vals <- unlist(lapply(seqs_split_as_list, length))
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
      # stop("Non DNA-alphabet character in the sequences: ", char_levels, "\n")
      warning("Non DNA-alphabet character in the sequences: ", char_levels, "\n")
      # seqs_split_as_list <- strsplit(givenSeqs, split = NULL)
      # seqs_split_as_matrix <- matrix(unlist(seqs_split_as_list), ncol = length(seqs_split_as_list[[1]]), byrow = T)
      # seqs_list_again <- apply(seqs_split_as_matrix, 1, function(x) if(!"N" %in% x) x)
      # seqs_matrix_again <- matrix(unlist(seqs_list_again), ncol = length(seqs_list_again[[1]]), byrow = T)
      # return (seqs_matrix_again)
      #
  } else {
      # All OK!
      message("Sequences OK, ", levels(length_vals)[1])
  }
}


#' @title
#' Given a FASTA file, prepares the data matrix with one-hot encoded features.
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
#' @return A matrix of sequences represented with one-hot-encoding. Dimensions
#' of the matrix: 4*(sequence length) x number of sequences.
#' @export
prepare_data_from_FASTA <- function(inputFastaFilename) {
  start <- Sys.time()
  if (file.exists(inputFastaFilename)) {
    givenSeqs <- seqinr::read.fasta(inputFastaFilename, seqtype = "DNA", as.string = TRUE)
  } else {
    stop("File not found, please check if it exists")
  }
  #
  givenSeqs <- toupper(givenSeqs)
  # print(givenSeqs[1])
  print(length(givenSeqs))
  # 
  # Keep only selected positions from all given sequences
  seqs_split_as_list <- strsplit(givenSeqs, split = NULL)
  seqs_split_as_matrix <- matrix(unlist(seqs_split_as_list), ncol = length(seqs_split_as_list[[1]]), byrow = T)
  givenSeqs <- apply(seqs_split_as_matrix[ , 56:145], 1, paste, collapse="")
  # 
  # 
  assert_attributes(givenSeqs)
  print(length(givenSeqs))
  #
  oheSeqs <- get_one_hot_encoded_seqs(givenSeqs)
  print(Sys.time() - start)
  return(t(oheSeqs))
}
