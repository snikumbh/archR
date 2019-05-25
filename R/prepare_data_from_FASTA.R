#' Title
#'
#' @param givenSeq
#'
#' @return
#'
#'
#' @examples
#'
one_hot_encode <- function(givenSeq){
  # Input: A DNA seq as a vector of caharacters (A/C/G/T)
  # Returns: A row matrix of size 4*seq_len
  dna_alphabet <- c("A", "C", "G", "T")
  seq_len <- length(givenSeq)
  one_hot_encoded <- matrix( rep(0, 4*seq_len), nrow = 1, byrow = T)
  # characters to match
  one_hot_encoded[ , 0*seq_len + which(givenSeq == "A")] <- 1
  one_hot_encoded[ , 1*seq_len + which(givenSeq == "C")] <- 1
  one_hot_encoded[ , 2*seq_len + which(givenSeq == "G")] <- 1
  one_hot_encoded[ , 3*seq_len + which(givenSeq == "T")] <- 1
  return (one_hot_encoded)
}


#' Title
#'
#' @param givenFastaSeqs
#'
#' @return
#'
#'
#' @examples
#'
get_one_hot_encoded_seqs <- function(givenFastaSeqs){
  # Input: list of FASTA seqs
  # Returns: A one-hot encoded feature matrix for the sequences
  seqs_split_as_list <- strsplit(givenFastaSeqs, split = NULL)
  encoded_seqs <- lapply(seqs_split_as_list, one_hot_encode)
  encoded_seqs <- do.call(rbind, encoded_seqs)
  return (encoded_seqs)
}


#' Title
#'
#' @param givenSeqs
#'
#' @return
#'
#'
#' @examples
#'
assert_attributes <- function(givenSeqs){
  # Check that all sequences are of same length
  length_vals <- as.factor(lapply(givenSeqs, length))
  if(length(levels(length_vals)) > 1){
    warning("Sequences are of different length\n")
    stop("Stopping")
  }
}


#' Title
#'
#' @param inputFastaFilename
#'
#' @return Returns one-hot-encoded sequences as matrix with dimensions #Features
#' x #Sequences
#' @export
#'
#' @examples
#'
prepare_data_from_FASTA <- function(inputFastaFilename){
  # Input: Complete FASTA filename with path
  #
  # Returns: data matrix
  start <- Sys.time()
  givenSeqs <- seqinr::read.fasta(inputFastaFilename, seqtype="DNA", as.string = TRUE)
  #
  assert_attributes(givenSeqs)
  #
  oheSeqs <- get_one_hot_encoded_seqs( toupper(givenSeqs) )
  print(Sys.time()-start)
  return (t(oheSeqs))
}



# testSeqs <- c("AACGTGACTA",
#               "ACCGATCGAT",
#               "GGCATCATGC",
#               "TCATCTAGAT"
#               )
# ohe_seqs <- one_hot_encoded_seqs(testSeqs)




