#' archR: A package for de novo discovery of different sequence 
#' architectures with characteristic combination(s) of sequence 
#' elements (motifs) or nucleotide compositions.
#'
#' Given a set of DNA sequences, \code{archR} enables unsupervised 
#' discovery of _de novo_ clusters with characteristic sequence 
#' architectures characterized by position-specific motifs or composition 
#' of stretches of nucleotides, e.g., CG-richness, etc.
#'
#' The archR package provides three categories of important functions:
#' related to data preparation and manipulation, performing non-negative 
#' matrix factorization, performing clustering, and visualization-related 
#' functions.
#'
#' @importFrom methods is
#'
#' @section Functions for data preparation and manipulation:
#' \itemize{
#' \item \code{\link{prepare_data_from_FASTA}}
#' \item \code{\link{get_one_hot_encoded_seqs}}
#' }
#'
#'
#' @section Functions for visualizations:
#' \itemize{
#' \item \code{\link{plot_ggseqlogo}}
#' \item \code{\link{plot_ggheatmap}}
#' \item \code{\link{viz_basis_vectors_as_heatmap}}
#' \item \code{\link{viz_basis_vectors_as_seqlogo}}
#' \item \code{\link{viz_basis_vectors_in_combined_heatmaps_seqlogos}}
#' \item \code{\link{plot_arch_for_clusters}}
#' plot_ggseqlogo_of_seqs
#' \item \code{\link{viz_seqs_as_acgt_mat_from_seqs}}
#' }
#'
#' @docType package
#' @name archR
NULL


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("archRconfig"))
