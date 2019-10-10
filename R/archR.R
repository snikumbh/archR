#' archR: A package for discovering different architectures of sequence
#' elements
#'
#' Given a set of DNA sequences, the \code{archR} package enables discovery of
#' clusters of sequence with characteristic sequence architectures.
#'
#'
#' The archR package provides three categories of important functions:
#' related to data preparation and manipulation, performing non-negative matrix
#' factorization, performing clustering, and visualization-related functions.
#'
#' @section Functions for data preparation and manipulation:
#' \itemize{
#' \item \code{\link{prepare_data_from_FASTA}}
#' \item \code{\link{make_sinuc_PWMs}}
#' \item \code{\link{sparsify_mat}}
#' }
#'
#' @section Functions related to performing NMF and model selection:
#' \itemize{
#' \item \code{\link{get_q2_using_py}}
#' \item \code{\link{compute_q2}}
#' \item \code{\link{cv_model_select_pyNMF}}
#' \item \code{\link{generate_folds}}
#' \item \code{\link{get_best_K}}
#' \item \code{\link{get_best_alpha}}
#' \item \code{\link{get_q2_aggregates_chosen_var}}
#' \item \code{\link{get_q2_threshold_by_K}}
#' }
#'
#' @section For related to clustering:
#' \itemize{
#' \item \code{\link{choose_clusters}}
#' \item \code{\link{get_clusters}}
#' }
#' @section Functions for visualizations:
#' \itemize{
#' \item \code{\link{plot_ggseqlogo}}
#' \item \code{\link{plot_ggheatmap}}
#' \item \code{\link{viz_all_factors_as_heatmap}}
#' \item \code{\link{viz_all_factors_as_seqlogo}}
#' \item \code{\link{viz_all_factors_in_combined_heatmaps_seqlogos}}
#' \item \code{\link{plot_arch_for_clusters}}
#' \item \code{\link{represent_matrix_of_acgt}}
#' \item \code{\link{plot_cv_K}}
#' \item \code{\link{plot_cv_Alpha}}
#' }
#'
#' @docType package
#' @name archR
NULL


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("archRconfig"))
