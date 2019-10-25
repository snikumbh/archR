## Function to check that the flags variable
## Expected to be:
## 1. not NULL
## 2. having the right fields and
## 3. all are LOGICAL values
.assert_archR_flags <- function(flags) {
    ## Expected fields
    expNames <-  c("debugFlag", "verboseFlag", "plotVerboseFlag", "timeFlag")
    if (is.null(flags)) {
        stop("'flags' are NULL")
    }
    if (!is.list(flags)) {
        stop("flags variable is not a list")
    } else {
        matchNames <- names(flags) %in% expNames
        if (is.null(names(flags)) && !all(matchNames)) {
            stop("Unexpected names or no elements in flags variable")
        } else {
            ## test all are set to logical values
            all_logical <- unlist(lapply(flags, is.logical))
            if (any(all_logical) == FALSE) {
                stop("Expected only LOGICAL values in flag variable,
                        found otherwise")
            }
        }
    }
}
## =============================================================================



## Function to check properties of samples matrix from NMF
## Expected to be
## - not NULL
## - A matrix
## - ncol > 0
## - nrows depends on number of factors of NMF
## - ? It has only 0s and 1s
.assert_archR_samplesMatrix <- function(samplesMatrix) {
    check_nrows <- 0
    if (is.null(samplesMatrix)) {
        stop("NULL value found, instead of a matrix")
    }
    if (!is.matrix(samplesMatrix)) {
        stop("Expected a matrix, found otherwise")
    } else {
        if (ncol(samplesMatrix) < 1) {
            stop("0 columns (sequences) in samplesMatrix")
        }
        ## assume sinuc
        if (nrow(samplesMatrix) == check_nrows) {
            stop("Check matrix, nrows == ", check_nrows)
        } else {
            # matElements <-
            # if () {
            #     stop("")
            # }
        }
    }

}
## =============================================================================


## Function to check properties of features matrix from NMF
## Expected to be
## - not NULL
## - A matrix
## - ncol > 0, nrow > 0, nrow %% 4 == 0
## - ncol depends on number of factors of NMF
## - ? It has only 0s and 1s
.assert_archR_featuresMatrix <- function(featuresMatrix) {
    check_ncols <- 0
    if (is.null(featuresMatrix)) {
        stop("NULL value found, instead of a matrix")
    }
    if (!is.matrix(featuresMatrix)) {
        stop("Expected a matrix, found otherwise")
    } else {
        if (any(is.na(featuresMatrix))) {
            stop("Factors have NA")
        }
        if (ncol(featuresMatrix) < 1) {
            stop("0 columns (sequences) in samplesMatrix")
        }
        if (ncol(featuresMatrix) == check_ncols) {
            stop("Check matrix, 'ncols' is: ", check_ncols)
        }
        if (all(featuresMatrix == 0)) {
            ## This will lead to an error if hopach is performed, hence throwing
            ## error
            stop("WARNING: All zeroes as factors")
        }
        if ((nrow(featuresMatrix) %% 4) != 0) {
            stop("#Rows in featuresMatrix not a multiple of 4")
        }
    }
}
## =============================================================================



## Function to independently check validity of minSeqs param in config
## Expected to be:
## 1. not NULL
## 2. numeric and > 0
.assert_archR_minSeqs_independent <- function(minSeqs_var) {
    if (is.null(minSeqs_var)) {
        stop("'minSeqs' is NULL")
    }
    if (!is.numeric(minSeqs_var)) {
        stop("'minSeqs' should be numeric and > 0")
    } else {
        if (minSeqs_var < 0) {
            stop("'misSeqs' should be > 0")
        }
    }
}
## =============================================================================



## Function to check validity of minSeqs w.r.t. given sequences
## Expected to be:
## 1. independently valid (see f .assert_archR_minSeqs_independent)
## 2. < #given sequences
##
.assert_archR_minSeqs_in_tandem <- function(minSeqs_var, given_seqs_size) {
    .assert_archR_minSeqs_independent(minSeqs_var)
    if (minSeqs_var > given_seqs_size) {
        stop("'minSeqs' is > number of input sequences.
            Typically, a small number")
    }
}
## =============================================================================



## Function to check validity of kFolds param in config
## Expected to be:
## 1. not NULL
## 2. numeric and > 0
.assert_archR_kFolds_independent <- function(kFolds_var) {
    if (is.null(kFolds_var)) {
        stop("'kFolds' is NULL")
    }
    if (!is.numeric(kFolds_var)) {
        stop("'kFolds' should be numeric and > 0")
    } else {
        if (kFolds_var < 1) {
            stop("'kFolds' should be > 0")
        }
    }
}
## =============================================================================



.assert_archR_kFolds_in_tandem <- function(kFolds_var, given_seqs_size) {
    if (is.null(kFolds_var)) {
        stop("'kFolds' is NULL")
    }
    if (!is.numeric(kFolds_var)) {
        stop("'kFolds' should be numeric and > 0")
    } else {
        if (kFolds_var < 1) {
            stop("'kFolds' should be > 0")
        }
        if (kFolds_var < 3) {
            stop("Set at least 3 cross-validation folds")
        } else if (kFolds_var > given_seqs_size) {
            stop("CV folds should be less than or equal to #sequences.\n
                Standard values: 3, 5, 10.")
        }
    }
}
## =============================================================================


## Function to check validity of parallelization-related param in config.
## parallelize and nCoresUse
## parallelize is expected to be:
## 1. not NULL
## 2. LOGICAL
## If parallelize is TRUE, nCoresUse is expected to be:
## 1. not NULL
## 2. numeric and > 0
.assert_archR_parallelization_params <- function(par_var, nCores_var) {
    if (is.null(par_var)) {
        stop("'parallelize' is NULL, expected LOGICAL")
    }
    if (!is.logical(par_var)) {
        stop("'parallelize' should be LOGICAL (TRUE/FALSE)")
    }
    ## nCoresUse, bother about it only when parallelize is TRUE
    if (par_var) {
        if (is.null(nCores_var)) {
            stop("'nCoresUse' is NULL")
        }
        if (!is.numeric(nCores_var) || nCores_var < 1) {
            stop("'nCoresUse' should be numeric and > 0 and < available #cores")
        }
        if (nCores_var > parallel::detectCores()) {
            stop("Specified more than available cores. Available cores: ",
                parallel::detectCores())
        }
    }
}
## =============================================================================



## Function to check properties of collection of chunks for next iteration
## Expected to be
## 1. not NULL
## 2. A list
## 3. all elements of non-zero/positive length
##
.assert_archR_OK_for_nextIteration <- function(nxtOuterChunksColl) {
    ## Check if it is NULL
    if (is.null(nxtOuterChunksColl)) {
        stop("Chunks for next iteration are NULL")
    }
    if (!is.list(nxtOuterChunksColl)) {
        stop("Chunks for next iteration expected as a list, found otherwise")
    }
    ## Check if any chunk in the collection is of zero-length?
    if (any(lapply(nxtOuterChunksColl, length) == 0)) {
        message("WARNING: Chunks for next iteration have a problem")
        stop(which(lapply(nxtOuterChunksColl, length) == 0),
                " have zero lengths")
    }

}
## =============================================================================



## Function to check properties of paramRanges used for model selection w/ NMF
## Expected to be
## 1. not NULL
## 2. List with elements having a fixed set of names
## 3. acceptable range for alpha: any thing non-negative
## 4. acceptable range for K_vals: any thing positive (> 0)
##
.assert_archR_model_selection_params <- function(param_ranges_var) {
    expNames <- c("alphaBase", "alphaPow", "k_vals")
    alphaBase <- param_ranges_var$alphaBase
    alphaPow <- param_ranges_var$alphaPow
    k_vals <- param_ranges_var$k_vals
    if (is.null(param_ranges_var)) {
        stop("'param_ranges' is NULL")
    }
    matchNames <- names(param_ranges_var) %in% expNames
    if (!all(matchNames)) {
        stop("Unexpected names of 'param_ranges' elements. Expected:", expNames)
    } else {
        if (is.numeric(k_vals) || is.numeric(alphaBase) ||
            is.numeric(alphaPow)) {
            ## all OK, check ranges now
            alphaVal <- alphaBase^alphaPow
            if (alphaVal < 0) {
                stop("Resulting alpha value is < 0. Check 'alphaBase' and
                    'alphaPow'")
            }
            if (!all(k_vals > 0)) {
                stop("'k_vals' should be > 0")
            }
        } else {
            stop("Either of 'k_vals', 'alphaBase' or 'alphaPow' is not a
                    numeric value")
        }
    }
}
## =============================================================================



## Function to check validity of nIterationsuse for NMF in config
## Expected to be:
## 1. not NULL
## 2. numeric and > 0
## 3. ?
.assert_archR_nIterations <- function(nIter_var) {
    if (is.null(nIter_var)) {
        stop("'nIterationsUse' is NULL")
    }
    if (!is.numeric(nIter_var)) {
        stop("'nIterationsUse' should be numeric and > 0")
    } else {
        if (nIter_var < 0) {
            stop("'nIterationsUse' should be > 0")
        }
    }
}
## =============================================================================




## Function to independently check validity of innerChunkSize
## Expected to be:
## 1. not NULL
## 2. numeric and >0
## When everything is satisfied, still whether it is valid depends on what is
## the total number of sequences
.assert_archR_innerChunkSize_independent <- function(innerChunkSize_var) {
    if (is.null(innerChunkSize_var)) {
        stop("'innerChunkSize' is NULL")
    }
    if (!is.numeric(innerChunkSize_var)) {
        stop("'innerChunkSize' should be numeric")
    } else {
        if (innerChunkSize_var <= 0) {
            stop("'innerChunkSize' should be > 0")
        }
    }
}
## =============================================================================



## Function to check validity of innerChunkSize w.r.t. given sequences
## Expected to be:
## 1. independently valid (see f .assert_archR_innerChunkSize_independent)
## 2. < #given sequences
##
.assert_archR_innerChunkSize_in_tandem <- function(innerChunkSize_var,
                                                given_seqs_size) {
    .assert_archR_innerChunkSize_independent(innerChunkSize_var)
    if (innerChunkSize_var > given_seqs_size) {
        stop("'innerChunkSize' should be <= number of input sequences")
    }
}
## =============================================================================



## Function to check properties of configuration variables
## Expected to be:
## 1. not NULL
## 1. A list with fixed set of names
## 2. Individual element assertions should pass
## 3.
.assert_archR_config <- function(config_var, given_seqs_size = NA) {
    expNames <- c("kFolds", "parallelize", "nCoresUse", "nIterationsUse",
                    "seedVal", "paramRanges", "innerChunkSize", "modSelLogFile",
                    "minSeqs", "flags")
    if (is.null(config_var)) {
        stop("'config' is NULL")
    }
    if (!is.list(config_var)) {
        stop("'config' is expected to be a list, found otherwise")
    } else {
        matchNames <- names(config_var) %in% expNames
        if (!all(matchNames)) {
            stop("Check names of elements in 'config'")
        }
        .assert_archR_kFolds_independent(config_var$kFolds)
        .assert_archR_innerChunkSize_independent(config_var$innerChunkSize)
        if (!is.na(given_seqs_size)) {
            .assert_archR_kFolds_in_tandem(config_var$kFolds, given_seqs_size)
            .assert_archR_innerChunkSize_in_tandem(config_var$innerChunkSize,
                                                    given_seqs_size)
            .assert_archR_minSeqs_in_tandem(config_var$minSeqs, given_seqs_size)
        }
        .assert_archR_minSeqs_independent(config_var$minSeqs)
        .assert_archR_model_selection_params(config_var$paramRanges)
        .assert_archR_parallelization_params(config_var$parallelize,
                                            config_var$nCoresUse)
        .assert_archR_nIterations(config_var$nIterationsUse)
        .assert_archR_flags(config_var$flags)
    }

}
## =============================================================================


## We need to write functions to sanity check the important variables
## throughout the archR algorithm. These are:
## -- globFactors Matrix
## -- globFactorsClustering Hopach Object
## -- globClustAssignments List
## -- outerChunk Matrix?
## -- outerChunksColl List
## -- innerChunksColl List
## -- intClustFactors Matrix
## -- clustFactors Matrix
## -- seqsClustLabels Vector
## -- collatedClustAssignments List
## --
##
## =============================================================================


## Function to check validity of NMFresult from .handle_chunk_w_NMF function
## Expected to be:
## 1. not NULL
## 2. a list w elements 'forGlobFactors' (a matrix) & 'forGlobClustAssignments'
## (a list)
## 3.
.assert_archR_NMFresult <- function(NMFresultObj) {
    if (is.null(NMFresultObj)) {
        stop("NMF result object is NULL")
    }
    expNames <- c("forGlobFactors", "forGlobClustAssignments")
    matchNames <- names(NMFresultObj) %in% expNames
    if (!all(matchNames)) {
        stop("Check names of elements in NMF result for inner chunk")
    }
    if (!is.list(NMFresultObj$forGlobFactors) &&
        !is.list(NMFresultObj$forGlobClustAssignments)) {
        stop("NMF result should hold information as list w/ elements per
                inner chunk")
    } else {
        .assert_archR_featuresMatrix(NMFresultObj$forGlobFactors)
        if (!is.list(NMFresultObj$forGlobClustAssignments)) {
            stop("'Global cluster assignments' in NMF result should be a list")
        } else {
            ## Check number of clusters == #factors (ncol(featuresMat))
            if (length(NMFresultObj$forGlobClustAssignments) !=
                ncol(NMFresultObj$forGlobFactors)) {
                    stop("In NMF result, #clusters != #factors")
            }
        }
    }
}
## =============================================================================


## Function to check validity of seqsClustLabels
## Expected to be:
## 1. not NULL
## 2. non-empty vector
## 3. all entries in the vector of same lengths
.assert_archR_seqsClustLabels <- function(seqsClustLabels) {
    if (is.null(seqsClustLabels)) {
        stop("Cluster labels for sequences NULL")
    }
    if (is.vector(seqsClustLabels) && length(seqsClustLabels) == 0) {
        stop("Expecting cluster labels for sequences as a non-empty vector")
    }
}

## Function to check validity of seqsClustLabels at the end of an iteration
## Expected to be:
## 1. general assertions passing
## 2. all entries in the vector of same lengths
.assert_archR_seqsClustLabels_at_end <- function(seqsClustLabels) {
    splitChar <- "-"
    one <- 1
    .assert_archR_seqsClustLabels(seqsClustLabels)
    all_lengths <- unlist(lapply(strsplit(seqsClustLabels,
                                            split = splitChar),
                                    length))
    check_length <- length(unique(all_lengths))
    if (check_length != one) {
        stop("Cluster labels for sequences have different lengths")
    }
}
## =============================================================================



## Function to check validity of hopach object for usability with archR
## Expected to be:
## 1. not NULL
## 2. a nested list w/ one element named 'clustering' which is also a list
## 3. The 'clustering' further has 3 elements: '$k', '$sizes' and '$order'
.assert_archR_hopachObj <- function(hopachObj, test_null = TRUE) {
    expNames <- c("clustering")
    expNames2 <- c("k", "sizes", "order")
    if (test_null) {
        if (is.null(hopachObj)) {
            stop("Hopach object is NULL")
        }
    }
    if (!is.list(hopachObj)) {
        stop("Check hopach object, need a list with element 'clustering'")
    } else if (is.list(hopachObj) && expNames %in% names(hopachObj)) {
        if (!all(expNames2 %in% names(hopachObj$clustering))) {
            stop("For archR, the element 'clustering' in hopach object should",
            " have elements named  'k', 'sizes' and 'order'")
        }
    }

}
## =============================================================================


## Function to check validity of lists in general
## Expected to be:
## 1. not NULL
## 2. a list
## 3. Not have any 0-length element
##
.assert_archR_list_properties <- function(listVar) {
    returnMessage <- "SAMARTH"
    if (is.null(listVar)) {
        returnMessage <- "NULL"
    }
    if (!is.list(listVar)) {
        returnMessage <- "Nlist"
    } else {
        check_lengths <- lapply(listVar, length)
        if (any(check_lengths == 0)) {
            returnMessage <- "0LengthEntry"
        }
    }
    return(returnMessage)
}
## =============================================================================



.assert_archR_globClustAssignments <- function(given_var) {
    returnMessage <- .assert_archR_list_properties(given_var)
    if (returnMessage == "NULL") {
        stop("Cluster assignments variable is NULL")
    }
    if (returnMessage == "Nlist") {
        stop("Cluster assignments variable is not a list")
    }
    if (returnMessage == "0LengthEntry") {
        stop("Cluster assignments variable has a 0-length entry")
    }
}
## =============================================================================


## Function to check consistency of nSeqs in clustLabels variable and in
## clustAssignments variable
## Expected to be:
## 1. Ensure list variable is OK and seqsClustLabels is OK
## 2. holding same number of sequences (the variable lengths)
##
.assert_archR_consistent_nSeqs_w_clusters <- function(seqsClustLabels,
                                                    clustAssignments) {

    .assert_archR_seqsClustLabels(seqsClustLabels)
    .assert_archR_globClustAssignments(clustAssignments)
    nSeqs_in_labels <- length(seqsClustLabels)
    nSeqs_in_assignments <- length(unlist(clustAssignments))
    if (!nSeqs_in_labels == nSeqs_in_assignments) {
        stop("Number of sequences in seqsClustLabels and clustAssignments not",
            " equal")
    }
}
## =============================================================================
