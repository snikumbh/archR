## NMF model selection functions II
## Other approaches for model selection
##
## Note on regularization: 
## We use dispersion for regularization together to Amari-type distance measure.
## The latter is a measure of stability of features, and the first of the 
## clusters.
##

.stability_model_select_pyNMF2 <- function(X,
                            param_ranges,
                            parallelDo = FALSE,
                            nCores = NA,
                            nIterations = 100,
                            bootstrap = TRUE,
                            returnBestK = TRUE,
                            ## TODO: Should tol stay as an argument?
                            tol = 10^-3,
                            bound = 10^-6,
                            flags = list(debugFlag = FALSE,
                                verboseFlag = TRUE, plotVerboseFlag = FALSE,
                                timeFlag = FALSE)
                            ){
    dbg <- flags$debugFlag
    vrbs <- flags$verboseFlag
    ##
    foo_scores <- expand.grid(list(kValue = param_ranges$k_vals,
                                    ScoreType = c("AmariTypeDistance"),
                                    nRuns = nIterations, Score = 100
                            ))
    ##
    .msg_pstr("Bound : ", bound, flg=dbg)
    ##
    bestK <- 1
    breakNow <- FALSE
    prev_amari <- NA
    for(kValue in param_ranges$k_vals){
        if(kValue > 1){
            .msg_pstr("Checking K =", kValue, flg=vrbs)
        }
        
        ## Get NMF results 
        resultList <- .perform_multiple_NMF_runs(X = X, kVal = kValue,
            alphaVal = 0, parallelDo = parallelDo,
            nCores = nCores, nRuns = nIterations,
            bootstrap = bootstrap)
        
        featMatList <- .get_feat_or_samp_matList(resultList, feat = TRUE)
        
        this_amari <- .get_amari_from_featMatList(featMatList)
        if(is.na(this_amari)){
            warning("NA in Amari-type distance computation")
        }
        .msg_pstr("AmariTypeDist : ", this_amari, flg=dbg)
        
        foo_scores[kValue, "Score"] <- this_amari
        
        
        if(this_amari > bound) breakNow <- TRUE
        ##
        check_errant <- FALSE
        if(kValue %% 5 == 0){
            check_errant <- .check_no_mag_change_fail_condition(foo_scores)
            if(check_errant) breakNow <- TRUE
        }
        ##
        if(breakNow){
            ## greater than bound, choose and break loop
            ## Here, if kValue = 1, bestK would be assigned 0. Avoid this.
            if(kValue > 1) bestK <- kValue - 1
            if(check_errant) bestK <- 1
            break
        }
        ##
        prev_amari <- this_amari
        
    }
    if(returnBestK) return(bestK)
    return(foo_scores)
}


.check_no_mag_change_fail_condition <- function(foo_scores){
    score_pows <- abs(floor(log10(foo_scores$Score)))
    score_pows <- score_pows[which(!is.nan(score_pows))]
    # print(score_pows)
    # diffs <- diff(score_pows)
    # print("Score_POWs")
    # print(diffs)
    # print(abs(diffs))
    if(all(score_pows >= 17)){
        return(TRUE)
    }
    return(FALSE)
}

.mag_change <- function(A, B){
    return(abs(floor(log10(A)) - floor(log10(B))))
}


.get_feat_or_samp_matList <- function(resultList, bootstrap, feat = TRUE){
    if(feat){
        featMatList <- lapply(resultList$nmf_result_list, get_features_matrix)
        return(featMatList)
    }
    sampMatList <- lapply(resultList$nmf_result_list, get_samples_matrix)
    
    if(bootstrap){
        sampMatListNew <-
            lapply(seq_len(length(sampMatList)), function(x){
                thisMat <- as.matrix(sampMatList[[x]])
                thisNR <- nrow(thisMat)
                thisNC <- ncol(thisMat)
                origOrdX <- matrix(rep(-100,thisNR*thisNC),
                    nrow = thisNR, ncol =  thisNC)
                origOrdX[,resultList$new_ord[[x]]] <- thisMat
                origOrdX
            })
        sampMatList <- sampMatListNew
    }
    return(sampMatList)
}

.tol_best_k <- function(kValue, this_amari, prev_amari, tol, verbose){
    magChange <- .mag_change(this_amari, prev_amari)
    .msg_pstr("Change is : ", magChange, flg=verbose)
    bestK <- NULL
    if(magChange >= abs(log10(tol))){
        ## detected fall based on tolerance, choose and break loop
        bestK <- kValue - 1
        .msg_pstr("magnitude change > tol, ", abs(log10(tol)), flg=verbose)
        .msg_pstr("This amariType Distance = ", this_amari, flg=verbose)
        .msg_pstr("Would have choosen bestK as : ", bestK, flg=verbose)
    }
    return(bestK)
}


# .get_coph_amari_disp_for_k <- function(X, kValue, parallelDo, nCores, nIterations, 
#     bootstrap#, amari = TRUE, coph = FALSE, disp = FALSE
#     ){
#     
#     resultList <- .perform_multiple_NMF_runs(X = X, kVal = kValue,
#         alphaVal = 0, parallelDo = parallelDo,
#         nCores = nCores, nRuns = nIterations,
#         bootstrap = bootstrap)
#     
#     # if(disp){
#         ## This is required when measures related to samplesMat are used.
#         ## Examples are dispersion and cophenetic correlation
#         ## 
#         # sampMatList <- lapply(resultList$nmf_result_list, function(x){
#         #     get_samples_matrix(x)
#         # })
#         ##
#         # disp_val <- .get_dispersion_from_sampleMatList(sampMatList)
#     # }
#     
#     # if(amari){
#         ## Required for Amari-type distance
#         featMatList <- lapply(resultList$nmf_result_list, function(x){
#             get_features_matrix(x)
#         })
#         amari_val <- .get_amari_from_featMatList(featMatList)
#     # }
#     # if(coph){
#         ## This is required when measures related to samplesMat are used.
#         ## Examples are dispersion and cophenetic correlation
#         ## 
#         sampMatList <- lapply(resultList$nmf_result_list, function(x){
#             get_samples_matrix(x)
#         })
#         if(bootstrap){
#             sampMatListNew <-
#                 lapply(seq_len(length(sampMatList)), function(x){
#                     thisMat <- as.matrix(sampMatList[[x]])
#                     thisNR <- nrow(thisMat)
#                     thisNC <- ncol(thisMat)
#                     origOrdX <- matrix(rep(-100,thisNR*thisNC),
#                         nrow = thisNR, ncol =  thisNC)
#                     origOrdX[,resultList$new_ord[[x]]] <- thisMat
#                     origOrdX
#                 })
#             sampMatList <- sampMatListNew
#         }
#         ##
#         coph_val <- .get_cophcor_from_sampleMatList(sampMatList)
#         disp_val <- .get_dispersion_from_sampleMatList(sampMatList)
#     # }
#     return(c(amari_val, coph_val, disp_val))
# }

# .get_dispersion_from_sampleMatList <- function(sampMatList){
#     ## Not required for Amari-type distance, but cophcor or dispersion
#     consensusMat <- getConsensusMat(sampMatList)
#     ##
#     return(NMF::dispersion(consensusMat))
# }


# .get_cophcor_from_sampleMatList <- function(sampMatList){
#     ## Not required for Amari-type distance, but cophcor or dispersion
#     consensusMat <- getConsensusMat(sampMatList)
#     ##
#     return(NMF::cophcor(consensusMat))
# }

.get_amari_from_featMatList <- function(featMatList){
    return(computeAmariDistances(featMatList))
}

# .get_cophcor_for_k <- function(X, kValue, parallelDo, nCores, nIterations, 
#     bootstrap){
#     resultList <- .perform_multiple_NMF_runs(X = X, kVal = kValue,
#         alphaVal = 0, parallelDo = parallelDo,
#         nCores = nCores, nRuns = nIterations,
#         bootstrap = bootstrap)
#     ## Required for Amari-type distance
#     # featMatList <- lapply(resultList$nmf_result_list, function(x){
#     #     get_features_matrix(x)
#     # })
#     ## This is required when measures related to samplesMat are used.
#     ## Examples are dispersion and cophenetic correlation
#     ## 
#     sampMatList <- lapply(resultList$nmf_result_list, function(x){
#         get_samples_matrix(x)
#     })
#     ##
#     if(bootstrap){
#         sampMatListNew <-
#             lapply(seq_len(length(sampMatList)), function(x){
#                 thisMat <- as.matrix(sampMatList[[x]])
#                 thisNR <- nrow(thisMat)
#                 thisNC <- ncol(thisMat)
#                 origOrdX <- matrix(rep(-100,thisNR*thisNC),
#                     nrow = thisNR, ncol =  thisNC)
#                 origOrdX[,resultList$new_ord[[x]]] <- thisMat
#                 origOrdX
#             })
#         sampMatList <- sampMatListNew
#     }
#     ## Not required for Amari-type distance, but cophcor or dispersion
#     consensusMat <- getConsensusMat(sampMatList)
#     ##
#     return(NMF::cophcor(consensusMat))
# }

# .get_amari_for_k <- function(X, kValue, parallelDo, nCores, nIterations, 
#                             bootstrap){
#     resultList <- .perform_multiple_NMF_runs(X = X, kVal = kValue,
#         alphaVal = 0, parallelDo = parallelDo,
#         nCores = nCores, nRuns = nIterations,
#         bootstrap = bootstrap)
#     ##
#     featMatList <- lapply(resultList$nmf_result_list, function(x){
#         get_features_matrix(x)
#     })
#     ## This is required when measures related to samplesMat are used.
#     ## Examples are dispersion and cophenetic correlation
#     ## 
#     # sampMatList <- lapply(resultList$nmf_result_list, function(x){
#     #     get_samples_matrix(x)
#     # })
#     # ##
#     # if(bootstrap){
#     #     sampMatListNew <-
#     #         lapply(seq_len(length(sampMatList)), function(x){
#     #             thisMat <- as.matrix(sampMatList[[x]])
#     #             thisNR <- nrow(thisMat)
#     #             thisNC <- ncol(thisMat)
#     #             origOrdX <- matrix(rep(-100,thisNR*thisNC),
#     #                 nrow = thisNR, ncol =  thisNC)
#     #             origOrdX[,resultList$new_ord[[x]]] <- thisMat
#     #             origOrdX
#     #         })
#     #     sampMatList <- sampMatListNew
#     # }
#     ## Not required for Amari-type distance
#     # consensusMat <- getConsensusMat(sampMatList)
#     ##
#     return(computeAmariDistances(featMatList))
# }

amariDistance <- function(matA, matB) {
    K <- dim(matA)[2]
    corrMat <- stats::cor(matA, matB)
    return(1 - (sum(apply(corrMat, 1, max)) + 
            sum(apply(corrMat, 2, max))) / (2 * K))
}


computeAmariDistances <- function(matrices){
    B <- length(matrices)
    distances.list <- unlist(lapply(seq_len(B - 1), function(b) {
        distances <- lapply(seq(from=b + 1, to=B, by=1), function(b.hat) {
            amariDistance(matrices[[b]], matrices[[b.hat]])
        })
    })
    )
    return(mean(distances.list))
}


getConsensusMat <- function(matrices){
    memberships <- lapply(matrices, function(x){
        apply(x, 2, which.max)
    })
    ## compute connectivity matrix per run
    ## finally, compute the consensus matrix (average of all connectivity 
    ## matrices)
    ## Compute the dispersion score (range: -1 to +1)
    connectivityMats <- lapply(seq_along(matrices), function(x){
        xMat <- matrices[[x]]
        xMemberships <- memberships[[x]]
        nSamples <- ncol(xMat)
        connMat <- matrix(rep(0, nSamples*nSamples), nrow = nSamples)
        for(i in seq_len(nrow(connMat))){
            relColIdx <- which(xMemberships == xMemberships[i])
            connMat[i,relColIdx] <- 1
        }
        connMat
    })
    ##
    consensusMat <- .getMeanOfListOfMatrices(connectivityMats)
}
