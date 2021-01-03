## NMF model selection functions II
## Other approaches for model selection


.stability_model_select_pyNMF2 <- function(X,
                                    param_ranges,
                                    parallelDo = FALSE,
                                    nCores = NA,
                                    nIterations = 100,
                                    bootstrap = TRUE,
                                    returnBestK = TRUE,
                                    tol = 10^-3,
                                    bound = 10^-6,
                                    flags = list(debugFlag = FALSE,
                                                verboseFlag = TRUE,
                                                plotFlag = FALSE,
                                                timeFlag = FALSE)
                                    ){
    ##
    sam_scores <- expand.grid(list(
        kValue = param_ranges$k_vals,
        ScoreType = c("AmariTypeDistance", "Dispersion"),
        nRuns = nIterations,
        Score = -0.01
    ))
    ##
    if(flags$debugFlag) message("Tolerance: ", tol, " & Bound : ", bound)
    ##
    bestK <- 1
    prev_amari <- NA
    for(kValue in param_ranges$k_vals){
        if(flags$verboseFlag) message("Checking K = ", kValue)
        resultList <- .perform_multiple_NMF_runs(X = X, kVal = kValue,
                                                    alphaVal = 0,
                                                    parallelDo = parallelDo,
                                                    nCores = nCores,
                                                    nRuns = nIterations,
                                                    bootstrap = bootstrap)
        featMatList <- lapply(resultList$nmf_result_list,
                                    function(x){
                                        get_features_matrix(x)
                                    })

        sampMatList <- lapply(resultList$nmf_result_list,
                                    function(x){
                                        get_samples_matrix(x)
                                    })
        ##
        if(bootstrap){
            sampMatListNew <- lapply(seq_len(length(sampMatList)),
                            function(x){
                                thisMat <- as.matrix(sampMatList[[x]])
                                thisNR <- nrow(thisMat)
                                thisNC <- ncol(thisMat)
                                origOrdX <- matrix(rep(-100,thisNR*thisNC),
                                                nrow = thisNR,
                                                ncol =  thisNC)
                                origOrdX[,resultList$new_ord[[x]]] <- thisMat
                                origOrdX
                            })
            sampMatList <- sampMatListNew

        }
        ##
        consensusMat <- getConsensusMat(sampMatList)
        ##
        this_disp <- NMF::dispersion(consensusMat)
        this_amari <- computeAmariDistances(featMatList)
        if(flags$debugFlag) message("AmariTypeDist : ", this_amari)

        sam_scores[sam_scores$kValue == kValue, "Score"] <-
            c(this_amari, this_disp)
        ##
        if(is.na(this_amari)){
            warning("NA in Amari-type distance computation")
        }
        ##
        if(this_amari == 0) {
            bestK <- 1
            break
        }
        ##
        if(this_amari > 0 && this_amari < bound){
            if(!is.na(prev_amari)){
                magnitudeChange <- abs(floor(log10(this_amari)) -
                                            floor(log10(prev_amari)))
                if(flags$debugFlag){
                    message("Change is : ", magnitudeChange)
                }
                if(magnitudeChange >= abs(log10(tol))){
                    ## detected fall based on tolerance, choose and break loop
                    bestK <- kValue - 1
                    if(flags$debugFlag) {
                        message("magnitude change higher than tolerance, ",
                            abs(log10(tol)))
                        message("This amariType Distance = ", this_amari)
                        message("Would have choosen bestK as : ", bestK)
                    }
                    #break
                }
            }
        }else{
            ## greater than bound, choose and break loop
            ## Bug: If this is reached at kValue = 1, bestK would be assigned 0.
            ## Avoid this.
            if(kValue > 1) bestK <- kValue - 1
            if(flags$debugFlag){
                message("This amariType Distance = ", this_amari)
            }
            if(flags$verboseFlag) {
                ## message("Choosing bestK as : ", bestK)
                ## message already provided by handle_chunk_w_NMF2
            }
            break
        }
        ##
        prev_amari <- this_amari
    }
    ##
    if(returnBestK) return(bestK)
    ##
    return(sam_scores)
}


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
