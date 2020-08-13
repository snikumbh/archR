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
    ## Use cophenetic correlation cofficient after multiple runs to select the
    ## best model
    ##
    sam_scores <- purrr::cross_df(list(
        kValue = param_ranges$k_vals,
        # ScoreType = c("CophCorrelation", "AmariTypeDistance", "Dispersion"),
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
            sampMatListNew <- lapply(1:length(sampMatList),
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
            c(this_amari,  this_disp)
        ##

        ##
        if(this_amari < bound){
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
                        message("magnitude change higher than tolerance which is ",
                            abs(log10(tol)))
                        message("This amariType Distance = ", this_amari)
                        message("Would have choosen bestK as : ", bestK)
                    }
                    #break
                }
            }
        }else{
            ## greater than bound, choose and break loop
            bestK <- kValue - 1
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
    # if(flags$debugFlag) {sam_scores %>% print(n = Inf)}
    if(returnBestK){
        ##
        return(bestK)
    }
    ##
    return(sam_scores)
}


amariDistance <- function(matrix.A, matrix.B) {
    K <- dim(matrix.A)[2]
    C <- cor(matrix.A, matrix.B)
    return(1 - (sum(apply(C, 1, max)) + sum(apply(C, 2, max))) / (2 * K))
}


computeAmariDistances <- function(matrices){
    B <- length(matrices)
    distances.list <- unlist(lapply(1:(B - 1), function(b) {
        distances <- lapply((b + 1):B, function(b.hat) {
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
    ## finally, compute the consensus matrix (average of all connectivity matrices)
    ## Compute the dispersion score (range: -1 to +1)
    connectivityMats <- lapply(seq_along(matrices), function(x){
        xMat <- matrices[[x]]
        xMemberships <- memberships[[x]]
        nSamples <- ncol(xMat)
        connMat <- matrix(rep(0, nSamples*nSamples), nrow = nSamples)
        for(i in 1:nrow(connMat)){
            relColIdx <- which(xMemberships == xMemberships[i])
            connMat[i,relColIdx] <- 1
        }
        connMat
    })
    ##
    consensusMat <- .getMeanOfListOfMatrices(connectivityMats)
}