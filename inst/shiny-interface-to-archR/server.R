
# Define server logic to process input ----
server <- function(input, output, session) {


getConfigDF <- function(namesList, valuesList, ifDefaultList){
    if(!is.null(ifDefaultList)){
        return(data.frame(
            ParameterName = namesList,
            ParameterValue = valuesList,
            IfDefault = unlist(ifDefaultList),
            stringsAsFactors = FALSE))
    }else{
        return(data.frame(
            ParameterName = namesList,
            ParameterValue = valuesList,
            stringsAsFactors = FALSE))
    }
}
##==========================================================================

setParallelFlag <- function(parallelFlag){
    parallelSetting <- FALSE
    if(!is.null(parallelFlag)){
        if(parallelFlag){
            parallelSetting <- TRUE
        }
    }
    return(parallelSetting)
}
##==========================================================================
setCheckpointFlag <- function(checkpointFlag){
    checkpointSetting <- FALSE
    if(!is.null(checkpointFlag)){
        if(checkpointFlag){
            checkpointSetting <- TRUE
        }
    }
    return(checkpointSetting)
}
##==========================================================================

setArchRConfig <- function(input){
    ##
    if(!is.null(input$parallelize)){
        parallelSetting <- setParallelFlag(input$parallelize)
    }
    ##
    if(!is.null(input$checkpoint)){
        checkpointSetting <- setCheckpointFlag(input$checkpoint)
    }

    fl <- flagsLogicalVector(input$flags)

    archRconfig <- archR::archRSetConfig(
        innerChunkSize = as.numeric(input$innerChunkSize),
        cvFolds = as.numeric(input$cvFolds),
        parallelize = parallelSetting,
        nCoresUse = as.numeric(input$nCoresUse),
        nIterationsUse = as.numeric(input$nIterationsUse),
        kMin = 1,
        kMax = as.numeric(input$kMax),
        modSelType = input$modSelType,
        tol = 10^as.numeric(input$tol),
        bound = 10^as.numeric(input$bound),
        checkpointing = checkpointSetting,
        modSelLogFile = "log.txt",
        flags = list(debugFlag = fl$debugFlag,
                     verboseFlag = fl$verboseFlag,
                     plotVerboseFlag = fl$plotVerboseFlag,
                     timeFlag = fl$timeFlag)
    )

    return(archRconfig)
}
##==========================================================================


getConfigDFParamNames <- function(){
    allParamNamesList <- c("Chunk size",
                           "# Factors (range)",
                           "Model selection by",
                           "Tolerance",
                           "Bound",
                           "Cross-validation folds",
                           "# Iterations",
                           "Parallelize?",
                           "# Cores",
                           "Checkpointing",
                           "# Flags")
    return(allParamNamesList)
}
##==========================================================================
##

defaultValues <- function(varname = NULL){
    defaultValues <- list("innerChunkSize" = 500,
                          "kMax" = 20,
                          "modSelType" = "stability",
                          "tolerance" = -3,
                          "bound" = -8,
                          "cvFolds" = 5,
                          "nIterationsUse" = 100)
    return(defaultValues[[varname]])
}

getConfigDFParamValues <- function(input, parallelSetting,
                                   checkpointSetting){
    allParamValuesList <- as.character(
        c(input$innerChunkSize,
          paste(1, input$kMax, sep = " to "),
          input$modSelType,
          10^as.numeric(input$tol),
          10^as.numeric(input$bound),
          input$cvFolds,
          input$nIterationsUse,
          parallelSetting,
          input$nCoresUse,
          checkpointSetting,
          paste(input$flags, collapse=", "))

    )
    ifDefaultList <- lapply(seq_len(11), function(i){
                            if(i <= 7) {""}
                            else{NA}
                        })
    names(ifDefaultList) <- getConfigDFParamNames()
    if(input$innerChunkSize == defaultValues("innerChunkSize")){
        ifDefaultList[[1]] <- "Y"
    }
    if(as.numeric(input$kMax) == defaultValues("kMax")){
        ifDefaultList[[2]] <- "Y"
    }
    if(input$modSelType == defaultValues("modSelType")){
        ifDefaultList[[3]] <- "Y"
    }
    if(as.numeric(input$tol) == defaultValues("tolerance")){
        ifDefaultList[[4]] <- "Y"
    }
    if(as.numeric(input$bound) == defaultValues("bound")){
        ifDefaultList[[5]] <- "Y"
    }
    if(as.numeric(input$cvFolds) == defaultValues("cvFolds")){
        ifDefaultList[[6]] <- "Y"
    }
    if(as.numeric(input$nIterationsUse) == defaultValues("nIterationsUse")){
        ifDefaultList[[7]] <- "Y"
    }
    # message(length(allParamValuesList))
    # message(length(ifDefaultList))
    stopifnot(length(ifDefaultList) == length(allParamValuesList))
    return(list(allParamValuesList, ifDefaultList))
}
##
setChooseIdx <- function(parallelSetting, modSelType, totalLength){
    stabilityIdx <- c(4,5)
    cvIdx <- c(6)
    parIdx <- c(9)
    chooseIdx <- rep(TRUE, totalLength)
    if(!parallelSetting) chooseIdx[parIdx] <- FALSE
    if(modSelType == "cv"){
        chooseIdx[stabilityIdx] <- FALSE
    }else{
        chooseIdx[cvIdx] <- FALSE
    }
    return(chooseIdx)
}
##==========================================================================

configSummaryTable <- reactive({

    req(input$performAction)

    if(is.null(input$performAction)){
        return(NULL)
    }else if(input$performAction == 'processNew'){
        ##
        parallelSetting <- setParallelFlag(input$parallelize)
        ##
        checkpointSetting <- setCheckpointFlag(input$checkpoint)
        ##
        # fl <- flagsLogicalVector(input$flags)
        ##
        allParamNamesList <- getConfigDFParamNames()
        temp <- getConfigDFParamValues(input, parallelSetting,
                                       checkpointSetting)
        allParamValuesList <- temp[[1]]
        ifDefaultList <- temp[[2]]
        chooseIdx <- setChooseIdx(parallelSetting, input$modSelType,
                                  length(allParamNamesList))
        ##
        df <- getConfigDF(allParamNamesList[chooseIdx],
                          allParamValuesList[chooseIdx],
                          ifDefaultList[chooseIdx])
        df
    }else{
        ##
        req(input$inputRDSResultFilename)
        ## TODO: Show error if not RDS file
        archRresult <- getResult(input$inputRDSResultFilename$datapath)
        ## if archRresult object is NULL, reading resuled in err/warning
        if(!is.null(archRresult)){
            allParamNamesList <- getConfigDFParamNames()
            ##
            allParamValuesList <- as.character(
                c(archRresult$config$innerChunkSize,
                  paste(min(archRresult$config$paramRanges$k_vals),
                        max(archRresult$config$paramRanges$k_vals),
                        sep = " to "),
                  archRresult$config$modSelType,
                  archRresult$config$tol,
                  archRresult$config$bound,
                  archRresult$config$kFolds,
                  archRresult$config$nIterationsUse,
                  archRresult$config$parallelize,
                  archRresult$config$nCoresUse,
                  archRresult$config$checkpointing,
                  paste(names(which(archRresult$config$flags == TRUE)),
                        collapse=", ")
                )
            )
            ##
            chooseIdx <- setChooseIdx(archRresult$config$parallelize,
                                      archRresult$config$modSelType,
                                      length(allParamNamesList))

            df <- getConfigDF(allParamNamesList[chooseIdx],
                              allParamValuesList[chooseIdx], NULL)
            return(df)
        }else{
            return(NULL)
        }
    }
})
##==========================================================================

flagsLogicalVector <- function(ipfl){
    setDebugFlag <- FALSE
    setVerboseFlag <- FALSE
    setTimeFlag <- FALSE
    setPlotFlag <- FALSE
    ##
    if(!is.null(ipfl) && length(ipfl) > 0){
        for(i in seq_along(ipfl)){
            if(ipfl[i] == "debugFlag"){
                setDebugFlag <- TRUE
            }else if(ipfl[i] == "verboseFlag"){
                setVerboseFlag <- TRUE
            }else if(ipfl[i] == "timeFlag"){
                setTimeFlag <- TRUE
            }else if(ipfl[i] == "plotVerboseFlag"){
                setPlotFlag <- TRUE
            }
        }
    }
    ##
    return(list(debugFlag = setDebugFlag,
                verboseFlag = setVerboseFlag,
                timeFlag = setTimeFlag,
                plotVerboseFlag = setPlotFlag))
}
##==========================================================================

getResult <- function(fname){
    expFormatStr <- "archR result file is expected to be an RDS file"
    archRresult <- NULL
    archRresult <- tryCatch({
        readRDS(fname)
    }, error = function(err){
        if(err$message == "unknown input format"){
            thisMsg <- paste("Wrong file format.",
                             expFormatStr)
            showModal(modalDialog(
                title = "Error",
                thisMsg,
                easyClose = TRUE
            ))
        }else{
            thisMsg <- paste("Unexpected error in reading file.",
                             expFormatStr)
            showModal(modalDialog(
                title = "Error",
                thisMsg,
                easyClose = TRUE
            ))
        }
    }, warning = function(warn){
        if(warn$message == "unknown input format"){
            thisMsg <- paste("Wrong file format.",
                             expFormatStr)
            showModal(modalDialog(
                title = "Error",
                thisMsg,
                easyClose = TRUE
            ))
        }else{
            thisMsg <- paste("Unexpected error in reading file.",
                             expFormatStr)
            showModal(modalDialog(
                title = "Error",
                thisMsg,
                easyClose = TRUE
            ))
        }
    }, finally = {
        ##
    })
    return(archRresult)
}
##==========================================================================

summarizeRDSinTable <- reactive({
    req(input$inputRDSResultFilename)
    archRresult <- getResult(input$inputRDSResultFilename$datapath)
    nFactors <- length(archRresult$seqsClustLabels)
    factorNames <- paste0("# Factors (", paste0("Iteration",
                                               c(1:nFactors), ")"))
    factorValues <- unlist(lapply(1:nFactors, function(n){
        archRresult$clustBasisVectors[[n]]$nBasisVectors
    }))

    statNames <- c("# Iterations", factorNames)
    statValues <- c(nFactors, factorValues)
    data.frame(Name = statNames,
               Values = statValues,
               stringsAsFactors = FALSE
    )
})
##==========================================================================

getFasta <- function(fname){
    expFormatStr <- "Input file is expected to be a valid FASTA file"
    raw_seq <- NULL
    raw_seq <- tryCatch({
        archR::prepare_data_from_FASTA(fname,
                                       sinuc_or_dinuc = "dinuc",
                                       rawSeq = TRUE)
    }, error = function(err){
        if(err$message == "unknown input format"){
            thisMsg <- paste("Wrong file format.",
                             expFormatStr)
            showModal(modalDialog(
                title = "Error",
                thisMsg,
                easyClose = TRUE
            ))
        }else{
            thisMsg <- paste("Unexpected error in reading file.",
                             expFormatStr)
            showModal(modalDialog(
                title = "Error",
                thisMsg,
                easyClose = TRUE
            ))
        }
    }, warning = function(warn){
        if(warn$message == "unknown input format"){
            thisMsg <- paste("Wrong file format.",
                             expFormatStr)
            showModal(modalDialog(
                title = "Error",
                thisMsg,
                easyClose = TRUE
            ))
        }else{
            thisMsg <- paste("Unexpected error in reading file.",
                             expFormatStr)
            showModal(modalDialog(
                title = "Error",
                thisMsg,
                easyClose = TRUE
            ))
        }
    }, finally = {
        ##
    })
    return(raw_seq)

}
##==========================================================================

summarizeFASTAinTable <- reactive({
    req(input$inputFastaFilename)
    fname <- input$inputFastaFilename$datapath
    raw_seq <- getFasta(fname)
    ## summary as a data.frame
    if(!is.null(raw_seq)){
        seqStatNames <- c("# Sequences", "Sequence length(s)", "Alphabet")
        nSeqs <- length(raw_seq)
        ## A holds sequence lengths/widths
        A <- unique(Biostrings::width(raw_seq))
        ## B holds summary of quantiles
        B <- summary(A)
        seqLens <- A # if(length(A) > 1){B[c(1,3,6)]}else{A}
        ## A holds the alphabetFreq matrix
        ## expected dimensions nSeqs x 5 when DNA/RNA w/ baseOnly=T
        A <- Biostrings::alphabetFrequency(raw_seq, baseOnly = TRUE)
        if(nSeqs == length(which(A[,5] == 0))){
            seqChars <- paste0(colnames(A)[1:4], collapse=",")
        }else{
            seqChars <- paste0(colnames(A), collapse=",")
        }

        seqStatValues <- c(nSeqs, paste(seqLens,collapse=", "), seqChars)
        data.frame(Name = seqStatNames,
                   Values = seqStatValues,
                   stringsAsFactors = FALSE
        )
    }
})
##==========================================================================


jobStatusTable <- function(sl_job, submitTorF){
    if(submitTorF){
        statusText <- archR_get_job_status(slr_job = sl_job)
        return(statusText$queue)
    }else{
        return(NULL)
    }
}
##==========================================================================


observeEvent(input$checkStatus, {
    ##
    output$statusTable <- renderTable({
        jobStatusTable(sl_job = global_sl_job, input$submitAsk)
    })

})
##==========================================================================

observeEvent(input$cancelJob, {
    ##
    rslurm::cancel_slurm(global_sl_job)
    ##TO-DO: Update job status?
    Sys.sleep(time = 2)
    output$statusTable <- renderTable({
        jobStatusTable(sl_job = global_sl_job, input$submitAsk)
    })
})
##==========================================================================




##
# hideTab(inputId = "archRTabs", target = "Result")
# hideTab(inputId = "archRTabs", target = "Response")
hideTab(inputId = "archRTabs", target = "Status")
hideTab(inputId = "archRTabs", target = "OutputLog")
# hideTab(inputId = "archRTabs", target = "OutputLog2")
##
observeEvent(input$callArchR, {
    # showTab(inputId = "archRTabs", target = "Response")
    showNotification(type = "message", duration = 15,
                     ui = "Head to the Response tab in the Summary box"
    )
    if(input$parallelize){
        if(input$useSlurm){
            if(input$submitAsk){
                showTab(inputId = "archRTabs", target = "Status")
            }
        }
    }
    ##
    updateTabsetPanel(session, inputId = "archRTabs", selected = "Response")
    withCallingHandlers({
        shinyjs::html("responseText", "")
        archRconfig <- setArchRConfig(input)
        ##
        if(is.null(input$inputFastaFilename)){
            # message("Select FASTA file")
            showModal(modalDialog(
                title = "Error",
                "Please provide a FASTA file",
                easyClose = TRUE
            ))
        }else{
        ##
        req(input$inputFastaFilename)

        fname <- input$inputFastaFilename$datapath
        tssSeqs_raw <- getFasta(fname)
        if(!is.null(tssSeqs_raw)){
        tssSeqs <- archR::prepare_data_from_FASTA(fname,
                                    sinuc_or_dinuc = "dinuc")
        nSeqs <- ncol(tssSeqs)
        ##
        positions <- seq(1, Biostrings::width(tssSeqs_raw[1]))
        ##
        if(archRconfig$parallelize){
            if(!is.null(input$useSlurm)){
                if(input$useSlurm){
                    job_dir <- archR::handle_dir_creation(
                        file.path(input$jobLocation,
                               paste0("_rslurm_", input$jobName)
                               )
                        )
                    # message(job_dir)
                    ## create dummy dir and file
                    dfname <- paste0(job_dir, "slurm_0.out")
                    write("", dfname)
                    # message("Out file written: ", dfname)
                    logfilename <- dfname
                    # message("Logfilename from dfname:", logfilename)
                    # logfilename2 <<- paste0(job_dir, "slurm_0.out")
                    # message("Logfilename2:", logfilename2)
                    # urmprint(Sys.time())
                    ####
                    temp <- unlist(strsplit(job_dir, split = "_rslurm_"))
                    manipulatedJobname <-
                         unlist(strsplit(temp[length(temp)], split="/"))[1]
                    # input$jobName
                    message("Job name: ", manipulatedJobname)
                    message("Slurm scripts destination: ", job_dir)
                    # message("Ready to submit job")
                    # Sys.sleep(time=60)
                    ## setup rslurm slurm_call
                    sl_job <-
                        archR_slurm_call(dir_path = input$jobLocation,
                        archR::archR,
                        list(config = archRconfig,
                            seqsMat = tssSeqs,
                            seqsRaw = tssSeqs_raw,
                            seqsPositions = positions,
                            thresholdItr = as.numeric(input$thresholdsItr),
                            oDir = input$resultLocation),
                       jobname = manipulatedJobname, #input$jobName,
                       r_template = system.file("templates",
                                                "archR_r_template.txt",
                                                package = "archR"),
                       ##
                       ## using custom template, so that we can call
                       ## load_all()
                       ## This is required when archR is not installed
                       ## (basically only during development)
                       pkgs = c("reticulate"),
                       slurm_options =
                           list("cpus-per-task" = as.numeric(input$nCoresUse),
                                "mem-per-cpu" = input$memPerCPU,
                                "nodes" = as.numeric(input$nodes)
                                ),
                       submit = input$submitAsk)
                    ##
                    if(!input$submitAsk){
                        message(paste0("Submission scripts output in",
                                       " directory:\n", job_dir))
                        # message(paste0("Submission scripts output in",
                        #                " directory _rslurm_", input$jobName))
                    }else{

                        statusText <- archR_get_job_status(sl_job)
                        # message(paste0("Submitted batch job ",
                        #                statusText$queue$JOBID))
                        ##
                        output$statusTable <- renderTable({
                            global_sl_job <<- sl_job
                            jobStatusTable(sl_job, input$submitAsk)
                        })
                        #### Reading slurm output log file ####

                        showTab(inputId = "archRTabs", target = "OutputLog")
                        ## From examples in shiny gallery
                        pollData <- reactivePoll(500, session,
                                        ## This function returns the time that
                                        ## the logfile was last modified
                                        checkFunc = function(){
                                         if (file.exists(logfilename)){
                                             file.info(logfilename)$mtime[1]
                                         }else{""}
                                        },
                                        ## This function returns the content of
                                        ## the logfile
                                        valueFunc = function(){
                                         readLines(logfilename)
                                        }
                        )

                        ##
                        output$logText <- renderText({
                            ## Read the text, we use a wellPanel to
                            ## restrict height
                            # logfilename <<- paste0(job_dir, "slurm_0.out")
                            text <- pollData()
                            paste(text, collapse = '\n')
                        })
                        #### Reading slurm output log file Ends ####
                        message("`Status' tab shows the slurm job status")
                        message("Slurm job output can be monitored in the `OutputLog' tab")
                    } ##ifElse "submit asked?" ends here
                }else{
                    ## Create cluster and run parallely directly
                    ## Have a short wrapper function do this outside shiny
                    message("Running archR in parallel")
                    archR_result <- archR::archR(config = archRconfig,
                                 seqsMat = tssSeqs,
                                 seqsRaw = tssSeqs_raw,
                                 seqsPositions = positions,
                                 thresholdItr = as.numeric(input$thresholdsItr),
                                 oDir = input$resultLocation
                                )
                }}
        }else{
            message("Running archR serially")
            ##
            archR_result <- archR::archR(config = archRconfig,
                            seqsMat = tssSeqs,
                            seqsRaw = tssSeqs_raw,
                            seqsPositions = positions,
                            thresholdItr = as.numeric(input$thresholdsItr),
                            oDir = input$resultLocation
                            )
            ##
        }
        } ## ifElse "error in reading provided FASTA file" ends here
        } ##ifElse "FASTA file not provided" ends here
    }, ## withCallingHandlers ends
    message = function(m) {
        shinyjs::html(id = "responseText", html = m$message, add = TRUE)
    },
    warning = function(m) {
        shinyjs::html(id = "responseText", html = m$message, add = TRUE)
    })
})
##==========================================================================

# output$pollText <- renderText({
#     # Read the text, and make it a consistent number of lines so
#     # that the output box doesn't grow in height.
#     text <- pollData()
#     length(text) <- 14
#     text[is.na(text)] <- ""
#     paste(text, collapse = '\n')
# })
##==========================================================================




# output$infoBoxProcessNew <- renderInfoBox({
#
#     shinydashboard::infoBox(title = "Help",
#             value = paste0("You can now run archR with default parameter settings",
#             "or customize them in the 'Parameter Settings' tab ",
#             "and then run archR"), icon = shiny::icon("hand-point-right"),
#             color = "purple")
# })

observeEvent(input$reset, {
    # str(input$resetSelect)
    if(any("all" == input$resetSelect)){
        resetAll <- TRUE
        # message("Updating  all")
    }
    if(resetAll || any("innerChunkSize" == input$resetSelect)){
        # message("Updating  chunk size")
        updateNumericInput(session, inputId = "innerChunkSize", value = 500)
    }
    if(resetAll || any("kMax" == input$resetSelect)){
        # message("Updating  kMax")
        updateNumericInput(session, inputId = "kMax", value = 20)
    }
    if(resetAll || any("modSelType" == input$resetSelect)){
        # message("Updating  modSelType")
        updateRadioButtons(session, inputId = "modSelType",
                            selected = "stability")
    }
    if(resetAll || any("tol" == input$resetSelect)){
        # message("Updating  tolerance")
        updateSliderInput(session, "tol", value = -3)
    }
    if(resetAll || any("bound" == input$resetSelect)){
        # message("Updating  bound")
        updateSliderInput(session, "bound", value = -8)
    }
    if(resetAll || any("cvFolds" == input$resetSelect)){
        # message("Updating  cv folds")
        updateSelectInput(session, inputId = "cvFolds", selected = 5)
    }
    if(resetAll || any("nIterationsUse" == input$resetSelect)){
        # message("Updating  nIterationsIse")
        updateNumericInput(session, inputId = "nIterationsUse",
                           value = 100)
    }
    if(resetAll || any("thresholdsItr" == input$resetSelect)){
        # message("Updating  thresholdsItr")
        updateNumericInput(session, inputId = "thresholdsItr",
                            value = 3)
    }
    if(resetAll || any("resultLocation" == input$resetSelect)){
        # message("Updating resultLocation")
        updateTextInput(session, inputId = "resultLocation",
                        value = "./archR_run")
    }
    # showNotification(type = "message", duration = 5,
    #                  ui = "All")
})


observeEvent(input$inputFastaFilename,{
    showNotification(type = "message", duration = 15,
         ui = paste0("Next: Set custom parameter values OR run archR with ",
                 "default parameter settings ",
                 "(See 'Parameter Settings' tab)")
    )
})


# #Return the summary of configuration for printing ----
output$configSummary <- renderTable({
    configSummaryTable()
})
##==========================================================================

# #Return the summary of configuration for printing ----
output$dataSummary <- renderTable({
    if(!is.null(input$performAction)){
        if(input$performAction == "processNew"){
            #cat("FASTA sequences read as a
            #DNAStringSet object for viewing:\n")
            summarizeFASTAinTable()
        }else{
            summarizeRDSinTable()
        }
    }
})
##==========================================================================


getMatImagePlot <- function(itr = NULL){
    req(input$inputRDSResultFilename)
    archRresult <- getResult(input$inputRDSResultFilename$datapath)
    if(!is.null(archRresult)){
        if(input$positionLabelsCheck){
            temp_labels <- seq(as.numeric(input$startAt),
                               as.numeric(input$endAt))
            if(input$unCheckZero){
                temp_labels <- setdiff(temp_labels, 0)
            }
            ## make checks
            ## 1. Length is right
            if(length(temp_labels) ==
               Biostrings::width(archRresult$rawSeqs[1])){
                pos_labels <- temp_labels
            }else{
                message("Warning: Length of position labels does not equal ",
                        "length of sequence. Reverting to default ",
                        "position labels")
                pos_labels <- seq(1,
                                  Biostrings::width(archRresult$rawSeqs[1]))
            }
        }else{
            pos_labels <- seq(1,
                              Biostrings::width(archRresult$rawSeqs[1]))
        }
        # pos_labels = seq(1, Biostrings::width(archRresult$rawSeqs[1]))
        if(itr == 0){
            message("[", Sys.time(), "] Plotting null sequences' image matrix")
            if(as.numeric(input$xtickFreqMAT) == 0){
                xt_val <- 1
            }else{
                xt_val <- as.numeric(input$xtickFreqMAT)
            }
            return(archR::viz_matrix_of_acgt_image(
                as.character(archRresult$rawSeqs),
                position_labels = pos_labels,
                savefilename = NULL,
                xt_freq = xt_val,
                yt_freq = as.numeric(input$ytickFreq)
                )
                )
        }else{
            seqsClustLabels <- archRresult$seqsClustLabels[[itr]]
            sorted_order <- sort(seqsClustLabels, index.return = TRUE)
            ##
            message("[", Sys.time(), "] Plotting sequence matrix, iteration ", itr)
            if(as.numeric(input$xtickFreqMAT) == 0){
                xt_val <- 1
            }else{
                xt_val <- as.numeric(input$xtickFreqMAT)
            }

            return(archR::viz_matrix_of_acgt_image(
                as.character(archRresult$rawSeqs[sorted_order$ix]),
                position_labels = pos_labels,
                savefilename = NULL,
                xt_freq = xt_val,
                yt_freq = as.numeric(input$ytickFreq)))
        }
    }
}
##==========================================================================

getArchPWMPlot <- function(itr = NULL, selectedRawSeqs = NULL,
                           plotTitle = NULL, forDownload = FALSE){
    req(input$inputRDSResultFilename)
    archRresult <- getResult(input$inputRDSResultFilename$datapath)
    if(!is.null(archRresult)){
        if(input$positionLabelsCheck){
            temp_labels <- seq(as.numeric(input$startAt),
                               as.numeric(input$endAt))
            if(input$unCheckZero){
                temp_labels <- setdiff(temp_labels, 0)
            }
            ## make checks
            ## 1. Length is right
            if(length(temp_labels) ==
               Biostrings::width(archRresult$rawSeqs[1])){
                pos_labels <- temp_labels
            }else{
                message("Warning: Length of position labels does not equal ",
                        "length of sequence. Reverting to default ",
                        "position labels")
                pos_labels <- seq(1,
                                  Biostrings::width(archRresult$rawSeqs[1]))
            }
        }else{
            pos_labels <- seq(1,
                              Biostrings::width(archRresult$rawSeqs[1]))
        }
        # pos_labels <- seq(1, Biostrings::width(archRresult$rawSeqs[1]))

        if(!is.null(itr) && itr == 0){
            message("[", Sys.time(), "] Plotting null PWM plot")
            if(as.numeric(input$xtickFreqPWM) == 0){
                xt_val <- 1
            }else{
                xt_val <- as.numeric(input$xtickFreqPWM)
            }
            return(archR::plot_ggseqlogo_of_seqs(
                seqs = archRresult$rawSeqs,
                position_labels = pos_labels,
                xt_freq = xt_val,
                title = paste("Sequence logo of all ",
                    length(archRresult$rawSeqs), " sequences" ))
            )
        }
        if(is.null(itr) && !is.null(selectedRawSeqs)){
            # message("SAMARTH -- plotting, aalo re")
            # seqsClustLabels <- archRresult$seqsClustLabels[[itr]]
            # sorted_order <- sort(seqsClustLabels, index.return = TRUE)
            # clusters_ord <- archR::get_seqs_clusters_in_a_list(seqsClustLabels)
            ##
            # message("Generating architectures...itr", itr)
            if(as.numeric(input$xtickFreqPWM) == 0){
                xt_val <- 1
            }else{
                xt_val <- as.numeric(input$xtickFreqPWM)
            }
            return(
                # archR::plot_arch_for_clusters_new(
                # archRresult$rawSeqs,
                # list_of_elements = clusters_ord,
                # position_labels = pos_labels,
                # xt_freq = xt_val,
                # PDFfname = NULL)
                archR::plot_ggseqlogo_of_seqs(
                    seqs = selectedRawSeqs,
                    position_labels = pos_labels,
                    xt_freq = xt_val,
                    title = plotTitle)
            )
        }
        if(forDownload){
            # message("SAMARTH -- plotting, aalo re")
            seqsClustLabels <- archRresult$seqsClustLabels[[itr]]
            sorted_order <- sort(seqsClustLabels, index.return = TRUE)
            clusters_ord <- archR::get_seqs_clusters_in_a_list(seqsClustLabels)
            ##
            # message("Generating architectures...itr", itr)
            if(as.numeric(input$xtickFreqPWM) == 0){
                xt_val <- 1
            }else{
                xt_val <- as.numeric(input$xtickFreqPWM)
            }
            return(
                archR::plot_arch_for_clusters_new(
                archRresult$rawSeqs,
                list_of_elements = clusters_ord,
                position_labels = pos_labels,
                xt_freq = xt_val,
                PDFfname = NULL)
                # archR::plot_ggseqlogo_of_seqs(
                #     seqs = selectedRawSeqs,
                #     position_labels = pos_labels,
                #     xt_freq = xt_val,
                #     title = plotTitle)
            )
        }
    }
}
##==========================================================================






## Adapted from https://gist.github.com/wch/5436415
get_plot_output_list <- function(rawSeqs, clust_list) {
    # Insert plot output objects the list
    nPlots <- length(clust_list)
    if(input$positionLabelsCheck){
        temp_labels <- seq(as.numeric(input$startAt),
                           as.numeric(input$endAt))
        if(input$unCheckZero){
            temp_labels <- setdiff(temp_labels, 0)
        }
        ## make checks
        ## 1. Length is right
        if(length(temp_labels) ==  Biostrings::width(rawSeqs[1])){
            pos_labels <- temp_labels
        }else{
            message("Warning: Length of position labels does not equal ",
                    "length of sequence. Reverting to default ",
                    "position labels")
            pos_labels <- seq(1, Biostrings::width(rawSeqs[1]))
        }
    }else{
        pos_labels <- seq(1, Biostrings::width(rawSeqs[1]))
    }
    ##
    cluster_lengths <- unlist(lapply(clust_list, length))
    cumsums_of_cluster_lengths <- cumsum(cluster_lengths)
    cluster_names <- sort(as.character(seq_along(clust_list)))
    ##
    # message(nPlots)

    ###############
    plot_output_list <- lapply(1:nPlots, function(my_j) {
        # message(my_j)
        ##
        if(my_j > 1){
            startN <- 1 + cumsums_of_cluster_lengths[my_j-1]
        }else{
            startN <- 1
        }
        endN <- cumsums_of_cluster_lengths[my_j]

        plot_title <- paste0("(", my_j , "/", nPlots, ") Arch `",
                             cluster_names[my_j], "': ",
                             cluster_lengths[my_j],
                             " sequences (",  startN, "-",  endN, ")" )
        message(plot_title)
        ##
        plotname <- paste("plot", my_j, sep="")
        plot_output_object <- plotOutput(plotname,
                                         height = "50px", width = "100%")
        plot_output_object <- renderPlot({

            if(as.numeric(input$xtickFreqPWM) == 0){
                xt_val <- 1
            }else{
                xt_val <- as.numeric(input$xtickFreqPWM)
            }
            # message(my_j, " here, samarth")
            withCallingHandlers({
                # shinyjs::html("responseText", "")
                message("[", Sys.time(), "] Generating architecture, ",
                        my_j , "/", nPlots)
                isolate(getArchPWMPlot(itr = NULL,
                                       selectedRawSeqs = rawSeqs[ clust_list[[my_j]] ],
                                       plotTitle = plot_title))
            }, ## withCallingHandlers ends
            message = function(m) {
                shinyjs::html(id = "responseText", html = m$message, add = TRUE)
            },
            warning = function(m) {
                shinyjs::html(id = "responseText", html = m$message, add = TRUE)
            })

            # isolate(
           # archR::plot_ggseqlogo_of_seqs(seqs = rawSeqs[ clust_list[[my_j]] ],
           #                   position_labels = pos_labels,
           #                   xt_freq = xt_val,
           #                   title = plot_title))

            ##
        }, height = 150)
    }) ## lapply ends


    do.call(tagList, plot_output_list) # needed to display properly.

    return(plot_output_list)
}
##==========================================================================

hideTab(inputId = "archRPlotTabs", target = "All sequences' PWM (null plot)")
hideTab(inputId = "archRPlotTabs", target = "ACGT sequence matrix (null plot)")
hideTab(inputId = "archRPlotTabs", target = "Cluster architectures (PWM) Iteration 1")
hideTab(inputId = "archRPlotTabs", target = "Cluster architectures (PWM) Iteration 2")
hideTab(inputId = "archRPlotTabs", target = "Cluster architectures (PWM) Iteration 3")
hideTab(inputId = "archRPlotTabs", target = "Cluster architectures (PWM) Iteration 4")
hideTab(inputId = "archRPlotTabs", target = "Cluster architectures (PWM) Iteration 5")
##
hideTab(inputId = "archRPlotTabs", target = "ACGT sequence matrix Iteration 1")
hideTab(inputId = "archRPlotTabs", target = "ACGT sequence matrix Iteration 2")
hideTab(inputId = "archRPlotTabs", target = "ACGT sequence matrix Iteration 3")
hideTab(inputId = "archRPlotTabs", target = "ACGT sequence matrix Iteration 4")
hideTab(inputId = "archRPlotTabs", target = "ACGT sequence matrix Iteration 5")




#### used with renderUI
observeEvent(input$plotACGTMatrix, {
    if(is.null(input$inputRDSResultFilename)){
        showModal(modalDialog(
            title = "Error",
            "Please choose an archR result file (RDS file)",
            easyClose = TRUE
        ))
    }else{
        req(input$inputRDSResultFilename)
        # # Show element once submit button is pressed
        updateTabsetPanel(session, inputId = "archRTabs",
                          selected = "Response")
        archRresult <- getResult(input$inputRDSResultFilename$datapath)
        max_itr <- length(archRresult$seqsClustLabels)
        nIters <- max_itr
        ##
        for(itrc in seq_len(nIters)){
            local({
            choose_itr <- itrc
            showTab(inputId = "archRPlotTabs",
                    target = paste0("ACGT sequence matrix Iteration ",
                                    choose_itr))
            #
            seqsClustLabels <- archRresult$seqsClustLabels[[choose_itr]]
            sorted_order <- sort(seqsClustLabels, index.return = TRUE)
            clusters_ord <- archR::get_seqs_clusters_in_a_list(seqsClustLabels)
            plotName <- paste0("plotMatImage",choose_itr)
            message("[", Sys.time(), "] Plotting ", plotName)
            output[[plotName]] <- renderPlot({
                withCallingHandlers({
                    # shinyjs::html("responseText", "")
                    isolate(getMatImagePlot(itr = choose_itr))
                }, ## withCallingHandlers ends
                message = function(m) {
                    shinyjs::html(id = "responseText", html = m$message, add = TRUE)
                },
                warning = function(m) {
                    shinyjs::html(id = "responseText", html = m$message, add = TRUE)
                })
            })
            })
        }
        # # Hide loading element when done
        # shinyjs::hideElement(id = 'loading')
        # output$plotMatImage <- renderPlot({
        #     isolate(getMatImagePlot(itr = choose_itr))
        # })
    }
})
##==========================================================================

## Observable event for plotting architecture PWMs
observeEvent(input$plotArchPWMs, {
    if(is.null(input$inputRDSResultFilename)){
        showModal(modalDialog(
            title = "Error",
            "Please choose an archR result file (RDS file)",
            easyClose = TRUE
        ))
    }else{
        req(input$inputRDSResultFilename)
        updateTabsetPanel(session, inputId = "archRTabs",
                          selected = "Response")
        archRresult <- getResult(input$inputRDSResultFilename$datapath)
        max_itr <- length(archRresult$seqsClustLabels)
        nIters <- max_itr
        ##
        for(itrc in seq_len(nIters)){
            local({
            choose_itr <- itrc
            showTab(inputId = "archRPlotTabs",
                    target = paste0("Cluster architectures (PWM) Iteration ",
                                    choose_itr))
            seqsClustLabels <- archRresult$seqsClustLabels[[choose_itr]]
            sorted_order <- sort(seqsClustLabels, index.return = TRUE)
            clusters_ord <- archR::get_seqs_clusters_in_a_list(seqsClustLabels)
            plotName <- paste0("plotPWM",choose_itr)
            output[[plotName]] <- renderUI({
                withCallingHandlers({
                    # shinyjs::html("responseText", "")
                    # message("=== Plotting architectures, iteration ",
                    #         choose_itr, "===")
                    isolate(get_plot_output_list(archRresult$rawSeqs,
                                         clusters_ord))
                    # message("=== Done ===")
                }, ## withCallingHandlers ends
                message = function(m) {
                    shinyjs::html(id = "responseText", html = m$message, add = TRUE)
                },
                warning = function(m) {
                    shinyjs::html(id = "responseText", html = m$message, add = TRUE)
                })
            })
            # output[[plotName]] <- renderUI({
            #     withCallingHandlers({
            #         shinyjs::html("responseText", "")
            #         isolate(get_plot_output_list(archRresult$rawSeqs,
            #                                  clusters_ord))
            #     }, ## withCallingHandlers ends
            #     message = function(m) {
            #         shinyjs::html(id = "responseText", html = m$message, add = TRUE)
            #     },
            #     warning = function(m) {
            #         shinyjs::html(id = "responseText", html = m$message, add = TRUE)
            #     })
            #
            # })
            })
        }


    }
})
##==========================================================================

## Observable event for plotting architecture PWMs
observeEvent(input$plotNullMat, {
    if(is.null(input$inputRDSResultFilename)){
        showModal(modalDialog(
            title = "Error",
            "Please choose an archR result file (RDS file)",
            easyClose = TRUE
        ))
    }else{
        updateTabsetPanel(session, inputId = "archRTabs",
                          selected = "Response")
        showTab(inputId = "archRPlotTabs",
                target = "ACGT sequence matrix (null plot)")
        output$nullPlotMatImage <- renderPlot(execOnResize = FALSE, {
            withCallingHandlers({
                # shinyjs::html("responseText", "")
                isolate(getMatImagePlot(itr = 0))
            }, ## withCallingHandlers ends
            message = function(m) {
                shinyjs::html(id = "responseText", html = m$message, add = TRUE)
            },
            warning = function(m) {
                shinyjs::html(id = "responseText", html = m$message, add = TRUE)
            })
        },
        # caheKeyExpr = { list(input$startAt, input$endAt)}
        )

    }
})
##==========================================================================

## Observable event for plotting null PWMs
observeEvent(input$plotNullPWM, {
    if(is.null(input$inputRDSResultFilename)){
        showModal(modalDialog(
            title = "Error",
            "Please choose an archR result file (RDS file)",
            easyClose = TRUE
        ))
    }else{
        updateTabsetPanel(session, inputId = "archRTabs",
                          selected = "Response")
        showTab(inputId = "archRPlotTabs",
                target = "All sequences' PWM (null plot)")
        output$nullPlotPWM <- renderPlot({
            withCallingHandlers({
                # shinyjs::html("responseText", "")
                isolate(getArchPWMPlot(itr = 0))
            }, ## withCallingHandlers ends
            message = function(m) {
                shinyjs::html(id = "responseText", html = m$message, add = TRUE)
            },
            warning = function(m) {
                shinyjs::html(id = "responseText", html = m$message, add = TRUE)
            })
        },
        # cacheKeyExpr = { list(input$xtickFreqPWM, input$positionLabelsCheck,
        #                       input$startAt, input$endAt, input$unCheckZero) }
        )
    }
})
##==========================================================================


# ## The plots tab panel is loaded dynamically at runtime.
# ## This is because the number of iterations in the provided RDS (result)
# ## file is known only at runtime.
# output$Tabpanels <- renderUI({
#
#     req(input$inputRDSResultFilename)
#     archRresult <- getResult(input$inputRDSResultFilename$datapath)
#
#     nIters <- length(archRresult$seqsClustLabels)
#
#     nTabs <- nIters
#     # nNullTabs <- 2
#     # nPWMTabs <- 3
#     # nMatTabs <- 3
#     # nullTabIds <- c(1,2)
#     # pwmTabIds <- seq_len(nIters) + length(nullTabIds)
#     # matTabIds <- seq_len(nIters) + length(pwmTabIds) + length(nullTabIds)
#     tabNames <- c("nullPWM", "nullACGTMatrix",
#                   paste0('PWMTab', seq_len(nTabs)),
#                   paste0('matrixTab', seq_len(nTabs))
#                   )
#     nullTabs <- vector("list", 2)
#
#     nullTabs[[1]] <- tabPanel("All sequences' PWM (null plot)",
#                               shinycssloaders::withSpinner(
#                                 plotOutput("nullPlotPWM", height = 150),
#                                 type = 4, proxy.height = "100px"
#                             )
#                     )
#     ##
#     nullTabs[[2]] <- tabPanel("ACGT sequence matrix (null plot)",
#                             wellPanel(id = "plotWellPanelInTab",
#                               style = "overflow-y:scroll; max-height: 900px",
#                               shinycssloaders::withSpinner(
#                                   plotOutput("nullPlotMatImage",
#                                              height = 950,
#                                              width = "70%"),
#                                   type = 4, proxy.height = "100px"
#                               )
#                             )
#
#                     )
#     ##
#     pwmTabs <- vector("list", nTabs)
#     ##
#     pwmTabs <- lapply(seq_len(nTabs), function(x){
#         tabPanel(paste0("Cluster architectures (PWM) Iteration ", x),
#                  wellPanel(id = paste0("plotWellPanelInTab", x),
#                            style = "overflow-y:scroll; max-height: 900px",
#                            shinycssloaders::withSpinner(
#                                 uiOutput(paste0("plotPWM",x)),
#                             type = 4, proxy.height = "100px"
#                            )
#                  )
#         )
#     })
#     ##
#     matTabs <- vector("list", nTabs)
#     ##
#     matTabs <- lapply(seq_len(nTabs), function(x){
#         tabPanel(paste0("ACGT sequence matrix Iteration ", x),
#                  wellPanel(id = paste0("plotWellPanelInTab", x),
#                            style = "overflow-y:scroll; max-height: 900px",
#                            # shinyjs::hidden(div(id = 'loading',
#                            shinycssloaders::withSpinner(
#                                plotOutput(paste0("plotMatImage", x),
#                                           height = 950,
#                                           width = "70%"),
#                                type = 4, proxy.height = "100px"
#                            )
#                            # )) ## div-id loading ends
#                  ),
#         )
#     })
#     ##
#     do.call(tabsetPanel, c(nullTabs, pwmTabs, matTabs))
#     # hideTab(inputId = "plotTabs", target = "nullPWM")
#     # hideTab(inputId = "plotTabs", target = "nullACGTMatrix")
#     ## Can we hide these tabs until button press?
#
# })
# ##==========================================================================

# output$nullPlotPWM <- renderPlot({
#     getArchPWMPlot(itr = 0)
# })
#
# output$nullPlotMatImage <- renderPlot({
#     getMatImagePlot(itr = 0)
# })

# output$plotPWM <- renderPlot({
#     getArchPWMPlot(itr = 0)
# })

# output$resultText3 <- renderText({
#     if(is.null(input$performAction))
#         return(NULL)
#     else if(input$performAction == "processNew"){
#         return(NULL)
#     }else{
#         "Ready to analyze existing archR result"
#     }
# })
##==========================================================================

plotInputNullPWM <- function() {
    getArchPWMPlot(itr = 0)
}

plotInputPWM <- function(itr) {
    getArchPWMPlot(itr = itr, forDownload = TRUE)
}


plotInputMat <- function(itr) {
    getMatImagePlot(itr = itr)
}




output$downloadPlotNullMat <- downloadHandler(
    filename = function() { paste(input$nullMatFilename, '.png', sep='') },
    content = function(file) {
        png(file, width = 450, height = 900, units = "px")
        print(plotInputMat(itr = 0))
        dev.off()
    }
)

output$downloadPlotMat1 <- downloadHandler(
    filename = function() { paste(input$Mat1Filename, '.png', sep='') },
    content = function(file) {
        png(file, width = 450, height = 900, units = "px")
        print(plotInputMat(itr = 1))
        dev.off()
    }
)

output$downloadPlotMat2 <- downloadHandler(
    filename = function() { paste(input$Mat2Filename, '.png', sep='') },
    content = function(file) {
        png(file, width = 450, height = 900, units = "px")
        print(plotInputMat(itr = 2))
        dev.off()
    }
)

output$downloadPlotMat3 <- downloadHandler(
    filename = function() { paste(input$Mat3Filename, '.png', sep='') },
    content = function(file) {
        png(file, width = 450, height = 900, units = "px")
        print(plotInputMat(itr = 3))
        dev.off()
    }
)

output$downloadPlotMat4 <- downloadHandler(
    filename = function() { paste(input$Mat4Filename, '.png', sep='') },
    content = function(file) {
        png(file, width = 450, height = 900, units = "px")
        print(plotInputMat(itr = 4))
        dev.off()
    }
)

output$downloadPlotMat5 <- downloadHandler(
    filename = function() { paste(input$Mat5Filename, '.png', sep='') },
    content = function(file) {
        png(file, width = 450, height = 900, units = "px")
        print(plotInputMat(itr = 5))
        dev.off()
    }
)


output$downloadPlotNullPWM <- downloadHandler(
    filename = function() { paste(input$nullPWMFilename, '.pdf', sep='') },
    content = function(file) {
        pdf(file, width = 11, height = 2)
        print(plotInputNullPWM())
        dev.off()
    }
)

output$downloadPlotPWM1 <- downloadHandler(
    filename = function() { paste(input$PWM1Filename, '.pdf', sep='') },
    content = function(file) {
        pdf(file, width = 11, height = 2)
        print(plotInputPWM(itr = 1))
        dev.off()
    }
)

output$downloadPlotPWM2 <- downloadHandler(
    filename = function() { paste(input$PWM2Filename, '.pdf', sep='') },
    content = function(file) {
        pdf(file, width = 11, height = 2)
        print(plotInputPWM(itr = 2))
        dev.off()
    }
)

output$downloadPlotPWM3 <- downloadHandler(
    filename = function() { paste(input$PWM3Filename, '.pdf', sep='') },
    content = function(file) {
        pdf(file, width = 11, height = 2)
        print(plotInputPWM(itr = 3))
        dev.off()
    }
)

output$downloadPlotPWM4 <- downloadHandler(
    filename = function() { paste(input$PWM4Filename, '.pdf', sep='') },
    content = function(file) {
        pdf(file, width = 11, height = 2)
        print(plotInputPWM(itr = 4))
        dev.off()
    }
)

output$downloadPlotPWM5 <- downloadHandler(
    filename = function() { paste(input$PWM5Filename, '.pdf', sep='') },
    content = function(file) {
        pdf(file, width = 11, height = 2)
        print(plotInputPWM(itr = 5))
        dev.off()
    }
)






session$onSessionEnded(function() {
    stopApp()
})

}
