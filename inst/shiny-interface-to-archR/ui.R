
## File upload size limit changed to 25MB
options(shiny.maxRequestSize = 25*1024^2)



alignRight <- function (x) {
    tags$div(style="float:right", x)
}

header <- shinydashboard::dashboardHeader(title = paste0("archR v", packageVersion("archR")))

sidebar <- shinydashboard::dashboardSidebar(width = "200px",
    shinydashboard::sidebarMenu(
        id = "performAction",
        shinydashboard::menuItem("Process New Data", tabName = "processNew"),
        shinydashboard::menuItem("Analyze Existing archR Result",
                 tabName = "analyzeExisting")
    )
)


body <- shinydashboard::dashboardBody(
    shinyjs::useShinyjs(),
    ##
    conditionalPanel(
        condition = "input.performAction == 'processNew'",
        titlePanel("Process New Data with archR")
    ),
    conditionalPanel(
        condition = "input.performAction == 'analyzeExisting'",
        titlePanel("Analyze Existing archR Result")
    ),
    fluidRow(
        column(width = 3,
        shinydashboard::box(title = "Input", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, collapsed = FALSE, width = NULL,#4,
            # div(style = 'overflow-y: scroll'),
            conditionalPanel(
                condition = "input.performAction == 'analyzeExisting'",
                h3("Upload data"),
                fileInput(inputId = "inputRDSResultFilename",
                          label = "Choose archR result (RDS file)",
                          multiple = FALSE, accept = c("text/plain"),
                          width = NULL, buttonLabel = "Browse",
                          placeholder = "No file selected"),
            ),
            #### ProcessNew conditional panel ####
            conditionalPanel(
                condition = "input.performAction == 'processNew'",
                ## Get input to process new data from FASTA file
                # tags$hr(style="border-color: black;"),
                h3("Upload data"),
                ## Get input FASTA file
                fileInput(inputId = "inputFastaFilename",
                          label = "Choose FASTA file",
                          multiple = FALSE, accept = NULL,
                          width = NULL, buttonLabel = "Browse",
                          placeholder = "No file selected"),
                # Include clarifying text ----
                helpText("Note: You can now run archR with default parameter settings",
                         "or customize them in the 'Parameter Settings' tab ",
                         "and then run archR"),
                ## action button to call archR
                actionButton("callArchR", "Run archR"),
                ##
            ) ## conditionalPanel for processNew ends
        ), ## input box ends
        ), ## first column ends
        column(width = 4,
        #### Settings box ####
        shinydashboard::box(title = "Parameter Settings", status = "primary",
           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
           width = NULL, #4,
           # background = "maroon",
           conditionalPanel(
               condition = "input.performAction == 'analyzeExisting'",
               ## tabsetPanels
               tabsetPanel(id = "plotParamsTabs", type = "tabs",
               tabPanel("Architecture (PWM) Plots",
                        h4(tags$b("Cluster architectures as PWMs")), br(),
                        ## Params
                        selectInput(inputId = "xtickFreqPWM",
                                    label = "X-axis tick frequency:",
                                    choices = c(seq(1,4,by=1),
                                                seq(5,50,length.out = 10)),
                                    selected = 5, multiple = FALSE,
                                    width = "150px", size = NULL
                                    ),
                        checkboxInput(inputId = "positionLabelsCheck",
                                      label = "Custom sequence position labels",
                                      value = FALSE, width = NULL),
                        conditionalPanel(
                            condition = "input.positionLabelsCheck == 1",
                            numericInput("startAt", "Start at:",
                                         value = 0, width = "150px",
                                         min = NA, max = NA, step = NA),
                            numericInput("endAt", "End at:",
                                         value = 0, width = "150px",
                                         min = NA, max = NA, step = NA),
                            helpText(
                            "Note: Use the checkbox below to include/exclude zero.",
                            "This is ignored unless the given sequence of ",
                            "labels contain a zero."
                            ),
                            checkboxInput(inputId = "unCheckZero",
                                          label = "Exclude zero",
                                          value = FALSE, width = NULL),
                        ),
                        #########
                        ##
                        actionButton("plotArchPWMs", "Plot architecture PWMs"),
                        alignRight(actionButton("plotNullPWM", "Plot null PWM")),
                        ##
                        #############
                        #############
               ),
               tabPanel("ACGT Matrix Plots",
                        h4(tags$b("sequence matrix")), br(),
                        ## Params
                        selectInput(inputId = "xtickFreqMAT",
                                    label = "X-axis tick frequency:",
                                    choices = c(seq(1,4,by=1),
                                                seq(5,50,length.out = 10)),
                                    selected = 5, multiple = FALSE,
                                    width = "150px", size = NULL
                        ),
                        selectInput(inputId = "ytickFreq",
                                    label = "Y-axis tick frequency:",
                                    choices = c(seq(100,1000,length.out = 19)),
                                    selected = 200, multiple = FALSE,
                                    width = "150px", size = NULL
                        ),
                        ##
                        actionButton("plotACGTMatrix", "Plot sequence matrices"),
                        ##
                        alignRight(actionButton("plotNullMat", "Plot null matrix")),
               )),

               ## wellpanel solution
               # wellPanel(id = "paramWellPanelAnalyzeExisting",
               #           style = "overflow-y:scroll; max-height: 500px",
               # ## Plot cluster architecturess (PWMs from sequences)
               # ## action button to plot architecture sequence logos
               # ## =======
               # # tags$hr(style="border-color: black;"),
               #
               # h4(tags$b("Cluster architectures as PWMs")), br(),
               # ## Params
               # numericInput("xtickFreqPWM", "X-axis tick frequency:",
               #              value = 5, width = "150px",
               #              min = 0, max = 50, step = 5),
               # ##
               # actionButton("plotArchPWMs", "Plot architectures as PWMs"),
               # ##
               # actionButton("plotNullPWM", "Generate null PWM"),
               # tags$hr(style="border-color: black;"),
               # ##
               # ## Plot ACGT matrix of clustered sequences
               # ## =======
               # # tags$hr(style="border-color: black;"),
               # h4(tags$b("ACGT sequence matrix")), br(),
               # ## Params
               # numericInput("xtickFreqMAT", "X-axis tick frequency:",
               #              value = 5, width = "150px",
               #              min = 0, max = 50, step = 5),
               # numericInput("ytickFreq", "Y-axis tick frequency:",
               #              value = 200, width = "150px",
               #              min = 100, max = 1000, step = 50),
               # ##
               # actionButton("plotACGTMatrix", "Plot ACGT sequence matrices"),
               # ##
               # actionButton("plotNullMat", "Generate null matrix"),
               # ), ## paramWellPanelAnalyzeExisting
           ),
           #### ProcessNew conditional panel ####
           conditionalPanel(
               condition = "input.performAction == 'processNew'",
               wellPanel(id = "paramWellPanelProcessNew",
                         style = "overflow-y:scroll; max-height: 500px;background: white",
               # h3("Set parameters"),
               ## Set parameters
               ## Ability to reset to default?
               ## Maybe, just selected inputs such as just the numeric ones
               ## Tolerance, Bound
               ##
               # alignRight(

               # ),
               ## innerChunkSize Numeric.
               sliderInput(inputId = "innerChunkSize",
                           label = "Chunk size:",
                           min = 100, max = 1000, value = 500, step = 100,
                           round = FALSE, ticks = TRUE,
                           animate = FALSE, width = NULL, sep = ","),
               ## kMin Numeric.
               ## kMax Numeric.
               sliderInput(inputId = "kMax", label = "Maximum number of factors:",
                           min = 10, max = 100, value = 20, step = 5,
                           round = FALSE, ticks = TRUE,
                           animate = FALSE, width = NULL, sep = ","),
               ##
               ## modSelType Character. Radio button
               ##
               radioButtons(inputId = "modSelType",
                            label = "Model Selection Procedure:",
                            choices = c("Stability" = "stability",
                                        "Cross-validation (CV)" = "cv"),
                            inline = TRUE),
               ##
               conditionalPanel(
                   condition = "input.modSelType == 'cv'",
                   ## cvFolds Numeric.
                   selectInput(inputId = "cvFolds", label = "Number of CV folds:",
                               choices = c(3,5,10), selected = 5, multiple = FALSE,
                               selectize = TRUE, width = "100px", size = NULL),
               ),
               ##
               ##
               conditionalPanel(
                   condition = "input.modSelType == 'stability'",
                   ## tol Numeric.
                   sliderInput(inputId = "tol", label = "Tolerance:",
                               min = -8, max = -1, value = -3, step = 1,
                               round = FALSE, ticks = TRUE,
                               animate = FALSE, width = NULL, sep = ",", pre = "10^"),
                   ## bound Numeric.
                   sliderInput(inputId = "bound", label = "Bound:",
                               min = -12, max = -1, value = -8,
                               step = 1,
                               round = FALSE, ticks = TRUE,
                               animate = FALSE, width = NULL, sep = ",", pre = "10^"),
               ),
               ##
               ## nIterationsUse
               numericInput("nIterationsUse", "# Iterations for NMF:",
                            value = 100, min = 50, max = 1000, step = 50,
                            width = "150px"),
               ## Checkpoints
               radioButtons(inputId = "checkpoint", label = "Save checkpoints?",
                            choices = c("Yes" = "TRUE", "No" = "FALSE"),
                            selected = "TRUE",
                            inline = TRUE),
               ## alphaBase,alphaPow
               ## minSeqs
               ## modSelLogFile
               ##
               numericInput("thresholdsItr", "# Iterations for archR:",
                            value = 3, min = 1, max = 5, step = 1,
                            width = "150px"),
               ##
               checkboxGroupInput("flags", "Set flags for output display:",
                                  c("Debug" = "debugFlag",
                                    "Verbosity" = "verboseFlag",
                                    "Plotting" = "plotVerboseFlag",
                                    "Time" = "timeFlag"), inline = TRUE),
               ##
               textInput(inputId = "resultLocation",
                         label = "Location to save archR result:",
                         value = "./archR_run", width = "250px", placeholder = NULL),
               ##
               ## parallelize Logical.
               radioButtons(inputId = "parallelize", label = "Parallelize?",
                            choices = c("Yes" = "TRUE", "No" = "FALSE"),
                            selected = "FALSE",
                            inline = TRUE),
               tags$hr(style="border-color: black;"),
               ##
               conditionalPanel(
                   condition = "input.parallelize == 'TRUE'",
                   numericInput("nCoresUse", "# Cores:", value = 4,
                                width = "100px", min = 2, max = 48, step = 1),
                   ## slurm-based Logical.
                   radioButtons(inputId = "useSlurm", label = "Use Slurm?",
                                choices = c("Yes" = "TRUE", "No" = "FALSE"),
                                selected = "TRUE",
                                inline = TRUE),
                   ##
                   conditionalPanel(
                       condition = "input.useSlurm == 'TRUE'",
                       textInput(inputId = "jobName", label = "Job name:",
                                 value = "archR_app_job", width = "300px"),
                       textInput(inputId = "memPerCPU", label = "`--mem-per-cpu`:",
                                 value = "5G", width = "200px", placeholder = NULL),
                       numericInput(inputId = "nodes", label = "`--nodes`:",
                                    value = "1", width = "100px"),
                       textInput(inputId = "jobLocation",
                                 label = "Location to write archR job files:",
                                 value = ".", width = "250px", placeholder = NULL),
                       checkboxInput(inputId = "submitAsk", label = "Submit job?",
                                     value = TRUE, width = NULL),
                       helpText(
                           "Note: Currently, you can only submit one job at a time."),
                   ),
               ), ## parallelize conditionalPanel ends
               tags$hr(style="border-color: black;"),
                checkboxInput(inputId = "resetCheck", label = "Reset parameters to default values?",
                              value = FALSE, width = "300px"),
               conditionalPanel(
                   condition = "input.resetCheck == 1",
               selectInput(inputId = "resetSelect",
                           label = "Select parameter(s) to reset:",
                           choices = c("All" = "all",
                                       "Chunk size" = "innerChunkSize",
                                       "Maximum number of factors" = "kMax",
                                        "Model Selection Procedure" = "modSelType",
                                       "Tolerance (Stability-based model selection)" = "tol",
                                       "Bound (Stability-based model selection)" = "bound",
                                       "Number of CV folds ( (CV-based model selection)" = "cvFolds",
                                       "# Iterations for NMF" = "nIterationsUse",
                                       "# Iterations for archR" = "thresholdsItr",
                                       "Location to save archR result" = "resultLocation"
                           ),
                           selected = "all", multiple = TRUE,
                           selectize = TRUE, width = "300px", size = NULL),
               actionButton(inputId = "reset", label = "Reset")
               ), ## conditionalPanel for reset ends
               ) ## all params wellPanel ends
           ) ## conditionalPanel for processNew ends
        ) ## parameter settings box ends
        # #### runArchR button box ####
        # box(title = "Parameter settings", status = "primary",
        #     solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
        #     width = NULL,#4,
        # ) ## runArchR button box ends
        ), ## second column ends

        column(width = 5,
        #### Info tabBox ####
        ####
        shinydashboard::box(title = "Summary", status = "primary",
            solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
            width = NULL, # height = "440px", #4,
        # tabBox(title = "Summary", width = 8, height = "440px",
               ## keeping height of this tabBox fixed is important
               ## because a change in the tabBox height triggers a
               ## redrawing/adjustment of the Plots tabBox
               ##
               tabsetPanel(id = "archRTabs", type = "tabs",
                    tabPanel("Configuration",
                             conditionalPanel(
                                 condition = "input.performAction == 'processNew'",
                                 # Include clarifying text ----
                                 helpText("About this tab: View chosen archR ",
                                          "configuration."),
                             ),
                             conditionalPanel(
                                 condition = "input.performAction == 'analyzeExisting'",
                                 # Include clarifying text ----
                                 helpText("About this tab: View the ",
                                          "then configuration settings for ",
                                          "the archR result file provided.")
                             ),
                             br(),
                            tableOutput("configSummary")),
                    ##
                    tabPanel("Data",
                         conditionalPanel(
                             condition = "input.performAction == 'processNew'",
                             # Include clarifying text ----
                             helpText("About this tab: View a brief summary of ",
                                      "the FASTA file provided. ",
                                      "Includes #sequences, lengths of sequences, ",
                                      "and the alphabet over all the sequences."),
                         ),
                         conditionalPanel(
                             condition = "input.performAction == 'analyzeExisting'",
                             # Include clarifying text ----
                             helpText("About this tab: View a brief summary ",
                                      "of the archR result file provided. ",
                                       "Includes #iterations, and # clusters ",
                                       "per iteration.")
                         ),
                         br(),
                        tableOutput("dataSummary")),
                    ##
                    # tabPanel("Result",
                    #          verbatimTextOutput("resultText3")),
                    ##
                    tabPanel(title = "Response", value = "Response",
                             conditionalPanel(
                                 condition = "input.performAction == 'processNew'",
                             # Include clarifying text ----
                             helpText("About this tab: View R console output upon ",
                                      "running archR")
                             ),
                             conditionalPanel(
                                 condition = "input.performAction == 'analyzeExisting'",
                                 # Include clarifying text ----
                                 helpText("About this tab: View response/progress ",
                                          "from plotting functions.")
                             ),
                             br(),
                            wellPanel(id = "paramWellPanelOutputLog",
                                style =
                                    "overflow-y:scroll; max-height: 500px; background: white",
                                verbatimTextOutput("responseText")
                            )
                            ),
                    ##
                    tabPanel("Status",
                             # Include clarifying text ----
                             helpText("About this tab: View status/queue of ",
                                    "archR jobs submitted via the archR GUI. ",
                                     ),
                            br(),
                            tableOutput("statusTable"),
                            br(),
                            actionButton("checkStatus", "Check Job Status"),
                            alignRight(actionButton("cancelJob", "Cancel Job")),
                            br(), br(), br(),
                            # Include clarifying text ----
                            helpText(
                            "Note: Only showing jobs submitted via archR UI.",
                            "Currently, you can only submit one job at a time."),
                    ),
                    ## OutputLog 1
                    tabPanel("OutputLog",
                             # Include clarifying text ----
                             helpText("About this tab: View archR output ",
                                      "from the job submitted to Slurm ",
                                      "via the archR GUI. ",
                             ),
                             br(),
                    wellPanel(id = "paramWellPanelOutputLog1",
                        style = "overflow-y:scroll; max-height: 500px; background: white",
                        verbatimTextOutput("logText")
                     )
                    )
                    ##
                    ##
                    ## TO-DO:
                    ## Can we manage multiple slurm submits from within
                    ## the archR UI?
                    ## 1. For every slurm-submit, we add a outputLog tab.
                    ## 2. Name this tab by the submitted job ID?
                    ## 3. Say, we allow a maximum of 3 jobs
                    ## 4. We can statically add tabs to the tabsetPanel,
                    ##  but keep them hidden, and show them only when
                    ##  additional jobs get submitted.
                    ## 5. How will we pair the tab with the slurm output
                    ## file?
               ) ## tabsetPanel ends

        ), ## tabBox ends
        ), ## third column ends

    ), ## first fluidRow ends

    fluidRow(
        # #### Info tabBox ####
        # conditionalPanel(
        #     condition = "input.performAction == 'processNew'",
        #     shinydashboard::infoBoxOutput("infoBoxProcessNew")
        # ),
        #### Plots tabBox ####
        conditionalPanel(
            condition = "input.performAction == 'analyzeExisting'",
            shinydashboard::box(title = "Plots", status = "primary",
               solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
               width = 12,
               # uiOutput("Tabpanels")
                   # ##
               tabsetPanel(id = "archRPlotTabs", type = "tabs",
               tabPanel("All sequences' PWM (null plot)",
                        textInput("nullPWMFilename", label = "Filename",
                                  value = "null_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotNullPWM",
                                           label = "Download this plot"),
                        tags$hr(style="border-color: black;"),
                        # shinycssloaders::withSpinner(
                            plotOutput("nullPlotPWM", height = 150),
                            # type = 4, proxy.height = "100px"
                        # )
               ),
                ##
               tabPanel("ACGT sequence matrix (null plot)",
                        textInput("nullMatFilename", label = "Filename",
                                  value = "null_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotNullMat",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = "plotWellPanelInTab",
                                  style = "overflow-y:scroll; max-height: 900px; background: white",
                                  # shinycssloaders::withSpinner(
                                      plotOutput("nullPlotMatImage",
                                                 height = 950,
                                                 width = "70%"),
                                  #     type = 4, proxy.height = "100px"
                                  # )
                        )

               ),
               ##
               tabPanel("Cluster architectures (PWM) Iteration 1",
                        textInput("PWM1Filename", label = "Filename",
                                  value = "iteration1_arch_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotPWM1",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab1"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",
                                  # shinycssloaders::withSpinner(
                                  uiOutput("plotPWM1"),
                                  #     type = 4, proxy.height = "100px"
                                  # )
                        )
               ),
               tabPanel("Cluster architectures (PWM) Iteration 2",
                        textInput("PWM2Filename", label = "Filename",
                                  value = "iteration2_arch_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotPWM2",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab2"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",
                                  # shinycssloaders::withSpinner(
                                  uiOutput("plotPWM2"),
                                  #     type = 4, proxy.height = "100px"
                                  # )
                        )
               ),
               tabPanel("Cluster architectures (PWM) Iteration 3",
                        textInput("PWM3Filename", label = "Filename",
                                  value = "iteration3_arch_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotPWM3",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab3"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",
                                  # shinycssloaders::withSpinner(
                                  uiOutput("plotPWM3")
                                  #     type = 4, proxy.height = "100px"
                                  # )
                        )
               ),
               tabPanel("Cluster architectures (PWM) Iteration 4",
                        textInput("PWM4Filename", label = "Filename",
                                  value = "iteration4_arch_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotPWM4",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab4"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",
                                  # shinycssloaders::withSpinner(
                                  uiOutput("plotPWM4")
                                  #     type = 4, proxy.height = "100px"
                                  # )
                        )
               ),
               tabPanel("Cluster architectures (PWM) Iteration 5",
                        textInput("PWM5Filename", label = "Filename",
                                  value = "iteration5_arch_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotPWM5",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab5"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",
                                  # shinycssloaders::withSpinner(
                                  uiOutput("plotPWM5")
                                  #     type = 4, proxy.height = "100px"
                                  # )
                        )
               ),
               # ),
               ##
               tabPanel("ACGT sequence matrix Iteration 1",
                        textInput("Mat1Filename", label = "Filename",
                                  value = "iteration1_mat_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotMat1",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab2"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",
                                  # shinycssloaders::withSpinner(
                                      plotOutput("plotMatImage1",
                                                 height = 950,
                                                 width = "70%"),
                                  #     type = 4, proxy.height = "100px"
                                  # )

                        )
               ),

               tabPanel("ACGT sequence matrix Iteration 2",
                        textInput("Mat2Filename", label = "Filename",
                                  value = "iteration2_mat_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotMat2",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab3"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",

                                  # shinycssloaders::withSpinner(
                                      plotOutput("plotMatImage2",
                                      height = 950,
                                      width = "70%"),
                                #   type = 4, proxy.height = "100px"
                                # )

                        )
               ),
               # ##
               tabPanel("ACGT sequence matrix Iteration 3",
                        textInput("Mat3Filename", label = "Filename",
                                  value = "iteration3_mat_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotMat3",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab4"),
                                  style = "overflow-y:scroll; max-height: 900px",

                                  # shinycssloaders::withSpinner(
                                      plotOutput("plotMatImage3",
                                      height = 950,
                                      width = "70%"),
                                #   type = 4, proxy.height = "100px"
                                # )
                        )
                ),
               ##
               tabPanel("ACGT sequence matrix Iteration 4",
                        textInput("Mat4Filename", label = "Filename",
                                  value = "iteration4_mat_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotMat4",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab5"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",

                                  # shinycssloaders::withSpinner(
                                  plotOutput("plotMatImage4",
                                             height = 950,
                                             width = "70%"),
                                  #   type = 4, proxy.height = "100px"
                                  # )
                        )
               ),
               ##
               tabPanel("ACGT sequence matrix Iteration 5",
                        textInput("Mat5Filename", label = "Filename",
                                  value = "iteration5_mat_plot", width = "200px"),
                        downloadButton(outputId = "downloadPlotMat5",
                                       label = "Download the plot"),
                        tags$hr(style="border-color: black;"),
                        wellPanel(id = paste0("plotWellPanelInTab6"),
                                  style = "overflow-y:scroll; max-height: 900px; background: white",

                                  # shinycssloaders::withSpinner(
                                  plotOutput("plotMatImage5",
                                             height = 950,
                                             width = "70%"),
                                  #   type = 4, proxy.height = "100px"
                                  # )
                        )
               )
               ),## tabsetPanel Ends
               ##
               ##
            ), ## tabBox ends
        )
    ),

) ## dashboardBody ends






ui <- shinydashboard::dashboardPage(header, sidebar, body)



