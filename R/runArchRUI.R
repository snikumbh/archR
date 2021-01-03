#' @title Launch archR's graphical user interface
#' 
#' @description Launch archR's Shiny app, a graphical user interface for archR.
#' The app supports the following functionalities.
#' 
#' \describe{
#' \item{Input data}{The app enables reading input FASTA files and provides a 
#' summary noting the number of sequences and existence of any non-ACGT 
#' characters.}
#' 
#' \item{Configuration settings}{The app enables the user to set custom 
#' configuration.}
#' 
#' \item{Processing}{The user can choose to run archR in either serial or 
#' parallel mode. When parallel mode is selected, the user can further choose 
#' to run archR via the Slurm workload manager or otherwise. Jobs can be 
#' submitted and managed (queue status monitoring and canceling a submitted 
#' job) from within the app.}
#' 
#' \item{Viewing result}{The app can help visualize results from an archR run 
#' when provided with archR's result file (as RDS file). Specifically, sequence 
#' logos of architectures and images of clustered sequences can be viewed and 
#' saved to disk as PDF and PNG files respectively}
#' } 
#' 
#' @param lbrowser Logical. Set this to FALSE if you do not wish to launch the
#' app in the browser. Default is TRUE, i.e., the app will launch in the
#' browser. Note: Saving PDF plots when the app is run via the external viewer 
#' in RStudio could have problems.
#'
#' @return 
#' Shiny application
#'
#' @examples 
#' 
#' if(interactive()) {
#'   run_archR_UI(lbrowser = TRUE)
#' }
#'
#' @export
run_archR_UI <- function(lbrowser = TRUE) {
    appDir <- function() system.file("shiny-interface-to-archR", 
        package = "archR")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `archR`.",
                        call. = FALSE)
    }

    if (!requireNamespace(c("shiny", "shinyjs", "shinydashboard",
                            "shinycssloaders"), quietly = TRUE)) {
        stop("Please install the R package 'shiny' to run the application")
    } else {
        shiny::runApp(appDir, display.mode = "normal", 
                                    launch.browser = lbrowser)
    }

}
