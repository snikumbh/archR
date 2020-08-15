#' @title Launch archR user interface
#' @description Lauch archR's graphical user interface
#'
#' @param lbrowser Logical. Set this to FALSE if you do not wish to launch the
#' app in the browser. Default is TRUE, i.e., the app will launch in the
#' browser. Note: Saving plots when the app is run via the external viewer in
#' RStudio could have problems.
#'
#' @export
runArchRUI <- function(lbrowser = TRUE) {
    appDir <- system.file("shiny-interface-to-archR", package = "archR")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `archR`.",
             call. = FALSE)
    }

    if (requireNamespace(c("shiny", "shinyjs", "shinydashboard",
                           "shinycssloaders"), quietly = TRUE)) {
        shiny::runApp(appDir, display.mode = "normal", launch.browser = lbrowser)
    } else {
        message("Please install the R package shiny to run the application")
    }

}
