#' @title Launch archR user interface
#' @description Lauch archR's graphical user interface
#'
#' @export
runArchRUI <- function(lbrowser = TRUE) {
    appDir <- system.file("shiny-interface-to-archR", package = "archR")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `archR`.",
             call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal", launch.browser = lbrowser)
}
