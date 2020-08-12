
#' @title Launch archR user interface
#' @description Lauch archR's graphical user interface
#'
#' @export
runArchRUI <- function(lbrowser = TRUE) {
    appDir <- system.file("shiny-interface-to-archR", package = "seqarchR")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `seqarchR`.",
             call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal", launch.browser = lbrowser)
}
