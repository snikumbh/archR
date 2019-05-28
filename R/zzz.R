
sklearn <- NULL

.onLoad <- function(libname, pkgname){
  # use superassignment to update global reference to scipy
  # sklearn <<- reticulate::import("sklearn", delay_load = TRUE)
  # reticulate::source_python('inst/py/perform_nmf.py')
  # reticulate::source_python(system.file("extdata", "example_data.fa", package = "archeR", mustWork = TRUE))
  # reticulate::source_python('inst/python/perform_nmf.py', envir = globalenv())
}
