#' @docType package
#' @name SeuratMapper-pacakge
#' @rdname SeuratMapper-pacakge
#'
#' @aliases SeuratMapper
#'
"_PACKAGE"


GetDefaultArguments <- function(f, method = NULL) {
  UseMethod(generic = 'GetDefaultArguments', object = f)
}

#' @method GetDefaultArguments character
#'
GetDefaultArguments.character <- function(f, method = NULL) {
  ''
}

#' @method GetDefaultArguments function
#'
GetDefaultArguments.function <- function(f, method = NULL) {
  browser()
  return(GetDefaultArguments(
    f = as.character(x = substitute(expr = f)),
    method = method
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  op <- options()
  # TODO: replace this
  options(shiny.maxRequestSize = 100 * (1024 ^ 2))
}
