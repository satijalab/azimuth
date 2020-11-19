#' @include zzz.R
#' @include seurat.R
#' @include helpers.R
#' @include ui.R
#'
NULL

#' Launch the mapping app
#'
#' @param config Path to JSON-formatted configuration file specifying options;
#' for an example config file, see
#' \code{system.file("resources", "config.json", package = "Azimuth")}
#' @param ... Options to set, see \code{?`\link{Azimuth-package}`} for details
#' on \pkg{Azimuth}-provided options
#'
#' @section Specifying options:
#' R options can be provided as named arguments to \code{AzimuthApp} through
#' dots (...), set in a config file, or set globally. Arguments provided to
#' \code{AzimuthApp} through dots take precedence if the same option is provided
#' in a config file. Options provided through dots or a config file take
#' precedence if the same option was set globally.
#'
#' Options in the \code{\link[Azimuth:Azimuth-package]{Azimuth.app}} namespace
#' can be specified using a shorthand notation in both the config file and as
#' arguments to \code{AzimuthApp}. For example, the option
#' \code{Azimuth.app.reference} can be shortened to \code{reference} in the
#' config file or as an argument to \code{AzimuthApp}
#'
#' @return None, launches the mapping Shiny app
#'
#' @importFrom shiny runApp shinyApp
#' @importFrom withr with_options
#' @importFrom jsonlite read_json
#'
#' @export
#'
#' @seealso \code{\link{Azimuth-package}}
#'
#' @examples
#' if (interactive()) {
#'   AzimuthApp(system.file("resources", "config.json", package = "Azimuth"))
#' }
#'
AzimuthApp <- function(config = NULL, ...) {
  useShinyjs()
  # If multiple items have the same name in the named list, with_options sets
  # the option to the last entry with that name in the list. Therefore, putting
  # the config file options first, followed by options set in dots, followed by
  # hardcoded options, achieves the desired precedence.
  opts <- list()
  # Add options set through config file
  if (!is.null(x = config)) {
    opts <- c(opts, read_json(path = config, simplifyVector = TRUE))
  }
  # Add options set through named arguments in dots
  args <- list(...)
  if (length(x = args) && !is.null(x = names(x = args))) {
    # only add named elements
    opts <- c(opts, args[names(x = args) != ""])
  }
  # if any arguments from dots or config file have no "." character,
  # prepend the "Azimuth.app" namespace
  for (i in seq_along(along.with = opts)) {
    if (!grepl(pattern = '\\.', x = names(x = opts)[i])) {
      names(x = opts)[i] <- paste0('Azimuth.app.', names(x = opts)[i])
    }
  }
  # Add sensible defaults
  # Shiny doesn't set shiny.maxRequestSize on load
  if (!'shiny.maxRequestSize' %in% names(x = opts) && is.null(x = getOption(x = 'shiny.maxRequestSize'))) {
    opts$shiny.maxRequestSize <- 500 * (1024 ^ 2)
  }
  # Add pageLength to jQuery DataTables options
  opts$DT.options <- as.list(x = c(
    opts$DT.options,
    getOption(x = 'DT.options')
  ))
  if (!'pageLength' %in% names(x = opts$DT.options)) {
    opts$DT.options$pageLength <- 10L
  }
  # Set future.globals.maxSize; this is not user-configurable
  maxcells <- with_options(
    new = opts,
    code = getOption(x = 'Azimuth.app.max_cells')
  )
  opts$future.globals.maxSize <- maxcells * 320000
  # Launch the app
  with_options(
    new = opts,
    code = runApp(appDir = shinyApp(ui = AzimuthUI, server = AzimuthServer))
  )
  return(invisible(x = NULL))
}
