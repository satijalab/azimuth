#' @importFrom shiny NS fileInput
#'
NULL

# tenXMultiFileInput <- function(id, label = "10X Multi File Input") {
#   ns <- NS(id)
# }

tenxH5Input <- function(id, label = "10X H5 File Input") {
  ns <- NS(id)
  fileInput(inputId = ns("file"), label = label)
}

tenxTarballInput <- function(id, label = "10X Tarball Input") {
  ns <- NS(id)
}

tenxMultiFileInput <- function(id, label = "10X Multi-file Input") {
  ns <- NS(id)
}
