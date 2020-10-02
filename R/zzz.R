#' @docType package
#'
#' @section Package options:
#'
#' SeuratMapper uses the following options to control the behaviour of the app,
#' users can configure these with \code{\link[base]{options}}:
#'
#' \describe{
#' \item{\code{Azimuth.map.ncells}}{
#'   Minimum number of cells required to accept uploaded file.
#'   Defaults to \code{100}
#'  }
#'  \item{\code{Azimuth.map.ngenes}}{
#'   Minimum number of genes in common with reference to accept uploaded file.
#'   Defaults to \code{250}
#'  }
#'  \item{\code{Azimuth.map.nanchors}}{
#'   Minimum number of anchors that must be found to complete mapping.
#'   Defaults to \code{50}
#'  }
#'  \item{\code{Azimuth.map.pbcorthresh}}{
#'   Only proceed to mapping if query dataset meets or exceeds this threshold in
#'   pseudobulk correlation test.
#'  }
#'  \item{\code{Azimuth.de.mincells}}{
#'   Minimum number of cells per cluster for differential expression; defaults
#'   to \code{15}
#'  }
#'  \item{\code{Azimuth.de.digits}}{
#'   Number of digits to round differential expression table to; defaults to
#'   \code{3}
#'  }
#'  \item{\code{Azimuth.sct.ncells}, \code{Azimuth.sct.nfeats}}{
#'   Number of cells and features to use for
#'   \code{\link[Seurat]{SCTransform}}, respectively. Defaults to \code{2000}
#'   for each
#'  }
#' }
#'
#' @inheritSection AzimuthApp App options
#'
#'
#' @aliases SeuratMapper
#'
"_PACKAGE"

default.options <- list(
  Azimuth.de.digits = 3L,
  Azimuth.de.mincells = 15L,
  Azimuth.map.ncells = 100L,
  Azimuth.map.ngenes = 250L,
  Azimuth.map.nanchors = 50L,
  Azimuth.map.pbcorthresh = 0.75,
  Azimuth.sct.ncells = 2000L,
  Azimuth.sct.nfeats = 2000L
)

#' Attach dependent packages
#'
#' Attaches the following packages
#' \itemize{
#'  \item shinyBS
#' }
#'
#' @return Attaches the required packages and returns invisible \code{NULL}
#'
#' @keywords internal
#'
AttachDeps <- function() {
  deps <- c(
    'shinyBS'
  )
  for (d in deps) {
    if (!paste0('package:', d) %in% search()) {
      packageStartupMessage("Attaching ", d)
      attachNamespace(ns = d)
    }
  }
}

#' Returns a dataframe of the frequency or percentage of levels of category.2
#' (column) for object split by each level of category.1 (row)
#'
#' @param object a Seurat object
#' @param category.1 a metadata field in the object
#' @param category.2 another metadata field in the object
#' @param percentage if TRUE, returns percentages; otherwise, counts
#'
#' @importFrom Seurat FetchData
#'
#' @keywords internal
#'
CategoryTable <- function(
  object,
  category.1,
  category.2,
  percentage = FALSE
) {
  data <- FetchData(object = object, vars = c(category.1, category.2))
  data[, category.1] <- droplevels(x = factor(x = data[, category.1]))
  data[, category.1] <- factor(x = data[, category.1], levels = sort(x = levels(x = data[, category.1])))
  data[, category.2] <- droplevels(x = factor(x = data[, category.2]))
  data[, category.2] <- factor(x = data[, category.2], levels = sort(x = levels(x = data[, category.2])))
  tbl <- table(
    data[, category.1],
    data[, category.2],
    useNA = "ifany"
  )
  if (percentage) {
    tbl <- t(x = apply(
      X = tbl,
      MARGIN = 1,
      FUN = function(x) {
        return(round(x = 100 * (x / sum(x)), digits = 1))
      }
    ))
    if (length(levels(data[, category.2])) == 1) {
      tbl <- t(tbl)
      colnames(x = tbl) <- levels(x = data[, category.2])
    }
  }
  return(as.data.frame.matrix(x = tbl))
}

#' Create an annoy index
#'
#' @note Function exists because it's not exported from \pkg{uwot}
#'
#' @param name,ndim Paramters
#'
#' @return An nn index object
#'
#' @importFrom methods new
#' @importFrom RcppAnnoy AnnoyAngular AnnoyManhattan AnnoyEuclidean AnnoyHamming
#'
#' @keywords internal
#'
CreateAnn <- function(name, ndim) {
  return(switch(
    EXPR = name,
    cosine = new(Class = AnnoyAngular, ndim),
    manhattan = new(Class = AnnoyManhattan, ndim),
    euclidean = new(Class = AnnoyEuclidean, ndim),
    hamming = new(Class = AnnoyHamming, ndim),
    stop("BUG: unknown Annoy metric '", name, "'")
  ))
}

#' Sanitize feature names for \code{\link[shiny]{selectInput}}
#'
#' \code{\link[shiny]{selectInput}} has some limitiations with biological
#' feature names. This function sanitizes feature names according to the
#' following rules:
#' \itemize{
#'  \item Names matching the regular expression \dQuote{\\.\\d+$} are
#'   \strong{removed}
#' }
#'
#' @param features A character vector of feature names
#'
#' @return \code{features}, but sanitized
#'
#' @keywords internal
#'
#' @seealso \code{\link[shiny]{selectInput}}
#'
FilterFeatures <- function(features) {
  return(sort(x = grep(
    pattern = '\\.\\d+$',
    x = features,
    value = TRUE,
    invert = TRUE
  )))
}

#' Return names of metadata columns in a Seurat object that have an
#' appropriate number of levels for plotting when converted to a factor
#'
#' @param object a Seurat object
#' @param min.levels minimum number of levels in a metadata factor to include
#' @param max.levels maximum number of levels in a metadata factor to include
#'
#' @keywords internal
#'
PlottableMetadataNames <- function(
  object,
  min.levels = 2,
  max.levels = 20
) {
  column.status <- sapply(
    X = object[[]],
    FUN = function(column) {
      length(x = levels(x = droplevels(x = as.factor(x = column)))) >= min.levels &&
        length(x = levels(x = droplevels(x = as.factor(x = column)))) <= max.levels
    }
  ) & (colnames(object[[]]) != "mapping.score") &
   (colnames(object[[]]) != "predicted.id")
  return(colnames(object[[]])[column.status])
}

#' Prepare differential expression results for rendering
#'
#' @param diff.exp A dataframe with differential expression results from
#' \code{\link[presto:wilcoxauc]{presto::wilcoxauc}}
#' @param groups.use Names of groups to filter \code{diff.exp} to; groups must
#' be found in \code{diff.exp$group}
#' @param n Number of feature to filter \code{diff.exp} to per group
#' @param logfc.thresh logFC threshold
#'
#' @return \code{diff.exp}, ordered by adjusted p-value, filtered to \code{n}
#' features per group in \code{group.use}
#'
#' @importFrom rlang %||%
#' @importFrom utils head
#'
#' @seealso \code{\link[presto]{wilcoxauc}}
#'
#' @keywords internal
#'
RenderDiffExp <- function(
  diff.exp,
  groups.use = NULL,
  n = 10L,
  logfc.thresh = 0L
) {
  cols.keep <- c('avgExpr', 'auc', 'padj', 'pct_in', 'pct_out')
  groups.use <- groups.use %||% unique(x = as.character(x = diff.exp$group))
  diff.exp <- lapply(
    X = groups.use,
    FUN = function(group) {
      group.de <- diff.exp[diff.exp$group == group, , drop = FALSE]
      group.de <- group.de[group.de$logFC > logfc.thresh, , drop = FALSE]
      group.de <- group.de[order(group.de$padj, -group.de$auc), , drop = FALSE]
      return(head(x = group.de, n = n))
    }
  )
  diff.exp <- do.call(what = 'rbind', diff.exp)
  rownames(x = diff.exp) <- make.unique(names = diff.exp$feature)
  diff.exp <- signif(
    x = diff.exp[, cols.keep, drop = FALSE],
    digits = getOption(
      x = "Azimuth.de.digits",
      default = default.options$Azimuth.de.digits
    )
  )
  return(diff.exp)
}

#' Manage tabs with \pkg{shinyjs}
#'
#' Quickly generate JavaScript IDs for Shiny tab panels. Also build a function
#' to hide Shiny tab panels in JavaScript using \pkg{shinyjs}
#'
#' @param id ID of a \code{\link[shiny]{tabsetPanel}}
#' @param values One or more values of a \code{\link[shiny]{tabPanel}}
#' (see the \code{value} parameter)
#' @param fxn Name of JavaScript call function
#'
#' @return \code{TabJSHide}: a string with a JavaScript function to hide a set
#' of tabs
#'
#' @name TabJS
#' @rdname TabJS
#'
#' @keywords internal
#'
#' @note These functions are designed to run custom JavaScript code using
#' \code{\link[shinyjs]{extendShinyJS}}; use of custom JavaScript code requires
#' the \pkg{V8} package. \pkg{V8} requires a local install of either the
#' \href{https://chromium.googlesource.com/v8/v8}{V8} JavaScript Engine or
#' \href{https://nodejs.org/en/}{Node.js}
#'
#' @seealso \code{\link[shinyjs]{extendShinyJS}} \code{\link[shiny]{tabPanel}}
#'
#' @examples
#' \dontrun{
#' TabJSHide('tabs', values = c('mapped', 'fexplorer'))
#' }
#'
TabJSHide <- function(id, values, fxn = 'hide') {
  keys <- paste0(
    '$(',
    sQuote(x = TabJSKey(id = id, values = values), q = FALSE),
    ').hide()',
    collapse = '; '
  )
  return(paste0(
    'shinyjs.',
    fxn,
    ' = function() {',
    keys,
    ';}'
  ))
}

#' @rdname TabJS
#'
#' @return \code{TabJSKey}: a string with the JavaScript ID for a given
#' set of tabs
#'
#' @examples
#' \dontrun{
#' TabJSKey('tabs', values = 'mapped')
#' }
#'

TabJSKey <- function(id, values) {
  return(paste0(
    "#",
    id,
    " li a[data-value=",
    dQuote(x = values, q = FALSE),
    "]"
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  # Attach deps
  AttachDeps()
  op <- options()
  # TODO: replace this
  options(shiny.maxRequestSize = 100 * (1024 ^ 2))
  # Set some default options
  toset <- !names(x = default.options) %in% names(x = op)
  if (any(toset)) {
    options(default.options[toset])
  }
}
