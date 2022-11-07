#' @docType package
#'
#' @section Package options:
#'
#' \pkg{Azimuth} uses the following options to control the behavior of the app.
#' Users can provide these as named arguments to \code{\link{AzimuthApp}}
#' through dots (...), specify these in the config file, or configure these with
#' \code{\link[base]{options}}.
#'
#' \subsection{App options}{
#'  The following options control app behavior
#'  \describe{
#'   \item{\code{Azimuth.app.default_adt}}{
#'    ADT to select by default in feature/violin plot
#'   }
#'   \item{\code{Azimuth.app.default_gene}}{
#'    Gene to select by default in feature/violin plot
#'   }
#'   \item{\code{Azimuth.app.default_metadata}}{
#'    Default metadata transferred from reference.
#'   }
#'   \item{\code{Azimuth.app.demodataset}}{
#'    Path to data file (in any Azimuth-supported format) to automatically load
#'    when the user clicks a button. The button is only available in the UI
#'    if this option is set
#'   }
#'   \item{\code{Azimuth.app.googlesheet}}{
#'    Google Sheet identifier (appropriate for use with
#'    \code{\link[googlesheets4:gs4_get]{googlesheets4:gs4_get}}) to write log
#'    records. Logging is only enabled if this and other \code{google*} options
#'    are set
#'   }
#'   \item{\code{Azimuth.app.googletoken}}{
#'    Path to directory containing Google Authentication token file.
#'    Logging is only enabled if this and other \code{google*} options are set
#'   }
#'   \item{\code{Azimuth.app.googletokenemail}}{
#'    Email address corresponding to the Google Authentication token file.
#'    Logging is only enabled if this and other \code{google*} options are set
#'   }
#'   \item{\code{Azimuth.app.max_cells}}{
#'    Maximum number of cells allowed to upload
#'   }
#'   \item{\code{Azimuth.app.metadata_notransfer}}{
#'    Metadata to annotate in reference but not transfer to query
#'   }
#'   \item{\code{Azimuth.app.mito}}{
#'    Regular expression pattern indicating mitochondrial features in query
#'    object
#'   }
#'   \item{\code{Azimuth.app.plotseed}}{
#'    Seed to shuffle colors for cell types
#'   }
#'   \item{\code{Azimuth.app.reference}}{
#'    URL or directory path to reference dataset; see
#'    \code{\link{LoadReference}} for more details
#'   }
#'   \item{\code{Azimuth.app.refuri}}{
#'    URL for publicly available reference dataset, used in the downloadable
#'    analysis script in case \code{Azimuth.app.reference} points to a directory
#'   }
#'   \item{\code{Azimuth.app.refdescriptor}}{
#'    Provide (as a string) the html to render the reference description on the
#'    welcome page
#'   }
#'   \item{\code{Azimuth.app.welcomebox}}{
#'    Provide (as a string) the code to render the box on the welcome page
#'    (quotes escaped). Example:
#'    ```
#'    box(
#'      h3(\"Header\"),
#'      \"body text\",
#'      a(\"link\", href=\"www.satijalab.org\", target=\"_blank\"),
#'      width = 12
#'    )
#'    ```
#'   }
#'   \item{\code{Azimuth.app.homologs}}{
#'    URL or path to file containing the human/mouse homolog table.
#'   }
#'   \item{\code{Azimuth.app.metatableheatmap}}{
#'    Display the meta.data table as a heatmap rather than in tabular form.
#'    defaults to FALSE.
#'   }
#'   \item{\code{Azimuth.app.overlayedreference}}{
#'    Display the mapped query on top of greyed out reference in the 'Cell
#'    Plots' tab.
#'    defaults to FALSE
#'   }
#'  }
#' }
#'
#' \subsection{Control options}{
#'  These options control mapping and analysis behavior
#'  \describe{
#'   \item{\code{Azimuth.map.ncells}}{
#'    Minimum number of cells required to accept uploaded file
#'    defaults to \code{100}
#'   }
#'   \item{\code{Azimuth.map.ngenes}}{
#'    Minimum number of genes in common with reference to accept uploaded file;
#'    defaults to \code{250}
#'   }
#'   \item{\code{Azimuth.map.nanchors}}{
#'    Minimum number of anchors that must be found to complete mapping.
#'    Defaults to \code{50}
#'   }
#'   \item{\code{Azimuth.map.panchorscolors}}{
#'    Configure the valuebox on the main page corresponding to the values for
#'    failure, warning, success for fraction of unique query cells that
#'    participate in anchor pairs. Failure corresponds to
#'    [0:\code{Azimuth.map.fracanchorscolors[1]}), warning to
#'    [\code{Azimuth.map.fracanchorscolors[1]}:\code{Azimuth.map.fracanchorscolors[2]}),
#'    and success is >= \code{Azimuth.map.fracanchorscolors[2]}.
#'    Defaults to \code{c(5, 15)}
#'   }
#'   \item{\code{Azimuth.map.postmapqccolors}}{
#'    Configure the valuebox on the main page corresponding to the values for
#'    failure, warning, success for the post mapping cluster based QC metric.
#'    Failure corresponds to [0:\code{Azimuth.map.postmapqc[1]}), warning to
#'    [\code{Azimuth.map.postmapqc[1]}:\code{Azimuth.map.postmapqc[2]}),
#'    and success is >= \code{Azimuth.map.postmapqc[2]}.
#'    Defaults to \code{c(0.15, 0.25)}
#'   }
#'   \item{\code{Azimuth.map.postmapqcds}}{
#'    Set the amount of query random downsampling to perform before computing
#'    the mapping QC metric.
#'    Defaults to \code{5000}
#'   }
#'   \item{\code{Azimuth.map.ntrees}}{
#'    Annoy (approximate nearest neighbor) n.trees parameter
#'    Defaults to \code{20}
#'   }
#'   \item{\code{Azimuth.map.ndims}}{
#'     Number of dimensions to use in FindTransferAnchors and TransferData
#'     Defaults to \code{50}
#'   }
#'   \item{\code{Azimuth.de.mincells}}{
#'    Minimum number of cells per cluster for differential expression; defaults
#'    to \code{15}
#'   }
#'   \item{\code{Azimuth.de.digits}}{
#'    Number of digits to round differential expression table to; defaults to
#'    \code{3}
#'   }
#'   \item{\code{Azimuth.sct.ncells}, \code{Azimuth.sct.nfeats}}{
#'    Number of cells and features to use for
#'    \code{\link[Seurat]{SCTransform}}, respectively. Defaults to \code{2000}
#'    for each
#'   }
#'  }
#' }
#'
#' \subsection{External options}{
#'  The following options are used by external dependencies that have an effect
#'  on \pkg{Azimuth}'s behavior. Refer to original package documentation for
#'  more details
#'  \describe{
#'   \item{\code{\link[shiny:shiny-options]{shiny.maxRequestSize}}}{
#'    User-configurable; used for controlling the maximum file size of uploaded
#'    datasets. Defaults to 500 Mb
#'   }
#'   \item{\code{\link[DT:datatable]{DT.options}}}{
#'   User-configurable; used for controlling biomarker table outputs.
#'   Defaults to setting \code{pageLength} to \code{10}
#'   }
#'   \item{\code{\link[future:future.options]{future.globals.maxSize}}}{
#'    \strong{Non-configurable}; used for parallelization. Defaults to
#'    \code{Azimuth.app.max_cells * 320000}
#'   }
#'  }
#' }
#'
#' @md
#'
"_PACKAGE"

app.title <- 'Azimuth'

default.options <- list(
  Azimuth.app.default_adt = "CD3-1",
  Azimuth.app.default_gene = "GNLY",
  Azimuth.app.default_metadata = NULL,
  Azimuth.app.demodataset = NULL,
  Azimuth.app.max_cells = 50000,
  Azimuth.app.mito = '^MT-',
  Azimuth.app.plotseed = NULL,
  Azimuth.app.reference = 'https://seurat.nygenome.org/references/pbmc',
  Azimuth.app.welcomebox = '',
  Azimuth.app.homologs = 'https://seurat.nygenome.org/azimuth/references/homologs.rds',
  Azimuth.app.refdescriptor = "",
  Azimuth.app.metatableheatmap = FALSE,
  Azimuth.app.overlayedreference = FALSE,
  Azimuth.app.ncount_min = NULL,
  Azimuth.app.ncount_max = NULL,
  Azimuth.app.nfeature_min = NULL,
  Azimuth.app.nfeature_max = NULL,
  Azimuth.app.pctmt_min = NULL,
  Azimuth.app.pctmt_max = NULL,
  Azimuth.de.digits = 3L,
  Azimuth.de.mincells = 15L,
  Azimuth.map.ncells = 100L,
  Azimuth.map.ngenes = 250L,
  Azimuth.map.nanchors = 50L,
  Azimuth.map.panchorscolors = c(5, 15),
  Azimuth.map.postmapqccolors = c(2, 3.75),
  Azimuth.map.postmapqcds = 5000L,
  Azimuth.map.ntrees = 20L,
  Azimuth.map.ndims = 50L,
  Azimuth.sct.ncells = 2000L,
  Azimuth.sct.nfeats = 2000L
)

qc.ids <- c(
  'map',
  'num.ncountmin',
  'num.ncountmax',
  'num.nfeaturemin',
  'num.nfeaturemax',
  'minmt',
  'maxmt',
  'check.qcscale',
  'check.qcpoints'
)

selectize.opts <- list(
  maxOptions = 1000L,
  maxItems = 1L
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
#' @importFrom SeuratObject FetchData
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

#' Format Time Differences
#'
#' @param dt A \code{\link[base]{difftime}} object
#'
#' @return The time difference in a nice string
#'
#' @keywords internal
#'
#' @seealso \code{\link[base:difftime]{base::difftime}}
#'
FormatDiffTime <- function(dt) {
  if (!inherits(x = dt, what = 'difftime')) {
    stop("'df' must be a difftime object")
  }
  dtfmt <- ifelse(
    test = dt < 60,
    yes = 'in %S seconds',
    no = 'in %M minutes %S seconds'
  )
  return(gsub(
    pattern = ' 0',
    replacement = ' ',
    x = format(x = .POSIXct(xx = dt), format = dtfmt),
    fixed = TRUE
  ))
}

#' Get Azimuth's CSS file
#'
#' Helper function to pull the location of Azimuth's CSS file
#'
#' @return The path to Azimuth's CSS file
#'
#' @keywords internal
#'
GetCSS <- function() {
  css <- system.file('www', 'azimuth.css', package = 'Azimuth')
}

#' Return names of metadata columns in a Seurat object that have an
#' appropriate number of levels for plotting when converted to a factor
#'
#' @param object a Seurat object
#' @param exceptions vector of metadata names to explicitly allow
#' @param min.levels minimum number of levels in a metadata factor to include
#' @param max.levels maximum number of levels in a metadata factor to include
#'
#' @keywords internal
#'
PlottableMetadataNames <- function(
  object,
  exceptions,
  min.levels = 2,
  max.levels = 20
) {
  column.status <- sapply(
    X = object[[]],
    FUN = function(column) {
      length(x = levels(x = droplevels(x = as.factor(x = column)))) >= min.levels &&
        length(x = levels(x = droplevels(x = as.factor(x = column)))) <= max.levels
    }
  ) & ! (grepl(pattern = ".score$", x = colnames(x = object[[]]))) |
    (grepl(pattern = "^predicted.", x = colnames(x = object[[]])) &
    ! (grepl(pattern = ".score$", x = colnames(x = object[[]])))) |
    colnames(x = object[[]]) %in% exceptions
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
  # cols.keep <- c('logFC', 'auc', 'padj', 'pct_in', 'pct_out')
  print("Rendering differential expression")
  cols.keep <- c('auc', 'padj', 'pct_in', 'pct_out')
  if (is.null(diff.exp)){
    print("Differential Expression is empty ")
  }
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
  # the things might be characters so they cant be ordered? 
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


#' Prepare differential expression motif results for rendering
#'
#' @param diff.exp A dataframe with differential expression results from 
#' \code{\link[Seurat:FindAllMarkers]{Seurat::FindAllMarkers}}
#' @param groups.use Names of groups to filter \code{diff.exp} to; groups must
#' be found in \code{diff.exp$cluster}
#' @param n Number of feature to filter \code{diff.exp} to per group
#'
#' @return \code{diff.exp}, ordered by adjusted p-value, filtered to \code{n}
#' features per group in \code{group.use}
#'
#' @importFrom rlang %||%
#' @importFrom utils head
#'
#' @seealso \code{\link[Seurat]{FindAllMarkers}}
#'
#' @keywords internal
#'
RenderDiffMotifExp <- function(
    diff.exp,
    groups.use = NULL,
    n = 10L
    #logfc.thresh = 0L
) {
  # cols.keep <- c('logFC', 'auc', 'padj', 'pct_in', 'pct_out')
  print("Rendering differential motifexpression")
  cols.keep <- c('avg_diff', 'p_val_adj', 'pct.1', 'pct.2')
  if (is.null(diff.exp)){
    print("Differential Expression is empty ")
  }
  print(head(diff.exp))
  groups.use <- groups.use %||% unique(x = as.character(x = diff.exp$cluster))
  diff.exp <- lapply(
    X = groups.use,
    FUN = function(cluster) {
      cluster.de <- diff.exp[diff.exp$cluster == cluster, , drop = FALSE]
      #cluster.de <- cluster.de[cluster.de$logFC > logfc.thresh, , drop = FALSE]
      cluster.de <- cluster.de[order(cluster.de$p_val_adj, -cluster.de$avg_diff), , drop = FALSE]
      return(head(x = cluster.de, n = n))
    }
  )
  # the things might be characters so they cant be ordered? 
  diff.exp <- do.call(what = 'rbind', diff.exp)
  rownames(x = diff.exp) <- make.unique(names = diff.exp$gene)
  diff.exp.num <- signif(
    x = diff.exp[, cols.keep, drop = FALSE],
    digits = getOption(
      x = "Azimuth.de.digits",
      default = default.options$Azimuth.de.digits
    )
  )
  diff.exp.num$motif_id <- diff.exp$motif_id
  print("FINALIZED DIFF EXP")
  print(head(diff.exp))
  return(diff.exp.num)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  # Attach dependencies
  AttachDeps()
  op <- options()
  # Set default options
  toset <- !names(x = default.options) %in% names(x = op)
  if (any(toset)) {
    options(default.options[toset])
  }
}
