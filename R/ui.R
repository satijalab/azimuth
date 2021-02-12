#' @include zzz.R
#'
#' @importFrom DT DTOutput
#' @importFrom htmltools div h3 h4 HTML includeCSS p tagList tags
#' @importFrom shinyjs disabled useShinyjs
#' @importFrom shiny actionButton checkboxInput downloadButton fileInput
#' fluidRow htmlOutput icon numericInput plotOutput radioButtons selectizeInput
#' tableOutput textOutput uiOutput verbatimTextOutput
#' @importFrom shinyBS bsButton bsPopover
#' @importFrom shinydashboard box dashboardBody dashboardHeader dashboardSidebar
#' dashboardPage menuItem sidebarMenu sidebarMenuOutput tabItem tabItems
#' valueBoxOutput
#'
NULL

AzimuthUI <- tagList(
  useShinyjs(),
  includeCSS(path = GetCSS()),
  dashboardPage(
    dashboardHeader(title = app.title),
    dashboardSidebar(
      fileInput(
        inputId = 'file',
        label = p(
          'File Upload',
          bsButton(
            'q1',
            label = '',
            icon = icon(name = 'question'),
            style = 'info',
            size = 'extra-small'
          )
        ),
        accept = c('.h5', '.h5ad', '.h5seurat', '.rds')
      ),
      bsPopover(
        id = 'q1',
        title = 'Supported file types',
        content = paste(
          '10x Genomics H5',
          'Seurat object (RDS)',
          'H5AD',
          'H5Seurat',
          'Matrix/matrix/data.frame (RDS)',
          sep = '; '
        ),
        placement = 'right',
        trigger = 'focus',
        options = list(container = 'body')
      ),
      actionButton(inputId = 'triggerdemo', label = 'Load demo dataset'),
      htmlOutput(outputId = 'message', inline = FALSE),
      sidebarMenu(
        menuItem(
          text = 'Welcome',
          tabName = 'tab_welcome'
        ),
        sidebarMenuOutput(outputId = 'menu1'),
        sidebarMenuOutput(outputId = 'menu2')
      ),
      htmlOutput(outputId = 'containerid', inline = FALSE)
    ),
    dashboardBody(
      tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }'))),
      tabItems(
        # Welcome tab
        tabItem(
          tabName = 'tab_welcome',
          htmlOutput(outputId = 'welcomebox')
        ),
        # Preprocessing + QC Tab
        tabItem(
          tabName = 'tab_preproc',
          fluidRow(
            box(
              title = p(
                'QC Filters',
                bsButton(
                  inputId = 'q2',
                  label = '',
                  icon = icon(name = 'question'),
                  style = 'info',
                  size = 'extra-small'
                )
              ),
              div(
                id = 'ncount',
                numericInput(
                  inputId = 'num.ncountmin',
                  label = NULL,
                  value = 0,
                  width = '90%'
                ),
                numericInput(
                  inputId = 'num.ncountmax',
                  label = NULL,
                  value = 0,
                  width = '90%'
                )
              ),
              div(
                id = 'nfeature',
                numericInput(
                  inputId = 'num.nfeaturemin',
                  label = NULL,
                  value = 0,
                  width = '90%'
                ),
                numericInput(
                  inputId = 'num.nfeaturemax',
                  label = NULL,
                  value = 0,
                  width = '90%'
                )
              ),
              div(
                id = 'pctmt',
                numericInput(
                  inputId = 'minmt',
                  label = NULL,
                  value = 0,
                  width = '90%'
                ),
                numericInput(
                  'maxmt',
                  label = NULL,
                  value = 0,
                  width = '90%'
                )
              ),
              textOutput(outputId = 'text.cellsremain'),
              div(
                id = 'xferopts',
                h4(
                  'Transfer Options',
                  bsButton(
                    inputId = 'xferinput',
                    label = '',
                    icon = icon(name = 'question'),
                    style = 'info',
                    size = 'extra-small'
                  )
                ),
                bsPopover(
                  id = 'xferinput',
                  title = 'Transfer Options',
                  content = 'Select the meta.data fields to transfer from the reference',
                  placement = 'right',
                  trigger = 'focus',
                  options = list(container = 'body')
                ),
                selectizeInput(
                  inputId = 'metadataxfer',
                  label = 'Reference Metadata to Transfer',
                  choices = '',
                  multiple = TRUE
                )
              ),
              disabled(actionButton(
                inputId = 'map',
                label = 'Map cells to reference'
              )),
              width = 4
            ),
            bsPopover(
              id = 'q2',
              title = 'QC Filters',
              content = paste(
                'Select a minimum and maximum value for nCount (number of molecules)',
                'nFeature (number of genes expressed)',
                'and mitochondrial percentage (if applicable)',
                sep = ', '
              ),
              placement = 'right',
              trigger = 'focus',
              options = list(container = 'body')
            ),
            box(
              checkboxInput(inputId = 'check.qcscale', label = 'Log-scale Y-axis'),
              checkboxInput(inputId = 'check.qcpoints', label = 'Hide points'),
              plotOutput(outputId = 'plot.qc'),
              tableOutput(outputId = 'table.qc'),
              width = 8
            )
          ),
          fluidRow(
            valueBoxOutput(outputId = 'valuebox.upload', width = 3),
            valueBoxOutput(outputId = 'valuebox.preproc', width = 3),
            valueBoxOutput(outputId = 'valuebox.mapped', width = 3)
          ),
        ),
        # Cell tab
        tabItem(
          tabName = 'tab_cell',
          box(
            title = 'Reference',
            checkboxInput(inputId = 'labels', label = 'Show labels'),
            selectizeInput(
              inputId = 'metacolor.ref',
              label = 'Metadata to color by',
              choices = '',
              multiple = TRUE,
            ),
            plotOutput(outputId = 'refdim'),
            width = 12
          ),
          box(
            title = 'Query',
            selectizeInput(
              inputId = 'metacolor.query',
              label = 'Metadata to color by',
              choices = '',
              multiple = TRUE,
            ),
            plotOutput(outputId = 'objdim'),
            width = 12
          ),
          box(
            title = 'Metadata table',
            div(
              style = 'display: inline-block; vertical-align: top; width: 25%',
              selectizeInput(
                inputId = 'metarow',
                label = 'Table rows',
                choices = ''
              )
            ),
            div(
              style = 'display: inline-block; vertical-align: top; width: 25%',
              selectizeInput(
                inputId = 'metacol',
                label = 'Table columns',
                choices = ''
              )
            ),
            div(
              style = 'display: inline-block; vertical-align: top; width: 50%',
              radioButtons(
                inputId = 'radio.pct',
                label = NULL,
                choices = c('Percentage','Frequency'),
                inline = TRUE
              )
            ),
            tableOutput(outputId = 'table.metadata'),
            width = 12
          )
        ),
        # Feature tab
        tabItem(
          tags$head(tags$style(HTML(".selectize-dropdown .optgroup-header { font-weight: bold; font-size: 13px; color: black; background: #f6f6f6}"))),
          tabName = 'tab_feature',
          box(
            title = 'Feature Plots',
            div(
              id = 'featureinput',
              class = 'thirds',
              selectizeInput(
                inputId = 'feature',
                label = 'Feature',
                choices = ''
              )
            ),
            div(
              id = 'imputedinput',
              class = 'thirds',
              selectizeInput(
                inputId = 'adtfeature',
                label = 'Imputed protein',
                choices = ''
              )
            ),
            div(
              id = 'continput',
              class = 'thirds',
              selectizeInput(
                inputId = 'metadata.cont',
                label = 'Prediction Scores and Metadata',
                choices = ''
              )
            ),
            plotOutput(outputId = 'edim'),
            selectizeInput(
              inputId = 'metagroup',
              label = 'Metadata to group by',
              choices = '',
              width = '25%'
            ),
            checkboxInput(inputId = 'check.featpoints', label = 'Hide points'),
            plotOutput(outputId = 'evln'),
            width = 12
          ),
          box(
            title = p(
              'Predicted cell type cluster biomarkers',
              bsButton(
                inputId = 'q3',
                label = '',
                icon = icon(name = 'question'),
                style = 'info',
                size = 'extra-small'
              )
            ),
            bsPopover(
              id = 'q3',
              title = 'Biomarkers Table',
              content = paste(
                'Only available for clusters with at least 15 cells.',
                paste(
                  # 'logFC: log fold-change between cells in the cluster specified and other cells',
                  'auc: area under ROC',
                  'padj: Benjamini-Hochberg adjusted p value',
                  'pct_in: percent of cells in the cluster with nonzero feature value',
                  'pct_out: percent of cells out of the cluster with nonzero feature value',
                  sep = '; '
                )
              ),
              placement = 'right',
              trigger = 'focus',
              options = list(container = 'body')
            ),
            div(
              id = 'markerclustersgroupinput',
              class = 'halves',
              selectizeInput(
                inputId = 'markerclustersgroup',
                label = 'Metadata group',
                choices = ''
              )
            ),
            div(
              id = 'markerclustersgroupinput',
              class = 'halves',
              selectizeInput(
                inputId = 'markerclusters',
                label = 'Predicted cell type',
                choices = ''
              )
            ),
            div(
              id = 'biotable',
              class = 'halves',
              h3('RNA biomarkers'),
              DTOutput(outputId = 'biomarkers')
            ),
            div(
              id = 'imputedtable',
              class = 'halves',
              uiOutput(outputId = 'imputedlabel'),
              DTOutput(outputId = 'adtbio')
            ),
            width = 12
          )
        ),
        # Downloads tab
        tabItem(
          tabName = 'tab_download',
          div(
            id = 'scriptdl',
            box(
              title = 'Analysis script template ',
              downloadButton(
                outputId = 'dlscript',
                label = 'Download'
              ),
              width = 6
            )
          ),
          div(
            id = 'umapdl',
            box(
              title = 'UMAP (Seurat Reduction RDS)',
              verbatimTextOutput(outputId = 'text.dlumap'),
              downloadButton(
                outputId = 'dlumap',
                label = 'Download'
              ),
              width = 6
            )
          ),
          div(
            id = 'imputeddl',
            box(
              title = 'Imputed protein (Seurat Assay RDS)',
              verbatimTextOutput(outputId = 'text.dladt'),
              downloadButton(
                outputId = 'dladt',
                label = 'Download'
              ),
              width = 6
            )
          ),
          div(
            id = 'predictionsdl',
            box(
              title = 'Predicted cell types and scores (TSV)',
              verbatimTextOutput(outputId = 'text.dlpred'),
              downloadButton(
                outputId = 'dlpred',
                label = 'Download'
              ),
              width = 6
            )
          )
        )
      )
    )
  )
)
