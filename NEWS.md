# Azimuth 0.5.0 (2023-11-04)

## Changes
- Azimuth ATAC
  - Functionality added to `RunAzimuth()` ([vignette here](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html#scatac-seq-queries)), to `AzimuthApp`, and to the [Azimuth Website](https://azimuth.hubmapconsortium.org/) to support sc/snATAC-seq queries. 
- Compatibility with Seurat v5

# Azimuth 0.4.6 (2022-08-23)

## Changes
- Compatibility fix for uwot >=v0.1.13 due to changes in UMAP model parameters

# Azimuth 0.4.5 (2022-04-27)

## Changes
- `RunAzimuth()` added to support local annotation.
- `LoadH5AD()` assumes matrix is CSR if metadata to resolve matrix type is unavailable.
- Bug fix for anndata files containing cells with 0 counts.
- Options to set default filtering thresholds (e.g. `Azimuth.app.ncount_max`)

# Azimuth 0.4.4 (2022-03-08)

## Changes
- `LoadH5AD()` now supports H5AD files containing CSC matrices at X or X/raw. Previously only CSR matrices were readable.

# Azimuth 0.4.3 (2021-08-18)

## Changes
- Bug fix for cell renaming when interactively subsetting
- Bug fix for references created with fewer than 50 dimensions

# Azimuth 0.4.2 (2021-07-28)

## Changes
- Demo dataset is now optional
- Add default hosted homolog file location
- Bug fix when query cell name(s) overlap with reference names
- Bug fix for Windows path specification for the reference directory

# Azimuth 0.4.1 (2021-06-01)

## Added
- Ability to download all Azimuth results on Downloads tab 

# Azimuth 0.4.0 (2021-05-04)

## Added
- Option to have `meta.data` in reference that will display on hover but not be transferable via `Azimuth.app.metadata_notransfer`
- Option to display `meta.data` table in heatmap form with `Azimuth.app.metatableheatmap`
- Option to switch default `DimPlot` display in the cell plots tab to show the query overlayed on the reference with `Azimuth.app.overlayedreference`
- Expose the `dims` parameter to `FindTransferAnchors` and `TransferData` as an option, `Azimuth.map.ndims`
- `ConvertGeneNames` function to convert human/mouse, ensemble gene names. Location of conversion table file given by `Azimuth.app.homologs`.
- Option to have multiple demo buttons. This can be configured via the `Azimuth.app.demodatasets` option.

# Azimuth 0.3.2 (2021-03-23)

## Added
- `LoadFileInput()` and `LoadReference()` are now exported functions

## Changes
- Fix to prevent error when no demo dataset is provided
- Bug fixes in template analysis script generation
- Fix neighbor graph used for clustering in `ClusterPreservationScore()`

# Azimuth 0.3.1 (2021-03-12)

## Changes
- Loading query datasets now filters out cells with zero counts and features with zero counts.
