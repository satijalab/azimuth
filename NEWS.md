# Azimuth 0.4.1 (2021-05-14)

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
