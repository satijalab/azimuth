


bioc_pkgs <- c('BSgenome.Hsapiens.UCSC.hg38', 'glmGamPoi', 'GenomeInfoDb', 'GenomicRanges', 'TFBSTools', 'JASPAR2020', 'EnsDb.Hsapiens.v86', 'IRanges', 'Rsamtools')
BiocManager::install(bioc_pkgs,ask=F, update = F, 'lib = /usr/local/lib/R/site-library')
