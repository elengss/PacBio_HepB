# PacBio_HepB

Packages used
1. Samtools (http://www.htslib.org/)
2. PBMM2 (https://github.com/PacificBiosciences/pbmm2)
3. R Version 4.1.2 with Bioconductor (https://www.bioconductor.org/)

Bioconductor packages:
 [1] Rsubread_2.8.0              edgeR_3.36.0
 [3] limma_3.50.0                GenomicAlignments_1.30.0
 [5] Rsamtools_2.10.0            Biostrings_2.62.0
 [7] XVector_0.34.0              SummarizedExperiment_1.24.0
 [9] Biobase_2.54.0              MatrixGenerics_1.6.0
[11] matrixStats_0.61.0          GenomicRanges_1.46.0
[13] GenomeInfoDb_1.30.0         IRanges_2.28.0
[15] S4Vectors_0.32.0            BiocGenerics_0.40.0
[16] clusterProfiler_4.2.2

They can be installed and loaded with these commands in R
bio_pkgs <- c("Rsubread","edgeR","limma","GenomicAlignments","Rsamtools","Biostrings","XVector","SummarizedExperiment","Biobase","MatrixGenerics","matrixStats","GenomicRanges","GenomeInfoDb","IRanges","S4Vectors","BiocGenerics","clusterProfiler")
# install:
BiocManager::install(bio_pkgs)
# load all at once
invisible(lapply(pkgs, function(x) library(x, character.only=TRUE)))
