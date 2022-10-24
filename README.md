# PacBio_HepB

This is a repository of the scripts used to prepare the paper. <br/>
It aligns PacBio Isoseq reads from Hepatitis B infected/transfected hepatocyte cell lines and performs differential gene/transcript expression and chimera detection. <br/>
The directories are as follows <br/>
Reference - stores the hepatitis B fasta and GTF files. Human files are large and can be downloaded from the Gencode website (further details in ReadMe in directory) <br/>
Data - this is where the raw fastq files are stored<br/>
Align_and_Quantify - this aligns the reads to viral and human genomes and quantifies them<br/>
Differential - this performs differential expression and gene set enrichment analysis<br/>
Chimera - this extracts soft clipped reads, remaps them and then then remaps them to the human genome<br/>

Packages used<br/>
1. R Version 4.1.2 with Bioconductor (https://www.bioconductor.org/)<br/>
2. Samtools (http://www.htslib.org/)<br/>
3. PBMM2 (https://github.com/PacificBiosciences/pbmm2)<br/>

Samtools and pbmm2 are both downloadable as conda packages:

conda install -c bioconda samtools<br/>
conda install -c bioconda pbmm2<br/>

More details on the Conda and Anaconda package manager system can be found here : https://anaconda.org

Bioconductor packages:<br/>
Rsubread_2.8.0<br/>
edgeR_3.36.0<br/>
limma_3.50.0<br/>
GenomicAlignments_1.30.0<br/>
Rsamtools_2.10.0<br/>
Biostrings_2.62.0<br/>
XVector_0.34.0<br/>
SummarizedExperiment_1.24.0<br/>
Biobase_2.54.0<br/>
MatrixGenerics_1.6.0<br/>
matrixStats_0.61.0<br/>
GenomicRanges_1.46.0<br/>
GenomeInfoDb_1.30.0<br/>
IRanges_2.28.0<br/>
S4Vectors_0.32.0<br/>
BiocGenerics_0.40.0<br/>
clusterProfiler_4.2.2<br/>

They can be installed and loaded with this command in R : <br/><br/>
bio_pkgs <- c("Rsubread","edgeR","limma","GenomicAlignments","Rsamtools","Biostrings","XVector","SummarizedExperiment","Biobase","MatrixGenerics","matrixStats","GenomicRanges","GenomeInfoDb","IRanges","S4Vectors","BiocGenerics","clusterProfiler")<br/>
BiocManager::install(bio_pkgs)<br/><br/>
Cut and paste tha above command into R or Rstudio to install all the required packages.
