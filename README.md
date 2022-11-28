# PacBio_HepB

This is a repository of the scripts used to prepare the paper. <br/>
It aligns PacBio Isoseq reads from Hepatitis B infected/transfected hepatocyte cell lines and performs differential gene/transcript expression and chimera detection. <br/>
The directories are as follows <br/>
Reference - stores the hepatitis B fasta and GTF files. Human files are large and can be downloaded from the Gencode website (further details in ReadMe in directory) <br/>
Data - this is where the raw fastq files are stored<br/>
Align_and_Quantify - this aligns the reads to viral and human genomes and quantifies them<br/>
Differential - this performs differential expression and gene set enrichment analysis<br/>
Chimera - this extracts soft clipped reads, remaps them and then then remaps them to the human genome<br/></br>
The shell scripts set a variable called "home" to define the directory paths used.</br>
There is a small utility "Findhome.sh" which can find the value for that variable on the system running these scripts.</br>
On most systems scripts are run by adding "./" in front of the script name, for example: ./Findhome.sh</br>
Many Linux systems require the script be declared as executable using the "chmod" command, for example: chmod +x Findhome.sh</br><hr>

Packages used<br/>
1. R Version 4.1.2 with Bioconductor (https://www.bioconductor.org/)<br/>
2. Samtools (http://www.htslib.org/)<br/>
3. PBMM2 (https://github.com/PacificBiosciences/pbmm2)<br/>

Samtools and pbmm2 are both downloadable as conda packages:

conda install -c bioconda samtools<br/>
conda install -c bioconda pbmm2<br/><br/>

On some systems the conda package manager has to be loaded before running the shell scripts, on our system the command is :</br>
conda activate XXXX</br>
where "XXXX" is the name of the conda environment</br></br>

More details on the Conda and Anaconda package manager system can be found here : https://anaconda.org<hr>

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

There are currently some problems running edgeR on M1 and M2 Macs. EdgeR requires fortran libraries to run. This requires a new installation of the gcc package which has these libraries (these are in addition to the compilers installed with Xtools).
