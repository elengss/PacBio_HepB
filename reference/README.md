<h3>This folder contains the reference files for the HBV genotype D clone (GenBank ID NC_003977.2) used in the experiments.</h3>

<p>The D3L fasta file used as our referent is 3581 bp and runs from the first TSS for preCore (preC-L, 1525) to the TATAAA Poly Adenylation Site (PAS, 1918-1923).</br> 
The coordinates used in the GTF file are adjusted accordingly. The original GTF file of spliced transcripts is from the following paper :</br>
Quantitative analysis of the splice variants expressed by the major hepatitis B virus genotypes</br>
C. S. Lim, V. Sozzi, M. Littlejohn, L. K. W. Yuen, N. Warner, B. Betz-Stablein, et al.</br>
Microb Genom 2021  DOI: 10.1099/mgen.0.000492. https://www.ncbi.nlm.nih.gov/pubmed/33439114v</br></br>
To it we have added the six major canonical transcripts, as identified here :</br>
Single-Nucleotide Resolution Mapping of Hepatitis B Virus Promoters in Infected Human Livers and Hepatocellular Carcinoma</br>
K. Altinel, K. Hashimoto, Y. Wei, C. Neuveut, I. Gupta, A. M. Suzuki, et al.</br>
J Virol 2016 <u>90</u>:10811-10822 DOI: 10.1128/JVI.01625-16 https://www.ncbi.nlm.nih.gov/pubmed/27681123</br>
</p>

The sequence of the plasmid backbone (pUC57.fa) used for the HBV1.3 construct is also included.</br>

The Ensemble database annotation file for the human transcriptome (ENST_ENSG_GN), used in the differential expression analysis 
(differential.R) and for identification of chimeric transcripts (quantify_human_chimare.R) is also included.

Human fasta and GTF files should be downloaded from here :

https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.p13.genome.fa.gz
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz

The downloaded .gz archive files should be decompressed (gunzip utility) before use to give the files GRCh38.p13.genome.fa and gencode.v40.chr_patch_hapl_scaff.annotation.gtf
