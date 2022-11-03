#!/bin/bash
### Alignment script using pbmm2
###
### "BMRC" is the name of the server cluster we used for analysis.
### other systems will require local versions of these initiation calls
### Local information for BMRC cluster :
echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 
### BMRC specific calls to initiate pbmm2
### load conda
module load Anaconda3/2022.05
###
### These next 3 steps only need to be done once
### 1) generate a local conda environment for our programs to reside in
### conda create -n PBMM python=3
### 2) install the pbmm2 program into our local environment
### conda install -c bioconda pbmm2
### set up bash shell to use conda
### This last command may be requested by the script, but is already set
### 3) conda init bash
### Activate the conda environment where pbmm2 was installed
conda activate PBMM
###
### define local paths for reference and data folders
home="/gpfs3/well/mckeating/projects/Pacbio-Runs-June-HBV/Home"
echo "Home directory is :"$home
data_dir="$home/data"
ref_dir="$home/reference"
echo "Reference fasta files & gtf files in : $ref_dir"
echo ""
### MAP TO HEP B D3L
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015; do 
echo aligning "$tmp" for HBV
date +”%H:%M:%S”
pbmm2 align --sort --preset ISOSEQ $ref_dir/D3L.fa $data_dir/demultiplex.dT_bc"$tmp"RC_PB_3p--dT_PB_5p.hifi_reads.fastq.gz $data_dir/"$tmp"_D3L.bam; done

module load samtools
### prepare files for R
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015; do
echo making and extracting sam file "$tmp" for HBV
samtools view $data_dir/"$tmp"_D3L.bam -o $data_dir/"$tmp"_D3L.noheader.sam
awk '{print $6}' $data_dir/"$tmp"_D3L.noheader.sam > cigar
awk '{print $5}' $data_dir/"$tmp"_D3L.noheader.sam > mapq
awk '{print $4}' $data_dir/"$tmp"_D3L.noheader.sam > pos
awk '{print $2}' $data_dir/"$tmp"_D3L.noheader.sam > bit
awk '{print $1}' $data_dir/"$tmp"_D3L.noheader.sam > query
paste query mapq bit pos cigar > $data_dir/"$tmp".D3L.short; done
### tidy up files
rm query bit pos mapq cigar
echo " "
date +”%H:%M:%S”

### MAP TO HUMAN
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015; do pbmm2 align --sort --preset ISOSEQ $ref_dir/GRCh38.p13.genome.fa $data_dir/demultiplex.dT_bc"$tmp"RC_PB_3p--dT_PB_5p.hifi_reads.fastq.gz $data_dir/"$tmp"_prec_iso_fq_human.bam; done
### BMRC specific calls to load samtools
