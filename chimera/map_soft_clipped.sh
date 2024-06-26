#!/bin/bash
### Soft clipping of reads script using pbmm2 and samtools
###
### "BMRC" is the name of the server cluster we used for analysis.
### other systems will require local versions of these initiation calls
### This is the call in terminal on our system (slurm interface) :
### sbatch --mem=0 -p short map_soft_clipped.sh
### The --mem=0 parameter sets memory to the maximum available 
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
eval "$(conda shell.bash hook)"
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
### use the utility "Findhome.sh" to obtain the value for this variable.
home="/gpfs3/well/mckeating/projects/Pacbio-Runs-June-HBV/Home"
echo "Home directory is : "$home
data_dir="$home/data"
ref_dir="$home/reference"
echo "Reference fasta files & gtf files in : "$ref_dir
echo ""
echo "the fastq positive and negative data file names generated by the R script "
echo "extract_soft_clipped.R"
echo " "
echo "As before, the barcodes are in the variable tmp"
echo " "
### BMRC specific calls to load samtools
module load samtools
### concatenate human, viral and plasmid reads into a single fasta file
cat $ref_dir/GRCh38.p13.genome.fa $ref_dir/D3L.fa $ref_dir/pUC57.fa > $ref_dir/humanviralplasmid.fa

for tmp in 1002 1003 1010 1012
do
echo Extracting fastq clipping files and aligning "$tmp"
pbmm2 align --sort --preset ISOSEQ $ref_dir/humanviralplasmid.fa $data_dir/"$tmp".positive.start.fastq $data_dir/"$tmp".positive.start.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ $ref_dir/humanviralplasmid.fa $data_dir/"$tmp".positive.end.fastq $data_dir/"$tmp".positive.end.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ $ref_dir/humanviralplasmid.fa $data_dir/"$tmp".negative.start.fastq $data_dir/"$tmp".negative.start.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ $ref_dir/humanviralplasmid.fa $data_dir/"$tmp".negative.end.fastq $data_dir/"$tmp".negative.end.humanviralplasmid.bam
done


module load samtools
for tmp in 1002 1003 1010 1012
do
echo Removing short (<100 bp) sequences and making sam files of "$tmp"
samtools view -h $data_dir/"$tmp".positive.start.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.bam
samtools view -h $data_dir/"$tmp".positive.end.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.bam
samtools view -h $data_dir/"$tmp".negative.start.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.bam
samtools view -h $data_dir/"$tmp".negative.end.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.bam
done

for tmp in 1002 1003 1010 1012
do
echo Indexing bam files of $data_dir/"$tmp" 
samtools index $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.bam
samtools index $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.bam
samtools index $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.bam
samtools index $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.bam
done

for tmp in 1002 1003 1010 1012
do
echo Filtering plasmid vector sequences in bam files of "$tmp" 
samtools view -h $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.bam pUC57 > $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.plasmid.bam
samtools view -h $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.bam pUC57 > $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.plasmid.bam
samtools view -h $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.bam pUC57 > $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.plasmid.bam
samtools view -h $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.bam pUC57 > $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.plasmid.bam
done

for tmp in 1002 1003 1010 1012
do
echo Removing low quality reads in bam files of "$tmp" 
samtools view -bq 30 $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.plasmid.bam > $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.plasmid.bam > $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.plasmid.bam > $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.plasmid.bam > $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.plasmid30.bam
done

for tmp in 1002 1003 1010 1012
do
samtools view -o $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.plasmid30.sam $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -o $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.plasmid30.sam $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.plasmid30.bam
samtools view -o $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.plasmid30.sam $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -o $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.plasmid30.sam $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.plasmid30.bam
done
### Remove counter files if pre-existing...
rm >> posstartlengths posendlengths negstartlengths negendlengths
for tmp in 1002 1003 1010 1012
do
echo Counting sam file linenumbers for $data_dir/"$tmp" 
wc -l $data_dir/"$tmp".positive.start.humanviralplasmid.filt100.plasmid30.sam >> posstartlengths
wc -l $data_dir/"$tmp".positive.end.humanviralplasmid.filt100.plasmid30.sam >> posendlengths
wc -l $data_dir/"$tmp".negative.start.humanviralplasmid.filt100.plasmid30.sam >> negstartlengths
wc -l $data_dir/"$tmp".negative.end.humanviralplasmid.filt100.plasmid30.sam >> negendlengths
done
echo " "
date +”%H:%M:%S”
echo " "
echo "Script finished."
