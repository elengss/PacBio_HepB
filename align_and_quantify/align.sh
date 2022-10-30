### Generate mapping file from D3L.fa
echo "STARTED AT" 
date +”%H:%M:%S”
echo " "
ECHO "PROCESSING FOR HBV READS"
echo "Time for a coffee break"
### MAP TO HEP B D3L
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015; do 
echo aligning "$tmp" for HBV
### unpack the .gz archive
gunzip /data/demultiplex.dT_bc"$tmp"RC_PB_3p--dT_PB_5p.hifi_reads.fastq.gz
### Count the lines in the file and store in the fle "fastqlines"
#### fastqlines is a file that contains a list of the number of lines in each fastq file, which can be divided by 4 to get number of reads for normalization in downstream steps
wc -l /data/demultiplex.dT_bc"$tmp"RC_PB_3p--dT_PB_5p.hifi_reads.fastq >> fastqlines
### align sequences to D3L.fa referent genome :
pbmm2 align --sort --preset ISOSEQ /reference/D3L.fa /data/demultiplex.dT_bc"$tmp"RC_PB_3p--dT_PB_5p.hifi_reads.fastq "$tmp"_D3L.bam; done
echo " "
echo "Files processed, HBV aignments complete, number of lines in each fastq file are :"
cat fastqlines
echo " "
### prepare HBV files for R
echo "EXTRACTING HBV DATA FROM BAM FILES FOR R"
echo "These files are the input for the quantify_viral.R script."
### prepare files for R
### File for recording line numbers in each file (# sequences is 1/4 of number of lines)
echo "Length   Filename" > linenumbers
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015; do
echo processing "$tmp"
samtools view "$tmp"_D3L.bam -o "$tmp"_D3L.noheader.sam
wc -l "$tmp"_D3L.noheader.sam >> linenumbers
awk '{print $6}' "$tmp"_D3L.noheader.sam > cigar
awk '{print $5}' "$tmp"_D3L.noheader.sam > mapq
awk '{print $4}' "$tmp"_D3L.noheader.sam > pos
awk '{print $2}' "$tmp"_D3L.noheader.sam > bit
awk '{print $1}' "$tmp"_D3L.noheader.sam > query
paste query mapq bit pos cigar > "$tmp".D3L.short; done
echo 'R files generated'
### tidy up files
rm query bit pos mapq cigar
echo " "
echo "Files processed, number of lines in each SAM file are :"
cat linenumbers
echo " "
date +”%H:%M:%S”
echo " "
echo "PROCESSING FOR HUMAN READS, this is the slow bit!"
echo "More like lunch than a coffee break"
echo " "
### MAP TO HUMAN
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015; do 
echo aligning "$tmp" for Human
pbmm2 align --sort --preset ISOSEQ /reference/GRCh38.p13.genome.fa /data/demultiplex.dT_bc"$tmp"RC_PB_3p--dT_PB_5p.hifi_reads.fastq.gz "$tmp"_prec_iso_fq_human.bam; done
echo " "
echo "SCRIPT FINISHED"
echo " "
