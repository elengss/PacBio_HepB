### Generate mapping file from D3L.fa
echo "STARTED AT" 
date +”%H:%M:%S”
echo " "
ECHO "PROCESSING FOR HBV READS"
echo "Time for a coffee break"
### MAP TO HEP B D3L
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015; do 
echo aligning "$tmp" for HBV
pbmm2 align --sort --preset ISOSEQ /reference/D3L.fa /data/demultiplex.dT_bc"$tmp"RC_PB_3p--dT_PB_5p.hifi_reads.fastq.gz "$tmp"_D3L.bam; done

echo " "
### prepare HBV files for R
echo "EXTRACTING HBV DATA FROM BAM FILES FOR R"
echo "These files are the input for the quantify_viral.R script."
### prepare files for R
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015; do
echo processing "$tmp"
samtools view "$tmp"_D3L.bam -o "$tmp"_D3L.noheader.sam
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

