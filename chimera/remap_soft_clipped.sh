### Align script for all 3 reference file types (human, viral and plasmid).
echo "Aligning reads and extracting chimeric regions."
echo " "
echo "STARTED AT" 
date +”%H:%M:%S”
echo "concatenate referent fasta files "
echo "cat /reference/GRCh38.p13.genome.fa /reference/D3L.fa /reference/pUC57.fa > humanviralplasmid.fa"
cat /reference/GRCh38.p13.genome.fa /reference/D3L.fa /reference/pUC57.fa > humanviralplasmid.fa

echo "Time for (another) coffee break"

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
echo aligning clips of "$tmp"
pbmm2 align --sort --preset ISOSEQ humanviralplasmid.fa "$tmp".positive.start.fastq "$tmp".positive.start.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ humanviralplasmid.fa "$tmp".positive.end.fastq "$tmp".positive.end.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ humanviralplasmid.fa "$tmp".negative.start.fastq "$tmp".negative.start.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ humanviralplasmid.fa "$tmp".negative.end.fastq "$tmp".negative.end.humanviralplasmid.bam
done


module load samtools
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
echo generating bam files for "$tmp"
samtools view -h "$tmp".positive.start.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > "$tmp".positive.start.humanviralplasmid.filt100.bam
samtools view -h "$tmp".positive.end.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > "$tmp".positive.end.humanviralplasmid.filt100.bam
samtools view -h "$tmp".negative.start.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > "$tmp".negative.start.humanviralplasmid.filt100.bam
samtools view -h "$tmp".negative.end.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > "$tmp".negative.end.humanviralplasmid.filt100.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
echo indexing bam files for "$tmp"
samtools index "$tmp".positive.start.humanviralplasmid.filt100.bam
samtools index "$tmp".positive.end.humanviralplasmid.filt100.bam
samtools index "$tmp".negative.start.humanviralplasmid.filt100.bam
samtools index "$tmp".negative.end.humanviralplasmid.filt100.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
echo filtering pUC57 (plasmid vector) files for "$tmp"
samtools view -h "$tmp".positive.start.humanviralplasmid.filt100.bam pUC57 > "$tmp".positive.start.humanviralplasmid.filt100.plasmid.bam
samtools view -h "$tmp".positive.end.humanviralplasmid.filt100.bam pUC57 > "$tmp".positive.end.humanviralplasmid.filt100.plasmid.bam
samtools view -h "$tmp".negative.start.humanviralplasmid.filt100.bam pUC57 > "$tmp".negative.start.humanviralplasmid.filt100.plasmid.bam
samtools view -h "$tmp".negative.end.humanviralplasmid.filt100.bam pUC57 > "$tmp".negative.end.humanviralplasmid.filt100.plasmid.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
echo filtering bq scores < 30 for "$tmp"
samtools view -bq 30 "$tmp".positive.start.humanviralplasmid.filt100.plasmid.bam > "$tmp".positive.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 "$tmp".positive.end.humanviralplasmid.filt100.plasmid.bam > "$tmp".positive.end.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 "$tmp".negative.start.humanviralplasmid.filt100.plasmid.bam > "$tmp".negative.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 "$tmp".negative.end.humanviralplasmid.filt100.plasmid.bam > "$tmp".negative.end.humanviralplasmid.filt100.plasmid30.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
echo converting bam to sam format for "$tmp"
samtools view -o "$tmp".positive.start.humanviralplasmid.filt100.plasmid30.sam "$tmp".positive.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -o "$tmp".positive.end.humanviralplasmid.filt100.plasmid30.sam "$tmp".positive.end.humanviralplasmid.filt100.plasmid30.bam
samtools view -o "$tmp".negative.start.humanviralplasmid.filt100.plasmid30.sam "$tmp".negative.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -o "$tmp".negative.end.humanviralplasmid.filt100.plasmid30.sam "$tmp".negative.end.humanviralplasmid.filt100.plasmid30.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
echo counting sam file lines for "$tmp"
wc -l "$tmp".positive.start.humanviralplasmid.filt100.plasmid30.sam
wc -l "$tmp".positive.end.humanviralplasmid.filt100.plasmid30.sam 
wc -l "$tmp".negative.start.humanviralplasmid.filt100.plasmid30.sam
wc -l "$tmp".negative.end.humanviralplasmid.filt100.plasmid30.sam
done
