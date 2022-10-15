#### the reference file is a concatenation of human, viral and plasmid fasta files
cat humanviral.fa pUC57.fa > humanviralplasmid.fa

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
pbmm2 align --sort --preset ISOSEQ /well/ansari/users/yem086/hepb/new/humanviralplasmid.fa "$tmp".positive.start.fastq "$tmp".positive.start.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ /well/ansari/users/yem086/hepb/new/humanviralplasmid.fa "$tmp".positive.end.fastq "$tmp".positive.end.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ /well/ansari/users/yem086/hepb/new/humanviralplasmid.fa "$tmp".negative.start.fastq "$tmp".negative.start.humanviralplasmid.bam
pbmm2 align --sort --preset ISOSEQ /well/ansari/users/yem086/hepb/new/humanviralplasmid.fa "$tmp".negative.end.fastq "$tmp".negative.end.humanviralplasmid.bam
done


module load samtools
for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
samtools view -h "$tmp".positive.start.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > "$tmp".positive.start.humanviralplasmid.filt100.bam
samtools view -h "$tmp".positive.end.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > "$tmp".positive.end.humanviralplasmid.filt100.bam
samtools view -h "$tmp".negative.start.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > "$tmp".negative.start.humanviralplasmid.filt100.bam
samtools view -h "$tmp".negative.end.humanviralplasmid.bam | awk 'length($10) > 100 || $1 ~ /^@/' | samtools view -bS - > "$tmp".negative.end.humanviralplasmid.filt100.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
samtools index "$tmp".positive.start.humanviralplasmid.filt100.bam
samtools index "$tmp".positive.end.humanviralplasmid.filt100.bam
samtools index "$tmp".negative.start.humanviralplasmid.filt100.bam
samtools index "$tmp".negative.end.humanviralplasmid.filt100.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
samtools view -h "$tmp".positive.start.humanviralplasmid.filt100.bam pUC57 > "$tmp".positive.start.humanviralplasmid.filt100.plasmid.bam
samtools view -h "$tmp".positive.end.humanviralplasmid.filt100.bam pUC57 > "$tmp".positive.end.humanviralplasmid.filt100.plasmid.bam
samtools view -h "$tmp".negative.start.humanviralplasmid.filt100.bam pUC57 > "$tmp".negative.start.humanviralplasmid.filt100.plasmid.bam
samtools view -h "$tmp".negative.end.humanviralplasmid.filt100.bam pUC57 > "$tmp".negative.end.humanviralplasmid.filt100.plasmid.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
samtools view -bq 30 "$tmp".positive.start.humanviralplasmid.filt100.plasmid.bam > "$tmp".positive.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 "$tmp".positive.end.humanviralplasmid.filt100.plasmid.bam > "$tmp".positive.end.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 "$tmp".negative.start.humanviralplasmid.filt100.plasmid.bam > "$tmp".negative.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -bq 30 "$tmp".negative.end.humanviralplasmid.filt100.plasmid.bam > "$tmp".negative.end.humanviralplasmid.filt100.plasmid30.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
samtools view -o "$tmp".positive.start.humanviralplasmid.filt100.plasmid30.sam "$tmp".positive.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -o "$tmp".positive.end.humanviralplasmid.filt100.plasmid30.sam "$tmp".positive.end.humanviralplasmid.filt100.plasmid30.bam
samtools view -o "$tmp".negative.start.humanviralplasmid.filt100.plasmid30.sam "$tmp".negative.start.humanviralplasmid.filt100.plasmid30.bam
samtools view -o "$tmp".negative.end.humanviralplasmid.filt100.plasmid30.sam "$tmp".negative.end.humanviralplasmid.filt100.plasmid30.bam
done

for tmp in 1002 1003 1004 1009 1010 1012 1013 1015
do
wc -l "$tmp".positive.start.humanviralplasmid.filt100.plasmid30.sam
wc -l "$tmp".positive.end.humanviralplasmid.filt100.plasmid30.sam 
wc -l "$tmp".negative.start.humanviralplasmid.filt100.plasmid30.sam
wc -l "$tmp".negative.end.humanviralplasmid.filt100.plasmid30.sam
done
