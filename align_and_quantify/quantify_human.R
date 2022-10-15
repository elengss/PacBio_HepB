
##### map human reads to human transcript gtf file
##### produces a list of mapped transcripts for each sample, and also outputs files with suffix featureCounts which describe which reads map to each transcript
listmycounts<-list()
c("1002","1003","1004","1009","1010","1012","1013","1015")->all
library(Rsubread)
for(i in 1:length(all)){
listmycounts[[i]]<-featureCounts(paste(all[i],"_prec_iso_fq_human.bam",sep=""),annot.ext="/reference/gencode.v40.chr_patch_hapl_scaff.annotation.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=FALSE,allowMultiOverlap=F,GTF.featureType=c("transcript"),GTF.attrType=c("gene_id"),isLongRead=T,useMetaFeatures=T,minMQS=30,countMultiMappingReads=F,reportReads="CORE")
print(i)
}
