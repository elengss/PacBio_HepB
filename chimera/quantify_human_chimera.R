### Identify chimeric reads
###
### load required library
library(Rsubread)
c("1002","1003","1010","1012")->all
list()->startpos
list()->endpos
list()->startneg
list()->endneg
for(i in 1:length(all)){
  paste("data/",all[i],".positive.start.humanviralplasmid.filt100.bam",sep="")->x
  mycounts<-featureCounts(x,annot.ext="reference/gencode.v42.chr_patch_hapl_scaff.annotation.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=FALSE,allowMultiOverlap=F,GTF.featureType=c("transcript"),GTF.attrType=c("gene_id"),primaryOnly=T,minMQS=50,isLongRead=T,useMetaFeatures=T,countMultiMappingReads=F,reportReads="CORE")
  mycounts->startpos[[i]]
  paste("data/",all[i],".positive.end.humanviralplasmid.filt100.bam",sep="")->x
  mycounts<-featureCounts(x,annot.ext="reference/gencode.v42.chr_patch_hapl_scaff.annotation.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=FALSE,allowMultiOverlap=F,GTF.featureType=c("transcript"),GTF.attrType=c("gene_id"),primaryOnly=T,minMQS=50,isLongRead=T,useMetaFeatures=T,countMultiMappingReads=F,reportReads="CORE")
  mycounts->endpos[[i]]
  paste("data/",all[i],".negative.start.humanviralplasmid.filt100.bam",sep="")->x
  mycounts<-featureCounts(x,annot.ext="reference/gencode.v42.chr_patch_hapl_scaff.annotation.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=FALSE,allowMultiOverlap=F,GTF.featureType=c("transcript"),GTF.attrType=c("gene_id"),primaryOnly=T,minMQS=50,isLongRead=T,useMetaFeatures=T,countMultiMappingReads=F,reportReads="CORE")
  mycounts->startneg[[i]]
  paste("data/",all[i],".negative.end.humanviralplasmid.filt100.bam",sep="")->x
  mycounts<-featureCounts(x,annot.ext="reference/gencode.v42.chr_patch_hapl_scaff.annotation.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=FALSE,allowMultiOverlap=F,GTF.featureType=c("transcript"),GTF.attrType=c("gene_id"),primaryOnly=T,minMQS=50,isLongRead=T,useMetaFeatures=T,countMultiMappingReads=F,reportReads="CORE")
  mycounts->endneg[[i]]}

save(startpos,file="startpos_incplasmid.Rd")
save(endpos,file="endpos_incplasmid.Rd")
save(startneg,file="startneg_incplasmid.Rd")
save(endneg,file="endneg.Rd_incplasmid")


matrix(nrow=length(startpos[[1]][[1]]),ncol=8)->mat
for(i in 1:8){
matrix(nrow=length(startpos[[1]][[1]]),ncol=4)->mat2
startpos[[i]][[1]]->mat2[,1]
endpos[[i]][[1]]->mat2[,2]
startneg[[i]][[1]]->mat2[,3]
endneg[[i]][[1]]->mat2[,4]
rowSums(mat2)->mat[,i]
}
rownames(startpos[[1]][[1]])->rownames(mat)
all->colnames(mat)
mat[which(rowSums(mat)>0),]->matnot0
read.table("reference/ENST_ENSG_GN",sep = " ",header=F)->gn
gn[match(rownames(matnot0),gn$V3),2]->gene
rownames(matnot0)<-gene
matnot0[order(rowSums(matnot0),decreasing=T),]->matnew
write.table(matnew,file="chimera_inc_plasmid.csv",col.names = NA,row.names = TRUE,sep = ",")
print("Script completed")
