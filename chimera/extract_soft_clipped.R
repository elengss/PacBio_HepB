
extractsoftclip<-function(inbam,outfq){
aln <- scanBam(inbam)
aln[[1]]$cigar -> cigar
cigarRangesAlongQuerySpace(cigar, ops="S") -> cigarir

#### separate the pos and negative strand and 5' 3' end
aln[[1]]$cigar[which(aln[[1]]$strand=="+")]->positivecigar
aln[[1]]$qname[which(aln[[1]]$strand=="+")]->positiveqname
names(positivecigar)<-positiveqname
cigarRangesAlongQuerySpace(positivecigar, ops="S",drop.empty.ranges=T) -> cigarir
startclip<-list()
endclip<-list()
for(i in 1:length(cigarir)){
cigarir[[i]][which(end(cigarir[[i]])==aln[[1]]$qwidth[i])]->endclip[[i]]
cigarir[[i]][which(start(cigarir[[i]])==1)]->startclip[[i]]
print(i)
}

##### do positive start
aln[[1]]$seq[which(aln[[1]]$strand=="+")]->positiveseq
mylist<-list()
for(i in 1:length(startclip)){
extractAt(positiveseq[i],startclip[[i]])->mylist[[i]]
print(i)
}

test2<-vector()
for(i in 1:length(mylist)){
length(mylist[[i]][[1]])->test2[i]
print(i)
}

aln[[1]]$qname[which(aln[[1]]$strand=="+")]->positiveqname
positiveqname[which(test2>0)]-> positiveqname2
mylist[which(test2>0)]->mylist2
LargeDNAStringSet <- DNAStringSet(sapply(sapply(mylist2, `[[`, 1),"[[",1))
positiveqname2->names(LargeDNAStringSet)
writeXStringSet(LargeDNAStringSet, paste(outfq,".positive.start.fastq",sep=""), append=FALSE,compress=FALSE, compression_level=NA, format="fastq")

### do positive end
aln[[1]]$seq[which(aln[[1]]$strand=="+")]->positiveseq
mylist<-list()
for(i in 1:length(endclip)){
extractAt(positiveseq[i],endclip[[i]])->mylist[[i]]
print(i)
}

test2<-vector()
for(i in 1:length(mylist)){
length(mylist[[i]][[1]])->test2[i]
print(i)
}
aln[[1]]$qname[which(aln[[1]]$strand=="+")]->positiveqname
positiveqname[which(test2>0)]-> positiveqname2
mylist[which(test2>0)]->mylist2
LargeDNAStringSet <- DNAStringSet(sapply(sapply(mylist2, `[[`, 1),"[[",1))
positiveqname2->names(LargeDNAStringSet)
writeXStringSet(LargeDNAStringSet, paste(outfq,".positive.end.fastq",sep=""), append=FALSE,compress=FALSE, compression_level=NA, format="fastq")

#### do negative start

aln[[1]]$cigar[which(aln[[1]]$strand=="-")]->negativecigar
aln[[1]]$qname[which(aln[[1]]$strand=="-")]->negativeqname
names(negativecigar)<-negativeqname
cigarRangesAlongQuerySpace(negativecigar, ops="S",drop.empty.ranges=T) -> cigarir
startclip<-list()
endclip<-list()
for(i in 1:length(cigarir)){
cigarir[[i]][which(end(cigarir[[i]])==aln[[1]]$qwidth[i])]->endclip[[i]]
cigarir[[i]][which(start(cigarir[[i]])==1)]->startclip[[i]]
print(i)
}
aln[[1]]$seq[which(aln[[1]]$strand=="-")]->negativeseq
mylist<-list()
for(i in 1:length(startclip)){
extractAt(negativeseq[i],startclip[[i]])->mylist[[i]]
print(i)
}
test2<-vector()
for(i in 1:length(mylist)){
length(mylist[[i]][[1]])->test2[i]
print(i)
}
aln[[1]]$qname[which(aln[[1]]$strand=="-")]->negativeqname
negativeqname[which(test2>0)]-> negativeqname2
mylist[which(test2>0)]->mylist2
LargeDNAStringSet <- DNAStringSet(sapply(sapply(mylist2, `[[`, 1),"[[",1))
negativeqname2->names(LargeDNAStringSet)
writeXStringSet(LargeDNAStringSet, paste(outfq,".negative.start.fastq",sep=""), append=FALSE,compress=FALSE, compression_level=NA, format="fastq")


### do negative end
aln[[1]]$seq[which(aln[[1]]$strand=="-")]->negativeseq
mylist<-list()
for(i in 1:length(endclip)){
extractAt(negativeseq[i],endclip[[i]])->mylist[[i]]
print(i)
}
test2<-vector()
for(i in 1:length(mylist)){
length(mylist[[i]][[1]])->test2[i]
print(i)
}
aln[[1]]$qname[which(aln[[1]]$strand=="-")]->negativeqname
negativeqname[which(test2>0)]-> negativeqname2
mylist[which(test2>0)]->mylist2
LargeDNAStringSet <- DNAStringSet(sapply(sapply(mylist2, `[[`, 1),"[[",1))
negativeqname2->names(LargeDNAStringSet)
writeXStringSet(LargeDNAStringSet, paste(outfq,".negative.end.fastq",sep=""), append=FALSE,compress=FALSE, compression_level=NA, format="fastq")}
