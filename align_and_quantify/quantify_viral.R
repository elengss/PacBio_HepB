read.table("/reference/D3L-everything+tss_unique.gtf",fill=T)->d3
d3[(nrow(d3)-6):nrow(d3),]->d3n

### this tests whether the intron junctions are within 12 BP of ref - can be modified accordingly
compv<-function(x,y){
return(ifelse(abs(x-y)<12,0,1))}

c("1002","1003","1004","1009","1010","1012","1013","1015")->all
allfiles<-list()
for(i in 1:length(all)){
data.frame(read.table(paste(all[i],".D3L.short",sep="")))->x
colnames(x)<-c("query","mapq","bit","pos","cigar")
x->allfiles[[i]]} 

### this processes the reads with no splice - ie. no 'N' in the cigar
list()->allcant
list()->alltestvec
for(i in 1:length(all)){
allfiles[[i]][which(allfiles[[i]][,"mapq"]>30),]->allfileh ### filters reads below MAPQ of 30
allfileh[grep("N",allfileh[,"cigar"],invert=T),]->allfileh2
allfileh2[,"cigar"]->cigar
allfileh2[,"pos"]->pos
end<-vector()
for(j in 1:length(cigar)){
explodeCigarOpLengths(cigar[j])->x
explodeCigarOps(cigar[j])->y
toMatch <- c("S","I")
matches <- grep(paste(toMatch,collapse="|"), y[[1]], invert=T)
sum(x[[1]][matches])->x2
x2->end[j]}
end+pos->truend

testvec1<-matrix(nrow=length(pos),ncol=nrow(d3n))
testvec2<-matrix(nrow=length(pos),ncol=nrow(d3n))
for(w in 1:length(pos)){
for(k in 1:nrow(d3n)){
compv(d3n$V4[k],pos[w])->testvec1[w,k]}}
vector()->cant
for(v in 1:7){
table(testvec1[,v])[which(names(table(testvec1[,v]))==0)]->cant[v]}
testvec1->alltestvec[[i]]
cant->allcant[[i]]
print(i)
}



#### fastqlines is a file that contains a list of the number of lines in each fastq file (obtained using the shell command “wc -l < filename”).
for i in 1002 1003 1004 1009 1010 1012 1013 1015
do
wc -l $i

read.table("/reference/fastqlines",header=F)->fqlines

#### Each sequence takes up 4 lines in the file, therefore the number of sequences in a file is 1/4 of the total number of lines present :

numseqs <- fqlines/4  

#### (sorry, I still can’t get my head round fqlines/4 -> numseqs).

 

Followed by :

                count=allcant[[i]]/numseqs[i,1])

in the relevant place.

 

Although it’s more efficient as you wrote it, I think a file entry like this :

File         # Lines

1002      1154712

1003      760564


#### fastqover4 is a file that contains a vector of library sizes to be normalised by
read.table("fastqover4",header=F)->fq
pdf("canonical_transcripts_D3L.start_only.norm_by_reads.pdf")
for(i in 1:length(all)){
df <- data.frame(transcript=d3n$V10,
                count=allcant[[i]]/fq[i,1])
print(p<-ggplot(data=df, aes(x=transcript, y=count)) +
  geom_bar(stat="identity")+labs(title=all[i]))}
dev.off()

matrix(nrow=7,ncol=8) -> mat
for(i in 1:length(all)){
allcant[[i]]->mat[,i]}
rownames(mat)<-d3n$V10
colnames(mat)<-all
write.csv(mat,file="can.csv")
#### this writes out a file containing canonical transcripts called can.csv

######################### Non-canonical transcripts

library(GenomicAlignments)
compv<-function(x,y){
return(ifelse(abs(x-y)<12,0,1))}

load("intronnew.Rd")

c("1002","1003","1004","1009","1010","1012","1013","1015")->all
allfiles<-list()
for(i in 1:length(all)){
data.frame(read.table(paste(all[i],".D3L.short",sep="")))->x
colnames(x)<-c("query","mapq","bit","pos","cigar")
x->allfiles[[i]]} 

unique(intronnew)->intronnewu
intronnewu[which(is.na(intronnewu[,1])==F),]->intronnewu
allinst<-list()
testintron2<-list()
for(i in 1:length(all)){
allfiles[[i]][which(allfiles[[i]][,"mapq"]>30),]->allfileh ### filters reads below MAPQ of 30
allfileh[grep("N",allfileh[,"cigar"]),]->allfilesplice
allfilesplice[,"cigar"]->cigar
allfilesplice[,"pos"]->pos

end<-vector()
for(j in 1:length(cigar)){
explodeCigarOpLengths(cigar[j])->x
explodeCigarOps(cigar[j])->y
toMatch <- c("S","I")
matches <- grep(paste(toMatch,collapse="|"), y[[1]], invert=T)
sum(x[[1]][matches])->x2
x2->end[j]}
end+pos->truend

allintron<-list()
countn<-vector()
library(stringr)
for(w in 1:length(cigar)){
str_count(cigar[w], "N")->countn[w]
explodeCigarOpLengths(cigar[w])->x
explodeCigarOps(cigar[w])->y
toMatch <- c("S","I")
matches <- grep(paste(toMatch,collapse="|"), y[[1]], invert=T)
x[[1]][matches]->xnew
y[[1]][matches]->ynew
list()->intron
grep("N",ynew)->ni
idx<-seq(1,length(xnew))
if(length(ni)>0){
for(k in 1:length(ni)){
as.numeric(as.character(pos[w]))+sum(xnew[which(idx<ni[k])])->ns
ns+xnew[ni[k]]->ne
c(ns,ne)->intron[[k]]}
intron->allintron[[w]]}}

firstins<-vector()
firstine<-vector()
secondins<-vector()
secondine<-vector()
thirdins<-vector()
thirdine<-vector()
for(q in 1:length(allintron)){
allintron[[q]][[1]][1]->firstins[q]
allintron[[q]][[1]][2]->firstine[q]
if(length(allintron[[q]])>1){
allintron[[q]][[2]][1]->secondins[q]
allintron[[q]][[2]][2]->secondine[q]}
if(length(allintron[[q]])>2){
allintron[[q]][[3]][1]->thirdins[q]
allintron[[q]][[3]][2]->thirdine[q]}
}

paste(firstins,firstine,secondins,secondine,thirdins,thirdine,sep="_")->allins
allins->allinst[[i]]

###########

testvecstart<-matrix(nrow=length(allintron),ncol=nrow(intronnewu))
testvecend<-matrix(nrow=length(allintron),ncol=nrow(intronnewu))
testvecin1s<-matrix(nrow=length(allintron),ncol=nrow(intronnewu))
testvecin1e<-matrix(nrow=length(allintron),ncol=nrow(intronnewu))
testvecin2s<-matrix(nrow=length(allintron),ncol=nrow(intronnewu))
testvecin2e<-matrix(nrow=length(allintron),ncol=nrow(intronnewu))
testvecin3s<-matrix(nrow=length(allintron),ncol=nrow(intronnewu))
testvecin3e<-matrix(nrow=length(allintron),ncol=nrow(intronnewu))
for(tmp1 in 1:length(allintron)){
for(tmp2 in 1:nrow(intronnewu)){
compv(intronnewu[,"start"][tmp2],pos[tmp1])->testvecstart[tmp1,tmp2]
compv(intronnewu[,"end"][tmp2],truend[tmp1])->testvecend[tmp1,tmp2]
compv(intronnewu[,"ins1"][tmp2],firstins[tmp1])->testvecin1s[tmp1,tmp2]
compv(intronnewu[,"ine1"][tmp2],firstine[tmp1])->testvecin1e[tmp1,tmp2]
#### only tests if there are 2 or more introns in the reference gtf 
if(is.na(intronnewu[,"ins2"][tmp2])==F){
compv(intronnewu[,"ins2"][tmp2],secondins[tmp1])->testvecin2s[tmp1,tmp2]
compv(intronnewu[,"ine2"][tmp2],secondine[tmp1])->testvecin2e[tmp1,tmp2]}
if(is.na(intronnewu[,"ins3"][tmp2])==F){
compv(intronnewu[,"ins3"][tmp2],thirdins[tmp1])->testvecin3s[tmp1,tmp2]
compv(intronnewu[,"ine3"][tmp2],thirdine[tmp1])->testvecin3e[tmp1,tmp2]}
}}

testvecin1s[which(countn==1),]->testvecin1s1
testvecin1e[which(countn==1),]->testvecin1e1

testvecin1s[which(countn==2),]->testvecin1s2
testvecin1e[which(countn==2),]->testvecin1e2
testvecin2s[which(countn==2),]->testvecin2s2
testvecin2e[which(countn==2),]->testvecin2e2

testvecin1s[which(countn==3),]->testvecin1s3
testvecin1e[which(countn==3),]->testvecin1e3
testvecin2s[which(countn==3),]->testvecin2s3
testvecin2e[which(countn==3),]->testvecin2e3
testvecin3s[which(countn==3),]->testvecin3s3
testvecin3e[which(countn==3),]->testvecin3e3

### check number of junctions

vector()->countj
for(t in 1:nrow(intronnewu)){
length(which(is.na(intronnewu[t,])==T))->countj[t]}
intronnewu[which(countj==8),]->in1
intronnewu[which(countj==6),]->in2
intronnewu[which(countj==4),]->in3


testvecin1s1[,which(countj==8)]->testvecin1s1s
testvecin1e1[,which(countj==8)]->testvecin1e1s

testvecin1s2[,which(countj==6)]->testvecin1s2s
testvecin1e2[,which(countj==6)]->testvecin1e2s
testvecin2s2[,which(countj==6)]->testvecin2s2s
testvecin2e2[,which(countj==6)]->testvecin2e2s

testvecin1s3[,which(countj==4)]->testvecin1s3s
testvecin1e3[,which(countj==4)]->testvecin1e3s
testvecin2s3[,which(countj==4)]->testvecin2s3s
testvecin2e3[,which(countj==4)]->testvecin2e3s
testvecin3s3[,which(countj==4)]->testvecin3s3s
testvecin3e3[,which(countj==4)]->testvecin3e3s



matrix(ncol=nrow(in1),nrow=nrow(testvecin1s1))->alltestvecstrict1
for(g in 1:nrow(in1)){
testvecin1s1s[,g]+testvecin1e1s[,g]->alltestvecstrict1[,g]}

matrix(ncol=nrow(in2),nrow=nrow(testvecin1s2))->alltestvecstrict2
for(g in 1:nrow(in2)){
testvecin1s2s[,g]+testvecin1e2s[,g]+testvecin2s2s[,g]+testvecin2e2s[,g]->alltestvecstrict2[,g]}

matrix(ncol=nrow(in3),nrow=nrow(testvecin1s3))->alltestvecstrict3
for(g in 1:nrow(in3)){
testvecin1s3s[,g]+testvecin1e3s[,g]+testvecin2s3s[,g]+testvecin2e3s[,g]+testvecin3s3s[,g]+testvecin3e3s[,g]->alltestvecstrict3[,g]
}

countin1<-vector()
countin2<-vector()
countin3<-vector()
for(r in 1:nrow(in1)){
length(which(alltestvecstrict1[,r]==0))->countin1[r]}
for(r in 1:nrow(in2)){
length(which(alltestvecstrict2[,r]==0))->countin2[r]}
for(r in 1:nrow(in3)){
length(which(alltestvecstrict3[,r]==0))->countin3[r]}
names(countin1)<-rownames(in1)
names(countin2)<-rownames(in2)
names(countin3)<-rownames(in3)

c(countin1,countin2,countin3)->countallin
countallin->testintron2[[i]]
print(i)
}

save(allinst,file="allinst.Rd")
save(testintron2,file="testintron2.Rd")

read.table("fastqover4",header=F) -> fq

matstrict<-matrix(ncol=8,nrow=length(testintron2[[1]]))
for(i in 1:ncol(matstrict)){
testintron2[[i]]->matstrict[,i]}
rownames(matstrict)<-names(testintron2[[1]])
colnames(matstrict)<-all
matstrict[which(rowSums(matstrict)!=0),]->matstrict2
#### this writes out a file containing non-canonical transcripts called non-can.csv
write.csv(matstrict2,file="non-can.csv")


##### this writes out intron junctions separately for forward and reverse strands


load("allinst.Rd")
allinst->allinslist

c("1002","1003","1004","1009","1010","1012","1013","1015")->all
allfiles<-list()
for(i in 1:length(all)){
data.frame(read.table(paste(all[i],".D3L.short",sep="")))->x
colnames(x)<-c("query","mapq","bit","pos","cigar")
x->allfiles[[i]]} 

vector()->allsplice
for(i in 1:length(allinslist)){
c(allsplice,names(table(allinslist[[i]])))->allsplice}

unique(allsplice)->allsplicenames
matrix(nrow=length(allsplicenames),ncol=8) ->bigmats
rownames(bigmats)<-allsplicenames
colnames(bigmats)<-all
for(i in 1:8){
allfiles[[i]][which(allfiles[[i]][,"mapq"]>30),]->allfileh ### filters reads below MAPQ of 30
allfileh[grep("N",allfileh[,"cigar"]),]->allfilesplice
allinslist[[i]]->temp
temp[which(allfilesplice$bit==0)]->temp2
table(temp2)->tmp
tmp[match(rownames(bigmats),names(tmp))]->bigmats[,i]}
rowSums(bigmats,na.rm=T)->bmr
bigmats[order(bmr,decreasing=T),]->bigmatso
bmr[order(bmr,decreasing=T)]->bmro1
bigmatso[which(bmro1>50),]->bigmats50forw
bigmatso[which(bmro1>10),]->bigmats10forw
bigmatso[which(bmro1>1),]->bigmats1forw
write.table(bigmats1forw,file="splice_forw.txt")


unique(allsplice)->allsplicenames
matrix(nrow=length(allsplicenames),ncol=8) ->bigmatn
rownames(bigmatn)<-allsplicenames
colnames(bigmatn)<-all
for(i in 1:8){
allfiles[[i]][which(allfiles[[i]][,"mapq"]>30),]->allfileh ### filters reads below MAPQ of 30
allfileh[grep("N",allfileh[,"cigar"]),]->allfilesplice
allinslist[[i]]->temp
temp[which(allfilesplice$bit==16)]->temp2
table(temp2)->tmp
tmp[match(rownames(bigmatn),names(tmp))]->bigmatn[,i]}
rowSums(bigmatn,na.rm=T)->bmr
bmr[order(bmr,decreasing=T)]->bmro2
bigmatn[order(bmr,decreasing=T),]->bigmatno
bigmatno[which(bmro2>50),]->bigmatn50reve
bigmatno[which(bmro2>10),]->bigmatn10reve
bigmatno[which(bmro2>1),]->bigmatn1reve
write.table(bigmatn1reve,file="splice_reve.txt")

