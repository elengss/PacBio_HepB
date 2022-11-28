##### 
##### perform differential analysis on human and viral reads
##### 
##### Other comparisons just require a minor re-write of the code here
##### All functions used are in the edgeR and Rsubread packages.
#####
library(edgeR)
library(Rsubread)
####
##### map human reads to human transcript gtf file using files generated in "align-script.sh"
##### produces a list of mapped transcripts for each sample, and also outputs files with suffix featureCounts which describe which reads map to each transcript
listmycounts<-list()
c("1002","1003","1004","1009","1010","1012","1013","1015")->all
for(i in 1:length(all)){
  listmycounts[[i]]<-featureCounts(paste("data/",all[i],"_prec_iso_fq_human.bam",sep=""),annot.ext="reference/gencode.v42.chr_patch_hapl_scaff.annotation.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=FALSE,allowMultiOverlap=F,GTF.featureType=c("transcript"),GTF.attrType=c("gene_id"),isLongRead=T,useMetaFeatures=T,minMQS=30,countMultiMappingReads=F,reportReads="CORE")
  print(i)
}
matrix(nrow=nrow(listmycounts[[1]][[1]]),ncol=8)->counts
for(i in 1:8){
  listmycounts[[i]][[1]]->counts[,i]}
rownames(listmycounts[[i]][[1]])->rownames(counts)
colnames(counts)<-all

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

#### transfection vs infection

group<-c("T","T","T","T","I","I","I","I")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupI - groupT, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
read.table("reference/ENST_ENSG_GN",header=F)->gn
gn[match(rownames(top.table),gn[,3]),2]->gene
cbind(gene,top.table)->TvI
write.csv(TvI,file="TvI.csv")

############### All infection

counts[,1:4]->counts2
d0 <- DGEList(counts2)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("siSCR","siSCR","siCTCF","siCTCF")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupsiSCR - groupsiCTCF, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
read.table("reference/ENST_ENSG_GN",header=F)->gn
gn[match(rownames(top.table),gn[,3]),2]->gene
cbind(gene,top.table)->TvI
write.csv(TvI,file="siSCRvssiCTCF.csv")

######### All transfection


counts[,5:8]->counts2
d0 <- DGEList(counts2)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("WT","WT","BS","BS")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupWT - groupBS, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
read.table("reference/ENST_ENSG_GN",header=F)->gn
gn[match(rownames(top.table),gn[,3]),2]->gene
cbind(gene,top.table)->TvI
write.csv(TvI,file="WTvsBS.csv")

######################## Viral transcripts

read.csv("can.csv")->can
can[,2:9]->can2
can[,1]->rownames(can2)
read.csv("non-can.csv")->noncan
noncan[,2:9]->noncan2
noncan[,1]->rownames(noncan2)
rbind(can2,noncan2)->counts
counts[,5:8]->counts2
d0 <- DGEList(counts2)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

########### wild type vs BS

group<-c("WT","WT","BS","BS")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupWT - groupBS, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="WTvsBS_viral.csv")

### siRNA vs si CTCF
counts[,1:4]->counts2
d0 <- DGEList(counts2)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("siSCR","siSCR","siCTCF","siCTCF")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupsiSCR - groupsiCTCF, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="siSCRvssiCTCF_viral.csv")

### transfection vs infection
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d0->d
###d <- d0[-drop,] ### drop doesn't work here as none are drop - it corrupts the object

group<-c("T","T","T","T","I","I","I","I")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupI - groupT, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="TvI_viral.csv")

