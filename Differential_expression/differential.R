##### 
##### perform differential analysis on human and viral reads
##### 
##### Other comparisons just require a minor re-write of the code here
##### All functions used are in the edgeR package.
#####
library(edgeR)
####
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
read.table("/well/ansari/users/yem086/ENST_ENSG_GN",header=F)->gn
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
read.table("/well/ansari/users/yem086/ENST_ENSG_GN",header=F)->gn
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
read.table("/well/ansari/users/yem086/ENST_ENSG_GN",header=F)->gn
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
read.table("/well/ansari/users/yem086/ENST_ENSG_GN",header=F)->gn
gn[match(rownames(top.table),gn[,3]),2]->gene
cbind(gene,top.table)->TvI
write.csv(TvI,file="WTvsBS_viral.csv")

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

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupsiSCR - groupsiCTCF, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
read.table("/well/ansari/users/yem086/ENST_ENSG_GN",header=F)->gn
gn[match(rownames(top.table),gn[,3]),2]->gene
cbind(gene,top.table)->TvI
write.csv(TvI,file="siSCRvssiCTCF_viral.csv")

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

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupI - groupT, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
read.table("/well/ansari/users/yem086/ENST_ENSG_GN",header=F)->gn
gn[match(rownames(top.table),gn[,3]),2]->gene
cbind(gene,top.table)->TvI
write.csv(TvI,file="TvI_viral.csv")
