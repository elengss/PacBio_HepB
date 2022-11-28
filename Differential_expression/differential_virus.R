##### load viral read counts and perform simple DGE analysis
##### comparing viral read frequencies between all pairs of groups
##### Other comparisons just require a minor re-write of the code here
##### All functions used are in the edgeR package.
#####
library(edgeR)
####
### Just a reminder of the filename barcodes (not that we use them here)
### "1002","1003","1004","1009","1010","1012","1013","1015"
###
### These 8 data files are structured in 4 replicate pairs :
### "A","A","B","B","C","C","D","D"

######################## import Viral transcript numbers
####
#### Read in viral counts previously saved by the "quantify_viral.R" script
read.csv("can.csv")->can
read.csv("non-can.csv")->noncan
can[,2:9]->can2
can[,1]->rownames(can2)
noncan[,2:9]->noncan2
noncan[,1]->rownames(noncan2)
rbind(can2,noncan2)->counts
###
### The following could be rolled up into a loop function, but this layout is 
### easier to edit and re-purpose.
###
### A vs B
counts1 <- counts[,1:4]
d0 <- DGEList(counts1)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("A","A","B","B")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupA - groupB, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="AvsB_viral.csv")

### A vs C
counts2 <- cbind(counts[,1:2],counts[,5:6])
d0 <- DGEList(counts2)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("A","A","C","C")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupA - groupC, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="AvsC_viral.csv")

### A vs D
counts3 <- cbind(counts[,1:2],counts[,7:8])
d0 <- DGEList(counts3)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("A","A","D","D")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupA - groupD, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="AvsD_viral.csv")

### B vs C
counts4 <- counts[,3:6]
d0 <- DGEList(counts4)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("B","B","C","C")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupB - groupC, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="BvsC_viral.csv")

### A vs D
counts5 <- cbind(counts[,3:4],counts[,7:8])
d0 <- DGEList(counts5)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("B","B","D","D")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupB - groupD, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="BvsD_viral.csv")

### C vs D
counts6 <- counts[,5:8]
d0 <- DGEList(counts6)
d0 <- calcNormFactors(d0)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group<-c("C","C","D","D")
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupC - groupD, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
write.csv(top.table,file="CvsD_viral.csv")
####
print("Finished")
