######## reads in gtf file and outputs matrix of start/end coordinates and intron junctions - saved in an R object called intronnew.Rd
read.table("/data/D3L-everything+tss_unique.gtf",fill=T)->d3
d3[(nrow(d3)-6):nrow(d3),]->d3n ### this gtf file is specific to HepB - the last 6 lines are removed as they as canonical transcripts with no introns
d3[1:(nrow(d3)-7),]->d3
d3[which(d3$V3=="exon"),]->exonname
exonname$V13->allname
unique(allname)->an
intronjunc<-matrix(nrow=length(an),ncol=10) ### this assumes maximum of 6 exon (5 introns)
transtart<-vector()
tranend<-vector()
for(i in 1:length(an)){
exonname[which(exonname$V13==an[i]),]->tmp
if(nrow(tmp)>1)
for(k in 1:(nrow(tmp)-1)){
tmp[k,"V5"]->intronjunc[i,(k-1)*2+1]
tmp[k+1,"V4"]->intronjunc[i,(k-1)*2+2]
tmp[1,"V4"]->transtart[i]
tmp[nrow(tmp),"V5"]->tranend[i]
}}
cbind(transtart,tranend,intronjunc)->intronnew
rownames(intronnew)<-an
colnames(intronnew)<-c("start","end","ins1","ine1","ins2","ine2","ins3","ine3","ins4","ine4","ins5","ine5")
save(intronnew,file="intronnew.Rd")
