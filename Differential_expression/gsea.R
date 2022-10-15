
read.csv("TvI.csv")->tvi
tvi[which(tvi$adj.P.Val<0.0000001),]->tvi
log(tvi$B)->original_gene_list
tvi$gene->names(original_gene_list)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)


gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

require(DOSE)
pdf("gut0_dotplot.pdf")
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()


gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

require(DOSE)
pdf("gut0_dotplot.pdf")
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()
