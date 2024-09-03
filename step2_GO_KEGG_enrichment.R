rm(list = ls())
options(stringsAsFactors = F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)


# load DEG data
DEG_edgeR <- read.csv("DEG_edgeR_all_reflid.xls",
                      sep="\t")

DEG <- DEG_edgeR[DEG_edgeR$regulated!="normal",1]

head(DEG)

# extract all the gene
gene_all <- DEG_edgeR$GeneID #20408

head(gene_all)

length(gene_all);length(DEG)

#### get gene entzezid from org.Hs.eg.db
keytypes(org.Hs.eg.db)


allID <- bitr(gene_all, fromType = "ENSEMBL", 
              toType = c( "ENTREZID" ), 
              OrgDb = org.Hs.eg.db ) 
head(allID)


degID <- bitr(DEG, fromType = "ENSEMBL", 
              toType = c( "ENTREZID" ), 
              OrgDb = org.Hs.eg.db ) 
head(degID)


# KEGG enrich analysis----
enrich <- enrichKEGG(gene = degID[,2],
                     organism='hsa',
                     universe=allID[,2],
                     pvalueCutoff=1,
                     qvalueCutoff=1)



# calculate the enrichment value
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) 
  as.numeric(x[1])/as.numeric(x[2])))
head(GeneRatio)

BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) 
  as.numeric(x[1])/as.numeric(x[2])  ))
head(BgRatio)

enrich_factor <- GeneRatio/BgRatio

out <- data.frame(enrich$ID,
                  enrich$Description,
                  enrich$GeneRatio,
                  enrich$BgRatio,
                  round(enrich_factor,2),
                  enrich$pvalue,
                  enrich$qvalue,
                  enrich$geneID,
                  enrich$Count)

colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor",
                   "pvalue","qvalue","geneID","count")

# save kegg enrichment result
write.csv(out,"kegg_enrichment_pathways.csv",
          row.names = F,quote = F)

out_sig0.05 <- out[out$qvalue<0.05,]

# barplot
bar <- barplot(enrich,showCategory=20,title="KEGG Pathway",
               colorBy="p.adjust")
bar

dev.new()

# save kegg barplot 
pdf(file = "kegg_bar_plot.pdf",width = 8,height = 6)
print(bar)
dev.off()



# =============================== GO enrichment analysis=====================
enrich <- enrichGO(gene =degID[,2],OrgDb='org.Hs.eg.db',
                   ont="BP",universe=allID[,2],pvalueCutoff=1,qvalueCutoff=1)

# calculate enrich factor
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) 
  as.numeric(x[1])/as.numeric(x[2])))

BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) 
  as.numeric(x[1])/as.numeric(x[2])))

enrich_factor <- GeneRatio/BgRatio

out <- data.frame(enrich$ID,
                  enrich$Description,
                  enrich$GeneRatio,
                  enrich$BgRatio,
                  round(enrich_factor,2),
                  enrich$pvalue,
                  enrich$qvalue,
                  enrich$geneID)

colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor",
                   "pvalue","qvalue","geneID")

write.csv(out,"enrich_GO-edgr.csv",row.names = F,quote = F)

out_sig0.05 <- out[out$qvalue<0.05,]


# barplot
bar <- barplot(enrich,showCategory=20,title="Biological Pathway",
               colorBy="p.adjust")
bar

# save go barplot
pdf(file = "BP_bar_plot-top20-edgr0731.pdf",width = 6,height = 6)
print(bar)
dev.off()


