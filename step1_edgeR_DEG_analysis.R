rm(list = ls())
options(stringsAsFactors = F)

# load package
library(edgeR)
library(ggplot2)

# load filter reads count and group data
exprSet <- read.csv("rawcounts_filtered_20408.csv", row.names = 1)
group_data <- read.csv("group_data.csv")


exprSet[1:4,1:4]
group_list <- group_data$infection
table(group_list)

# give group_list factor levels
group_list <- factor(group_list,levels = c("Noninfection","Infection"))
head(group_list)

# establish linear model
design <- model.matrix(~0+factor(group_list))
rownames(design) <- colnames(exprSet)
colnames(design) <- levels(factor(group_list))
head(design)

# establish edgeR DGEList object
DEG <- DGEList(counts=exprSet, 
               group=factor(group_list))

DEG <- calcNormFactors(DEG)

DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

fit <- glmFit(DEG, design)

lrt <- glmLRT(fit, contrast=c(-1,1)) 

DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

# idetifiy DEGs
fc_cutoff <- 2 
pvalue <- 0.05

DEG_edgeR$regulated <- "normal"

loc_up <- intersect(which(DEG_edgeR$logFC>log2(fc_cutoff)),
                    which(DEG_edgeR$PValue<pvalue))
loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))


DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"

table(DEG_edgeR$regulated) 

DEG_edgeR$GeneID <- rownames(DEG_edgeR)


# identify the symbol of DEGs。ref_id are the ensemble id and symbol id table 
# from gif reference annotation files
ref_id <- read.table("ensemble_and_symbol_id.txt", sep="\t",header = F, 
                      col.names = c("ens","coding","SYMBOL"))

DEG_edgeR_symbol <- merge(DEG_edgeR,ref_id[,c("ens","SYMBOL")], by.x = "GeneID",
                    by.y = "ens",all.x = TRUE)



# or using clusterprofiele
library(clusterProfiler)
id2symbol <- bitr(rownames(DEG_edgeR), 
                  fromType = "ENSEMBL", 
                  toType = "SYMBOL", 
                  OrgDb = org.Hs.eg.db)
head(id2symbol)
DEG_edgeR_symbol <- merge(id2symbol,DEG_edgeR,
                          by.x="ENSEMBL",by.y="GeneID",all.y=T)
head(DEG_edgeR_symbol)
dim(DEG_edgeR_symbol)

#-------------------------------------------------------------------------------
library(tidyverse)
DEG_edgeR_symbol_Sig <- filter(DEG_edgeR_symbol,regulated!="normal")

table(DEG_edgeR_symbol_Sig$regulated) # 共776行,76,669

