rm(list = ls())
options(stringsAsFactors = F)
set.seed(123)

# 加载包
library(edgeR)
library(ggplot2)
library(vegan)
library(FactoMineR)
library(factoextra)
library(tidyverse)

# ===============Read gene expression matrix information and metadata ===========
lname <- load(file = "Step03_infection_array.Rdata")
lname
dat_infect[1:5,1:5]


#======================step 1 PCA analysis======================================
# bacter vs nonbacteria 
group_bac <- groupdata_ad %>% 
  filter(!bacteria %in% c("Bacteria-coinfection", "Noninfection"))
dat_bac <- dat_infect[group_bac$sample_id,]
table(group_bac$bacteria)

# draw pca figure
dat_bac_infect <- as.data.frame(dat_bac)
dat_bac_pca <- PCA(dat_bac_infect, graph = FALSE)
class(dat_bac_pca)
dim(dat_bac_infect)

# statistic analysis
permanova_result <- adonis(dat_bac_infect ~ group_bac$bacteria, method = "bray")
print(permanova_result$aov.tab) 


# get P value 和 R2
p_value <- permanova_result$aov.tab$`Pr(>F)`[1]
r_squared <- permanova_result$aov.tab$R2[1]

# draw Pca figure with pvalue
p2 <- fviz_pca_ind(dat_bac_pca,
                   geom.ind = "point", # 只显示点，不显示文字
                   col.ind = group_bac$bacteria, # 用不同颜色表示分组
                   # palette = "rainbow",
                   palette = c("#00AFBB", "#E7B800"),
                   addEllipses = T, # 是否圈起来
                   legend.title = "Groups") + 
  theme_bw() +
  annotate("text", x = Inf, y = Inf, 
           label = paste("P-value:", round(p_value, 3), "\nR²:", round(r_squared, 3)),
           hjust = 1.1, vjust = 2, size = 5, color = "black")

print(p2)

#ggsave("advs_revision_result/bacteria_vs_nonbacteria_PCA.pdf")



# ============================step 2 edgr abnalysis============================
# load read count value
filter_count <- read.csv("flter_count.txt",sep="\t")
dim(filter_count) 
filter_count[1:5,1:5]

# remove non-bacteria infection
group_bacno <- groupdata_ad %>% 
  filter(!bacteria %in% c("Bacteria-coinfection", "Nonbacteria"))

group_list <- group_bacno$infection
table(group_list) 

group_list <- factor(group_list,levels = c("Noninfection","Infection"))
head(group_list)

design <- model.matrix(~0+factor(group_list))
rownames(design) <- colnames(exprSet)
colnames(design) <- levels(factor(group_list))
head(design)

#  make edgeR DGEList
DEG <- DGEList(counts=exprSet, 
               group=factor(group_list))

DEG <- calcNormFactors(DEG)

DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

# fit the model
fit <- glmFit(DEG, design)
lrt <- glmLRT(fit, contrast=c(-1,1)) 

DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

DEG_edgeR_symbol <- merge(id2symbol,DEG_edgeR,
                          by.x="ENSEMBL",by.y="GeneID",all.y=T)
DEG_edgeR_symbol_Sig <- filter(DEG_edgeR_symbol,regulated!="normal")
