setwd("C:/Users/Yanara Bernal/Desktop")

BiocManager::install("clusterProfiler", force = TRUE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("DOSE", force = TRUE)
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(vctrs)
library(dplyr)



BiocManager::install("DBI", force = TRUE)


BiocManager::install("AnnotationDbi", force = TRUE)





BiocManager::install("org.Hs.eg.db", character.only = TRUE, force = TRUE)
library(org.Hs.eg.db)
install.packages("rlang")
library(rlang)
install.packages("cli")
library(cli)
library(DOSE)
install.packages("ggridges")
library(ggridges)

# reading in data from deseq2
df <- read.csv("C:/Users/Yanara Bernal/Desktop/Finales DES/parpi_vep_rmsk_cosmic_test_cosmic_52818_fc_fdr.tsv", sep="")
df <- read.csv("C:/Users/Yanara Bernal/Desktop/Finales DES/antra_vep_rmsk_cosmic_test_cosmic_54567_fc_fdr.tsv", sep="")
df <- read.csv("C:/Users/Yanara Bernal/Desktop/Finales DES/alky_vep_rmsk_cosmic_test_cosmic_56425_fc_fdr.tsv", sep="")

df <- subset(df, differential != "non_differential") 

df <- subset(df, SYMBOL != "-")
# we want the log2 fold change 
original_gene_list <- df$fold_change_log

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 1000, 
             minGSSize = 3, 
             maxGSSize = 100, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")


########
write.table(gse, "gse_fclog_parpi.tsv", row.names = FALSE)

dotplot_parpi <- dotplot(gse)
dotplot_parpi2 <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ridgeplot_parpi3 <- ridgeplot(gse) + labs(x = "enrichment distribution")
dotplot_parpi4 <- dotplot(gse, showCategory=20)




########
write.table(gse, "gse_fclog_antra.tsv", row.names = FALSE)

dotplot_antra <- dotplot(gse)
dotplot_antra2 <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ridgeplot_antra3 <- ridgeplot(gse) + labs(x = "enrichment distribution")
dotplot_antra4 <- dotplot(gse, showCategory=20)
########
write.table(gse, "gse_fclog_alky.tsv", row.names = FALSE)

dotplot_alky <- dotplot(gse)
dotplot_alky2 <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ridgeplot_alky3 <- ridgeplot(gse) + labs(x = "enrichment distribution")
dotplot_alky4 <- dotplot(gse, showCategory=20)



library(patchwork)

(dotplot_parpi | dotplot_antra |dotplot_alky)
(dotplot_parpi2 | dotplot_antra2 |dotplot_alky2)
(ridgeplot_parpi3 | ridgeplot_antra3 |ridgeplot_alky3)
(dotplot_parpi4 | dotplot_antra4 |dotplot_alky4)


