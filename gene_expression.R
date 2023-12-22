library(dplyr)

setwd("/Users/macbook/Desktop/PAPER1")
parpi <- read.csv("~/Desktop/PAPER1/parpi/output/data_parpi_all_52818.tsv", sep="")
antra <- read.csv("~/Desktop/PAPER1/antra/output/data_antra_all_54567.tsv", sep="")
alky <- read.csv("~/Desktop/PAPER1/alky/output/data_alky_all_56425.tsv", sep="")


CCLE_expression <- read.csv("~/Desktop/PAPER1/gene_expression/CCLE_expression.csv")
sample_info <- read.csv("~/Desktop/PAPER1/gene_expression/sample_info.csv")
sample_info <- filter(sample_info, lineage == "breast") 


#antra

sample_info <- filter(sample_info, stripped_cell_line_name == "AU565" 
               | stripped_cell_line_name == "BT474"
               | stripped_cell_line_name == "BT483"
               | stripped_cell_line_name == "BT549"
               | stripped_cell_line_name == "CAL51"
               | stripped_cell_line_name == "CAL120"
               | stripped_cell_line_name == "CAL148"
               | stripped_cell_line_name == "CAMA1"
               | stripped_cell_line_name == "EFM192A"
               | stripped_cell_line_name == "HCC38"
               | stripped_cell_line_name == "HCC70"
               | stripped_cell_line_name == "HCC202"
               | stripped_cell_line_name == "HCC1419"
               | stripped_cell_line_name == "HCC1428"
               | stripped_cell_line_name == "HCC1500"
               #| stripped_cell_line_name == "HCC1937"
               | stripped_cell_line_name == "HCC1954"
               | stripped_cell_line_name == "HCC2218"
               | stripped_cell_line_name == "JIMT1"
               #| stripped_cell_line_name == "MCF7"
               | stripped_cell_line_name == "MDAMB175VII"
               #| stripped_cell_line_name == "MDA-MB-231"
               #| stripped_cell_line_name == "MDA-MB-436"
               | stripped_cell_line_name == "MDAMB468"
               | stripped_cell_line_name == "UACC812" 
               | stripped_cell_line_name == "ZR7530")


ID <- (sample_info$DepMap_ID)

CCLE_expression$expression <- ifelse(CCLE_expression$X %in% ID, "yes", "no") 
CCLE_expression <- filter(CCLE_expression, expression == "yes")
sample_info1 <- select(sample_info, DepMap_ID, stripped_cell_line_name )
colnames(sample_info1) <- c("X", "cell_line_name")
CCLE_expression <- merge(CCLE_expression, sample_info1, by = "X")



colnames(CCLE_expression) <- gsub("\\..*", "", colnames(CCLE_expression))



############ gene cross #############################

alky <- filter(alky, differential != "non_differential" & SYMBOL != "-")
antra <- filter(antra, differential != "non_differential" & SYMBOL != "-")
parpi <- filter(parpi, differential != "non_differential" & SYMBOL != "-")

parpi <- select(parpi, SYMBOL)
parpi_unique <- unique(parpi$SYMBOL)

antra <- select(antra, SYMBOL)
antra_unique <- unique(antra$SYMBOL)

alky <- select(alky, SYMBOL)
alky_unique <- unique(alky$SYMBOL)

gene <- merge(parpi, antra, all = TRUE)
gene <- merge(gene, alky, all = TRUE)


gene1 <- unique(gene$SYMBOL)


######## fin cross gene ##########
library(dplyr)
#gene_high <- as.data.frame(gene_high)

data <- CCLE_expression[, colnames(CCLE_expression) %in% gene1]
data <- cbind(data, cell_line_name = CCLE_expression$cell_line_name)

#data <- CCLE_expression[, colnames(CCLE_expression) %in% gene_low]

#CCLE_expression1$PARPi <- as.factor(CCLE_expression1$PARPi)





library(stringr)
library(dplyr)
library(tidyr)
require(reshape)
library(viridisLite)
library(viridis)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(car)
library(purrr)
# Load required R packages
library(highcharter)
# Set highcharter options
options(highcharter.theme = hc_theme_smpl(tooltip = list(valueDecimals = 2)))




CCLE <- CCLE_expression1[-c(1 , 6856)]
CCLE <- colnames(CCLE)
CCLE <- as.character(CCLE)


CCLE1 <- CCLE_expression1[-c(1 )]
CCLE1 <- colnames(CCLE1)
CCLE1 <- as.character(CCLE1)


########## antra high ##########


library(dplyr)
library(tidyr)
#data_processed <- CCLE_expression1 %>%
#  select(all_of(CCLE1)) %>%
#  pivot_longer(cols = all_of(CCLE), names_to = "Gene", values_to = "Gene_expression")




########### test ##############
library(ggplot2)


id_cell_lines <- read.delim("~/Desktop/PAPER1/id_cell_lines.txt")
colnames(id_cell_lines) <- c("cell_line_name",
                             "PARPi", "Anthracyclines", "AlkylatingAgents")
id_cell_lines <- select(id_cell_lines, cell_line_name, Anthracyclines)
CCLE_expression1 <- merge(data, id_cell_lines, by = "cell_line_name")



CCLE <- CCLE_expression1[-c(1 , 6856)]
CCLE <- colnames(CCLE)
CCLE <- as.character(CCLE)


CCLE1 <- CCLE_expression1[-c(1 )]
CCLE1 <- colnames(CCLE1)
CCLE1 <- as.character(CCLE1)

data_processed <- CCLE_expression1 %>%
  select(all_of(CCLE1), Anthracyclines) %>%
  pivot_longer(cols = -Anthracyclines, names_to = "Gene", values_to = "Gene_expression")




data_processed[data_processed == ""] <- NA


data_processed$Anthracyclines<-as.factor(data_processed$Anthracyclines) 
data_processed<-data_processed[!is.na(data_processed$Anthracyclines),]


stat.test <- data_processed %>%
  group_by(Gene) %>%
  wilcox_test(Gene_expression ~ Anthracyclines) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")
stat.test
write.table(stat.test, "wilcox_antra_dic.txt", row.names = FALSE)


stat.test <- data_processed %>%
  group_by(Gene) %>%
  t_test(Gene_expression ~ Anthracyclines) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")
stat.test
write.table(stat.test, "tt_antra_dic.txt", row.names = FALSE)


####### PARPI #####################

id_cell_lines <- read.delim("~/Desktop/PAPER1/id_cell_lines.txt")
colnames(id_cell_lines) <- c("cell_line_name",
                             "PARPi", "Anthracyclines", "AlkylatingAgents")
id_cell_lines <- select(id_cell_lines, cell_line_name, PARPi)
CCLE_expression1 <- merge(data, id_cell_lines, by = "cell_line_name")

data_processed <- CCLE_expression1 %>%
  select(all_of(CCLE1), PARPi) %>%
  pivot_longer(cols = -PARPi, names_to = "Gene", values_to = "Gene_expression")




data_processed[data_processed == ""] <- NA
data_processed$PARPi<-as.factor(data_processed$PARPi) 
data_processed<-data_processed[!is.na(data_processed$PARPi),]


stat.test <- data_processed %>%
  group_by(Gene) %>%
  wilcox_test(Gene_expression ~ PARPi) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")
stat.test

write.table(stat.test, "wilcox_parpi_dic.txt", row.names = FALSE)

stat.test <- data_processed %>%
  group_by(Gene) %>%
  t_test(Gene_expression ~ PARPi) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

write.table(stat.test, "tt_parpi_dic.txt", row.names = FALSE)
#####


id_cell_lines <- read.delim("~/Desktop/PAPER1/id_cell_lines.txt")
colnames(id_cell_lines) <- c("cell_line_name",
                             "PARPi", "Anthracyclines", "AlkylatingAgents")
id_cell_lines <- select(id_cell_lines, cell_line_name, AlkylatingAgents)
CCLE_expression1 <- merge(data, id_cell_lines, by = "cell_line_name")



CCLE <- CCLE_expression1[-c(1 , 6856)]
CCLE <- colnames(CCLE)
CCLE <- as.character(CCLE)

CCLE1 <- CCLE_expression1[-c(1 )]
CCLE1 <- colnames(CCLE1)
CCLE1 <- as.character(CCLE1)

data_processed <- CCLE_expression1 %>%
  select(all_of(CCLE1), AlkylatingAgents) %>%
  pivot_longer(cols = -AlkylatingAgents, names_to = "Gene", values_to = "Gene_expression")


data_processed[data_processed == ""] <- NA


data_processed$AlkylatingAgents<-as.factor(data_processed$AlkylatingAgents) 
data_processed<-data_processed[!is.na(data_processed$AlkylatingAgents),]


stat.test <- data_processed %>%
  group_by(Gene) %>%
  wilcox_test(Gene_expression ~ AlkylatingAgents) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")
stat.test
write.table(stat.test, "wilcox_alky_dic.txt", row.names = FALSE)


stat.test <- data_processed %>%
  group_by(Gene) %>%
  t_test(Gene_expression ~ AlkylatingAgents) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")
stat.test

write.table(stat.test, "tt_alky_dic.txt", row.names = FALSE)



######


stat_sig <- filter(stat.test, p.signif != "ns")
stat_sig1 <- c(stat_sig$Gene)
data_processed$sig <- ifelse(data_processed$Gene %in% stat_sig1, "yes", "no")
data_processed1 <- filter(data_processed, sig == "yes")


# Crear un box plot
bxp <- ggboxplot(
  data_processed1, x = "Gene", y = "Gene_expression", 
  color = "Anthracyclines", palette = c("#00AFBB", "#E7B800"), add = "jitter"
) +
  facet_grid(Anthracyclines ~ ., scales = "free_x")
# agregar p-values a  los box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Gene", dodge = 0.5)
p3 <- bxp + stat_pvalue_manual(
  stat.test, label = "{p}", y.position = 12,
  tip.length = 0.001
) 
p4 <- p3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("p4_antra_high.pdf", width = 10, height = 6)  # Cambia el nombre y la extensión según tu preferencia
# Tu código para crear el gráfico (por ejemplo, bxp <- ggplot() + ...)

# Dibuja el gráfico en el archivo PDF
print(p4)

# Cierra el dispositivo PDF
dev.off()




