library(dplyr)
setwd("/PARPi/")


#load data 19Q3:hg19 from https://depmap.org/portal/download/all/?releasename=DepMap+Public+19Q3&filename=CCLE_mutations.csv
ccle <- read.csv("~/PARPi/input/CCLE_mutations_19Q3.csv")

# load data from https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=sample_info.csv 
sample_info <- read.csv("~/PARPi/input/sample_info.csv")
data <- merge(ccle, sample_info)



data <- filter(data, lineage == "breast")
data <- filter(data, cell_line_name == "AU565" 
               | cell_line_name == "BT-474"
               | cell_line_name == "BT-483"
               | cell_line_name == "BT-549"
               | cell_line_name == "CAL-51"
               | cell_line_name == "CAL-120"
               | cell_line_name == "CAL-148"
               | cell_line_name == "CAMA-1"
               | cell_line_name == "HCC38"
               | cell_line_name == "HCC70"
               | cell_line_name == "HCC1419"
               | cell_line_name == "HCC1428"
               | cell_line_name == "HCC1500"
               | cell_line_name == "HCC1954"
               | cell_line_name == "JIMT-1"
               | cell_line_name == "MDA-MB-175-VII"
               | cell_line_name == "MDA-MB-468"
               | cell_line_name == "UACC-812" 
               | cell_line_name == "ZR-75-30")

data <- filter(data, Reference_Allele == "A" | Reference_Allele == "T")
data <- filter(data, Reference_Allele == "A" & Strand == "+" 
               | Reference_Allele == "T" & Strand == "-" )
data <- filter(data, Variant_Type == "SNP")

data$allele <- paste(data$Reference_Allele, data$Tumor_Seq_Allele1, sep = "/")

DNA <- filter(data, allele == "A/G" & Strand == "+" | allele == "T/C" & Strand == "-")
library(dplyr)
DNA$up <- paste(DNA$Chromosome, DNA$Start_position, sep = "_")
DNA$Uploaded_variation <- paste(DNA$up, DNA$allele, sep = "_")
DNA$up <- NULL



DNA <- c(DNA$Uploaded_variation)




####### load master table data ##############################

data <- read.delim2("~/PARPi/input/CCLE24_master_table_PARPi_groupCounts.tsv")


library(tidyr)
#install.packages("stringr")
library("stringr")

data <- filter(data, HS_count > 4 & LS_count >5)
data$ensembl <- ifelse(data$Reference == "A", "A/G", "T/C")
data$up <- paste(data$Region, data$Position, sep = "_")
data$Uploaded_variation <- paste(data$up, data$ensembl, sep = "_")
data$up <- NULL

#I remove from my data those variants reported in DNA, according to chromosome, position, and strand

data$dna_status <- ifelse(data$Uploaded_variation %in% DNA, "yes", "no")

#table(data$dna_status)
#no   yes 
#52818   225 

data <- filter(data, dna_status != "yes")


######## preparacion para test #################
data <- data %>% 
  select("Uploaded_variation", "ensembl", contains("bam"))


colnames(data) <- c("Uploaded_variation","ensembl",
                    "AU565","CAL148", "CAL51","CAMA1","HCC1954","HCC70",  "MDAMB468",
                    "UACC812","BT474","BT483", "BT549","CAL120" , "HCC1419",
                    "HCC1428", "HCC1500", "HCC38","MDAMB175VII", "ZR7530")


data <- as.data.frame(lapply(data, function(y) gsub("\\[", "", y)))
data <- as.data.frame(lapply(data, function(y) gsub("\\]", "", y)))



data[ , 21:24] <- str_split_fixed(data$AU565, ",", 4)
data[ , 25:28] <- str_split_fixed(data$CAL148, ",", 4)
data[ , 29:32] <- str_split_fixed(data$CAL51, ",", 4)
data[ , 33:36] <- str_split_fixed(data$CAMA1, ",", 4)
data[ , 37:40] <- str_split_fixed(data$HCC1954, ",", 4)
data[ , 41:44] <- str_split_fixed(data$HCC70, ",", 4)
data[ , 45:48] <- str_split_fixed(data$MDAMB468, ",", 4)
data[ , 49:52] <- str_split_fixed(data$UACC812, ",", 4)
data[ , 53:56] <- str_split_fixed(data$BT474, ",", 4)
data[ , 57:60] <- str_split_fixed(data$BT483, ",", 4)
data[ , 61:64] <- str_split_fixed(data$BT549, ",", 4)
data[ , 65:68] <- str_split_fixed(data$CAL120, ",", 4)
data[ , 69:72] <- str_split_fixed(data$HCC1419, ",", 4)
data[ , 73:76] <- str_split_fixed(data$HCC1428, ",", 4)
data[ , 77:80] <- str_split_fixed(data$HCC1500, ",", 4)
data[ , 81:84] <- str_split_fixed(data$HCC38, ",", 4)
data[ , 85:88] <- str_split_fixed(data$MDAMB175VII, ",", 4)
data[ , 89:92] <- str_split_fixed(data$ZR7530, ",", 4)

data <- subset(data, select = -c(3:20))


colnames(data) <- c("Uploaded_variation", "strand", 
                    "AU565_A",  "AU565_C",  "AU565_G",  "AU565_T",
                    "CAL148_A",  "CAL148_C",  "CAL148_G",  "CAL148_T",
                    "CAL51_A",  "CAL51_C",  "CAL51_G",  "CAL51_T",
                    "CAMA1_A", "CAMA1_C", "CAMA1_G", "CAMA1_T",
                    "HCC1954_A", "HCC1954_C", "HCC1954_G", "HCC1954_T",
                    "HCC70_A", "HCC70_C", "HCC70_G", "HCC70_T",
                    "MDAMB468_A",  "MDAMB468_C",  "MDAMB468_G",  "MDAMB468_T",
                    "UACC812_A", "UACC812_C", "UACC812_G", "UACC812_T",
                    "BT474_A", "BT474_C", "BT474_G", "BT474_T",
                    "BT483_A", "BT483_C", "BT483_G", "BT483_T",
                    "BT549_A", "BT549_C", "BT549_G", "BT549_T",
                    "CAL120_A", "CAL120_C", "CAL120_G", "CAL120_T",
                    "HCC1419_A",  "HCC1419_C",  "HCC1419_G",  "HCC1419_T",
                    "HCC1428_A", "HCC1428_C", "HCC1428_G", "HCC1428_T",
                    "HCC1500_A", "HCC1500_C", "HCC1500_G", "HCC1500_T", 
                    "HCC38_A", "HCC38_C", "HCC38_G", "HCC38_T",
                    "MDAMB175VII_A", "MDAMB175VII_C", "MDAMB175VII_G", "MDAMB175VII_T",
                    "ZR7530_A", "ZR7530_C", "ZR7530_G", "ZR7530_T")


cell_lines <- c( "AU565", "CAL148", "CAL51", "CAMA1", 
                 "HCC1954", "HCC70", "MDAMB468", "UACC812", "BT474", "BT483", "BT549", "CAL120", "HCC1419", "HCC1428", "HCC1500", 
                "HCC38", "MDAMB175VII", "ZR7530")

# Iterating over each cell line to create reads of the reference and alternate allele.
for (cell_line in cell_lines) {
  ref_col <- paste0(cell_line, "_ref")
  alt_col <- paste0(cell_line, "_alt")
  data[[ref_col]] <- ifelse(data$strand == "A/G", data[[paste0(cell_line, "_A")]], data[[paste0(cell_line, "_T")]])
  data[[alt_col]] <- ifelse(data$strand == "A/G", data[[paste0(cell_line, "_G")]], data[[paste0(cell_line, "_C")]])
}


data <- subset(data, select = -c(3:26))
data <- subset(data, select = -c(3:22))
data <- subset(data, select = -c(3:4))
data <- subset(data, select = -c(3:28))


data[data == ""] <- NA
###############################################################################
#### Statistical inference of differential RNA-editing sites ##################
###############from RNA-sequencing data by hierarchical modeling###############
################### doi: 10.1093/bioinformatics/btaa066 #######################
###############################################################################



# I apply a statistical test
### (A) I create the groups ###################################################
Low_sensitivity <-  c("BT474", "BT483", "BT549",
                      "CAL120", "HCC1419", "HCC1428",
                      "HCC1500",  "HCC38", "MDAMB175VII",
                      "ZR7530")

create.matrix <- function(site){
  ref <- c()
  alt <- c()
  for(i in seq(1,36, by=2)){
    ref.tmp <- site[[i]]
    alt.tmp <- site[[i+1]]
    ref <- append(ref, ref.tmp)
    alt <- append(alt, alt.tmp)
  }
  
  
  mat.editions <- rbind(alt, ref)
  colnames(mat.editions) <- c("AU565",  "CAL148", "CAL51",
                              "CAMA1", "HCC1954", "HCC70", "MDAMB468",
                              "UACC812", "BT474", "BT483", "BT549",
                              "CAL120", "HCC1419", "HCC1428",
                              "HCC1500",  "HCC38", "MDAMB175VII",
                              "ZR7530")
  #mat.editions[is.na(mat.editions)] <- 0
  #mat.editions[na.omit(mat.editions)]
  
  #  mat.ed.group <- rbind(colSums(mat.editions[1:3,]), colSums(mat.editions[4:6,]))
  return(mat.editions)
}


###(B)############ [ AG ]  CARGAR REDIT #####################
get_maximm_likelihood_parameters_beta_binomial_unimodal = function(data){
  #data is matrix. First row are G reads. Second row are A reads
  N = data[1,]+data[2,]
  y = data[1,]
  LL = function(par){
    alpha = par[1]
    beta = par[2]
    if(alpha < 1 | beta < 1 | alpha > 1e6 | beta > 1e6){ #boundaries of the optimization
      return(NA)
    }
    log_likelihood = -sum( lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) - lgamma(alpha+beta+N) + lgamma(alpha+y) + lgamma(beta+N-y) )
    return(log_likelihood)
  }
  parameters = optim(par=c(1,1), LL)
  log_likelihood = -1 * parameters$value
  all_parameters = parameters$par; names(all_parameters) = c('alpha','beta')
  converged = parameters$convergence
  output_object = list()
  output_object$output_parameters = all_parameters
  output_object$converged = converged
  output_object$log_likelihood = log_likelihood
  return(output_object)
}

#a generic function for log likelihood ratio test
run_likelihood_ratio_test = function(L0,La1, La2){
  #L0 is maximum log likelihood of null hypothesis
  #La1 is maximum log likelihood of alternative hypothesis for group1
  #La2 is max log likelihood of alternative hypothesis for group2
  chi_square_stat = -2*(L0 - (La1 + La2))
  df = 2
  p_value = pchisq(chi_square_stat,df=df,lower.tail=FALSE)
  return(p_value)
}

validate_input_data = function(data,groups){ #validates user input data
  #make sure no data is NA
  if(any(is.na(data))){
    stop("data cannot have NA values. You can probably substitute zero for missing (NA) values")
  }
  if(any(is.na(groups))){
    stop("groups cannot have NA values")
  }
  #make sure data is a matrix
  if(!is.matrix(data)){
    stop("data argument must be a matrix") #try on data = data.frame(apple=c(1,2),pears=c(2,3))
  }
  #make sure data is a 2xn matrix 
  if(nrow(data) != 2){
    stop("data argument must be a 2xn matrix") #try on data = matrix(c(1,2,3,2,3,4),nrow=3)
  }
  #make sure data is numeric, without decimals
  if(!all(is.numeric(data))){ 
    stop("elements in data must be integers") #try on data = matrix(c('1',2,3,4,5,6),nrow=2)
  }
  if(!all(data %% 1 == 0)){
    stop("elements in data must be integers") #try on data = matrix(c(0.2,2,1,7),nrow=2)
  }
  #make sure groups is a vector
  if(!is.vector(groups)){
    stop("groups must be a vector") #test on groups = matrix(c(1,2,3,4))
  }
  if(!is.character(groups)){
    stop("groups must be a chracter vector") #test on groups = as.factor(c('a','b','a','b')). groups=c(1,2,3,4)
  }
  #make sure number of groups and columns in data match 
  if(ncol(data) != length(groups)){
    stop("length(groups) must equal ncol(data)") 
  }
  #make sure there are only two groups
  if(length(unique(groups)) != 2){
    stop("must have exactly two groups") #test on groups=c('a','b','c','a'); data= matrix(c(1,2,3,4,5,6,7,8),nrow=2)
  }
}

#splits data by experimental condition
split_data_by_group = function(data,groups){
  type_groups = unique(groups)
  data1 = data[,groups==type_groups[1],drop=FALSE]
  data2 = data[,groups==type_groups[2],drop=FALSE]
  group1 = type_groups[1]
  group2 = type_groups[2]
  return(list(data1=data1,data2=data2,group1 = group1,group2=group2))
}

#the actual function for running the REDIT-LLR
REDIT_LLR = function(data,groups){
  validate_input_data(data,groups)
  split_data = split_data_by_group(data,groups)
  data1 = split_data$data1                                                                                               
  data2 = split_data$data2
  group1 = split_data$group1
  group2 = split_data$group2
  model_fit0 = get_maximm_likelihood_parameters_beta_binomial_unimodal(data=data)
  model_fit1 = get_maximm_likelihood_parameters_beta_binomial_unimodal(data=data1)
  model_fit2 = get_maximm_likelihood_parameters_beta_binomial_unimodal(data=data2)
  parameters0 = model_fit0$output_parameters
  parameters1 = model_fit1$output_parameters
  parameters2 = model_fit2$output_parameters
  L0 = model_fit0$log_likelihood
  La1 = model_fit1$log_likelihood
  La2 = model_fit2$log_likelihood
  p_value = run_likelihood_ratio_test(L0=L0,La1=La1,La2=La2)
  output_object = list()
  output_object$data = data
  output_object$groups = groups
  output_object[[paste0('mle.for.group.',group1)]] = parameters1
  output_object[[paste0('mle.for.group.',group2)]] = parameters2
  output_object[[paste0('mle.for.null.model')]] = parameters0
  output_object[[paste0('log.likelihood.for.group.',group1)]] = La1
  output_object[[paste0('log.likelihood.for.group.',group2)]] = La2
  output_object[[paste0('log.likelihood.for.null')]] = L0
  output_object[['p.value']] = p_value
  return(output_object)
}




####(C)########### [ AG ]  APLICAR #####################

source("REDIT_LLR.R")
out.df <- data.frame()

#max_rows <- nrow(prueba)


# Preallocate the data frame with empty rows
#out.df <- data.frame(matrix(NA, ncol = 47, nrow = max_rows))

for(i in 1:nrow(data)){
  mat.ed <- create.matrix(data[i,3:38])
  mat.ed <- mat.ed[ , colSums(is.na(mat.ed))==0]
  group <- ifelse(colnames(mat.ed) %in% Low_sensitivity, "Low", "High")
  #group <-as.numeric(as.factor(group)) 
  group <-as.factor(group)
  group <- as.vector(group) 
  mat.ed <- matrix(as.numeric(mat.ed), ncol = length(colnames(mat.ed)))
  test.tmp <- REDIT_LLR(data=mat.ed, groups=group)
  out.df <- rbind(out.df, data.frame(data[i,], pvalue=test.tmp$p.value,
                                     mle_group_resistance_alpha=test.tmp$mle.for.group.Low[[1]],
                                     mle_group_resistance_beta=test.tmp$mle.for.group.Low[[2]],
                                     mle_group_sensitive_alpha=test.tmp$mle.for.group.High[[1]],
                                     mle_group_sensitive_beta=test.tmp$mle.for.group.High[[2]],
                                     mle_null_model=test.tmp$mle.for.null.model,
                                     log_l_resistance=test.tmp$log.likelihood.for.group.Low,
                                     log_l_sensitive=test.tmp$log.likelihood.for.group.High,
                                     log_l_null=test.tmp$log.likelihood.for.null))
}



data <- out.df[!duplicated(out.df$Uploaded_variation), ]


data$EL_R <- data$mle_group_resistance_alpha/(data$mle_group_resistance_alpha + data$mle_group_resistance_beta) 
data$EL_S <- data$mle_group_sensitive_alpha/(data$mle_group_sensitive_alpha + data$mle_group_sensitive_beta) 
data$delta <- data$EL_R - data$EL_S
data$fdr <- p.adjust(data$pvalue, method="BH")
data$class <- ifelse(data$delta >= 0, "Low_sensitivity", "High_sensitivity")
data$fold_change = ((data$EL_R/data$EL_S))
data$fold_change_log = (log(data$fold_change, 2))



data$differential <- ifelse(data$pvalue < 0.01 & data$fdr < 0.1 & data$fold_change_log >= 2.5, 
                            "low_sensitivity",
                            ifelse(data$pvalue < 0.01 & data$fdr < 0.1 & data$fold_change_log <= - 2.5, 
                                   "high_sensitivity", 
                                   "non_differential"))


table(data$differential)


#high_sensitivity  low_sensitivity non_differential 
#5277             6996            40545 

data2 <- filter(data, differential != "non_differential")
#12,273

output_vep_parpi <- read.delim("~/Desktop/PAPER1/output/output_vep_parpi.txt")

data <- merge(data, output_vep_parpi, by = "Uploaded_variation", all.x = TRUE)

library(data.table)
rmsk <- read.delim("~/Desktop/PAPER1/input/rmsk.txt.gz", header=FALSE)

rmsk$start = (rmsk$V7 - 150)
rmsk$end = (rmsk$V8 + 150)
library(data.table)
rmsk$V6<-gsub("chr","", rmsk$V6)
# convert your datasets in data.table
setDT(data) 
setDT(rmsk)

#V6 chr
#V7 lim min
#V8 lim max
## merge de "data" con "rmsk" , segun condicionales de on = [igual chromosona, y POSITION de "data" entre V7 y V8 de "rmsk"
merge <- data[rmsk, on = .(Region == V6 , Position >= start,               
                           Position <= end), nomatch = NA,
              .(Uploaded_variation,V6,start,end,V10,V11, V12, V13)] ## selecion de columnas para el nuevo dataframe "merge"
merge1 <- merge[!is.na(merge$Uploaded_variation),]
merge2 <- merge1[!duplicated(merge1), ]

data2 <- merge(x = data, y = merge2, all.x = TRUE, no.dups = TRUE)

nrow(data2[duplicated(data2$Uploaded_variation), ])
data <- data2[!duplicated(data2$Uploaded_variation), ]

##### cosmic data #############################################################

cosmic <- read.delim("~/udd.cl/Ricardo Armisen - Yanara/BMC_Journal of Translational Medicine/ADDITIONAL_ANALYSIS/Census_allTue Jul 18 14_56_48 2023.tsv")
data$cosmic_gene <- ifelse(data$SYMBOL %in% cosmic$Gene.Symbol, "si", "no")

data <- merge(data, cosmic,
              by.x=c("SYMBOL"),
              by.y=c("Gene.Symbol"), all.x = TRUE)


