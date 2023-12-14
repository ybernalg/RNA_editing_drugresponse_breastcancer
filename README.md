# RNA_editing_drugresponse_breastcancer

## Data Filters for the identification of differential A>I(G) RNA-edited sites (DES)
Base counts per site with allele AG or TC.

Apply filters for the identification of RNA edited sites (Master_table_CCLE, see code):

> a. (Coverage at position >= 10)   
> b. (Variant == "AG" OR Variant == "TC")   
> c. (Coverage of alternative allele >= 5)   
> d. (Frequency of alternative allele >= 0.05)   

Merge cell line data into three groups (PARPi, Anthracyclines, Alkylating agents).


## DNA:RNA pair for to exclude variant in DNA  

Exclude all sites reported in DepMap (using version 19Q3:hg19): https://depmap.org/portal/download/all/?releasename=DepMap+Public+19Q3
Downloaded sample_info.csv & CCLE_mutations.csv

> a. Load data from DepMapPortal (https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=sample_info.csv)     
> b. Load data 19Q3:hg19 from DepMapPortal (https://depmap.org/portal/download/all/?releasename=DepMap+Public+19Q3&filename=CCLE_mutations.csv)     


### PARP inhibitors 

DNA data
> a. 1,239,235 total   
> b. 35,520 lineage == "breast"
> c. 14,397 selected cell lines   
> d. 2,949 Reference_Allele == "A" | Reference_Allele == "T"   
> e. 1,450 Reference_Allele == "A" & Strand == "+" | Reference_Allele == "T" & Strand == "-"    
> f. 1,394 Variant_Type == "SNP"   
> g. 630 allele == "A/G" & Strand == "+" | allele == "T/C" & Strand == "-"   
> h. 225 match with 53,043 (position)    

RNA data from Master table
> 53,043 A>I(G) RNA edited sites - 225 variants in DNA = 52,818 RNA-edited sites   

### Anthracyclines Group:

DNA data
> a. 1,239,235 total
> b. 35,520 lineage == breast
> c. 9,386 selected cell lines
> d. 1,703 Reference_Allele == "A" | Reference_Allele == "T"
> e. 840 Reference_Allele == "A" & Strand == "+" | Reference_Allele == "T" & Strand == "-"
> f. 795 Variant_Type == "SNP"
> g. 508 allele == "A/G" & Strand == "+" | allele == "T/C" & Strand == "-"
> h. 225 match with 53,043 (position)

RNA data from Master table


### Alkylating agents Group:

DNA data
> a. 1,239,235 total   
> b. 35,520 lineage == breast   
> c. 10,751 selected cell lines    
> d. 2,119 Reference_Allele == "A" | Reference_Allele == "T"    
> e. 1,039 Reference_Allele == "A" & Strand == "+" | Reference_Allele == "T" & Strand == "-"   
> f. 988 Variant_Type == "SNP"    
> g. 639 allele == "A/G" & Strand == "+" | allele == "T/C" & Strand == "-"    
> h. 237 match with 53,043 (position)    

RNA data from Master table

# Apply REDITs test (https://github.com/gxiaolab/REDITs) 
From the output, we created this variables for the selection of sites differentially edited 
```
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

```
# Characterization of A>I(G) RNA-edited sites

> Annotate with VEP hg19 from Variant effect Predictor (https://grch37.ensembl.org/Homo_sapiens/Tools/VEP)
> Cross-reference with COSMIC GENE CONSENSUS
> RepeatMasker from 





# Summary of sites identified 

| summary           |     PARPi         | Anthracyclines |  Alkylating agents |
|-------------------------|-------------------|----------------|--------------------|
| n group                 |    L= 10 ; H= 8  |   L= 8 ; H=7   |     L=10 ; H=9     |
|   Todos                 |      53,043       |     54,760     |  52,455            |
|  DNA out                |      52,818       |     54,567     | 56,425             |
|  diff p-value           |      17,679       |     14,852     | 18,792             |
|  diff p-value/FDR       |      16,708       |    13,425      | 17,730             |
|  diff p-value/FDR/FC    |      12,541       |   10,015       | 13,158             |
|  diff L                 |     6,996         |           |              |
|  diff H                 |     5,277         |          |             |
| liftover to hg38          |     12,273       |         9,920       | 13,016               |





