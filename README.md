# RNA_editing_drugresponse_breastcancer

# I. Differential A>I(G) RNA-edited sites (DES) from Breast Cancer cell lines    

## Data Filters for the identification of differential A>I(G) RNA-edited sites (DES)   
Base counts per site with allele AG or TC.

Apply filters for the identification of RNA edited sites (Master_table_CCLE, see code):

> a. (Coverage at position >= 10)   
> b. (Variant == "AG" OR Variant == "TC")   
> c. (Coverage of alternative allele >= 5)   
> d. (Frequency of alternative allele >= 0.05)   

Merge cell line data into three groups (PARPi, Anthracyclines, Alkylating agents).


### DNA:RNA pair for to exclude variant in DNA  

Exclude all sites reported in DepMap (using version 19Q3:hg19) (https://depmap.org/portal/download/all/?releasename=DepMap+Public+19Q3)    
Downloaded sample_info.csv & CCLE_mutations.csv    

> a. Load data from DepMapPortal (https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=sample_info.csv)          
> b. Load data 19Q3:hg19 from DepMapPortal (https://depmap.org/portal/download/all/?releasename=DepMap+Public+19Q3&filename=CCLE_mutations.csv)           


#### PARP inhibitors 

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

#### Anthracyclines Group:

DNA data
> a. 1,239,235 total   
> b. 35,520 lineage == breast   
> c. 9,386 selected cell lines   
> d. 1,703 Reference_Allele == "A" | Reference_Allele == "T"   
> e. 840 Reference_Allele == "A" & Strand == "+" | Reference_Allele == "T" & Strand == "-"   
> f. 795 Variant_Type == "SNP"   
> g. 508 allele == "A/G" & Strand == "+" | allele == "T/C" & Strand == "-"   
> h. 193 match with 53,043 (position)   

RNA data from Master table
> 54,760 A>I(G) RNA edited sites - 193 variants in DNA = 54,567 RNA-edited sites   

#### Alkylating agents Group:

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

### Apply REDITs test (https://github.com/gxiaolab/REDITs) 
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
### Characterization of A>I(G) RNA-edited sites

> a. Annotate and predict the functional consequences with VEP hg19 from Variant effect Predictor (https://grch37.ensembl.org/Homo_sapiens/Tools/VEP): our results by filtering for one selected consequence per variant allele.   
> b. Cross-reference with COSMIC GENE CONSENSUS (https://cancer.sanger.ac.uk/census)   
> c. RepeatMasker version hg19 from (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz) 





### Summary of sites identified 

| summary           |     PARPi         | Anthracyclines |  Alkylating agents |
|-------------------------|-------------------|----------------|--------------------|
| n group                 |    L= 10 ; H= 8  |   L= 8 ; H=7   |     L=10 ; H=9     |
|  All A>I(G)                  |      53,043       |     54,760     |  52,455            |
|  DNA out  A/G SNP         |      52,818       |     54,567     | 56,662             |
|  diff p-value/FDR/FC    |      12,541       |   10,015       | 13,384             |
|  diff L                 |     6,996         |       5,320    |      6,152     |
|  diff H                 |     5,277         |       4,695   |     7,232       |
| liftover to hg38          |     12,273       |              |            |


# II. A>I(G) RNA-edited sites in women with breast cancer      

## Identification of RNA editing level    
## Machine learning algorithms using A>I(G) RNA edited sites and clinical data     


