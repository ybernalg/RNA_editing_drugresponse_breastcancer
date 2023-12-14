# RNA_editing_drugresponse_breastcancer

Data Filters
Base counts per site with allele AG or TC.

Apply filters:

a. (Coverage at position >= 10)
b. (Variant == "AG" OR Variant == "TC")
c. (Coverage of alternative allele >= 5)
d. (Frequency of alternative allele >= 0.05)

Merge cell line data into three groups (PARPi, Anthracyclines, Alkylating agents).

Apply REDITs test (https://github.com/gxiaolab/REDITs) and decide whether to exclude NAs or replace them with 0:
a. Method 1 ==> Exclude NAs, group sizes change at each editing site.
b. Method 2 ==> Replace NAs with zeros.

On September 7th, Method 1 is selected.

Exclude all editing sites reported in DepMap (using version 19Q3:hg19 since 22Q4:hg38): https://depmap.org/portal/download/all/?releasename=DepMap+Public+19Q3
Downloaded sample_info.csv & CCLE_mutations.csv
a. 1,239,235 total
b. 35,520 lineage == breast
c. 14,397 selected cell lines
d. 2,949 Reference_Allele == "A" | Reference_Allele == "T"
e. 1,450 Reference_Allele == "A" & Strand == "+" | Reference_Allele == "T" & Strand == "-"
f. 1,394 Variant_Type == "SNP"
g. 630 allele == "A/G" & Strand == "+" | allele == "T/C" & Strand == "-"
h. 225 match with 53,043 (position)

B. Result 53,043 - 225 = 52,818 editing sites

Annotate with VEP hg19.

Select significant sites:
A. 17,679 p-value < 0.05 from REDITs test.
B. 16,708 p-value < 0.05 & FDR <= 0.1 with BH correction.
C. 5996 Overedited in low sensitivity + 6545 Overedited in high sensitivity p-value < 0.05 & FDR <= 0.1 & [>=2.5 <= -2.5 log fold change].

Cross-reference with COSMIC GENE CONSENSUS:
A. All: no (n= 50121) yes (n=2650)
B. 692 cosmic_gene == "yes" & expression != "Stable"

Anthra Group:

Downloaded sample_info.csv & CCLE_mutations.csv
a. 1,239,235 total
b. 35,520 lineage == breast
c. 9,386 selected cell lines
d. 1,703 Reference_Allele == "A" | Reference_Allele == "T"
e. 840 Reference_Allele == "A" & Strand == "+" | Reference_Allele == "T" & Strand == "-"
f. 795 Variant_Type == "SNP"
g. 508 allele == "A/G" & Strand == "+" | allele == "T/C" & Strand == "-"
h. 225 match with 53,043 (position)

Alkylating Group:

Downloaded sample_info.csv & CCLE_mutations.csv
a. 1,239,235 total
b. 35,520 lineage == breast
c. 10,751 selected cell lines
d. 2,119 Reference_Allele == "A" | Reference_Allele == "T"
e. 1,039 Reference_Allele == "A" & Strand == "+" | Reference_Allele == "T" & Strand == "-"
f. 988 Variant_Type == "SNP"
g. 639 allele == "A/G" & Strand == "+" | allele == "T/C" & Strand == "-"
h. 237 match with 53,043 (position)





