# miRSM 2.0_case_study
The input of the case study of miRSM 2.0 includes single-cell miRNA-mRNA co-sequencing data, and prior information of miRNA targets.

Script.R: Code for reproducibility of the analysis in the case study

GSE114071_NW_scsmRNA_K562_norm_log2.gct: miRNA expression data

GSE114071_NW_half_cell_K562_RNAseq_processed_data.gct: mRNA expression data

miRTarBase_v9.0.csv: Putative miRNA-mRNA interactions in miRTarBase v9.0

TarBase_v8.0.csv: Putative miRNA-mRNA interactions in TarBase v8.0

CML_genes.csv: A list of CML-related genes (including miRNAs and mRNAs)

# Troubleshooting
## Problem 1:
Some large required dependencies (e.g. org.Hs.eg.db, reactome.db) are not successfully installed when installing the miRSM R package.

## Potential solution: 
Install the miRSM R package under a good internet speed or download the source package of the large required dependencies to be locally installed. 

## Problem 2:
The non-negative matrix factorization (NMF) method fails to identify gene modules.

## Potential solution: 
The NMF method is limited to non-negative matrix. Please check if the input RNA-seq data is a non-negative matrix. If the input RNA-seq data is not a non-negative matrix, the data is recommended to remove negative values. If the input RNA-seq data is a non-negative matrix, it is recommended to add a very small constant to the data e.g. 1.0E-06.

## Problem 3:
The sensitivity canonical correlation (SCC) method fails to infer miRNA sponge modules.

## Potential solution: 
There is a problem computing SVD (singular value decomposition) for the SCC method. When the input RNA-seq data contains many missing values or zero values, the data is recommended to be imputed. Please refer to the corresponding Data preparation for details.

Contact: zjp@dali.edu.cn
