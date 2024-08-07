######################################################################################################################################################################################
########################################################################### miRSM 2.0 application in K562 dataset ####################################################################
######################################################################################################################################################################################

## Function for computing the average expression values of duplicate genes
# Exp_scRNA: Gene expression values of miRNAs or mRNAs in single cells, rows are cells and 
# columns are miRNAs or mRNAs
# Output: temp is single cell expression data without duplicate genes
Averg_Duplicate <- function(Exp_scRNA){    
      uniqueNameList <- unique(colnames(Exp_scRNA))
      noOfgenes <- length(uniqueNameList)
      temp <- matrix(0, nrow = nrow(Exp_scRNA), ncol = noOfgenes)
      colnames(temp) <- uniqueNameList
      rownames(temp) <- rownames(Exp_scRNA)
      for(c in 1:noOfgenes){
          GeneList <- which(colnames(Exp_scRNA) == colnames(temp)[c])
          for(r in 1:nrow(temp)) {
              temp[r, c] <- mean(as.numeric(Exp_scRNA[r, GeneList]))  
          }
      }
    return(temp)
}

########################################################################### Step 1: Installation of miRSM 2.0 ####################################################################
## Install and load miRSM 2.0 R package
if (!require("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") 
BiocManager::install("miRSM")
suppressMessages(library(miRSM))
packageVersion("miRSM")

################################################################################ Step 2: Data preparation #########################################################################
## The single-cell miRNA-mRNA co-sequencing data (GSE114071_NW_scsmRNA_K562_norm_log2.gct.gz for miRNA expression data, 
## and GSE114071_NW_half_cell_K562_RNAseq_processed_data.gct.gz for mRNA expression data) is from the Gene Expression Omnibus 
## (GEO, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114071). The obtained miRNA and mRNA expression matrix data is further transposed.
## Please paste all downloaded files into a single folder (set the folder as the directory of R environment).
miR_data <- read.table("GSE114071_NW_scsmRNA_K562_norm_log2.gct", header = TRUE, sep='\t')
mR_data <-  read.table("GSE114071_NW_half_cell_K562_RNAseq_processed_data.gct", header = TRUE, skip = 2, sep='\t')
miR_names <- miR_data[, 1]
mR_names <- mR_data[, 1]
miR_data <- miR_data[, -seq(2)]
mR_data <- mR_data[, -seq(2)]
miR_data[is.na(miR_data)] <- 0
mR_data[is.na(mR_data)] <- 0
miR_data_transposed <- t(miR_data)
mR_data_transposed <- t(mR_data)

## To ensure the consistency of the downloaded miRNA and mRNA expression data, the downloaded miRNA expression data is firstly processed with transformation. 
## Moreover, the unmatched K562 cells between miRNA and mRNA expression data are removed. Finally, the matched miRNA and mRNA expression data are further normalized.
miR_data_transposed <- 2^(miR_data_transposed)
miR_scRNA_matched <- miR_data_transposed[seq(19), ] 
mR_scRNA_matched <- mR_data_transposed[-c(6, 21, 22), ]
miR_scRNA_norm <- log2(miR_scRNA_matched + 1)
mR_scRNA_norm <- log2(mR_scRNA_matched + 1)

## Preprocess miRNA and mRNA expression data to SummarizedExperiment format using the SummarizedExperiment R package.
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
colData_miR <- DataFrame(row.names = miR_names)
colData_mR <- DataFrame(row.names = mR_names)
miR_scRNA_norm <- SummarizedExperiment(assays=list(miR_scRNA_norm = miR_scRNA_norm), colData = colData_miR)
mR_scRNA_norm <- SummarizedExperiment(assays=list(mR_scRNA_norm = mR_scRNA_norm), colData = colData_mR)

## Compute the average expression values of duplicate genes
miR_scRNA_average <- Averg_Duplicate(assay(miR_scRNA_norm))
mR_scRNA_average <- Averg_Duplicate(assay(mR_scRNA_norm))

## Remove genes with constant expression values in all cells
miR_scRNA_sd <- unlist(lapply(seq(dim(miR_scRNA_average)[2]), function(i) sd(miR_scRNA_average[, i])))
mR_scRNA_sd <- unlist(lapply(seq(dim(mR_scRNA_average)[2]), function(i) sd(mR_scRNA_average[, i])))
miR_scRNA_filter <- miR_scRNA_average[, which(miR_scRNA_sd > 0)]
mR_scRNA_filter <- mR_scRNA_average[, which(mR_scRNA_sd > 0)]

## Reserve high expression and variable miRNAs and mRNAs
miR_me <- unlist(lapply(seq(ncol(miR_scRNA_filter)), function(i) mean(miR_scRNA_filter[, i])))
mR_me <- unlist(lapply(seq(ncol(mR_scRNA_filter)), function(i) mean(mR_scRNA_filter[, i])))
miR_scRNA <- miR_scRNA_filter[, which(miR_me > median(miR_me))]
mR_scRNA <- mR_scRNA_filter[, which(mR_me > median(mR_me))]
miR_scRNA <- SummarizedExperiment(assays=list(miR_scRNA = miR_scRNA))
mR_scRNA <- SummarizedExperiment(assays=list(mR_scRNA =mR_scRNA))

## The prior information of miRNA targets is from miRTarBase v9.0 (https://mirtarbase.cuhk.edu.cn/) and TarBase v8.0
## (https://dianalab.e-ce.uth.gr/html/diana/web/index.php?r=tarbasev8). The downloaded miRNA targets of miRTarBase v9.0 
## and TarBase v8.0 in human can be obtained at https://github.com/zhangjunpeng411/miRSM_2.0_case_study (miRTarBase_v9.0.csv for 
## miRNA targets in miRTarBase v9.0, and TarBase_v8.0.csv for mRNA targets in TarBase v8.0). The input prior information of 
## miRSM 2.0 should also be SummarizedExperiment objects.
miRTarBase_v9.0 <- read.csv("miRTarBase_v9.0.csv", header = TRUE)
TarBase_v8.0 <- read.csv("TarBase_v8.0.csv", header = TRUE)
miRTarget <- unique(rbind(miRTarBase_v9.0, TarBase_v8.0))
miRTarget <- SummarizedExperiment(assays=list(miRTarget=miRTarget))

############################################################################## Step 3: Identifying gene modules #######################################################################
## Identifying gene modules at multi-sample and single-sample levels with GFA
set.seed(123)
modulegenes_GFA_all <- module_GFA(mR_scRNA)
length(modulegenes_GFA_all) == 10
nsamples <- nrow(mR_scRNA)
modulegenes_GFA_exceptk <- lapply(seq(nsamples), function(i) module_GFA(mR_scRNA[-i, ]))
lapply(seq(modulegenes_GFA_exceptk), function(i) length(modulegenes_GFA_exceptk[[i]]) == 9)

########################################################################## Step 4: Inferring miRNA sponge modules #######################################################################
## Inferring miRNA sponge modules at multi-sample level with sensitivity distance correlation
miRSM_GFA_SDC_all <- miRSM(miR_scRNA, 
                           mR_scRNA, 
			   NULL, 
			   miRTarget,
                           modulegenes_GFA_all,
			   method = "SDC")

## Inferring miRNA sponge modules at single-sample level with statistical perturbation strategy
miRSM_GFA_SDC_exceptk <- lapply(seq(nsamples), function(i) 
                          miRSM(miR_scRNA[-i, ], 
                          mR_scRNA[-i, ], 
			  NULL, 
                          miRTarget, 
			  modulegenes_GFA_exceptk[[i]],			    
                          method = "SDC"))

Modulegenes_all <- miRSM_GFA_SDC_all[[2]]
Modulegenes_exceptk <- lapply(seq(nsamples), function(i) 
                        miRSM_GFA_SDC_exceptk[[i]][[2]])
 
miRSM_GFA_SS <- miRSM_SS(Modulegenes_all, Modulegenes_exceptk)
Modulegenes_k <- lapply(seq(nsamples), function(i) 
                        miRSM_GFA_SS[[i]])

# Visualization of the number of miRNA sponge modules specific to each half K562 cell 
install.packages("ggplot2")
library(ggplot2)
Sample1 <- rownames(miR_scRNA)
Number1 <- unlist(lapply(seq(miRSM_GFA_SS), function(i) length(miRSM_GFA_SS[[i]])))
df1 <- data.frame(Sample = Sample1, Number = Number1)
ggplot(data = df1, mapping = aes(x = Sample, y = Number, fill = Sample)) + 
       geom_bar(stat = 'identity') +
       geom_text(mapping = aes(label = Number), hjust = 1.2) +       
       theme(axis.text.x = element_text(angle = 90), legend.position='none') +
       coord_flip() +
       labs(y="Number of miRNA sponge modules")

################################################################### Step 5: Evaluating module group similarity ##############################################################################
## Calculating sample-sample similarity with miRNA sponge module group similarity
Module_sim <- matrix(, nrow = nsamples, ncol = nsamples)
for (i in seq(nsamples)){
  for (j in seq(nsamples)){    
    Module_sim[i, j] <- module_group_sim(Modulegenes_k[[i]], Modulegenes_k[[j]])
  }
}

# Heatmap of the similarity between 19 half K562 cells in miRNA sponge modules
colnames(Module_sim) <- rownames(miR_scRNA)
rownames(Module_sim) <- rownames(miR_scRNA)
install.packages("corrplot")
library(corrplot)
corrplot(Module_sim, method = "shade", type = "upper", is.corr = FALSE, diag = FALSE, tl.col = "black")

######################################################### Step 6: Performaing Modular analysis of miRNA sponge modules #######################################################################
## Modular analysis of miRNA sponge modules at multi-sample and single-sample levels
# Functional analysis
miRSM_GFA_SDC_all_FEA <- module_FA(Modulegenes_all, Analysis.type = 'FEA')
miRSM_GFA_SDC_all_DEA <- module_FA(Modulegenes_all, Analysis.type = 'DEA')

# CML enrichment analysis. We get a list of CML (Chronic Myelogenous Leukemia) related genes (including miRNAs and mRNAs) from HMDD v4.017 and DisGeNET v7.018. The downloaded CML genes can be 
# obtained at https://github.com/zhangjunpeng411/miRSM_2.0_case_study (CML_genes.csv for CML genes in HMDD v4.0 and DisGeNET v7.0).
CML_genes <- read.csv("CML_genes.csv", header = FALSE)
CML_genes <- SummarizedExperiment(assays=list(CML_genes = CML_genes))
miRSM_GFA_SDC_all_pvalue <- module_CEA(mR_scRNA, NULL, CML_genes, Modulegenes_all)
miRSM_GFA_SS_pvalue <- lapply(seq(nsamples), function(i)
                              module_CEA(mR_scRNA, NULL, CML_genes, Modulegenes_k[[i]]))

# Validation analysis
Groundtruthcsv_high <- system.file("extdata", "Groundtruth_high.csv", package="miRSM")
Groundtruth_high <- read.csv(Groundtruthcsv_high, header = TRUE, sep = ",")
miRSM_GFA_SDC_all_Validate_high <- module_Validate(Modulegenes_all, Groundtruth_high)
miRSM_GFA_SS_Validate_high <- lapply(seq(nsamples), function(i)
                                       module_Validate(Modulegenes_k[[1]], Groundtruth_high))

# Co-expression analysis
miRSM_GFA_SDC_all_Coexpress <-  module_Coexpress(mR_scRNA, NULL, Modulegenes_all, 
                                                   resample = 100, method = "mean", 
						   test.method = "t.test")
miRSM_GFA_SS_Coexpress <-  lapply(seq(nsamples), function(i)
                                    module_Coexpress(mR_scRNA, NULL, Modulegenes_k[[i]], 
				    resample = 100, method = "mean", 
				    test.method = "t.test"))

# Visualization of the number of statistically co-expressed miRNA sponge modules specific to each half K562 cell 
Sample2 <- rownames(miR_scRNA)
Number2 <- unlist(lapply(seq(miRSM_GFA_SS_Coexpress), function(i) length(which(miRSM_GFA_SS_Coexpress[[i]][[3]] < 0.05))))

df2 <- data.frame(Sample = Sample2, Number = Number2)
ggplot(data = df2, mapping = aes(x = Sample, y = Number, fill = Sample)) + 
       geom_point(aes(x = Sample, y = Number, color= Sample), size = 5) +
       geom_bar(aes(fill = Sample), stat = "identity", width = 0.2) +
       geom_text(mapping = aes(label = Number), hjust = -0.4) +       
       theme(axis.text.x = element_text(angle = 90), legend.position='none') +
       coord_flip() +
       ylim(0, 18)  +
       labs(y="Number of statistically co-expressed modules")

# Distribution analysis of sharing miRNAs
miRSM_GFA_SDC_all_share_miRs <-  share_miRs(miR_scRNA, miRTarget, Modulegenes_all)
miRSM_GFA_SS_share_miRs <-  lapply(seq(nsamples), function(i)
                                     share_miRs(miR_scRNA[-i, ], miRTarget, Modulegenes_k[[i]]))
miRSM_GFA_SDC_all_miRdistribute <- module_miRdistribute(miRSM_GFA_SDC_all_share_miRs)
miRSM_GFA_SS_miRdistribute <- lapply(seq(nsamples), function(i)
                                       module_miRdistribute(miRSM_GFA_SS_share_miRs[[i]]))

# Visualization of the number of sharing miRNAs existing in the miRNA sponge modules specific to each half K562 cell 
Sample3 <- rownames(miR_scRNA)
Number3 <- unlist(lapply(seq(miRSM_GFA_SS_miRdistribute), function(i) length(which(as.numeric(miRSM_GFA_SS_miRdistribute[[i]][, 3]) == max(as.numeric(miRSM_GFA_SS_miRdistribute[[i]][, 3]))))))

df3 <- data.frame(Sample = Sample3, Number = Number3)
ggplot(data = df3, mapping = aes(x = Sample, y = Number, group = 1, color = Sample)) + 
       #geom_point(aes(x = Sample, y = Number, color= Sample), size = 5) +
       #geom_bar(aes(fill = Sample), stat = "identity", width = 0.2) +
       geom_line() +
       geom_point() +
       geom_text(mapping = aes(label = Number), hjust = -0.4) +       
       theme(axis.text.x = element_text(angle = 90), legend.position='none') +
       coord_flip() +
       ylim(0, 70)   +
       labs(y="Number of sharing miRNAs")

# Predicting miRNA-target interactions
miRSM_GFA_SDC_all_miRtarget <- module_miRtarget(miRSM_GFA_SDC_all_share_miRs, Modulegenes_all)
miRSM_GFA_SS_miRtarget <- lapply(seq(nsamples), function(i)
                                   module_miRtarget(miRSM_GFA_SS_share_miRs[[i]], Modulegenes_k[[i]]))

miRtarget_unique_number <- nrow(unique(do.call(rbind, miRSM_GFA_SDC_all_miRtarget)))
miRtarget_SS_unique_number <- unlist(lapply(seq(nsamples), function(i) nrow(unique(do.call(rbind, miRSM_GFA_SS_miRtarget[[i]])))))

# Visualization of the number of predicted miRNA-mRNA interactions specific to each half K562 cell 
Sample4 <- rownames(miR_scRNA)
df4 <- data.frame(Sample = Sample4, Number = miRtarget_SS_unique_number)
ggplot(data = df4, mapping = aes(x = Sample, y = Number, fill = Sample)) + 
       geom_bar(stat = 'identity') +
       geom_text(mapping = aes(label = Number), hjust = 1.2) +       
       theme(axis.text.x = element_text(angle = 90), legend.position='none') +
       labs(y = "Number of miRNA-mRNA interactions") +
       scale_fill_manual(values = heat.colors(19)) +
       coord_flip()

# Predicting miRNA sponge interactions
miRSM_GFA_SDC_all_miRsponge <- module_miRsponge(Modulegenes_all)
miRSM_GFA_SS_miRsponge <- lapply(seq(nsamples), function(i)
                                   module_miRsponge(Modulegenes_k[[i]]))

miRsponge_unique_number <- nrow(unique(do.call(rbind, miRSM_GFA_SDC_all_miRsponge)))
miRsponge_SS_unique_number <- unlist(lapply(seq(nsamples), function(i) nrow(unique(do.call(rbind, miRSM_GFA_SS_miRsponge[[i]])))))

# Visualization of the number of predicted miRNA sponge interactions specific to each half K562 cell 
Sample5 <- rownames(miR_scRNA)
df5 <- data.frame(Sample = Sample5, Number = miRsponge_SS_unique_number)
ggplot(data = df5, mapping = aes(x = Sample, y = Number, fill = Sample)) + 
       geom_bar(stat = 'identity') +
       geom_text(mapping = aes(label = Number), hjust = 1.2) +       
       theme(axis.text.x = element_text(angle = 90), legend.position='none') +
       labs(y = "Number of miRNA sponge interactions") +
       scale_fill_manual(values = terrain.colors(19)) +
       coord_flip()

save.image("K562_scRNA.RData")
