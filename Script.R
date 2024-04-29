## Load required dataset (https://github.com/zhangjunpeng411/miRSM2.0_case_study) and R packages
load("K562_scRNA.RData")
library(miRSM)

## Identifying gene modules at multi-sample and single-sample levels with GFA
nsamples <- nrow(mR_scRNA)
modulegenes_GFA_all <- module_GFA(mR_scRNA)
modulegenes_GFA_exceptk <- lapply(seq(nsamples), function(i) module_GFA(mR_scRNA[-i, ]))

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

## Modular analysis of miRNA sponge modules at multi-sample and single-sample levels
# Functional analysis
miRSM_GFA_SDC_all_FEA <- module_FA(Modulegenes_all, Analysis.type = 'FEA')
miRSM_GFA_SDC_all_DEA <- module_FA(Modulegenes_all, Analysis.type = 'DEA')

# CML enrichment analysis
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

# Distribution analysis of sharing miRNAs
miRSM_GFA_SDC_all_share_miRs <-  share_miRs(miR_scRNA, miRTarget, Modulegenes_all)
miRSM_GFA_SS_share_miRs <-  lapply(seq(nsamples), function(i)
                                     share_miRs(miR_scRNA[-i, ], miRTarget, Modulegenes_k[[i]]))
miRSM_GFA_SDC_all_miRdistribute <- module_miRdistribute(miRSM_GFA_SDC_all_share_miRs)
miRSM_GFA_SS_miRdistribute <- lapply(seq(nsamples), function(i)
                                       module_miRdistribute(miRSM_GFA_SS_share_miRs[[i]]))

# Predicting miRNA-target interactions
miRSM_GFA_SDC_all_miRtarget <- module_miRtarget(miRSM_GFA_SDC_all_share_miRs, Modulegenes_all)
miRSM_GFA_SS_miRtarget <- lapply(seq(nsamples), function(i)
                                   module_miRtarget(miRSM_GFA_SS_share_miRs[[i]], Modulegenes_k[[i]]))

# Predicting miRNA sponge interactions
miRSM_GFA_SDC_all_miRsponge <- module_miRsponge(Modulegenes_all)
miRSM_GFA_SS_miRsponge <- lapply(seq(nsamples), function(i)
                                   module_miRsponge(Modulegenes_k[[i]]))

# Calculating sample-sample similarity with miRNA sponge module group similarity
Module_sim <- matrix(, nrow = nsamples, ncol = nsamples)
for (i in seq(nsamples)){
  for (j in seq(nsamples)){    
    Module_sim[i, j] <- module_group_sim(Modulegenes_k[[i]], Modulegenes_k[[j]])
  }
}

# Identifying differential miRNA sponge modules
miRSM_diff <- matrix(, nrow = nsamples, ncol = nsamples)
for (i in seq(nsamples)){
    for (j in seq(nsamples)){    
        miRSM_diff[i, j] <- length(diff_module(Modulegenes_k[[i]], Modulegenes_k[[j]]))
    }
}
diag(miRSM_diff) <- 0

save.image("K562_scRNA.RData")
