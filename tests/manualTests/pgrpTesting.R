# This must be done locally because getParticipantData() is user-specific and takes a long time
testPgrp <- function(dt, groupId){
  print(dt)
  pgrp_T <- con$getParticipantData(group = groupId, dataType = dt, original_view = T)
  pgrp_F <- con$getParticipantData(group = groupId, dataType = dt, original_view = F)
  orig_T <- con$getDataset(dt, original_view = T)
  orig_F <- con$getDataset(dt, original_view = F)
  
  res <- list()
  res$pgrp_T_cols <- colnames(pgrp_T)
  res$pgrp_F_cols <- colnames(pgrp_F)
  res$orig_T_cols <- colnames(orig_T)
  res$orig_F_cols <- colnames(orig_F)
  
  res$view_T <- all.equal(colnames(pgrp_T), colnames(orig_T))
  res$view_F <- all.equal(colnames(pgrp_F), colnames(orig_F))
  
  print(paste0("view T: ", res$view_T))
  print(paste0("view F: ", res$view_F))
  
  return(res)
}

library(ImmuneSpaceR)
con <- CreateConnection("")
dataTypes <- c(	"neut_ab_titer",
                "fcs_sample_files",
                "fcs_analyzed_result",
                "demographics",
                "mbaa",
                "pcr",
                "fcs_control_files",
                "elisa",
                "gene_expression_files",
                "hai",
                "hla_typing",
                "elispot",
                "cohort_membership")

chk <- lapply(dataTypes, testPgrp, groupId = 155)

