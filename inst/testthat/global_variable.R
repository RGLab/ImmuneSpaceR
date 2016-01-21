# 
#  Global variables
#
library(testthat)
library(ImmuneSpaceR)
library(Biobase)

demo <- c("age_reported", "gender", "race")
idCols <- c("participant_id")
common_cols <- c(demo, idCols)


haiCols     <- data.frame(name = c("value_reported", "virus_strain"),
                          type = c("numeric", "character"))
nabCols     <- data.frame(name = c("value_reported", "virus_strain"),
                          type = c("numeric", "character"))
hlaCols     <- data.frame(name = c("allele_1", "allele_2", "locus_name"),
                          type = c("character", "character", "character"))
elisaCols   <- data.frame(name = c("value_reported", "analyte"),
                          type = c("numeric", "character"))
elispotCols <- data.frame(name = c("spot_number_reported", "analyte"),
                          type = c("numeric", "character"))
pcrCols     <- data.frame(name = c("value_reported", "entrez_gene_id"),
                          type = c("numeric", "character"))
gefCols     <- data.frame(name = c("file_info_name", "cohort"),
                          type = c("character", "character"))
mbaaCols    <- data.frame(name = c("analyte_name", "concentration_value"))
farCols     <- data.frame(name = c(#"population_cell_number", # declared as VARCHAR(500)
  "population_definition_reported"),
                          type = c(#"numeric",
                            "character"))
fcsCols     <- data.frame(name = c("file_info_name"),
                          type = c("character"))
fccCols     <- data.frame(name = c("sample_file", "control_file"),
                          type = c("character", "character"))

test_dataset <- function(con, name, common_cols, specif_cols){
  x <- con$getDataset(name, reload = TRUE)
  # Dataset is of the right class and not empty
  expect_is(x, "data.table")
  expect_more_than(nrow(x), 0)
  # All required columns are here
  expect_true(all(common_cols %in% colnames(x)))
  expect_true(all(specif_cols$name %in% colnames(x)))
  # Important columns have no NAs
  expect_false(all(is.na(x[, specif_cols$name, with = FALSE])))
  # Important columns have the right class
  expect_true(all(sapply(x[, specif_cols$name,
                           with = FALSE], class) == specif_cols$type))
}


test_EM <- function(con, name, cohort = NULL){
  EM <- con$getGEMatrix(name)
  EMsum <- con$getGEMatrix(name, summary = TRUE)
  # EM have the right class and aren't empty
  expect_is(EM, "ExpressionSet")
  expect_is(EMsum, "ExpressionSet")
  expect_more_than(ncol(EM), 0)
  expect_more_than(nrow(EM), 0)
  expect_more_than(ncol(EMsum), 0)
  expect_more_than(nrow(EMsum), 0)  
  
  # In summary, no gene is NA
  expect_false(any(is.na(fData(EMsum)$gene_symbol)))
  
  # Download by cohort name
#  con$clear_cache()
#   con$getGEMatrix(cohort = cohort)
#   con$getGEMatrix(cohort = cohort, summary = TRUE)
}

