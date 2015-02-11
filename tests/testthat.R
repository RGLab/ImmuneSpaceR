# TODO: EMs, flow datsets (waiting for immport update)
# Potential tests: getDataset: colFilter, colSelect
#                  test that caching is behaving as expected
library(ImmuneSpaceR)
library(Biobase)
library(testthat)

# 
#  Global variables
#
demo <- c("age_reported", "gender", "race")
idCols <- c("subject_accession")
common_cols <- c(demo, idCols)


haiCols     <- data.frame(name = c("value_reported", "virus_strain"),
                          type = c("numeric", "character"))
nabCols     <- data.frame(name = c("value_reported", "virus_strain"),
                          type = c("numeric", "character"))
hlaCols     <- data.frame(name = c("allele_1", "allele_2", "locus_name"),
                          type = c("character", "character", "character"))
#                           type = c("numeric", "numeric", "character"))
elisaCols   <- data.frame(name = c("value_reported", "analyte"),
                          type = c("numeric", "character"))
elispotCols <- data.frame(name = c("spot_number_reported", "analyte"),
                          type = c("numeric", "character"))
pcrCols     <- data.frame(name = c("value_reported", "entrez_gene_id"),
                          type = c("numeric", "character"))
gefCols     <- data.frame(name = c("file_info_name", "arm_name"),
                          type = c("character", "character"))
mbaaCols    <- data.frame(name = c("analyte_name", "concentration_value"))
farCols     <- data.frame(name = c(#"population_cell_number",
  "population_definition_reported"),
                          type = c(#"numeric",
                            "character"))


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
#   con$clear_cache()
#   con$getGEMatrix(cohort = cohort)
#   con$getGEMatrix(cohort = cohort, summary = TRUE)
}



# TESTS
sdy269 <- CreateConnection("SDY269")
sdy180 <- CreateConnection("SDY180")
sdy28 <- CreateConnection("SDY28")
# datasets
test_that("get_hai", {
  test_dataset(sdy269, "hai", common_cols, specif_cols = haiCols)
})
test_that("get_elisa", {
  test_dataset(sdy269, "elisa", common_cols, specif_cols = elisaCols)
})
test_that("get_elispot", {
  test_dataset(sdy269, "elispot", common_cols, specif_cols = elispotCols)
})
test_that("get_pcr", {
  test_dataset(sdy269, "pcr", common_cols, specif_cols = pcrCols)
})
test_that("get_gene_expression_files", {
  test_dataset(sdy269, "gene_expression_files", common_cols, specif_cols = gefCols)
})
test_that("fcs_analyzed_result", {
  test_dataset(sdy269, "fcs_analyzed_result", common_cols, specif_cols = farCols)
})
test_that("get_neut_ab_titer", {
  test_dataset(sdy180, "neut_ab_titer", common_cols, specif_cols = nabCols)
})
test_that("get_hla_typing", {
  test_dataset(sdy28, "hla_typing", common_cols, specif_cols = hlaCols)
})
# expression matrices
test_that("get_TIV2008", {
  test_EM(sdy269, "TIV_2008", cohort = "TIV Group 2008")
})
test_that("get_multiple", {
  test_EM(sdy269, c("TIV_2008", "LAIV_2008"))
})
