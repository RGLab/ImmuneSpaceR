
# Source depdencies --------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")

# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
sdy180 <- CreateConnection("SDY180", verbose = TRUE)

# Helper -------------------------------------------------------

test_dataset <- function(con, name, common_cols , specif_cols){
  x <- con$getDataset(name, reload = TRUE)
  
  # Dataset is of the right class and not empty
  expect_is(x, "data.table")
  expect_gt(nrow(x), 0)
  
  # All required columns are here
  expect_true(all(common_cols %in% colnames(x)))
  expect_true(all(specif_cols$name %in% colnames(x)))
  
  # Important columns have no NAs
  expect_false(all(is.na(x[, specif_cols$name, with = FALSE])))
  
  # Important columns have the right class
  expect_true(all(sapply(x[, specif_cols$name,
                           with = FALSE], class) == specif_cols$type))
}

# Dataset Tests ------------------------------------------------
context("getDatasets")

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

test_that("get_neut_ab_titer", {
  test_dataset(sdy180, "neut_ab_titer", common_cols, specif_cols = nabCols)
})

# test_that("fcs_analyzed_result", {
#  test_dataset(sdy269, "fcs_analyzed_result", common_cols, specif_cols = farCols)
# })

#test_that("get_fcs_sample_files", {
#  test_dataset(sdy180, "fcs_sample_files", common_cols, specif_cols = fcsCols)
#})
#test_that("get_fcs_control_files", {
#  test_dataset(sdy180, "fcs_control_files", common_cols, specif_cols = fccCols)
#})
#test_that("get_hla_typing", {
#  test_dataset(sdy28, "hla_typing", common_cols, specif_cols = hlaCols)
#})

# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


