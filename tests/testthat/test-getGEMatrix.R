context("ISCon$getGEMatrix()")

# CreateConnection ------------------------------------------
sdy <- CreateConnection("")
IS1 <- CreateConnection("IS1")


# Helper Function -------------------------------------------
phenoCols <- c(
  "participant_id", "study_time_collected", "study_time_collected_unit",
  "cohort", "cohort_type", "biosample_accession", "exposure_material_reported",
  "exposure_process_preferred"
)

test_EM <- function(EM, summary) {
  expect_is(EM, "ExpressionSet")
  expect_gt(ncol(Biobase::exprs(EM)), 0)
  expect_gt(nrow(Biobase::exprs(EM)), 0)
  if (summary == TRUE) {
    # In summary, no gene is NA
    expect_false(any(is.na(fData(EM)$gene_symbol)))
  }
}

test_PD <- function(PD, phenoCols) {
  expect_is(PD, "data.frame")
  expect_gt(nrow(PD), 0)
  expect_true(all(phenoCols %in% colnames(PD)))
}


# Main Tests ------------------------------------------------
test_that("gets TIV_2008 eSet non-summary", {
  EM <- sdy$getGEMatrix("SDY269_PBMC_TIV_Geo", outputType = "normalized")
  test_EM(EM, summary = FALSE)
})

# tests general raw output
test_that("gets TIV_2008 eSet raw", {
  EM <- sdy$getGEMatrix("SDY269_PBMC_TIV_Geo", outputType = "raw")
  test_EM(EM, summary = FALSE)
})

# ensures that constructExpressionSet is working ok
test_that("gets TIV_young eSet raw", {
  EM <- sdy$getGEMatrix("SDY56_PBMC_Young", outputType = "raw")
  test_EM(EM, summary = FALSE)
})

test_that("gets TIV_2008 eSet summary", {
  EM <- sdy$getGEMatrix("SDY269_PBMC_TIV_Geo", outputType = "summary", annotation = "latest")
  test_EM(EM, summary = TRUE)
})

test_that("get_multiple matrices non-summary", {
  EM <- sdy$getGEMatrix(c("SDY269_PBMC_TIV_Geo","SDY269_PBMC_LAIV_Geo"), outputType = "normalized")
  test_EM(EM, summary = FALSE)
})

test_that("get_multiple matrices summary", {
  EM <- sdy$getGEMatrix(c("SDY269_PBMC_TIV_Geo","SDY269_PBMC_LAIV_Geo"), outputType = "summary", annotation = "latest")
  test_EM(EM, summary = TRUE)
})

test_that("get_multiple matrices summary without cache error", {
  EM <- sdy$getGEMatrix(c("SDY269_PBMC_TIV_Geo","SDY269_PBMC_LAIV_Geo"), outputType = "summary", annotation = "latest")
  test_EM(EM, summary = TRUE)
})

test_that("get_multiple matrices summary with reload", {
  EM <- sdy$getGEMatrix(c("SDY269_PBMC_TIV_Geo","SDY269_PBMC_LAIV_Geo"),
                        outputType = "summary",
                        annotation = "latest",
                        reload = TRUE)
  test_EM(EM, summary = TRUE)
})

# Use specific tests here to ensure the IS1 report will load correctly
test_that("get ImmSig Study - SDY212 with correct anno and summary", {
  mats <- IS1$cache$GE_matrices$name[ grep("SDY212", IS1$cache$GE_matrices$name) ]
  EM <- IS1$getGEMatrix(mats,
                        outputType = "raw",
                        annotation = "ImmSig")
  test_EM(EM, summary = FALSE)
  expect_true("BS694717.1" %in% colnames(Biobase::exprs(EM)))
  expect_true("BS694717.1" %in% Biobase::sampleNames(EM))
  expect_true(all.equal(dim(Biobase::exprs(EM)), c(48771, 92)))
})

test_that("get ImmSig Study - SDY67 fixing 'X' for 'FeatureId'", {
  mats <- IS1$cache$GE_matrices$name[ grep("SDY67", IS1$cache$GE_matrices$name) ]
  EM <- IS1$getGEMatrix(mats,
                        outputType = "raw",
                        annotation = "ImmSig")
  test_EM(EM, summary = FALSE)
})

test_that("check pheno data.frame", {
  EM <- sdy$getGEMatrix("SDY269_PBMC_TIV_Geo", outputType = "normalized")
  PD <- Biobase::pData(EM)
  test_PD(PD, phenoCols)
})
