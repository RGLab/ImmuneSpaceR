context("ISCon$getGEMatrix()")

# Helper Functions -------------------------------------------
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
    expect_false(any(is.na(Biobase::fData(EM)$gene_symbol)))
  }
}

test_PD <- function(PD, phenoCols) {
  expect_is(PD, "data.frame")
  expect_gt(nrow(PD), 0)
  expect_true(all(phenoCols %in% colnames(PD)))
}


# Main Tests ------------------------------------------------
test_that("gets combined summary eset by default if at study level", {
  EM <- SDY269$getGEMatrix()
  test_EM(EM, summary = TRUE)
})

test_that("gets TIV_2008 eSet non-summary", {
  EM <- ALL$getGEMatrix(
    matrixName = "SDY269_PBMC_TIV_Geo",
    outputType = "normalized"
  )
  test_EM(EM, summary = FALSE)
})

# tests general raw output
test_that("gets TIV_2008 eSet raw", {
  EM <- ALL$getGEMatrix(
    matrixName = "SDY269_PBMC_TIV_Geo",
    outputType = "raw"
  )
  test_EM(EM, summary = FALSE)
})

# ensures that constructExpressionSet is working ok
test_that("gets TIV_young eSet raw", {
  EM <- ALL$getGEMatrix(
    matrixName = "SDY56_PBMC_Young",
    outputType = "raw"
  )
  test_EM(EM, summary = FALSE)
})

test_that("handles multiple samples per subject * timepoint combination", {
  warningMsg <- tryCatch(
    ALL$getGEMatrix(
      matrixName = "SDY1328_WholeBlood_HealthyAdults_geo",
      outputType = "raw",
      verbose = TRUE
    ),
    warning = function(w) return(w)
  )

  expect_true(grepl("Averaging the expression values", warningMsg))

  # should return matrix from cache
  EM <- ALL$getGEMatrix(
    matrixName = "SDY1328_WholeBlood_HealthyAdults_geo",
    outputType = "raw"
  )
  test_EM(EM, summary = FALSE)
  expect_equal(length(colnames(EM)), 164)
  # Note: orig matrix has 169 Biosamples, but 5 are assumed to be technical replicates
})

test_that("gets TIV_2008 eSet summary", {
  EM <- ALL$getGEMatrix(
    matrixName = "SDY269_PBMC_TIV_Geo",
    outputType = "summary",
    annotation = "latest"
  )
  test_EM(EM, summary = TRUE)
})

test_that("get_multiple matrices non-summary", {
  EM <- ALL$getGEMatrix(
    matrixName = c("SDY269_PBMC_TIV_Geo", "SDY269_PBMC_LAIV_Geo"),
    outputType = "normalized"
  )
  test_EM(EM, summary = FALSE)
})

test_that("get_multiple matrices summary", {
  EM <- ALL$getGEMatrix(
    matrixName = c("SDY269_PBMC_TIV_Geo", "SDY269_PBMC_LAIV_Geo"),
    outputType = "summary",
    annotation = "latest"
  )
  test_EM(EM, summary = TRUE)
})

test_that("get_multiple matrices summary without cache error", {
  EM <- ALL$getGEMatrix(
    matrixName = c("SDY269_PBMC_TIV_Geo", "SDY269_PBMC_LAIV_Geo"),
    outputType = "summary",
    annotation = "latest"
  )
  test_EM(EM, summary = TRUE)
})

test_that("get_multiple matrices summary with reload", {
  EM <- ALL$getGEMatrix(
    matrixName = c("SDY269_PBMC_TIV_Geo", "SDY269_PBMC_LAIV_Geo"),
    outputType = "summary",
    annotation = "latest",
    reload = TRUE
  )
  test_EM(EM, summary = TRUE)
})

test_that("get multiple matrices summary from different studies", {
  EM <- ALL$getGEMatrix(
    matrixName = c("SDY269_PBMC_TIV_Geo", "SDY180_WholeBlood_Grp2Pneunomax23_Geo"),
    outputType = "summary",
    annotation = "latest"
  )
  test_EM(EM, summary = TRUE)
})

test_that("gets SDY300 eset without error in .constructExpressionSet", {
  EM <- ALL$getGEMatrix(
    matrixName = "SDY300_otherCell_dcMonoFlu2011",
  )
  test_EM(EM, summary = TRUE)
})

test_that("loading from cache works correctly", {

  # Should load both matrices from cache
  expect_message(ALL$getGEMatrix(
    matrixName = c("SDY269_PBMC_TIV_Geo", "SDY180_WholeBlood_Grp2Pneunomax23_Geo"),
    outputType = "summary",
    annotation = "latest"
  ),
  "(cache)|(Combining ExpressionSets)",
  all = TRUE
  )

  # summary with a different annotation should download a new matrix
  expect_message(ALL$getGEMatrix(
    matrixName = "SDY180_WholeBlood_Grp2Pneunomax23_Geo",
    outputType = "summary",
    annotation = "default"
  ),
  "(Reading|Downloading)|(Constructing)",
  all = TRUE
  )

  # Should load eset from cache
  expect_message(
    ALL$getGEMatrix(
      matrixName = "SDY269_PBMC_TIV_Geo",
      outputType = "normalized",
      annotation = "latest"
    ),
    "Returning SDY269_PBMC_TIV_Geo_normalized_latest_eset from cache"
  )

  # Should load eset from cache without error if verbose uses
  expect_message(
    ALL$getGEMatrix(
      matrixName = "SDY269_PBMC_TIV_Geo",
      outputType = "normalized",
      annotation = "latest",
      verbose = TRUE
    ),
    "Returning SDY269_PBMC_TIV_Geo_normalized_latest_eset from cache"
  )

  # Should load matrix from cache and construct a new expressionset
  expect_message(ALL$getGEMatrix(
    matrixName = "SDY269_PBMC_TIV_Geo",
    outputType = "normalized",
    annotation = "default"
  ),
  "(Returning normalized matrix from cache)|(Downloading Features)|(Constructing ExpressionSet)",
  all = TRUE
  )

  # Should download a new matrix, load annotation from cache, and construct a new expressionset
  expect_message(ALL$getGEMatrix(
    matrixName = "SDY56_PBMC_Young",
    outputType = "normalized"
  ),
  "(Reading|Downloading)|(Returning latest annotation from cache)|(Constructing ExpressionSet)",
  all = TRUE
  )
})

# Use specific tests here to ensure the IS1 report will load correctly
test_that("get ImmSig Study - SDY212 with correct anno and summary", {
  mats <- IS1$cache$GE_matrices$name[ grep("SDY212", IS1$cache$GE_matrices$name) ]
  EM <- IS1$getGEMatrix(
    matrixName = mats,
    outputType = "raw",
    annotation = "ImmSig"
  )
  test_EM(EM, summary = FALSE)
  expect_true("BS694717.1" %in% colnames(Biobase::exprs(EM)))
  expect_true("BS694717.1" %in% Biobase::sampleNames(EM))
  expect_true(all.equal(dim(Biobase::exprs(EM)), c(48770, 92)))
})

test_that("get ImmSig Study - SDY67 fixing 'X' for 'FeatureId'", {
  mats <- IS1$cache$GE_matrices$name[ grep("SDY67", IS1$cache$GE_matrices$name) ]
  EM <- IS1$getGEMatrix(
    matrixName = mats,
    outputType = "raw",
    annotation = "ImmSig"
  )
  test_EM(EM, summary = FALSE)
})

test_that("check pheno data.frame", {
  EM <- ALL$getGEMatrix(
    matrixName = "SDY269_PBMC_TIV_Geo",
    outputType = "normalized"
  )
  PD <- Biobase::pData(EM)
  test_PD(PD, phenoCols)
})
