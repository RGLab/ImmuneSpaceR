context("Data Integration")

# Data --------------------------------------------------
eset <- suppressMessages(SDY269$getGEMatrix())
pheno <- Biobase::pData(eset)
hai <- SDY269$getDataset("hai")


# Tests --------------------------------------------------------
test_that("participant_id in both eset and hai", {
  expect_true("participant_id" %in% colnames(hai) & "participant_id" %in% colnames(pheno))
})
