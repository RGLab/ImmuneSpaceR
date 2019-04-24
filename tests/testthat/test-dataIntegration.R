context("Data Integration")

# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269")
mats <- sdy269$cache$GE_matrices$name
eset <- suppressMessages(sdy269$getGEMatrix(mats))
pheno <- Biobase::pData(eset)
hai <- sdy269$getDataset("hai")


# Tests --------------------------------------------------------
test_that("participant_id in both eset and hai", {
  expect_true("participant_id" %in% colnames(hai) & "participant_id" %in% colnames(pheno))
})
