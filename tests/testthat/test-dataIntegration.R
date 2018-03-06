
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
library(Biobase)
sdy269 <- CreateConnection("SDY269")
mats <- sdy269$cache$GE_matrices$name
eset <- suppressMessages(sdy269$getGEMatrix(mats))
pheno <- pData(eset)
hai <- sdy269$getDataset("hai")

# Tests --------------------------------------------------------
context("Data Integration")

test_that("participant_id in both eset and hai", {
  expect_true( "participant_id" %in% colnames(hai) & "participant_id" %in% colnames(pheno) )
})

# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


