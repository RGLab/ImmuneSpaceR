context("getGEMatrix")

# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")

# CreateConnection ------------------------------------------
sdy <- CreateConnection("")
IS1 <- CreateConnection("IS1")

# Helper Function -------------------------------------------
test_EM <- function(EM, summary){
  expect_is(EM, "ExpressionSet")
  expect_gt(ncol(Biobase::exprs(EM)), 0)
  expect_gt(nrow(Biobase::exprs(EM)), 0)
  if(summary == T){
    # In summary, no gene is NA
    expect_false(any(is.na(fData(EM)$gene_symbol)))
  }
}

test_PD <- function(PD, phenoCols){
  expect_is(PD, "data.frame")
  expect_gt(nrow(PD), 0)
  expect_true(all(phenoCols %in% colnames(PD)))
}

# Main Tests ------------------------------------------------
test_that("gets TIV_2008 eSet non-summary", {
  EM <- sdy$getGEMatrix("SDY269_PBMC_TIV_Geo", outputType = "normalized")
  test_EM(EM, summary = F)
})

# tests general raw output
test_that("gets TIV_2008 eSet raw", {
  EM <- sdy$getGEMatrix("SDY269_PBMC_TIV_Geo", outputType = "raw")
  test_EM(EM, summary = F)
})

# ensures that constructExpressionSet is working ok
test_that("gets TIV_young eSet raw", {
  EM <- sdy$getGEMatrix("SDY56_PBMC_Young", outputType = "raw")
  test_EM(EM, summary = F)
})

test_that("gets TIV_2008 eSet summary", {
  EM <- sdy$getGEMatrix("SDY269_PBMC_TIV_Geo", outputType = "summary", annotation = "latest")
  test_EM(EM, summary = T)
})

test_that("get_multiple matrices non-summary", {
  EM <- sdy$getGEMatrix(c("SDY269_PBMC_TIV_Geo","SDY269_PBMC_LAIV_Geo"), outputType = "normalized")
  test_EM(EM, summary = F)
})

test_that("get_multiple matrices summary", {
  EM <- sdy$getGEMatrix(c("SDY269_PBMC_TIV_Geo","SDY269_PBMC_LAIV_Geo"), outputType = "summary", annotation = "latest")
  test_EM(EM, summary = T)
})

test_that("get_multiple matrices summary without cache error", {
  EM <- sdy$getGEMatrix(c("SDY269_PBMC_TIV_Geo","SDY269_PBMC_LAIV_Geo"), outputType = "summary", annotation = "latest")
  test_EM(EM, summary = T)
})

test_that("get_multiple matrices summary with reload", {
  EM <- sdy$getGEMatrix(c("SDY269_PBMC_TIV_Geo","SDY269_PBMC_LAIV_Geo"),
                        outputType = "summary",
                        annotation = "latest",
                        reload = TRUE)
  test_EM(EM, summary = T)
})

test_that("get ImmSig Study - SDY212 with correct anno and summary", {
  mats <- IS1$cache$GE_matrices$name[ grep("SDY212", IS1$cache$GE_matrices$name) ]
  EM <- IS1$getGEMatrix(mats,
                        outputType = "raw",
                        annotation = "ImmSig")
  test_EM(EM, summary = F)
})

test_that("check pheno data.frame", {
  EM <- sdy$getGEMatrix("SDY269_PBMC_TIV_Geo", outputType = "normalized")
  PD <- Biobase::pData(EM)

  test_PD(PD, phenoCols)
})

# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}
