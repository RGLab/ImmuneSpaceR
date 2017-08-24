
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")
library(Biobase)

# CreateConnection ------------------------------------------
sdy <- CreateConnection("SDY269", onTest = TRUE)
print(sdy)
# Helper Function -------------------------------------------

test_EM <- function(EM, summary){
  expect_is(EM, "ExpressionSet")
  expect_gt(ncol(exprs(EM)), 0)
  expect_gt(nrow(exprs(EM)), 0)
  if(summary == T){
    # In summary, no gene is NA
    expect_false(any(is.na(fData(EM)$gene_symbol)))
  }
}

# Main Tests ------------------------------------------------

context("getGEMatrix")

test_that("gets TIV_2008 eSet non-summary", {
  EM <- sdy$getGEMatrix("TIV_2008", outputType = "normalized")
  test_EM(EM, summary = F)
})

test_that("gets TIV_2008 eSet summary", {
  EM <- sdy$getGEMatrix("TIV_2008", outputType = "summary", annotation = "latest")
  test_EM(EM, summary = T) 
})

test_that("get_multiple matrices non-summary", {
  EM <- sdy$getGEMatrix(c("TIV_2008","LAIV_2008"), outputType = "normalized")
  test_EM(EM, summary = F)
})

test_that("get_multiple matrices summary", {
  EM <- sdy$getGEMatrix(c("TIV_2008","LAIV_2008"), outputType = "summary", annotation = "latest")
  test_EM(EM, summary = T)
})

test_that("get_multiple matrices summary without cache error", {
  EM <- sdy$getGEMatrix(c("TIV_2008","LAIV_2008"), outputType = "summary", annotation = "latest")
  test_EM(EM, summary = T)
})

# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}