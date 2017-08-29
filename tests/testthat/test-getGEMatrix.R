
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")

# CreateConnection ------------------------------------------
sdy <- CreateConnection("SDY269")

# Helper Function -------------------------------------------

test_EM <- function(EM, summary){
  expect_is(EM, "ExpressionSet")
  expect_gt(ncol(EM), 0)
  expect_gt(nrow(EM), 0)
  if(summary == T){
    # In summary, no gene is NA
    expect_false(any(is.na(fData(EM)$gene_symbol)))
  }
}

# Main Tests ------------------------------------------------

context("getGEMatrix")

test_that("gets TIV_2008 eSet non-summary", {
  EM <- sdy$getGEMatrix("TIV_2008", summary = F)
  test_EM(EM, summary = F)
})

test_that("gets TIV_2008 eSet summary", {
  EM <- sdy$getGEMatrix("TIV_2008", summary = T, currAnno = T)
  test_EM(EM, summary = T) 
})

test_that("get_multiple matrices non-summary", {
  EM <- sdy$getGEMatrix(c("TIV_2008","LAIV_2008"), summary = F)
  test_EM(EM, summary = F)
})

test_that("get_multiple matrices summary", {
  EM <- sdy$getGEMatrix(c("TIV_2008","LAIV_2008"), summary = T, currAnno = T)
  test_EM(EM, summary = T)
})

# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}