
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")

# CreateConnection ------------------------------------------
sdy <- CreateConnection("SDY269", verbose = TRUE)

# Helper Function -------------------------------------------

test_EM <- function(EM){
  expect_is(EM, "ExpressionSet")
  expect_gt(ncol(EM), 0)
  expect_gt(nrow(EM), 0)
}

# Main Tests ------------------------------------------------

context("getGEMatrix")

test_that("gets TIV_2008 eSet non-summary", {
  EM <- sdy$getGEMatrix("TIV_2008", summary = F)
  test_EM(EM)
})

test_that("gets TIV_2008 eSet summary", {
  EM <- sdy$getGEMatrix("TIV_2008", summary = T)
  test_EM(EM) 
})

test_that("get_multiple matrices non-summary", {
  EM <- sdy$getGEMatrix(c("TIV_2008","LAIV_2008"), summary = F)
  test_EM(EM)
})

test_that("get_multiple matrices summary", {
  EM <- sdy$getGEMatrix(c("TIV_2008","LAIV_2008"), summary = T)
  test_EM(EM)
})

# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}