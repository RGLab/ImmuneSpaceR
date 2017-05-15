
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
sdy67 <- CreateConnection("SDY67")
allsdy <- CreateConnection("")


# Helper Functions ---------------------------------------------
try_lgea <- function(con){
  tryCatch(
    con$listGEAnalysis(),
    warning = function(w) return(w),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
context("listGEAnalysis")

test_that("returns df of GE analysis for single study if present", {
  res <- try_lgea(sdy269)
  expect_true( dim(res)[1] > 0 )
})

test_that("fails gracefully if GE analysis not present", {
  res <- try_lgea(sdy67)
  expect_true( res$message == "Study does not have Gene Expression Analyses" )
})

test_that("returns df of GE analysis for all studies", {
  res <- try_lgea(allsdy)
  expect_true( dim(res)[1] > 0 )
})


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


