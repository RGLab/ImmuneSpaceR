
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")
library(Rlabkey)


# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
suppressMessages(sdy67 <- CreateConnection("SDY67"))
allsdy <- CreateConnection("")


# Helper Functions ---------------------------------------------
try_ggea <- function(con, ...){
  tryCatch(
    suppressMessages(con$getGEAnalysis(...)),
    warning = function(w) return(w),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
context("getGEAnalysis")

test_that("returns df of GE analysis for single study if present", {
  res <- try_ggea(sdy269)
  expect_true( dim(res)[1] > 0 )
})

test_that("returns df of GE analysis using cohort filter", {
  filt <- makeFilter(c("cohort","equals","TIV Group 2008"))
  res <- try_ggea(sdy269, colFilter = filt)
  expect_true( dim(res)[1] > 0 )
})

test_that("fails gracefully if GE analysis not present", {
  res <- try_ggea(sdy67)
  expect_true( res$message == "Study does not have Gene Expression Analyses" )
})

test_that("returns df of GE analysis for all studies", {
  res <- try_ggea(allsdy)
  expect_true( dim(res)[1] > 0 )
})


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


