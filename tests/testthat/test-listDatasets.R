
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
sdy67 <- suppressMessages(CreateConnection("SDY67", verbose = TRUE))


# Helper Functions ---------------------------------------------
try_ld <- function(con, ...){
  tryCatch(
    capture.output(con$listDatasets(...)),
    warning = function(w) return(w),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
context("listDatasets")

test_that("both datasets and EM returned without argument", {
  res <- try_ld(sdy269)
  expect_true( all(c("Datasets", "Expression Matrices") %in% res) )
})

test_that("strings other than datasets and EM return error", {
  res <- try_ld(sdy269, output = "my_fav_sdy")
  expect_true( res$message == "arguments other than 'datasets' and 'expression' not accepted")
})

test_that("argument of datasets returns only datasets", {
  res <- try_ld(sdy269, output = "datasets")
  expect_true( !("Expression Matrices" %in% res) )
})

test_that("argument of EM returns only EM", {
  res <- try_ld(sdy269, output = "expression")
  expect_true( !("Datasets" %in% res) )
})

test_that("argument of EM returns error when no EM present", {
  res <- try_ld(sdy67, output = "expression")
  expect_true( res$message == "No expression matrices available")
})


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


