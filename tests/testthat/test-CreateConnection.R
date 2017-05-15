
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Helper Functions ---------------------------------------------

try_con <- function(study){
  tryCatch(
    con <- CreateConnection(study),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
context("CreateConnection")

test_that("Study argument is not NULL", {
  res <- try_con(NULL)
  expect_true("study cannot be NULL" == res$message)
})

test_that("Study argument accepts only one study", {
  res <- try_con(c("SDY269","SDY180"))
  msg_list <- strsplit(res$message, split = " ")[[1]]
  expect_true("multiple" %in% msg_list)
})

test_that("all studies can loaded with empty string", {
  res <- try_con("")
  expect_true(res$study == "Studies")
})

test_that("Nonexistent study fails to load", {
  res <- try_con("SDY4000")
  msg_list <- strsplit(res$message, split = " ")[[1]]
  expect_true("HTTP" %in% msg_list)
})

test_that("Existing study can be loaded", {
  res <- try_con("SDY400")
  expect_true(res$study == "SDY400")
})


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


