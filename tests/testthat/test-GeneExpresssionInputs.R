
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
sdy67 <- suppressMessages(CreateConnection("SDY67", verbose = TRUE))

# Helper Functions ---------------------------------------------
try_gei <- function(con){
  tryCatch(
    con$GeneExpressionInputs(),
    warning = function(w) return(w),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
context("GeneExpressionInputs")

test_that("returns GE inputs df if study has inputs", {
  res <- try_gei(sdy269)
  expect_true( (dim(res)[1] > 0) & (dim(res)[2] > 0) )
})

test_that("returns error if study does not have inputs", {
  res <- try_gei(sdy67)
  expect_true( res$message == "Gene Expression Inputs not found for study.")
})


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


