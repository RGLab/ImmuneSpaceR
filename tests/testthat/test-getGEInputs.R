context("getGEInputs")

# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
sdy87 <- suppressMessages(CreateConnection("SDY87", verbose = TRUE))

# Helper Functions ---------------------------------------------
try_gei <- function(con){
  tryCatch(
    con$getGEInputs(),
    warning = function(w) return(w),
    error = function(e) return(e)
  )
}

# Tests --------------------------------------------------------
test_that("returns GE inputs df if study has inputs", {
  res <- try_gei(sdy269)
  expect_true( (dim(res)[1] > 0) & (dim(res)[2] > 0) )
})

test_that("returns error if study does not have inputs", {
  res <- try_gei(sdy87)
  expect_true( res$message == "Gene Expression Inputs not found for study.")
})

# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


