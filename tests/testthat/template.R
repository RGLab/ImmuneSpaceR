
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
sdy180 <- CreateConnection("SDY180", verbose = TRUE)


# Helper Functions ---------------------------------------------
my_func <- function(foo){
  return("bar")
}


# Tests --------------------------------------------------------
context("my context")

test_that("my test", {
  expect_true( x == "TRUE")
})


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


