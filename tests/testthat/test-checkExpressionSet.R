
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)


# Helper Functions ---------------------------------------------
test_ES <- function(sdy, outputType, fail = FALSE) {
  out <- sdy269$.__enclos_env__$private$.checkExpressionSet(outputType = outputType)
  expect_that(nrow(out), equals(1))
  expect_that(ncol(out), equals(8))
  expect_true(out$outputType == outputType)
  expect_true(all(as.vector(t(out[1,-1]))))
  }


# Tests --------------------------------------------------------
context(".checkExpressionSet")

test_that("gets summary expression matrix and tests", {
  test_ES(sdy269, "summary")
})

test_that("gets raw expression matrix and tests", {
  test_ES(sdy269, "raw")
})

# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


