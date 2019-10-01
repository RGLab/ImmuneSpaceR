context("ISCon$listGEAnalysis()")

# Connections --------------------------------------------------
SDY34 <- CreateConnection("SDY34")


# Helper Functions ---------------------------------------------
try_lgea <- function(con) {
  res <- tryCatch(con$listGEAnalysis(),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
test_that("returns df of GE analysis for single study if present", {
  res <- try_lgea(SDY269)
  expect_true(dim(res)[1] > 0)
})

test_that("fails gracefully if GE analysis not present", {
  res <- try_lgea(SDY34)
  expect_true(res$message == "Study does not have Gene Expression Analyses.")
})

test_that("returns df of GE analysis for all studies", {
  res <- try_lgea(ALL)
  expect_true(dim(res)[1] > 0)
})
