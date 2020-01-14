context("ISCon$getGEAnalysis()")

# Helper Functions ---------------------------------------------
try_ggea <- function(con, ...) {
  tryCatch(
    suppressMessages(con$getGEAnalysis(...)),
    warning = function(w) {
      return(w)
    },
    error = function(e) {
      return(e)
    }
  )
}


# Tests --------------------------------------------------------
test_that("returns df of GE analysis for single study if present", {
  res <- try_ggea(SDY269)
  expect_true(dim(res)[1] > 0)
})

test_that("returns df of GE analysis using cohort filter", {
  filt <- makeFilter(c("cohort", "equals", "TIV Group 2008_PBMC"))
  res <- try_ggea(SDY269, colFilter = filt)
  expect_true(dim(res)[1] > 0)
})

test_that("fails gracefully if GE analysis not present", {
  res <- try_ggea(SDY87)
  expect_true(res$message == "Gene Expression Analysis not found for study.")
})

test_that("returns df of GE analysis for all studies", {
  res <- try_ggea(ALL)
  expect_true(dim(res)[1] > 0)
})
