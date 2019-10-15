context("ISCon$getGEInputs()")


# Helper Functions ---------------------------------------------
try_gei <- function(con) {
  res <- tryCatch(
    con$getGEInputs(),
    warning = function(w) return(w),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
test_that("returns GE inputs df if study has inputs", {
  res <- try_gei(SDY269)
  expect_true((dim(res)[1] > 0) & (dim(res)[2] > 0))
})

test_that("returns error if study does not have inputs", {
  res <- try_gei(SDY87)
  expect_true(res$message == "Gene Expression Inputs not found for study.")
})
