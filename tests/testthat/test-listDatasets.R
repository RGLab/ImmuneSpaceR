context("ISCon$listDatasets()")


# Helper Functions ---------------------------------------------
try_ld <- function(con, ...) {
  tryCatch(
    capture.output(con$listDatasets(...)),
    warning = function(w) return(w),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
test_that("both datasets and EM returned without argument", {
  res <- try_ld(SDY269)
  expect_true(all(c("datasets", "Expression Matrices") %in% res))
})

test_that("strings other than datasets and EM return error", {
  res <- try_ld(SDY269, output = "my_fav_sdy")
  expect_true(res$message == "output other than datasets and expressions not allowed")
})

test_that("argument of datasets returns only datasets", {
  res <- try_ld(SDY269, output = "datasets")
  expect_true(!("Expression Matrices" %in% res))
})

test_that("argument of EM returns only EM", {
  res <- try_ld(SDY269, output = "expression")
  expect_true(!("Datasets" %in% res))
})

test_that("argument of EM returns error when no EM present", {
  res <- try_ld(SDY87, output = "expression")
  expect_true(res == "No Expression Matrices Available")
})
