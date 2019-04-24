context("CreateConnection()")

# Helper Functions ---------------------------------------------
try_con <- function(study) {
  tryCatch(
    con <- CreateConnection(study),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
test_that("Study argument is not NULL", {
  res <- try_con(NULL)
  expect_true("study cannot be NULL" == res$message)
})

test_that("Study argument accepts only one study", {
  res <- try_con(c("SDY269", "SDY180"))
  msg_list <- strsplit(res$message, split = " ")[[1]]
  expect_true("multiple" %in% msg_list)
})

test_that("all studies can loaded with empty string", {
  res <- try_con("")
  expect_true(res$study == "Studies")
})

test_that("Nonexistent study fails to load", {
  res <- try_con("SDY4000")
  errMsg <- "SDY4000 is not a valid study. \n Use `verbose = TRUE` to see list of valid studies."
  expect_true(errMsg == res$message)
})

test_that("Existing study can be loaded", {
  res <- try_con("SDY400")
  expect_true(res$study == "SDY400")
})

test_that("ImmuneSignatures study can be loaded", {
  res <- try_con("IS1")
  expect_equal(res$study, "IS1")
  expect_equal(res$config$labkey.url.path, "/HIPC/IS1")
})

test_that("Lyoplate study cannot be loaded", {
  res <- try_con("Lyoplate")
  errMsg <- "Lyoplate is not a valid study. \n Use `verbose = TRUE` to see list of valid studies."
  expect_true(errMsg == res$message)
})

test_that("Connection loads with cohort-type in GE_matrices", {
  con <- CreateConnection("SDY400")
  expect_true("cohort_type" %in% colnames(con$cache$GE_matrices))
})
