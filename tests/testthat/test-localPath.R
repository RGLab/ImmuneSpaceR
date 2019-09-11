context("private$.localStudyPath()")


# Helper Functions ---------------------------------------------
testPath <- function(labkey.url.base, link, pipelineRoot) {
  SDY269$config$labkey.url.base <- labkey.url.base
  .localStudyPath <- SDY269$.__enclos_env__$private$.localStudyPath
  res <- tryCatch(
    path <- .localStudyPath(link = link),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
test_that("individual dev machine", {
  res <- testPath(
    labkey.url.base = "http://local-server:8080/labkey/",
    link = "http://local-server:8080/labkey//_webdav/Studies/SDY269/@files/analysis/exprs_matrices/SDY269_PBMC_TIV_Geo.tsv.summary"
  )

  expect_true(res == "/share/files/Studies/SDY269/@files/analysis/exprs_matrices/SDY269_PBMC_TIV_Geo.tsv.summary")
})

test_that("server", {
  res <- testPath(
    labkey.url.base = "https://test.immunespace.org",
    link = "https://test.immunespace.org/_webdav/Studies/SDY269/@files/analysis/exprs_matrices/SDY269_PBMC_TIV_Geo.tsv.summary"
  )

  expect_true(res == "/share/files/Studies/SDY269/@files/analysis/exprs_matrices/SDY269_PBMC_TIV_Geo.tsv.summary")
})
