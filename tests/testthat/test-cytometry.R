context("cytometry")

# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
con <- CreateConnection("")


# Tests --------------------------------------------------------
test_that("listWorkspaces", {
  ws <- con$listWorkspaces()

  expect_is(ws, "data.table")
  expect_gt(nrow(ws), 0)
})

test_that("listGatingSets", {
  gs <- con$listGatingSets()

  expect_is(gs, "data.table")
  expect_gt(nrow(gs), 0)
})

test_that("summarizeCyto", {
  expect_output(con$summarizeCyto(), "FCS sample files")
})

# test_that("summarizeGatingSet", {
#
# })

# test_that("loadGatingSet", {
#
# })


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}

