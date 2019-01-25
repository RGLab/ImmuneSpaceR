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
})

test_that("listGatingSets", {
  gs <- con$listGatingSets()

  expect_is(gs, "data.table")
})

test_that("summarizeCytometryData", {
  expect_output(con$summarizeCytometryData(), "FCS sample files")
})

# test_that("loadGatingSet", {
#
# })


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}

