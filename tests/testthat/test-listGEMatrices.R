context("ISCon$listGEMatrices()")


# Variables --------------------------------------------------
em_list <- list()
em_list[["SDY269"]] <- c(
  "SDY269_PBMC_LAIV_Geo",
  "SDY269_PBMC_TIV_Geo"
)
em_list[["SDY404"]] <- c(
  "SDY404_PBMC_Young_Geo",
  "SDY404_PBMC_Older_Geo"
)


# Helper Functions ---------------------------------------------
chk_mats <- function(con, exp_mat_names) {
  expect_true(!is.null(con$config$labkey.url.base))
  mats <- con$listGEMatrices()
  expect_true(all(exp_mat_names %in% mats$name))
}


# Tests --------------------------------------------------------
test_that("gets correct matrices for SDY269", {
  chk_mats(SDY269, em_list[["SDY269"]])
})

test_that("gets correct matrices for SDY404", {
  chk_mats(SDY404, em_list[["SDY404"]])
})
