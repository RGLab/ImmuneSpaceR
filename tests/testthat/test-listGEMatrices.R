context("ISCon$listGEMatrices()")

# Variables --------------------------------------------------
em_list <- list()
em_list[["SDY269"]] <- c("SDY269_PBMC_LAIV_Geo", "SDY269_PBMC_TIV_Geo")
em_list[["SDY404"]] <- c("SDY404_PBMC_Young_Geo",
                         "SDY404_PBMC_Older_Geo")


# Helper Functions ---------------------------------------------
chk_mats <- function(sdy, exp_mat_names){
  con <- CreateConnection(sdy, verbose = TRUE)
  expect_true(!is.null(con$config$labkey.url.base))

  mats <- con$listGEMatrices()
  expect_true(all(exp_mat_names %in% mats$name))
}


# Tests --------------------------------------------------------
test_that("gets correct matrices for SDY269", {
  sdy <- "SDY269"
  chk_mats(sdy, em_list[[sdy]])
})

test_that("gets correct matrices for SDY404", {
  sdy <- "SDY404"
  chk_mats(sdy, em_list[[sdy]])
})
