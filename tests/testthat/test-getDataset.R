context("ISCon$getDataset()")


# Test getDataset() ------------------------------------------------------------
test_getDataset("SDY269", "demographics")
test_getDataset("SDY269", "elisa")
test_getDataset("SDY269", "elispot")
test_getDataset("SDY269", "fcs_analyzed_result")
test_getDataset("SDY180", "fcs_control_files")
test_getDataset("SDY180", "fcs_sample_files")
test_getDataset("SDY269", "gene_expression_files")
test_getDataset("SDY269", "hai")
test_getDataset("SDY28", "hla_typing")
test_getDataset("SDY180", "mbaa")
test_getDataset("SDY180", "neut_ab_titer")
test_getDataset("SDY269", "pcr")
# test_getDataset("IS1", "hai")

test_that("returns empty data.table for invalid dataset name", {
  # check if it returned an empty data frame
  data <- suppressWarnings(try(CONNECTIONS$SDY269$getDataset("fakeData", reload = TRUE)))
  expect_is(data, "data.table")
  expect_equal(nrow(data), 0)

  warning <- capture_warning(CONNECTIONS$SDY269$getDataset("fakeData", reload = TRUE))
  expect_match(warning$message, "not a valid dataset")
})
