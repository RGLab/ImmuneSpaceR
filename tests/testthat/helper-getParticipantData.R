# Define helper test functions -------------------------------------------------
test_getParticipantData <- function(dataset, original_view = FALSE) {
  test_that(paste("auto_test", dataset), {
    skip_if_not(Sys.getenv("ISR_login") == "readonly@rglab.org")

    data <- try(ALL$getParticipantData("auto_test", dataset, original_view = original_view))
    columns <- colnames(ALL$getDataset(dataset, original_view = original_view, maxRows = 1))

    expect_is(data, "data.table")
    expect_gt(nrow(data), 0)
    expect_named(data, columns)
  })
}
