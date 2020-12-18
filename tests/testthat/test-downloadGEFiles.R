context("ISCon$downloadGEFiles()")

# Helper Functions ---------------------------------------------
getFileList <- function(con) {
  gef <- con$getDataset("gene_expression_files")
  nms <- unique(gef$file_info_name)
  nms <- nms[!is.na(nms)]
  if (length(nms) > 5) {
    nms <- nms[1:5]
  }
  return(nms)
}

try_ggef <- function(con) {
  files <- getFileList(con)
  dmp <- tempdir()
  res <- tryCatch(
    con$downloadGEFiles(files = files, destdir = dmp),
    warning = function(w) {
      return(w)
    },
    error = function(e) {
      return(e)
    }
  )
}

# Main Tests ------------------------------------------------
test_that("gets files when files are present", {
  res <- try_ggef(con = SDY269)
  expect_true(all(res))
})

# TODO: IS1 once virtual study import is changed to pull rawfiles
