
# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
sdy67 <- suppressMessages(CreateConnection("SDY67"))
allsdy <- CreateConnection("")


# Helper Functions ---------------------------------------------
getFileList <- function(con){
  gef <- con$getDataset("gene_expression_files")
  nms <- unique(gef$name)
  if(length(nms) > 5){ nms <- nms[1:5] }
  return(nms)
}

try_ggef <- function(con){
  files <- getFileList(con)
  tryCatch(
    capture.output(con$getGEFiles(files = files[1:5]), destdir = destdir),
    warning = function(w) return(w),
    error = function(e) return(e)
  )
}


# Tests --------------------------------------------------------
context("getGEFiles")

destdir <- tempdir()
# Still a WIP - commenting out for Pull Request 4/28/17 EH
# test_that("returns df of GE analysis for single study if present", {
#   sink(file = file.path(destdir,"aux"))
#   res <- try_ggef(con)
#   sink(NULL)
#   expect_true( dim(res)[1] > 0 )
# })
# 
# test_that("returns df of GE analysis using cohort filter", {
#   filt <- makeFilter(c("cohort","equals","TIV Group 2008"))
#   res <- try_ggef(sdy269, colFilter = filt)
#   expect_true( dim(res)[1] > 0 )
# })
# 
# test_that("fails gracefully if GE analysis not present", {
#   res <- try_ggef(sdy67)
#   expect_true( res$message == "Study does not have Gene Expression Analyses" )
# })
# 
# test_that("returns df of GE analysis for all studies", {
#   res <- try_ggef(allsdy)
#   expect_true( dim(res)[1] > 0 )
# })


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


