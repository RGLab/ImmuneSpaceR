
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


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


