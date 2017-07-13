library(testthat)
library(ImmuneSpaceR)

if (Sys.info["user"] != "biocbuild") test_check("ImmuneSpaceR")
