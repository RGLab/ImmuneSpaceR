library(testthat)
library(ImmuneSpaceR)

if (any(file.exists("~/.netrc", "~/_netrc"))) test_check("ImmuneSpaceR")
