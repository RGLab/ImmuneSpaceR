library(testthat)
library(ImmuneSpaceR)

labkey.url.base="https://www.immunespace.org" 
labkey.url.path="Studies/SDY269"
labkey.email.user="wjiang2@fhcrc.org"

# internal unit tests
#expectRes <<- readRDS("~/rglab/workspace/ImmuneSpaceR/misc/expectRes.rds")
#test_file("../inst/tests/test-IScon.R")
test_file("../inst/tests/test_ISR.R")
