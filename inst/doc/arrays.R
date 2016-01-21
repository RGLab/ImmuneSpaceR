## ----knitr, echo = FALSE-------------------------------------------------
library(knitr)
opts_chunk$set(echo = TRUE)
opts_chunk$set(cache = FALSE)

## ----netrc_req, echo = FALSE---------------------------------------------
# This chunk is only useful for BioConductor checks and shouldn't affect any other setup
if(Sys.getenv("ISR_login") != ""  & Sys.getenv("ISR_pwd") != ""){
  netrc_file <- tempfile("ImmuneSpaceR_tmp_netrc")
  netrc_string <- paste("machine www.immunespace.org login", ISR_login, "password", ISR_pwd)
  write(x = netrc_string, file = netrc_file)
  labkey.netrc.file <- netrc_file
}

## ----CreateConection, cache=FALSE----------------------------------------
library(ImmuneSpaceR)
sdy269 <- CreateConnection("SDY269")
all <- CreateConnection("")

## ----listDatasets--------------------------------------------------------
sdy269$listDatasets()

## ----listDatasets-which--------------------------------------------------
all$listDatasets(which = "expression")

## ----getGEMatrix---------------------------------------------------------
TIV_2008 <- sdy269$getGEMatrix("TIV_2008")
TIV_2011 <- all$getGEMatrix(x = "TIV_2011")

## ----ExpressionSet-------------------------------------------------------
TIV_2008

## ----getGEMatrix-cohorts-------------------------------------------------
LAIV_2008 <- sdy269$getGEMatrix(cohort = "LAIV group 2008")

## ----summary-------------------------------------------------------------
TIV_2008_sum <- sdy269$getGEMatrix("TIV_2008", summary = TRUE)

## ----summary-print-------------------------------------------------------
TIV_2008_sum

## ----multi---------------------------------------------------------------
# Within a study
em269 <- sdy269$getGEMatrix(c("TIV_2008", "LAIV_2008"))
# Combining accross studies
TIV_seasons <- all$getGEMatrix(c("TIV_2008", "TIV_2011"), summary = TRUE)

## ----caching-dataset-----------------------------------------------------
names(sdy269$data_cache)

## ----caching-reload------------------------------------------------------
TIV_2008 <- sdy269$getGEMatrix("TIV_2008", reload = TRUE)

## ----caching-clear-------------------------------------------------------
sdy269$clear_cache()

## ----sessionInfo---------------------------------------------------------
sessionInfo()

