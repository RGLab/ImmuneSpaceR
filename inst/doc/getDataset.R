## ----knitr, echo = FALSE-------------------------------------------------
library(knitr)
opts_chunk$set(echo = TRUE)

## ----netrc_req, echo = FALSE---------------------------------------------
# This chunk is only useful for BioConductor checks and shouldn't affect any other setup
if(Sys.getenv("ISR_login") != "" & Sys.getenv("ISR_pwd") != ""){
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
all$listDatasets()

## ----getDataset----------------------------------------------------------
hai_269 <- sdy269$getDataset("hai")
hai_all <- all$getDataset("hai")

print(head(hai_269))

## ----makeFilters, cache = FALSE------------------------------------------
library(Rlabkey)
# Get participants under age of 30
young_filter <- makeFilter(c("age_reported", "LESS_THAN", 30))
# Get a specific list of two participants
pid_filter <- makeFilter(c("participantid", "IN", "SUB112841.269;SUB112834.269"))

## ----filters-------------------------------------------------------------
# HAI data for participants of study SDY269 under age of 30
hai_young <- sdy269$getDataset("hai", colFilter = young_filter)
# List of participants under age 30
demo_young <- all$getDataset("demographics", colFilter = young_filter)
# ELISPOT assay results for two participants
mbaa_pid2 <- all$getDataset("elispot", colFilter = pid_filter)

## ----views---------------------------------------------------------------
full_hai <- sdy269$getDataset("hai", original_view = TRUE)
print(colnames(full_hai))

## ----caching-dataset-----------------------------------------------------
pcr <- sdy269$getDataset("pcr")
names(sdy269$data_cache)

## ----caching-views-------------------------------------------------------
pcr_ori <- sdy269$getDataset("pcr", original_view = TRUE)
names(sdy269$data_cache)

## ----caching-reload------------------------------------------------------
hai_269 <- sdy269$getDataset("hai", reload = TRUE)

## ----caching-clear-------------------------------------------------------
sdy269$clear_cache()
names(sdy269$data_cache)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

