# test_file("test_framework")
source("global_variable.R")

# Set curlOptions for download 
ISR_login <- Sys.getenv("ISR_login")
ISR_pwd <- Sys.getenv("ISR_pwd")
if(ISR_login != ""  &  ISR_pwd != ""){
  netrc_file <- tempfile("ImmuneSpaceR_tmp_netrc")
  netrc_string <- paste("machine www.immunespace.org login", ISR_login, "password", ISR_pwd)
  write(x = netrc_string, file = netrc_file)
  labkey.netrc.file <- netrc_file
}

# Connections
sdy269 <- CreateConnection("SDY269", verbose = TRUE)
sdy180 <- CreateConnection("SDY180", verbose = TRUE)
#sdy28 <- CreateConnection("SDY28", verbose = TRUE)

# datasets
test_that("get_hai", {
  test_dataset(sdy269, "hai", common_cols, specif_cols = haiCols)
})
test_that("get_elisa", {
  test_dataset(sdy269, "elisa", common_cols, specif_cols = elisaCols)
})
test_that("get_elispot", {
  test_dataset(sdy269, "elispot", common_cols, specif_cols = elispotCols)
})
test_that("get_pcr", {
  test_dataset(sdy269, "pcr", common_cols, specif_cols = pcrCols)
})
test_that("get_gene_expression_files", {
  test_dataset(sdy269, "gene_expression_files", common_cols, specif_cols = gefCols)
})
#test_that("fcs_analyzed_result", {
#  test_dataset(sdy269, "fcs_analyzed_result", common_cols, specif_cols = farCols)
#})
test_that("get_neut_ab_titer", {
  test_dataset(sdy180, "neut_ab_titer", common_cols, specif_cols = nabCols)
})
#test_that("get_fcs_sample_files", {
#  test_dataset(sdy180, "fcs_sample_files", common_cols, specif_cols = fcsCols)
#})
#test_that("get_fcs_control_files", {
#  test_dataset(sdy180, "fcs_control_files", common_cols, specif_cols = fccCols)
#})
#test_that("get_hla_typing", {
#  test_dataset(sdy28, "hla_typing", common_cols, specif_cols = hlaCols)
#})
# expression matrices
test_that("get_TIV2008", {
  test_EM(sdy269, "TIV_2008", cohort = "TIV Group 2008")
})
test_that("get_multiple", {
  test_EM(sdy269, c("TIV_2008", "LAIV_2008"))
})

# cleanup
if(exists("netrc_file")){
  file.remove(netrc_file)
}
