context("getParticipantData")

# Source depdencies -------------------------------------------
source("global_variable.R")
source("global_dependencies.R")
source("set_curlOptions.R")


# Connections --------------------------------------------------
# allSdy <- CreateConnection("", onTest = T) # for debugging
allSdy <- CreateConnection("")

# Pgrp Id - for "readonly" and ehenrich users -----------------
# is composed of SDY80, SDY180, SDY269, and SDY422 to incl. all datatypes
groupId <- "auto_test" # both travis and bioc on test / prod

# Helper Functions ---------------------------------------------
testPgrp <- function(dt, groupId){
  pgrp_T <- allSdy$getParticipantData(group = groupId, dataType = dt, original_view = TRUE, maxRows = 1)
  pgrp_F <- allSdy$getParticipantData(group = groupId, dataType = dt, original_view = FALSE, maxRows = 1)
  orig_T <- allSdy$getDataset(dt, original_view = TRUE, maxRows = 1)
  orig_F <- allSdy$getDataset(dt, original_view = FALSE, maxRows = 1)

  res <- list()
  res$view_T <- all.equal(colnames(pgrp_T), colnames(orig_T))
  res$view_F <- all.equal(colnames(pgrp_F), colnames(orig_F))

  return(res)
}

# Tests --------------------------------------------------------
if (identical(Sys.getenv("ISR_login"), "readonly@rglab.org")) {
  test_that("", {
    groups <- try(con$listParticipantGroups())

    expect_is(groups, "data.table")
    expect_gt(nrow(groups), 0)
    expect_named(groups, c("group_id", "group_name", "created", "subjects", "studies"))
    expect_true(groupId %in% groups$group_name)
  })

  test_that("Pdata neut_ab_titer", {
    res <- testPgrp(dt = "neut_ab_titer", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata fcs_sample_files", {
    res <- testPgrp(dt = "fcs_sample_files", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata fcs_analyzed_result", {
    res <- testPgrp(dt = "fcs_analyzed_result", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata demographics", {
    res <- testPgrp(dt = "demographics", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata mbaa", {
    res <- testPgrp(dt = "mbaa", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata pcr", {
    res <- testPgrp(dt = "pcr", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata fcs_control_files", {
    res <- testPgrp(dt = "fcs_control_files", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata elisa", {
    res <- testPgrp(dt = "elisa", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata gene_expression_files", {
    res <- testPgrp(dt = "gene_expression_files", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata hai", {
    res <- testPgrp(dt = "hai", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata hla_typing", {
    res <- testPgrp(dt = "hla_typing", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata elispot", {
    res <- testPgrp(dt = "elispot", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })

  test_that("Pdata cohort_membership", {
    res <- testPgrp(dt = "cohort_membership", groupId = groupId )
    expect_true( res$view_T == "TRUE")
    expect_true( res$view_F == "TRUE")
  })
}


# cleanup ------------------------------------------------------
if(exists("netrc_file")){
  file.remove(netrc_file)
}


