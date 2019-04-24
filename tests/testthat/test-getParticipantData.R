context("ISCon$getParticipantData()")


# Test listParticipantGroups() -------------------------------------------------
test_that("listParticipantGroups() works", {
  skip_if_not(Sys.getenv("ISR_login") == "readonly@rglab.org")

  groups <- try(CONNECTIONS$ALL$listParticipantGroups())
  expect_is(groups, "data.table")
  expect_gt(nrow(groups), 0)
  expect_named(groups, c("group_id", "group_name", "created", "subjects", "studies"))
  expect_true("auto_test" %in% groups$group_name)
})


# Test getParticipantData() ----------------------------------------------------
test_getParticipantData("neut_ab_titer")
test_getParticipantData("fcs_sample_files")
test_getParticipantData("fcs_analyzed_result")
test_getParticipantData("demographics")
test_getParticipantData("mbaa")
test_getParticipantData("pcr")
test_getParticipantData("fcs_control_files")
test_getParticipantData("elisa")
test_getParticipantData("gene_expression_files")
test_getParticipantData("hai")
test_getParticipantData("hla_typing")
test_getParticipantData("elispot")
test_getParticipantData("cohort_membership")
