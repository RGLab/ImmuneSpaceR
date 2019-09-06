context("ISCon$getParticipantData()")


# Test helpers -----------
# Test getParticipantIdsFromGroup()
test_that("getParticipantIdsFromGroup() works correctly", {
  skip_if_not(Sys.getenv("ISR_login") == "readonly@rglab.org")
  ids <- CONNECTIONS$ALL$.__enclos_env__$private$.getParticipantIdsFromGroup("travis_test")
  expect_is(ids, "character")
  expect_gt(length(ids), 0)

  expect_error(
    CONNECTIONS$SDY180$.__enclos_env__$private$.getParticipantIdsFromGroup("travis_test"),
    "This method only works with connection to all studies"
  )
})

test_that("checkParticipantGroup() works correctly", {
  expect_equal(CONNECTIONS$ALL$.__enclos_env__$private$.checkParticipantGroup("auto_test"), "auto_test")
  expect_equal(CONNECTIONS$ALL$.__enclos_env__$private$.checkParticipantGroup(163), "auto_test")
  expect_error(
    CONNECTIONS$ALL$.__enclos_env__$private$.checkParticipantGroup("fake"),
    "'fake' is not in the set of `group_name`"
  )
})

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



test_that("listParticipantGEMatrices() works correctly", {
  skip_if_not(Sys.getenv("ISR_login") == "readonly@rglab.org")
  matrices <- CONNECTIONS$ALL$listParticipantGEMatrices("auto_test")
  expect_is(matrices, "data.table")
  expect_gt(nrow(matrices), 0)
  expect_lt(nrow(matrices), 100)

  matrices <- CONNECTIONS$ALL$listParticipantGEMatrices("gem_test")
  expect_is(matrices, "data.table")
  expect_gt(nrow(matrices), 0)
  expect_lt(nrow(matrices), 100)

  expect_error(
    CONNECTIONS$ALL$listParticipantGEMatrices("fake"),
    "'fake' is not in the set of `group_name`"
  )
  expect_error(
    CONNECTIONS$SDY28$listParticipantGEMatrices("fake"),
    "This method only works with connection to all studies"
  )
})

# Test getParticipantGEMatrix()
test_that("getParticipantGEMatrix() works correctly", {
  skip_if_not(Sys.getenv("ISR_login") == "readonly@rglab.org")

  expect_message(
    {
      EM <- CONNECTIONS$ALL$getParticipantGEMatrix("gem_test")
    },
    "4 matrices found for gem_test"
  )
  expect_is(EM, "ExpressionSet")
  expect_gt(nrow(Biobase::exprs(EM)), 0)

  ids <- CONNECTIONS$ALL$.__enclos_env__$private$.getParticipantIdsFromGroup("gem_test")
  expect_lte(length(unique(EM$participant_id)), length(ids))

  # In summary, no gene is NA
  expect_false(any(is.na(Biobase::fData(EM)$gene_symbol)))
})
