#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------

# List participant groups
ISCon$set(
  which = "public",
  name = "listParticipantGroups",
  value = function() {
    private$.assertAllStudyConnection()

    participantGroupApi <- paste0(
      self$config$labkey.url.base,
      "/participant-group",
      "/Studies",
      "/browseParticipantGroups.api?",
      "distinctCatgories=false&",
      "type=participantGroup&",
      "includeUnassigned=false&",
      "includeParticipantIds=false"
    )
    # execute via Rlabkey's standard GET function
    response <- Rlabkey:::labkey.get(participantGroupApi)

    # parse JSON response via rjson's fromJSON parsing function
    parsed <- fromJSON(response)

    # construct a data.table for each group
    groupsList <- lapply(parsed$groups, function(group) {
      data.table(
        GroupId = group$id,
        GroupName = group$label,
        Subjects = length(group$category$participantIds)
      )
    })

    # merge the list to data.table
    participantGroups <- rbindlist(groupsList)

    if (nrow(participantGroups) == 0) {
      warning(
        "No participant groups found for the current user",
        immediate. = TRUE
      )
    } else {
      # set order by id
      setorder(participantGroups, GroupId)
      setkey(participantGroups, GroupId)
    }

    participantGroups
  }
)


# Retrieve a dataset by participant group
ISCon$set(
  which = "public",
  name = "getParticipantData",
  value = function(group, dataType, original_view = FALSE, ...) {
    private$.assertAllStudyConnection()

    if (!(dataType %in% self$availableDatasets$Name)) {
      warning(
        "'", dataType, "' is not a valid data type. ",
        "Valid data types are: ", paste(self$availableDatasets$Name, collapse = ", "),
        immediate. = TRUE
      )
      return(data.table())
    }

    # Must rerun this to ensure valid groups are only for that user and are not
    # changed within cache.
    validGroups <- self$listParticipantGroups()

    if (is.numeric(group)) {
      col <- "GroupId"
      groupId <- group
    } else if (is.character(group)) {
      col <- "GroupName"
      groupId <- validGroups$GroupId[validGroups$GroupName == group]
    } else {
      stop(
        "`group` should be a number or string. ",
        "Try again with valid `group`.",
        "\n Call `listParticipantGroups()` to see the available groups."
      )
    }

    if (!(group %in% validGroups[[col]])) {
      stop(
        "'", group, "' is not in the set of `", col,
        "` created by the current user",
        " on ", self$config$labkey.url.base,
        "\n Call `listParticipantGroups()` to see the available groups."
      )
    }

    dt <- dataType # for brevity

    # Get assay data + participantGroup + demographic data
    # Handle special cases
    if (dt == "demographics") {
      varSelect <- "cohort_membership.name AS Cohort "
      varJoin <- paste0(" JOIN cohort_membership ",
                        "ON ",
                        dt, ".ParticipantId = cohort_membership.ParticipantId ")
    } else {
      varSelect <- paste0("demo.age_reported, ",
                          "demo.race, ",
                          "demo.gender, ")
      varJoin <- paste0(" JOIN immport.arm_or_cohort ",
                        "ON ",
                        dt, ".arm_accession = arm_or_cohort.arm_accession ")
      if (dt %in% c("gene_expression_files", "cohort_membership", "fcs_sample_files")) {
        varSelect <- paste0( varSelect,
                             "demo.age_unit, ",
                             "demo.age_event, ",
                             "demo.ethnicity, ",
                             "demo.species, ",
                             "demo.phenotype ")
      } else {
        varSelect <- paste0(varSelect, "arm_or_cohort.name AS arm_name ")
      }
    }

    sqlAssay <- paste0("SELECT ",
                       dt, ".*, ",
                       "pmap.GroupId AS groupId, ",
                       varSelect,
                       "FROM ",
                       dt,
                       " JOIN ParticipantGroupMap AS pmap ",
                       "ON ", dt, ".ParticipantId = pmap.ParticipantId",
                       " JOIN Demographics AS demo ",
                       "ON ", dt, ".ParticipantId = demo.ParticipantId",
                       varJoin,
                       " WHERE pmap.GroupId = ", as.character(groupId))

    assayData <- labkey.executeSql(
      baseUrl = self$config$labkey.url.base,
      folderPath = self$config$labkey.url.path,
      schemaName = "study",
      sql = sqlAssay,
      colNameOpt = "fieldname",
      ...
    ) # allow for params to be passed, e.g. maxRows

    # Want to match getDataset() results in terms of colnames / order
    defaultCols <- colnames(
      self$getDataset(
        x = dt,
        original_view = original_view,
        maxRows = 1
      )
    )

    # Some names from assayData do not match default cols and need to changed manually.
    colnames(assayData)[grep("ParticipantId", colnames(assayData))] <- "participant_id"

    if (!original_view) {
      if (dt == "demographics") {
        changeCol <- "Cohort"
      } else if (dt == "cohort_membership") {
        changeCol <- "name"
      } else {
        changeCol <- "arm_name"
      }

      colnames(assayData)[grep(changeCol, colnames(assayData))] <- "cohort"
    } else {
      if (dt %in% c("gene_expression_files", "fcs_control_files")) {
        colnames(assayData)[grep("arm_name", colnames(assayData))] <- "cohort"
      }
    }

    filtData <- assayData[, colnames(assayData) %in% defaultCols]
    filtData <- filtData[, order(match(colnames(filtData), defaultCols))]

    filtData
  }
)



# PRIVATE ----------------------------------------------------------------------
# Check if all study connection
ISCon$set(
  which = "private",
  name = ".assertAllStudyConnection",
  value = function() {
    if (!identical(self$config$labkey.url.path, "/Studies/")) {
      stop(
        "This method only works with connection to all studies. ",
        'Create a connection to all studies by `con <- CreateConnection("")`'
      )
    }
  }
)



# HELPER -----------------------------------------------------------------------
