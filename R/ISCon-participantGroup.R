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
  value = function(group, dataType, original_view = FALSE, reload = FALSE, colFilter = NULL, transformMethod = "none", ...) {
    private$.assertAllStudyConnection()

    # Must rerun this to ensure valid groups are only for that user and are not
    # changed within cache.
    validGroups <- self$listParticipantGroups()

    if (is.numeric(group)) {
      col <- "GroupId"
      groupName <- validGroups$GroupName[validGroups$GroupId == group]
    } else if (is.character(group)) {
      col <- "GroupName"
      groupName <- group
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

    colFilter <- rbind(
      colFilter,
      makeFilter(
        c(paste0("ParticipantID/", groupName), "EQUAL", groupName)
      )
    )

    self$getDataset(
      dataType,
      original_view = original_view,
      reload = reload,
      colFilter = colFilter,
      transformMethod = transformMethod,
      ...
    )
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
