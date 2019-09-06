#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------

#' @importFrom jsonlite fromJSON
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

    # parse JSON response via jsonlite's fromJSON parsing function
    parsed <- fromJSON(response, simplifyDataFrame = FALSE)

    # construct a data.table for each group
    groupsList <- lapply(parsed$groups, function(group) {
      data.table(
        group_id = group$id,
        group_name = group$label,
        created = as.Date(group$created),
        subjects = length(group$category$participantIds),
        studies = length(unique(gsub("SUB\\d+.", "", group$category$participantIds)))
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
      setorder(participantGroups, group_id)
      setkey(participantGroups, group_id)
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
    groupName <- private$.checkParticipantGroup(group)

    colFilter <- rbind(
      colFilter,
      makeFilter(
        c(paste0("ParticipantId/", groupName), "EQUAL", groupName)
      )
    )

    return(self$getDataset(
      dataType,
      original_view = original_view,
      reload = reload,
      colFilter = colFilter,
      transformMethod = transformMethod,
      ...
    ))
  }
)

ISCon$set(
  which = "public",
  name = "listParticipantGEMatrices",
  value = function(group, verbose = FALSE) {
    private$.assertAllStudyConnection()
    groupName <- private$.checkParticipantGroup(group)

    participantIds <- private$.getParticipantIdsFromGroup(groupName)
    matrices <- self$listGEMatrices(verbose = verbose, participantIds = participantIds)

    return(matrices)
  }
)


ISCon$set(
  which = "public",
  name = "getParticipantGEMatrix",
  value = function(group, outputType = "summary", annotation = "latest", reload = FALSE) {
    private$.assertAllStudyConnection()
    groupName <- private$.checkParticipantGroup(group)

    ids <- private$.getParticipantIdsFromGroup(groupName)
    matNames <- self$listParticipantGEMatrices(groupName)$name

    message(paste0(length(matNames), " matrices found for ", groupName))

    mat <- self$getGEMatrix(matNames,
      outputType = outputType,
      annotation = annotation,
      reload = reload
    )

    return(mat[, mat$participant_id %in% ids])
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

ISCon$set(
  which = "private",
  name = ".checkParticipantGroup",
  value = function(group) {

    # Must rerun this to ensure valid groups are only for that user and are not
    # changed within cache.
    validGroups <- self$listParticipantGroups()

    if (is.numeric(group)) {
      col <- "group_id"
      groupName <- validGroups$group_name[validGroups$group_id == group]
    } else if (is.character(group)) {
      col <- "group_name"
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
    return(groupName)
  }
)

# Return a vector of participantIDs for a group
ISCon$set(
  which = "private",
  name = ".getParticipantIdsFromGroup",
  value = function(group) {
    private$.assertAllStudyConnection()

    # ---- Get groups and associated participantIDs ----
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

    # parse JSON response via jsonlite's fromJSON parsing function
    parsed <- jsonlite::fromJSON(response, simplifyDataFrame = FALSE)

    # Transform parsed json into a data.table with a row for each group
    # and a column containing a vector of relevant subjectids
    groupsList <- lapply(parsed$groups, function(group) {
      data.table(
        group_id = group$id,
        group_name = group$label,
        subjects = list(group$category$participantIds)
      )
    })
    validGroups <- rbindlist(groupsList)


    # ---- Get groupName ----

    if (is.numeric(group)) {
      col <- "group_id"
      groupName <- validGroups$group_name[validGroups$group_id == group]
    } else if (is.character(group)) {
      col <- "group_name"
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

    # ---- Return participantIds ----

    return(validGroups[group_name == groupName, subjects][[1]])
  }
)



# HELPER -----------------------------------------------------------------------
