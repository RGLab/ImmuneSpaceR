#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------
ISCon$set(
  which = "public",
  name = "listWorkspaces",
  value = function(reload = FALSE) {
    if (is.null(self$cache$cyto$listWorkspaces) || reload) {
      self$cache$cyto$listWorkspaces <- data.table(
        suppressWarnings(labkey.selectRows(
          baseUrl = self$config$labkey.url.base,
          folderPath = self$config$labkey.url.path,
          schemaName = "study",
          queryName = "flow_workspace",
          colNameOpt = "fieldname"
        ))
      )
    }

    self$cache$cyto$listWorkspaces
  }
)


ISCon$set(
  which = "public",
  name = "listGatingSets",
  value = function(reload = FALSE) {
    if (is.null(self$cache$cyto$listGatingSets) || reload) {
      self$cache$cyto$listGatingSets <- tryCatch(
        data.table(suppressWarnings(labkey.selectRows(
          baseUrl = self$config$labkey.url.base,
          folderPath = self$config$labkey.url.path,
          schemaName = "cytometry_processing",
          queryName = "SelectedData",
          colNameOpt = "fieldname"
        ))),
        error = function(e) data.table()
      )
    }

    self$cache$cyto$listGatingSets
  }
)


ISCon$set(
  which = "public",
  name = "summarizeCyto",
  value = function() {
    # retrieve stats
    nSamples <- suppressWarnings(labkey.executeSql(
      folderPath = self$config$labkey.url.path,
      schemaName = "study",
      sql = "SELECT COUNT(DISTINCT file_info_name) as n,
                    file_info_purpose
             FROM fcs_sample_files
             GROUP BY file_info_purpose;",
      colNameOpt = "fieldname"
    ))
    nControls <- suppressWarnings(labkey.executeSql(
      folderPath = self$config$labkey.url.path,
      schemaName = "study",
      sql = "SELECT COUNT(DISTINCT control_file) as n FROM fcs_control_files;",
      colNameOpt = "fieldname"
    ))
    nAnalyzed <- suppressWarnings(labkey.executeSql(
      folderPath = self$config$labkey.url.path,
      schemaName = "study",
      sql = "SELECT COUNT(DISTINCT ParticipantID) as participants,
                    COUNT(DISTINCT study_time_collected) as timepoints,
                    COUNT(DISTINCT arm_accession.name) as cohorts,
                    COUNT(DISTINCT population_name_reported) as names
             FROM fcs_analyzed_result;",
      colNameOpt = "fieldname"
    ))
    nWorkspaces <- nrow(self$listWorkspaces())
    nGatingsets <- nrow(self$listGatingSets())

    # print summary
    cat("<CytometryDataSummary>\n")
    cat("  -", sum(nSamples$n), "FCS sample files")
    if (nrow(nSamples) > 1) {
      cat(" (")
      cat(
        paste(nSamples$n, gsub(" result", "", nSamples$file_info_purpose)),
        sep = ", "
      )
      cat(")")
    }
    cat("\n")
    cat("  -", nControls$n, "FCS control files\n")
    cat("  -", nWorkspaces, "workspace files\n")
    cat("  -", nGatingsets, "gating sets\n")
    if (nrow(nAnalyzed) > 0) {
      cat("  - flow cytometry analyzed results\n")
      cat("    -", nAnalyzed$participants, "participants\n")
      cat("    -", nAnalyzed$timepoints, "time points\n")
      cat("    -", nAnalyzed$cohorts, "cohorts\n")
      cat("    -", nAnalyzed$names, "reported population names\n")
    }

    invisible(NULL)
  }
)


#' @importFrom flowWorkspace pData sampleNames
#' @importFrom flowCore markernames
ISCon$set(
  which = "public",
  name = "summarizeGatingSet",
  value = function(gatingSet) {
    .assertRstudio()

    gsList <- self$listGatingSets()
    study <- gsList[gatingSet == gating_set, study][1]
    wsId <- gsList[gatingSet == gating_set, wsId][1]
    if (!gatingSet %in% gsList$gating_set) {
      stop("'", gatingSet, "' is not a valid gating set name.")
    }

    # load gating set (RDS only) and merge metadata
    gs <- readRDS(
      paste0(.buildGSPath(study, wsid, gatingSet), "/", gatingSet, ".rds")
    )
    gs <- private$.mergePD(gs, study)

    # summarize pData
    pd <- pData(gs)
    lapply(pd, function(x) {
      n <- length(unique(x))
      if (n > 50) {
        n
      } else if (is.numeric(x) && n > 10) {
        summary(x)
      } else {
        table(x, dnn = NULL)
      }
    })
  }
)


#' @importFrom flowWorkspace load_gs
ISCon$set(
  which = "public",
  name = "loadGatingSet",
  value = function(gatingSet) {
    .assertRstudio()

    gsList <- self$listGatingSets()
    study <- gsList[gatingSet == gating_set, study][1]
    wsId <- gsList[gatingSet == gating_set, wsid][1]
    if (!gatingSet %in% gsList$gating_set) {
      stop("'", gatingSet, "' is not a valid gating set name.")
    }

    # Load gating set and merge metadata
    gs <- load_gs(.buildGSPath(study, wsId, gatingSet))
    private$.mergePD(gs, study)
  }
)



# PRIVATE ----------------------------------------------------------------------
#' @importFrom flowWorkspace openWorkspace
ISCon$set(
  which = "private",
  name = ".openWorkspace",
  value = function(workspace) {
    .assertRstudio()
    if (grepl(".jo$", workspace)) {
      stop(".jo files cannot be processed currently.")
    }

    wsList <- self$listWorkspaces()
    study <- wsList[workspace == file_info_name, study_accession][1]
    if (!workspace %in% wsList$file_info_name) {
      stop("'", workspace, "' is not a valid workspace file name.")
    }

    # Open workspace
    openWorkspace(
      paste0(
        "/share/files/Studies/",
        study,
        "/@files/rawdata/flow_cytometry/",
        workspace
      )
    )
  }
)


ISCon$set(
  which = "private",
  name = ".downloadCytoData",
  value = function(study = self$study, outputDir = ".") {
    .assertRstudio()
    if (identical(study, "Studies")) {
      stop(
        "This is a cross study connection. ",
        "Please specify the study argument."
      )
    }

    file.copy(
      paste0(
        "/share/files/Studies/", study, "/@files/rawdata/",
        study, "_cytometry.tar.gz"
      ),
      to = outputDir
    )
  }
)


ISCon$set(
  which = "private",
  name = ".mergePD",
  value = function(gs, study) {
    names <- names(pData(gs))

    # merge panel
    if (!"panel" %in% names) {
      pData(gs)$panel <- unlist(lapply(gs@data@frames, function(x) {
        markers <- markernames(x)
        paste(markers[markers != "-"], collapse = "|")
      }))[rownames(pData(gs))]
    }

    # merge demograhpics
    if (!"participant_id" %in% names) {
      pData(gs)$sampleNames <- sampleNames(gs)
      pd <- pData(gs)
      fcsSampleFiles <- self$getDataset(
        "fcs_sample_files",
        colFilter = makeFilter(c("study_accession", "EQUAL", study))
      )
      pd <- merge(
        x = pd, y = fcsSampleFiles,
        by.x = "FILENAME", by.y = "file_info_name",
        all.x = TRUE
      )
      row.names(pd) <- pd$sampleNames
      pData(gs) <- pd
      pData(gs)$sampleNames <- NULL
    }

    gs
  }
)



# HELPER -----------------------------------------------------------------------
.isRstudioDocker <- function() {
  dir.exists("/share/files/Studies") & dir.exists("/home/rstudio/.rstudio")
}


.assertRstudio <- function() {
  if (!.isRstudioDocker()) {
    stop("You are not in the ImmuneSpace RStudio Session.")
  }
}


.buildGSPath <- function(study, wsId, gatingSet) {
  paste0(
    "/share/files/Studies/",
    study,
    "/@files/analysis/gating_set/",
    wsId,
    "/",
    gatingSet
  )
}
