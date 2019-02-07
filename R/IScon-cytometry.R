#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------
ISCon$set(
  which = "public",
  name = "listWorkspaces",
  value = function() {
    study <- ifelse(self$study == "Studies", "", self$study)
    data.table(
      suppressWarnings(labkey.selectRows(
        baseUrl = self$config$labkey.url.base,
        folderPath = "/Studies",
        schemaName = "immport",
        queryName = "ds_workspace",
        colNameOpt = "fieldname",
        parameters = paste0("$STUDY=", study)
      ))
    )
  }
)


ISCon$set(
  which = "public",
  name = "listGatingSets",
  value = function() {
    study <- ifelse(self$study == "Studies", "", self$study)
    data.table(
      suppressWarnings(labkey.selectRows(
        baseUrl = self$config$labkey.url.base,
        folderPath = self$config$labkey.url.path,
        schemaName = "assay.General.gatingset",
        queryName = "Data",
        colNameOpt = "fieldname",
        containerFilter = "CurrentAndSubfolders"
      ))
    )
  }
)


ISCon$set(
  which = "public",
  name = "summarizeCytometryData",
  value = function() {
    n_sample <- suppressWarnings(labkey.executeSql(
      folderPath = self$config$labkey.url.path,
      schemaName = "study",
      sql = "SELECT COUNT(DISTINCT file_info_name) as n, file_info_purpose FROM fcs_sample_files GROUP BY file_info_purpose;",
      colNameOpt = "fieldname"
    ))
    n_control <- suppressWarnings(labkey.executeSql(
      folderPath = self$config$labkey.url.path,
      schemaName = "study",
      sql = "SELECT COUNT(DISTINCT control_file) as n FROM fcs_control_files;",
      colNameOpt = "fieldname"
    ))
    n_analyzed <- suppressWarnings(labkey.executeSql(
      folderPath = self$config$labkey.url.path,
      schemaName = "study",
      sql = "SELECT COUNT(DISTINCT ParticipantID) as participants,
                    COUNT(DISTINCT study_time_collected) as timepoints,
                    COUNT(DISTINCT arm_accession.name) as cohorts,
                    COUNT(DISTINCT population_name_reported) as names
             FROM fcs_analyzed_result;",
      colNameOpt = "fieldname"
    ))
    n_workspaces <- nrow(self$listWorkspaces())
    n_gatingsets <- nrow(self$listGatingSets())

    cat("\n<CytometryDataSummary>\n")
    cat("  -", sum(n_sample$n), "FCS sample files")
    if (nrow(n_sample) > 1) {
      cat(" (")
      cat(paste(n_sample$n, gsub(" result", "", n_sample$file_info_purpose)), sep = ", ")
      cat(")")
    }
    cat("\n")
    cat("  -", n_control$n, "FCS control files\n")
    cat("  -", n_workspaces, "workspace files\n")
    cat("  -", n_gatingsets, "gating sets\n")
    if (nrow(n_analyzed) > 0) {
      cat("  - flow cytometry analyzed results\n")
      cat("    -", n_analyzed$participants, "participants\n")
      cat("    -", n_analyzed$timepoints, "time points\n")
      cat("    -", n_analyzed$cohorts, "cohorts\n")
      cat("    -", n_analyzed$names, "reported population names\n")
    }
  }
)


#' @importFrom flowWorkspace load_gs
ISCon$set(
  which = "public",
  name = "loadGatingSet",
  value = function(gatingSet) {
    assertRstudio()

    gs_list <- self$listGatingSets()
    study <- gs_list[gatingSet == gating_set, study]
    if (gatingSet %in% gs_list$gating_set) {
      load_gs(
        paste0("/share/files/Studies/", study[1], "/@files/analysis/gating_set/", gsub("SDY\\d+_", "", gatingSet))
      )
    } else {
      stop(gatingSet, " is not a valid gating set name.")
    }
  }
)



# PRIVATE ----------------------------------------------------------------------
ISCon$set(
  which = "private",
  name = ".downloadCytometryData",
  value = function(study = self$study, outputDir = ".") {
    assertRstudio()

    if (identical(study, "Studies")) {
      warning("This is a cross study connection. Please specify the study argument.")
    } else {
      file.copy(
        paste0("/share/files/Studies/", study, "/@files/rawdata/", study, "_cytometry.tar.gz"),
        to = outputDir
      )
    }
  }
)



# HELPER -----------------------------------------------------------------------
isRstudioDocker <- function() {
  dir.exists("/share/files/Studies") &
    dir.exists("/home/rstudio/.rstudio") &
    !is.null(Rlabkey:::ifApiKey())
}

assertRstudio <- function() {
  if (!isRstudioDocker()) {
    stop("You are not in the ImmuneSpace RStudio Session.")
  }
}
