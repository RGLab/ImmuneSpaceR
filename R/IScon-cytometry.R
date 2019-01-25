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
        parameters = paste0("$STUDY=", study)
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
    n_population <- suppressWarnings(labkey.executeSql(
      folderPath = self$config$labkey.url.path,
      schemaName = "study",
      sql = "SELECT COUNT(DISTINCT population_name_reported) as n FROM fcs_analyzed_result;",
      colNameOpt = "fieldname"
    ))

    cat(sum(n_sample$n), "FCS sample files")
    if (nrow(n_sample) > 1) {
      cat(" (")
      cat(paste(n_sample$n, n_sample$file_info_purpose), sep = ", ")
      cat(")")
    }
    cat("\n")
    cat(n_control$n, "FCS control files\n")
    cat(n_population$n, "reported population names\n")
    cat(nrow(self$listWorkspaces()), "workspace files\n")
    cat(nrow(self$listGatingSets()), "gating sets\n")
  }
)


#' @importFrom flowWorkspace load_gs
ISCon$set(
  which = "public",
  name = "loadGatingSet",
  value = function(gatingSetName) {
    assertRstudio()

    gs_list <- self$listGatingSets()
    study <- gs_list[gatingSetName == gating_set_file, study]
    if (gatingSetName %in% gs_list$gating_set_file) {
      load_gs(
        paste0("/share/files/Studies/", study[1], "/@files/rawdata/gating_set/", gatingSetName)
      )
    } else {
      stop(gatingSetName, " is not a valid gating set name.")
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
  assert_that(
    isRstudioDocker(),
    msg = "You are not in the ImmuneSpace RStudio Session."
  )
}
