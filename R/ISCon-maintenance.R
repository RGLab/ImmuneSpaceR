#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------



# PRIVATE ----------------------------------------------------------------------

# Returns a list of data frames where TRUE in file_exists column marks files
# that are accessible. This function is used for administrative purposes to
# check that the raw files are properly loaded and accessible to the users.
#' @importFrom rjson fromJSON
#' @importFrom parallel mclapply detectCores
ISCon$set(
  which = "private",
  name = ".checkRawFiles",
  value = function(what = c("gene_expression_files",
                            "fcs_sample_files",
                            "fcs_control_files",
                            "protocols",
                            "ge_matrices")) {
    ## HELPERS
    ..checkLinks <- function (dataset, folder) {
      res <- data.frame(
        file_info_name = NULL,
        study_accession = NULL,
        file_link = NULL,
        file_exists = NULL,
        stringsAsFactors = FALSE
      )

      if (dataset %in% self$availableDatasets$Name) {
        temp <- self$getDataset(dataset, original_view = TRUE)

        if (dataset == "fcs_control_files") {
          temp <- temp[, file_info_name := control_file]
          temp <- temp[, c("pid", "sid") := data.table::tstrsplit(participant_id, "\\.")]
          temp <- temp[, study_accession := paste0("SDY", sid)]
        }

        temp <- temp[!is.na(file_info_name)]
        temp <- unique(temp[, list(study_accession, file_info_name)])

        file_link <- paste0(
          self$config$labkey.url.base,
          "/_webdav/Studies/",
          temp$study_accession,
          "/%40files/rawdata/",
          folder,
          "/",
          sapply(temp$file_info_name, URLencode)
        )

        studies <- unique(temp$study_accession)
        folder_link <- paste0(
          self$config$labkey.url.base,
          "/_webdav/Studies/",
          studies,
          "/%40files/rawdata/",
          folder,
          "?method=JSON"
        )

        file_list <- unlist(
          mclapply(
            folder_link,
            private$.listISFiles,
            mc.cores = detectCores()
          )
        )

        file_exists <- temp$file_info_name %in% file_list

        res <- data.frame(
          study = temp$study_accession,
          file_link = file_link,
          file_exists = file_exists,
          stringsAsFactors = FALSE
        )

        print(
          paste0(
            sum(res$file_exists),
            "/",
            nrow(res),
            " ",
            dataset,
            " with valid links."
          )
        )
      }

      res
    }


    ## MAIN
    ret <- list()
    what <- tolower(what)

    if ("gene_expression_files" %in% what) {
      ret$gene_expression_files <- ..checkLinks(
        "gene_expression_files",
        "gene_expression"
      )
    }

    if ("fcs_sample_files" %in% what) {
      ret$fcs_sample_files <- ..checkLinks(
        "fcs_sample_files",
        "flow_cytometry"
      )
    }

    if ("fcs_control_files" %in% what) {
      ret$fcs_control_files <- ..checkLinks(
        "fcs_control_files",
        "flow_cytometry"
      )
    }

    if ("protocols" %in% what) {
      if (private$.isProject()) {
        folders_list <- labkey.getFolders(
          baseUrl = self$config$labkey.url.base,
          folderPath = "/Studies/"
        )
        folders <- folders_list[, 1]
        folders <- folders[!folders %in% c("SDY_template","Studies")]
      } else {
        folders <- basename(self$config$labkey.url.path)
      }

      file_link <- paste0(
        self$config$labkey.url.base,
        "/_webdav/Studies/",
        folders,
        "/%40files/protocols/",
        folders,
        "_protocol.zip"
      )

      file_exists <- unlist(
        mclapply(
          file_link,
          private$.checkUrl,
          mc.cores = detectCores()
        )
      )

      print(
        paste0(
          sum(file_exists),
          "/",
          length(file_exists),
          " protocols with valid links."
        )
      )

      ret$protocols <- data.frame(
        study = folders,
        file_link = file_link,
        file_exists = file_exists,
        stringsAsFactors = FALSE
      )
    }

    if ("ge_matrices" %in% what) {
      mx <- .getLKtbl(
        con = self,
        schema = "assay.ExpressionMatrix.matrix",
        query = "Runs",
        colNameOpt = "rname"
      )

      mxLinks <- paste0(
        self$config$labkey.url.base,
        "/_webdav/Studies/",
        mx$folder_name,
        "/@files/analysis/exprs_matrices/",
        mx$name,
        ".tsv"
      )

      file_exists <- unlist(
        mclapply(
          mxLinks,
          private$.checkUrl,
          mc.cores = detectCores()
        )
      )

      print(
        paste0(
          sum(file_exists),
          "/",
          length(file_exists),
          " ge_matrices with valid links."
        )
      )

      ret$ge_matrices <- data.frame(
        file_link = mxLinks,
        file_exists = file_exists,
        stringsAsFactors = FALSE
      )

    } else {
      ret$ge_matrices <- data.frame(
        file_link = NULL,
        file_exists = NULL,
        stringsAsFactors = FALSE
      )
    }

    ret
  }
)


# Remove gene expression matrices that do not correspond to a run currently on
# PROD/TEST in the query assay.ExpressionMatrix.matrix.Runs.
# NOTE: Important to change the labkey.url.base variable depending on PROD/TEST
# to ensure you are not deleting any incorrectly.
#' @importFrom httr DELETE HEAD http_type http_error
ISCon$set(
  which = "private",
  name = ".removeOrphanGEMs",
  value = function() {
    # Double check we are working at project level and on correct server!
    if (!private$.isProject()) {
      stop("Can only be run at project level")
    }
    chkBase <- readline(
      prompt = paste0(
        "You are working on ",
        self$config$labkey.url.base,
        ". Continue? [T / f] "
      )
    )
    if (!(chkBase %in% c("T", "t", ""))) {
      return("Operation Aborted.")
    }

    # get runs listed in the proper table
    runs <- data.table(
      labkey.selectRows(
        baseUrl = self$config$labkey.url.base,
        folderPath = self$config$labkey.url.path,
        schemaName = "assay.ExpressionMatrix.matrix",
        queryName = "Runs",
        showHidden = TRUE
      )
    )

    noRunPres <- private$.getNoRunPres(runs)

    # if files to-be-rm, confirm rm ok, attempt delete, check results and report
    if (length(noRunPres) == 0) {
      return("No orphans found.")
    } else {
      print(noRunPres)
    }
    ok2rm <- readline(prompt = "Ok to remove all files listed above? [Y / n] ")
    if (toupper(ok2rm) == "Y" | ok2rm == "") {
      for (i in 1:length(noRunPres)) {
        private$.curlDelete(
          noRunPres[i],
          names(noRunPres)[i]
        )
      }

      noRunPresPost <- private$.getNoRunPres(runs)
      if (length(noRunpresPost) == 0) {
        return("No orphans found after removal. Success!")
      } else {
        print("Problems Occurred. Remaining Files with No Runs Present")
        print(noRunPresPost)
        print("************")
      }
    } else {
      print("Operation aborted.")
    }
  }
)


# This method allows admin to check which studies are compliant with the
# following modules/data files (GEF, RAW, GEO, GEM, DE, GEE, IRP, GSEA)
ISCon$set(
  which = "private",
  name = ".checkStudyCompliance",
  value = function(filterNonGE = TRUE,
                   showAllCols = FALSE,
                   onlyShowNonCompliant = TRUE,
                   verbose = FALSE) {
    ## HELPERS
    # Validation based on modules being turned on
    ..getModSdys <- function(name) {
      url <- paste0("https://www.immunespace.org/immport/studies/containersformodule.api?name=", name)
      res <- unlist(lapply(fromJSON(Rlabkey:::labkey.get(url))[[1]], function(x) {x[["name"]]}))
      res <- .spSort(res[res != "Studies" & res != "SDY_template"])
    }

    ## MAIN
    # For labkey.executeSql calls
    baseUrl <- self$config$labkey.url.base

    # res table
    if (self$study != "Studies") {
      Sdys <- c(self$study) # single study
    } else {
      Sdys <- private$.getSdyVec()
      Sdys <- .spSort(Sdys[Sdys != "angelica_test"])
    }

    mods <- c("GEF", "RAW", "GEO", "GEM", "DE_implied", "GSEA_implied", "GEE_implied", "IRP_implied")

    compDF <- data.frame(
      matrix(
        nrow = length(Sdys),
        ncol = length(mods)
      ),
      row.names = Sdys
    )
    colnames(compDF) <- mods

    # GEM
    withGems <- unique(self$cache$GE_matrices$folder)
    compDF$GEM <- rownames(compDF) %in% withGems

    # RAW
    file_list <- private$.getGEFileNames(TRUE)
    file_list <- file_list[file_list != "NULL"]
    compDF$RAW <- rownames(compDF) %in% names(file_list)[file_list == TRUE]

    # GEF
    gef <- self$getDataset("gene_expression_files")
    compDF$GEF <- rownames(compDF) %in% .subidsToSdy(gef$participant_id)

    # GEO
    geoGef <- gef[!is.na(gef$geo_accession), ]
    compDF$GEO <- rownames(compDF) %in% .subidsToSdy(geoGef$participant_id)

    # GEE - Studies with subjects having both GEM and response data for any timepoints
    # NOTE: when GEE is changed to allow NAb, can uncomment nab / respSubs lines
    hai <- self$getDataset("hai")
    # nab <- self$getDataset("neut_ab_titer")
    # resp <- union(resp$participant_id, nab$participant_id)
    inputSmpls <- labkey.selectRows(
      baseUrl = baseUrl,
      folderPath = "/Studies",
      schemaName = "study",
      queryName = "HM_InputSamplesQuery",
      containerFilter = "CurrentAndSubfolders",
      colNameOpt = "rname"
    )

    inputSmpls$study <- gsub("^SUB\\d{6}\\.", "SDY", inputSmpls$participantid)

    exprResp <- merge(
      inputSmpls,
      hai,
      by.x = c("participantid", "study_time_collected"),
      by.y = c("participant_id", "study_time_collected")
    )

    compDF$GEE_implied <- rownames(compDF) %in% .subidsToSdy(unique(exprResp$participantid))

    # GSEA - studies with subjects having GEAR results (i.e. compared multiple GEM timepoints)
    gear <- sapply(withGems, FUN = function(sdy){
      res <- tryCatch(labkey.executeSql(baseUrl = baseUrl,
                                        folderPath = paste0("/Studies/", sdy),
                                        schemaName = "gene_expression",
                                        sql = "SELECT COUNT (*) FROM gene_expression_analysis_results"),
                      error = function(e){ return( NA ) })
      output <- res[[1]] > 0 & !is.na(res)
    })
    compDF$GSEA_implied <- rownames(compDF) %in% names(gear)[gear == TRUE]

    # Do gea results exist? are they complete?

    existGEA <- labkey.selectRows(
      baseUrl = baseUrl,
      folderPath = "/Studies/",
      schemaName = "gene_expression",
      queryName = "gene_expression_analysis",
      colNameOpt = "rname",
      showHidden = TRUE)

    containers <- labkey.selectRows(
      baseUrl = baseUrl,
      folderPath = "/Studies/",
      schemaName = "core",
      queryName = "Containers",
      containerFilter = "CurrentAndSubfolders",
      showHidden = TRUE
    )

    existGEA$sdy <- containers$`Display Name`[ match(existGEA$container, containers$`Entity Id`)]

    gea <- lapply(withGems, FUN = function(sdy) {
      impliedGEA <- data.table(labkey.selectRows(
        baseUrl = baseUrl,
        folderPath = paste0("/Studies/", sdy),
        schemaName = "assay.expressionMatrix.matrix",
        queryName = "inputSamples",
        colNameOpt = "rname",
        showHidden = TRUE))

      # Does each time point have at least three participants -- only check those that do for existing GEA
      impliedGEA[, key := paste(biosample_arm_name,
                                biosample_study_time_collected,
                                biosample_study_time_collected_unit,
                                sep = " ")]
      # remove 0 day and negative time points
      impliedGEA <- impliedGEA[impliedGEA$biosample_study_time_collected > 0,]

      # count subjects per key
      impliedGEA <- impliedGEA[ , .(pids = length(unique(biosample_participantid))), by = key ]

      # find coefs with more than 3 pids
      impliedGEA <- impliedGEA[pids > 3,]

      if (nrow(impliedGEA) > 0) {
        currGEA <- existGEA[ existGEA$sdy == sdy, ]
        currGEA$key = paste(currGEA$arm_name, currGEA$coefficient)
        diff <- setdiff(impliedGEA$key, currGEA$key) # what is in implied and NOT in current
        missing_data <- if(length(diff) == 0 ){ "no diff" }else{ paste(diff, collapse = "; ") }
        res <- c( nrow(impliedGEA) > 0, length(diff) == 0, missing_data)


      } else {
        res <- c(FALSE, FALSE, NA)
      }

      return(res)
    })

    gea <- data.frame(do.call(rbind, gea), stringsAsFactors = F)
    colnames(gea) <- c("DGEA_implied", "DGEA_actual","DGEA_missing")
    rownames(gea) <- withGems
    compDF <- merge(compDF, gea, by = 0, all = T)
    compDF$DGEA_implied[ is.na(compDF$DGEA_implied) ] <- FALSE
    compDF$DGEA_actual[ is.na(compDF$DGEA_actual) ] <- FALSE
    compDF$DGEA_actual <- as.logical(compDF$DGEA_actual)
    compDF$DGEA_implied <- as.logical(compDF$DGEA_implied)
    row.names(compDF) <- compDF$Row.names
    compDF <- compDF[,-1]

    # IRP - Studies with subjects from multiple cohorts with GEM data at both target
    # timepoint and baseline + response data
    nab <- self$getDataset("neut_ab_titer")
    resp <- rbind(nab, hai, fill = TRUE) # nab has a col that hai does not
    suppressPackageStartupMessages(library(dplyr, quietly = T))
    # NOTE: At least SDY180 has overlapping study_time_collected for both hours and days
    # so it is important to group by study_time_collected_unit as well. This is reflected
    # in IRP_timepoints_hai/nab.sql
    geCohortSubs <- inputSmpls %>%
      group_by(study, study_time_collected, study_time_collected_unit) %>%
      # need multiple cohorts per timePoint (from IRP.js)
      filter((length(unique(cohort)) > 1 ) == TRUE) %>%
      ungroup() %>%
      group_by(study, cohort, study_time_collected_unit) %>%
      # need baseline + later timePoint (from IRP.Rmd)
      filter((length(unique(study_time_collected)) > 1 ) == TRUE)

    geRespSubs <- geCohortSubs[ geCohortSubs$participantid %in% unique(resp$participant_id), ]
    compDF$IRP_implied <- rownames(compDF) %in% .subidsToSdy(geRespSubs$participantid)
    studyTimepoints <- geRespSubs %>%
      group_by(study) %>%
      summarize(timepoints = paste(unique(study_time_collected), collapse = ","))
    compDF$IrpTimepoints <- studyTimepoints$timepoints[ match(rownames(compDF), studyTimepoints$study) ]

    # DE
    deSets <- c(
      "Neutralizing antibody titer",
      "Enzyme-linked immunosorbent assay (ELISA)",
      "Enzyme-Linked ImmunoSpot (ELISPOT)",
      "Hemagglutination inhibition (HAI)",
      "Polymerisation chain reaction (PCR)",
      "Flow cytometry analyzed results",
      "Multiplex bead array asssay"
    )
    compDF$DE_implied <- sapply(rownames(compDF), FUN = function(sdy) {
      res <- suppressWarnings(
        tryCatch(
          labkey.executeSql(
            baseUrl = baseUrl,
            folderPath = paste0("/Studies/", sdy),
            schemaName = "study",
            sql = "SELECT Label FROM ISC_datasets"
          ),
          error = function(e) {return( NA )}
        )
      )
      ret <- any(res[[1]] %in% deSets)
    })

    compDF$DE_actual <- rownames(compDF) %in% ..getModSdys("DataExplorer")
    compDF$GEE_actual <- rownames(compDF) %in% ..getModSdys("GeneExpressionExplorer")
    compDF$GSEA_actual <- rownames(compDF) %in% ..getModSdys("GeneSetEnrichmentAnalysis")
    compDF$IRP_actual <- rownames(compDF) %in% ..getModSdys("ImmuneResponsePredictor")

    colOrder <- c(
      "RAW",
      "GEF",
      "GEO",
      "GEM",
      "DE_implied",
      "DE_actual",
      "GEE_implied",
      "GEE_actual",
      "DGEA_actual",
      "DGEA_implied",
      "DGEA_missing_data",
      "IRP_implied",
      "IRP_actual",
      "IrpTimepoints",
      "GSEA_implied",
      "GSEA_actual")

    rowOrder <- Sdys[order(gsub("([A-Z]+)([0-9]+)", "\\1", Sdys),
                             as.numeric(gsub("([A-Z]+)([0-9]+)", "\\2", Sdys)))]

    compDF <- compDF[order(match(row.names(compDF), rowOrder)) , order(match(colnames(compDF), colOrder))]

    # Filter out studies that don't have GE since this is basis for everything
    if (filterNonGE == TRUE) {
      compDF <- compDF[compDF$GEO == TRUE | compDF$GEF == TRUE, ]
    }

    # Subset to only show problematic studies
    if (onlyShowNonCompliant == TRUE) {
      redux <- compDF[ , grep("implied|actual", colnames(compDF))]
      mod_sub <- c("DE", "GEE", "IRP", "GSEA", "DGEA")
      compliant <- lapply(mod_sub, FUN = function(mod) {
        idx <- grepl(mod, names(redux))
        sub <- redux[,idx]
        compliant <- sub[,1] == sub[,2]
      })
      compliant <- do.call(cbind, compliant)
      row.names(compliant) <- row.names(redux)
      idx <- which(apply(compliant, 1, all))
      compDF <- compDF[-(idx),]
    }

    # Defaults to showing only the actual module status and the difference with the implied
    if (showAllCols == FALSE) {
      compDF <- compDF[, grep("act|implied", colnames(compDF))]
    }
    if (verbose == TRUE) {
      message("NOTE: \n Return objects have an actual column that was generated by a call to the module url and an implied column that \n was created by looking at the filesystem. \n Please note that the DGEA_actual column is generated by looking at the file system and checking to ensure that all possible timepoints \n have been analyzed.")
    }


    compDF
  }
)


# Get vector of study folders
ISCon$set(
  which = "private",
  name = ".getSdyVec",
  value = function() {
    studies <- labkey.getFolders(
      baseUrl = self$config$labkey.url.base,
      folderPath = "/Studies/"
    )[, 1]
    studies <- studies[!studies %in% c("SDY_template","Studies")]

    studies
  }
)


# Get names of files in a single folder from webdav link
ISCon$set(
  which = "private",
  name = ".listISFiles",
  value = function(link) {
    opts <- self$config$curlOptions
    opts$options$netrc <- 1L

    response <- NULL

    res <- GET(url = link, config = opts)
    if (!http_error(res)) {
      response_json <- httr::content(res)
      response <- unlist(lapply(response_json$files, function(x) x$text))
    }

    response
  }
)


# Check if the url exists (is accessible)
ISCon$set(
  which = "private",
  name = ".checkUrl",
  value = function(url) {
    opts <- self$config$curlOptions
    opts$options$netrc <- 1L

    res <- HEAD(url, config = opts)

    if (http_error(res)) {
      ret <- FALSE
    } else {
      if (http_type(res) == "application/json") {
        res <- GET(url, config = opts)
        cont <- httr::content(res)
        ret <- is.null(cont$exception)
      } else {
        ret <- TRUE
      }
    }

    ret
  }
)


# Generate named list of files in either rawdata or analysis/exprs_matrices folders
ISCon$set(
  which = "private",
  name = ".getGEFileNames",
  value = function(rawdata) {
    studies <- private$.getSdyVec()

    # check webdav folder for presence of rawdata
    file_list <- lapply(studies, FUN = function(sdy) {
      suffix <- ifelse(rawdata,
                       "/%40files/rawdata/gene_expression?method=JSON",
                       "/%40files/analysis/exprs_matrices?method=JSON")

      dirLink <-  paste0(self$config$labkey.url.base,
                         "/_webdav/Studies/",
                         sdy,
                         suffix)
      files <- private$.listISFiles(dirLink)

      if (rawdata) {
        if (!is.null(files)) {
          files <- length(files) > 0
        }
      }

      return(files)
    })

    names(file_list) <- studies

    return(file_list)
  }
)


# Get no run pres (?)
ISCon$set(
  which = "private",
  name = ".getNoRunPres",
  value = function(runs) {
    emFls <- private$.getGEFileNames(FALSE)
    emFls <- emFls[emFls != "NULL"]
    tmpNms <- rep(x = names(emFls), times = lengths(emFls))
    emFls <- unlist(emFls)
    names(emFls) <- tmpNms

    # get names from emFls for comparison
    emNms <- sapply(emFls, FUN = function(x) {
      return(strsplit(x, "\\.tsv")[[1]][1])
    })

    # check runs for file names - assume if not in runs$Name then ok to delete!
    noRunPres <- emNms[!(emNms %in% runs$Name)]
    noRunPres <- noRunPres[!duplicated(noRunPres)]

    noRunPres
  }
)


# Delete the files
ISCon$set(
  which = "private",
  name = ".curlDelete",
  value = function(baseNm, sdy) {
    opts <- self$config$curlOptions
    opts$options$netrc <- 1L

    tsv <-  paste0(
      self$config$labkey.url.base,
      "/_webdav/Studies/",
      sdy,
      "/%40files/analysis/exprs_matrices/",
      baseNm,
      ".tsv"
    )
    smry <- paste0(tsv, ".summary")

    tsvRes <- DELETE(url = tsv, config = opts)
    smryRes <- DELETE(url = smry, config = opts)

    list(tsv = tsvRes, summary = smryRes)
  }
)



# HELPER -----------------------------------------------------------------------

# Sort studies by number
.spSort <- function(vec) {
  if (length(vec) > 0) {
    vec <- sort(as.numeric(gsub("SDY", "", vec)))
    vec <- paste0("SDY", as.character(vec))
  } else {
    vec <- NULL
  }

  vec
}


# Get SDY IDs from subids
.subidsToSdy <- function(subids) {
  sdys <- unique(gsub("^SUB.+", NA, unlist(strsplit(subids, split = "\\."))))
  sdys <- sdys[!is.na(sdys)]
  sdys <- paste0("SDY", sdys)

  sdys
}
