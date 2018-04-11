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


# This function outputs a list of lists. There are six possible sub-lists >
# 1. $gemAndRaw = study has gene expression matrix flat file and raw data
# 2. $gemNoRaw = study has gene expression matrix flat file but no raw data (unlikely)
# 3. $rawNoGem = study has raw data but no gene expression flat file, which is likely
#                in the case that no annotation pkg is availble in bioconductor or
#                it is RNAseq and had trouble being processed.
# 4. $gefNoGem = con$getDataset("gene_expression_files") reports raw data available, but
#                no gene expression matrix flat file has been generated. Similar to rawNoGem.
# 5. $gefNoRaw = con$getDataset() reports files being available, but no rawdata found. This
#                may be due to files being in GEO, but not having been downloaded to ImmPort
#                and ImmuneSpace.
# 6. $rawNoGef = This would be unexpected. Rawdata present on server, but no gene expression
#                files found by con$getDataset().
ISCon$set(
  which = "private",
  name = ".getSdysWithoutGEMs",
  value = function() {
    res <- list()

    # get list of matrices and determine which sdys they represent
    gems <- self$cache$GE_matrices
    withGems <- unique(gems$folder)

    file_list <- private$.getGEFileNames(TRUE)
    file_list <- file_list[file_list != "NULL"]
    emptyFolders <- names(file_list)[file_list == FALSE]
    withRawData <- names(file_list)[file_list == TRUE]

    # Compare lists
    res$gemAndRaw <- .spSort(intersect(withRawData, withGems))
    res$gemNoRaw <- .spSort(setdiff(withGems, withRawData))
    res$rawNoGem <- .spSort(setdiff(withRawData, withGems))

    # Check which studies without gems have gef in IS
    ge <- self$getDataset("gene_expression_files")
    geNms <- unique(ge$participant_id)
    gefSdys <- unique(sapply(geNms, FUN = function(x) {
      res <- strsplit(x, ".", fixed = TRUE)
      return(res[[1]][2])
    }))
    gefSdys <- paste0("SDY", gefSdys)

    res$gefNoGem <- .spSort(gefSdys[!(gefSdys %in% withGems)])
    res$gefNoRaw <- .spSort(setdiff(res$gefNoGem, res$rawNoGem))
    res$rawNoGef <- .spSort(setdiff(res$rawNoGem, res$gefNoGem))

    res
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
                   onlyShowNonCompliant = TRUE) {
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
      Sdys <- .spSort(private$.getSdyVec())
    }

    mods <- c("GEF", "RAW", "GEO", "GEM", "DE_implied", "GEE_implied", "IRP_implied")
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
      containerFilter = "CurrentAndSubfolders"
    )
    exprResp <- merge(
      inputSmpls,
      hai,
      by.x = c("Participant Id", "Study Time Collected"),
      by.y = c("participant_id", "study_time_collected")
    )
    compDF$GEE_implied <- rownames(compDF) %in% .subidsToSdy(unique(exprResp$`Participant Id`))

    # GSEA - studies with subjects having GEAR results (i.e. compared multiple GEM timepoints)
    # Do gea results exist? are they complete?

    gear <- sapply(withGems, FUN = function(sdy) {
      impliedGEA <- data.table(
        labkey.selectRows(
          baseUrl = baseUrl,
          folderPath = paste0("/Studies/", sdy),
          schemaName = "assay.ExpressionMatrix.matrix",
          queryName = "InputSamples",
          colNameOpt = "rname"))

      # Does each time point have at least three participants -- only check those that do for existing GEA
      if (nrow(impliedGEA) > 0) {
        q1 <- quote(biosample_arm_name)
        q2 <- quote(biosample_study_time_collected)
        q3 <- quote(biosample_study_time_collected_unit)
        impliedGEA[, key := paste(eval(q1), eval(q2), eval(q3), sep = " ")]
      }
      # count keys
      impliedGEA <- impliedGEA[, list(nr.appearances=length(key)),
                                  by = list(unique.values = key)]
      # remove zero day coefficient
      impliedGEA <- impliedGEA[!grepl("0 Days", impliedGEA$unique.values),]
      # only keep keys that have > 3 apperances (need at least three for GEA)
      impliedGEA <- impliedGEA[impliedGEA$nr.appearances > 3, ]

      # Find existing GEAR comparisons
      existGEA <- tryCatch(
        data.table(
          labkey.selectRows(
            baseUrl = baseUrl,
            folderPath = paste0("/Studies/", sdy),
            schemaName = "gene_expression",
            queryName = "gene_expression_analysis",
            colNameOpt = "rname")),
        error = function(e) {data.frame(analysis_acession = character(),
                                        expression_matrix = character())
        })

      if (nrow(existGEA) > 0) {
        q1 <- quote(arm_name)
        q2 <- quote(coefficient)
        existGEA[, key := paste(eval(q1), eval(q2), sep = " ")]
      }

      # Is there gene expression analysis?
      output_GEAR <- nrow(impliedGEA) != 0 & nrow(existGEA) != 0
      # if there is GEA is it complete?
      if (output_GEAR) {
        output_GEARcomp <- length(setdiff(impliedGEA$unique.values, existGEA$key)) == 0
      } else{
        output_GEARcomp <- NA
      }

      # if not complete what time points are missing
      if (!is.na(output_GEARcomp) & !output_GEARcomp) {
        missingdat <- paste(setdiff(impliedGEA$unique.values, existGEA$key), collapse = ", ")
      } else {
        missingdat <- NA
      }
      output <- c(output_GEAR, output_GEARcomp, missingdat)
    })
    gear <- data.frame(t(gear))
    colnames(gear) <- c("GSEA_implied", "GSEA_complete", "missing_GSEA_data")
    compDF <- merge(compDF, gear, by = 0, all = T)
    compDF[is.na(compDF$GSEA_implied),"GSEA_implied"] <- FALSE
    row.names(compDF) <- compDF$Row.names
    compDF <- compDF[,-1]

    # IRP - Studies with subjects from multiple cohorts with GEM data at both target
    # timepoint and baseline + response data
    nab <- self$getDataset("neut_ab_titer")
    resp <- rbind(nab, hai, fill = TRUE) # nab has a col that hai does not
    inputSmpls$study <- gsub("^SUB\\d{6}\\.", "SDY", inputSmpls$`Participant Id`)
    library(dplyr)
    # NOTE: At least SDY180 has overlapping study_time_collected for both hours and days
    # so it is important to group by study_time_collected_unit as well. This is reflected
    # in IRP_timepoints_hai/nab.sql
    geCohortSubs <- inputSmpls %>%
      group_by(study, `Study Time Collected`, `Study Time Collected Unit`) %>%
      # need multiple cohorts per timePoint (from IRP.js)
      filter((length(unique(Cohort)) > 1 ) == TRUE) %>%
      ungroup() %>%
      group_by(study, Cohort, `Study Time Collected Unit`) %>%
      # need baseline + later timePoint (from IRP.Rmd)
      filter((length(unique(`Study Time Collected`)) > 1 ) == TRUE)

    geRespSubs <- geCohortSubs[ geCohortSubs$`Participant Id` %in% unique(resp$participant_id), ]
    compDF$IRP_implied <- rownames(compDF) %in% .subidsToSdy(geRespSubs$`Participant Id`)
    studyTimepoints <- geRespSubs %>%
      group_by(study) %>%
      summarize(timepoints = paste(unique(`Study Time Collected`), collapse = ","))
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

    # Differences between implied and actual
    compDF$DE_diff <- compDF$DE_implied == compDF$DE_actual
    compDF$GEE_diff <- compDF$GEE_implied == compDF$GEE_actual
    compDF$GSEA_diff <- compDF$GSEA_implied == compDF$GSEA_actual
    compDF$IRP_diff <- compDF$IRP_implied == compDF$IRP_actual

    colOrder <- c(
      "RAW",
      "GEF",
      "GEO",
      "GEM",
      "DE_implied",
      "DE_actual",
      "DE_diff",
      "GEE_implied",
      "GEE_actual",
      "GEE_diff",
      "GSEA_implied",
      "GSEA_actual",
      "GSEA_diff",
      "GSEA_complete",
      "missing_GSEA_data",
      "IRP_implied",
      "IRP_actual",
      "IRP_diff"
    )

    compDF <- compDF[ , order(match(colnames(compDF), colOrder))]

    # Filter out studies that don't have GE since this is basis for everything
    if (filterNonGE == TRUE) {
      compDF <- compDF[compDF$GEO == TRUE | compDF$GEF == TRUE, ]
    }

    # Subset to only show problematic studies
    if (onlyShowNonCompliant == TRUE) {
      redux <- compDF[ , grep("diff", colnames(compDF))]
      idx <- which(apply(redux, 1, all))
      compDF <- compDF[-(idx), ]
    }

    # Defaults to showing only the actual module status and the difference with the implied
    if (showAllCols == FALSE) {
      compDF <- compDF[, grep("diff|act|time|comp", colnames(compDF))]
      message("NOTE: \n Return object has columns marked with 'diff' that indicate a difference between the actual column that was generated by a call to the module url and the believed column that \n was created by looking at the filesystem.  The 'act' column is included as a reference.")
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
