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
                            "gene_expression_matrices"),
                   mc.cores = 1) {
    ## HELPERS
    ..messageResults <- function(dataset, file_exists) {
      message(
        paste0(
          sum(file_exists),
          "/",
          length(file_exists),
          " ",
          dataset,
          " with valid links."
        )
      )
    }

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
            mc.cores = mc.cores
          )
        )

        file_exists <- temp$file_info_name %in% file_list

        res <- data.frame(
          file_info_name = temp$file_info_name,
          study_accession = temp$study_accession,
          file_link = file_link,
          file_exists = file_exists,
          stringsAsFactors = FALSE
        )

        ..messageResults(dataset, res$file_exists)
      }

      res
    }


    ## MAIN

    startTimeTotal <- Sys.time()

    ret <- list()
    what <- tolower(what)

    if ("gene_expression_files" %in% what) {
      startTime <- Sys.time()

      ret$gene_expression_files <- ..checkLinks(
        "gene_expression_files",
        "gene_expression"
      )

      endTime <- Sys.time()
      print( endTime - startTime )
    }

    if ("fcs_sample_files" %in% what) {
      startTime <- Sys.time()

      ret$fcs_sample_files <- ..checkLinks(
        "fcs_sample_files",
        "flow_cytometry"
      )

      endTime <- Sys.time()
      print( endTime - startTime )
    }

    if ("fcs_control_files" %in% what) {
      startTime <- Sys.time()

      ret$fcs_control_files <- ..checkLinks(
        "fcs_control_files",
        "flow_cytometry"
      )

      endTime <- Sys.time()
      print( endTime - startTime )
    }

    if ("protocols" %in% what) {
      startTime <- Sys.time()

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
          mc.cores = mc.cores
        )
      )

      ..messageResults("protocols", file_exists)

      ret$protocols <- data.frame(
        file_info_name = paste0(folders, "_protocol.zip"),
        study_accession = folders,
        file_link = file_link,
        file_exists = file_exists,
        stringsAsFactors = FALSE
      )

      endTime <- Sys.time()
      print( endTime - startTime )
    }

    if ("gene_expression_matrices" %in% what) {
      startTime <- Sys.time()

      suppressWarnings(
        mx <- .getLKtbl(
          con = self,
          schema = "assay.ExpressionMatrix.matrix",
          query = "Runs",
          colNameOpt = "rname"
        )
      )

      if (nrow(mx) > 0) {
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
            mc.cores = mc.cores
          )
        )

        ..messageResults("gene_expression_matrices", file_exists)

        ret$gene_expression_matrices <- data.frame(
          file_info_name = paste0(mx$name, ".tsv"),
          study_accession = mx$folder_name,
          file_link = mxLinks,
          file_exists = file_exists,
          stringsAsFactors = FALSE
        )
      } else {
        ret$gene_expression_matrices <- data.frame(
          file_info_name = NULL,
          study_accession = NULL,
          file_link = NULL,
          file_exists = NULL,
          stringsAsFactors = FALSE
        )
      }

      endTime <- Sys.time()
      print( endTime - startTime )
    }

    endTimeTotal <- Sys.time()
    print( '===========' )
    print( 'TOTAL TIME:' )
    print( endTimeTotal - startTimeTotal )

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
# following modules/data files (GEF, RAW, GEO, GEM, DE, GEE, IRP, GSEA, DR)
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
      url <-  url <- paste0(baseUrl, "/immport/studies/containersformodule.api?name=", name)
      res <- unlist(lapply(rjson::fromJSON(Rlabkey:::labkey.get(url))[[1]], function(x) {x[["name"]]}))
      res <- .spSort(res[grepl("SDY[0-9]+", res)])
    }

    ## MAIN
    # For labkey.executeSql calls
    baseUrl <- self$config$labkey.url.base

    # res table
    if (self$study != "Studies") {
      Sdys <- c(self$study) # single study
    } else {
      Sdys <- private$.getSdyVec()
      Sdys <- .spSort(Sdys)
    }

    mods <- c("GEF",
              "RAW",
              "GEO",
              "GEM_implied",
              "GEM_actual",
              "DE_implied",
              "GSEA_implied",
              "GEE_implied",
              "IRP_implied")

    compDF <- data.frame(
      matrix(
        nrow = length(Sdys),
        ncol = length(mods)
      ),
      row.names = Sdys
    )
    colnames(compDF) <- mods

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

    # GEM possible
    compDF$GEM_implied <- compDF$RAW == T | compDF$GEO == T

    # GEM actual
    withGems <- unique(self$cache$GE_matrices$folder)
    compDF$GEM_actual <- rownames(compDF) %in% withGems

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

    inputSmpls$study <- gsub("SUB[^>]+\\.", "SDY",inputSmpls$participantid)

    exprResp <- merge(
      inputSmpls,
      hai,
      by.x = c("participantid", "study_time_collected"),
      by.y = c("participant_id", "study_time_collected")
    )

    compDF$GEE_implied <- rownames(compDF) %in% .subidsToSdy(unique(exprResp$participantid))

    # GSEA - studies with subjects having GEAR results for multiple non-baseline timepoints
    gearSql <- "SELECT DISTINCT analysis_accession.coefficient FROM gene_expression_analysis_results"
    gear <- sapply(withGems, FUN = function(sdy){
      res <- suppressWarnings(
               tryCatch(
                 labkey.executeSql(baseUrl = baseUrl,
                                        folderPath = paste0("/Studies/", sdy),
                                        schemaName = "gene_expression",
                                        sql = gearSql),
                  error = function(e){ return( NA ) }
                 )
               )
      output <- !is.na(res) && nrow(res) > 1
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

    gea <- suppressWarnings(lapply(withGems, FUN = function(sdy) {
      impliedGEA <- data.table(labkey.selectRows(
        baseUrl = baseUrl,
        folderPath = paste0("/Studies/", sdy),
        schemaName = "assay.expressionMatrix.matrix",
        queryName = "inputSamples",
        colNameOpt = "rname",
        showHidden = TRUE))

      # ---- summarize to have same form and info as currGEA ----

      # 1. Remove all arm_name * study_time_collected with less than 4 replicates
      # otherwise predictive modeling cannot work
      impliedGEA[, subs := unique(length(biosample_participantid)), by = .(biosample_arm_name, biosample_study_time_collected)]
      impliedGEA <- impliedGEA[ subs > 3 ]

      # 2. Check for baseline within each arm_name and then filter out baseline
      impliedGEA[, baseline := any(biosample_study_time_collected <= 0), by = .(biosample_arm_name) ]
      impliedGEA <- impliedGEA[ baseline == TRUE ]
      impliedGEA <- impliedGEA[ biosample_study_time_collected > 0 ]

      # 3. Generate key
      impliedGEA[, key := paste(biosample_arm_name, biosample_study_time_collected, biosample_study_time_collected_unit)]

      # 4. Summarize by arm_name * study_time_collected for number of subs and key
      smryGEA <- impliedGEA[ , list(key = unique(key), subs = unique(subs)), by = .(biosample_arm_name, biosample_study_time_collected)]

      # -------------------------------------------

      # Get current GEA and compare
      currGEA <- existGEA[ existGEA$sdy == sdy, ]
      currGEA$key <- paste(currGEA$arm_name, currGEA$coefficient)

      if (nrow(smryGEA) > 0) {
        diff <- sort(setdiff(smryGEA$key, currGEA$key)) # In implied and NOT in current
        missing_data <- if(length(diff) == 0 ){ "no diff" }else{ paste(diff, collapse = "; ") }
      } else {
        missing_data <- NA
      }

      res <- c( nrow(smryGEA) > 0, missing_data)
    }))

    gea <- data.frame(do.call(rbind, gea), stringsAsFactors = FALSE)
    colnames(gea) <- c("DGEA_implied", "DGEA_missing")
    rownames(gea) <- withGems
    compDF <- merge(compDF, gea, by = 0, all = TRUE)
    compDF$DGEA_implied[ is.na(compDF$DGEA_implied) ] <- FALSE
    compDF$DGEA_implied <- as.logical(compDF$DGEA_implied)
    row.names(compDF) <- compDF$Row.names
    compDF <- compDF[,-1]

    # IRP - Studies with subjects from multiple cohorts with GEM data at both target
    # timepoint and baseline + response data that has baseline and a later timepoint.
    # GEM and response later timepoints do NOT need to be the same!
    nab <- self$getDataset("neut_ab_titer")
    resp <- rbind(nab, hai, fill = TRUE) # nab has a col that hai does not
    resp <- resp[, list(study_time_collected, study_time_collected_unit,
                    response = value_preferred/mean(value_preferred[study_time_collected <= 0]), na.rm = T),
             by = "virus,participant_id"]
    resp <- resp[ !is.na(response) ]

    # NOTE: At least SDY180 has overlapping study_time_collected for both hours and days
    # so it is important to group by study_time_collected_unit as well. This is reflected
    # in IRP_timepoints_hai/nab.sql.
    inputSmpls <- data.table(inputSmpls)
    geCohortSubs <- inputSmpls[ participantid %in% resp$participant_id ]
    geCohortSubs <- geCohortSubs[ , .SD[length(unique(cohort)) > 1], by = .(study, study_time_collected, study_time_collected_unit)]
    geCohortSubs <- geCohortSubs[, .SD[length(unique(study_time_collected)) > 1 & 0 %in% unique(study_time_collected)], by = .(study, cohort, study_time_collected_unit)]
    compDF$IRP_implied <- rownames(compDF) %in% unique(geCohortSubs$study)

    studyTimepoints <- geCohortSubs[ , list(timepoints = paste(sort(unique(study_time_collected)),
                                                               collapse = ",")),
                                     by = .(study)]
    compDF$IrpTimepoints <- studyTimepoints$timepoints[ match(rownames(compDF), studyTimepoints$study) ]

    compDF$DE_actual <- rownames(compDF) %in% ..getModSdys("DataExplorer")
    compDF$GEE_actual <- rownames(compDF) %in% ..getModSdys("GeneExpressionExplorer")
    compDF$GSEA_actual <- rownames(compDF) %in% ..getModSdys("GeneSetEnrichmentAnalysis")
    compDF$IRP_actual <- rownames(compDF) %in% ..getModSdys("ImmuneResponsePredictor")
    compDF$DGEA_actual <- rownames(compDF) %in% ..getModSdys("DifferentialExpressionAnalysis")

    # DE - b/c ISC_study_datasets cannot provide gene_expression info, we use
    # compDF$DGEA_actual as a proxy since it pulls the current GEA query, which should
    # have the same info as DGEA_filteredGEAR ( what the con$plot() uses via con$getGEAnalysis() )
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
      ret <- any(res[[1]] %in% deSets) | compDF$DGEA_actual[rownames(compDF) == sdy] == TRUE
    })

    # ---- dimension reduction (DR) ----
    # NOTE:  Dimension reduction is pointless or impossible to run on a study where
    # NOTE:  there is not enough data to come up with a matrix where both dimensions
    # NOTE:  are 3 or greater.
    # NOTE:  For now, it seems that simply checking to see if there are any assay/subject
    # NOTE:  combinations where number of subjects and features are both greater than 3
    # NOTE:  is a good enough estimate of whether or not
    # NOTE:  dimension reduction makes sense or is possible without doing too much of
    # NOTE:  the checks and filtering that already happens in the module. It may mark as false
    # NOTE:  some studies where dimension reduction might be possible when including
    # NOTE:  multiple timepoints or assays. To fix this, we could further group
    # NOTE:  by timepoint or assay to get an idea of what the dimensions would be.

    # Get current
    compDF$DR_actual <- rownames(compDF) %in% ..getModSdys("DimensionReduction")

    # get dimension reduction assay info
    dimRedux_assay_data <- labkey.selectRows(
      baseUrl = baseUrl,
      folderPath = "/Studies",
      schemaName = "study",
      queryName = "DimRedux_assay_data",
      containerFilter = "CurrentAndSubfolders",
      colNameOpt = "rname"
    )

    # Add a column for study
    dimRedux_assay_data$study <- gsub("SUB[^>]+\\.", "SDY",dimRedux_assay_data$participantid)
    setDT(dimRedux_assay_data)

    # > dimRedux_assay_data
    # participantid                name          label timepoint features  study
    # 1:  SUB102424.28               elisa          ELISA    1 Days       10  SDY28
    # 2:  SUB102425.28               elisa          ELISA    1 Days       24  SDY28
    # 3:  SUB102426.28               elisa          ELISA    1 Days       18  SDY28
    # 4:  SUB102427.28               elisa          ELISA    1 Days       12  SDY28
    # 5:  SUB102428.28 fcs_analyzed_result Flow Cytometry    1 Days     1140  SDY28
    # ---
    # 28442:  SUB86640.241             elispot        ELISPOT   14 Days        1 SDY241
    # 28443:  SUB86641.241               elisa          ELISA   14 Days        6 SDY241
    # 28444:  SUB86641.241             elispot        ELISPOT   14 Days        2 SDY241
    # 28445:  SUB86642.241               elisa          ELISA   14 Days        6 SDY241
    # 28446:  SUB86642.241             elispot        ELISPOT   14 Days        1 SDY241
    #

    # Group by study, timepoint, and assay, and get the number of subjects and features for that
    # assay
    dimensionInfo <- dimRedux_assay_data[, .(subjectCount = length(unique(participantid)),
                                             featureCount = min(features)),
                                         by = c("study", "timepoint", "name")]

    # Are there any lines where subject count and feature count are both greater than three?
    dimensionInfo[, dimMinMet := subjectCount >= 3 & featureCount >= 3]
    DRPossible <- dimensionInfo[, .(dimMinMet = any(dimMinMet)), by = "study"]

    # Add to compDF
    compDF$DR_implied <- rownames(compDF) %in% DRPossible[dimMinMet == TRUE, study]

    colOrder <- c(
      "RAW",
      "GEF",
      "GEO",
      "GEM_implied",
      "GEM_actual",
      "DE_implied",
      "DE_actual",
      "GEE_implied",
      "GEE_actual",
      "DGEA_implied",
      "DGEA_actual",
      "DGEA_missing",
      "IRP_implied",
      "IRP_actual",
      "IrpTimepoints",
      "GSEA_implied",
      "GSEA_actual",
      "DR_implied",
      "DR_actual")

    rowOrder <- Sdys[order(gsub("([A-Z]+)([0-9]+)", "\\1", Sdys),
                             as.numeric(gsub("([A-Z]+)([0-9]+)", "\\2", Sdys)))]

    compDF <- compDF[order(match(row.names(compDF), rowOrder)) , order(match(colnames(compDF), colOrder))]

    # Filter out studies that don't have GE since this is basis for everything
    if (filterNonGE == TRUE) {
      compDF <- compDF[compDF$GEO == TRUE | compDF$GEF == TRUE, ]
    }
        # Subset to only show problematic studies
    if (onlyShowNonCompliant == TRUE) {
      redux <- compDF[ , grep("implied|actual|missing", colnames(compDF))]
      mod_sub <- c("DE", "GEE", "IRP", "GSEA", "DGEA")
      compliant <- lapply(mod_sub, FUN = function(mod) {
        idx <- grepl(mod, names(redux))
        sub <- redux[,idx]
        if (mod == "DGEA") {
          imp_vs_act <- sub[,1] == sub[,2]
          missing_dat <- is.na(sub[,3]) == T | sub[,3] == "no diff"
          compliant <- imp_vs_act == missing_dat
        } else {compliant <- sub[,1] == sub[,2]}
      })
      compliant[[6]] <- compDF$GEM_implied == compDF$GEM_actual
      compliant <- do.call(cbind, compliant)
      row.names(compliant) <- row.names(redux)
      idx <- which(apply(compliant, 1, all))
      compDF <- compDF[-(idx),]
    }

    # Defaults to showing only the actual module status and the difference with the implied
    if (showAllCols == FALSE) {
      compDF <- compDF[, grep("act|implied$", colnames(compDF))]
    }

    if (verbose == TRUE) {
      message("NOTE: \n Return objects have an actual column that was generated by a call to the module url and an implied column that \n was created by looking at the filesystem.")
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
    studies <- studies[grepl('SDY[0-9]+', studies)]

    studies
  }
)


# Get names of files in a single folder from webdav link
ISCon$set(
  which = "private",
  name = ".listISFiles",
  value = function(link) {
    response <- NULL
    res <- tryCatch(
      Rlabkey:::labkey.get(link),
      warning = function(w) return(w),
      error = function(e) return(NULL)
    )
    if (!is.null(res)) {
      tmp <- rjson::fromJSON(res)
      response <- sapply(tmp$files, function(x){ return(x$text) }) # basename only
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
          files <- files[ grep("\\.(tsv|csv|cel|txt)$", files, ignore.case = T) ]
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


# getGEMatrix test function
ISCon$set(
  which = "private",
  name = ".checkExpressionSet",
  value = function(allMatrices = FALSE,
                   ...) {
  # Grab expression set
  mat_names <- ifelse(allMatrices == FALSE, self$cache$GE_matrices$name[[1]], self$cache$GE_matrices$name)
  es <- self$getGEMatrix(mat_names, ...)
  opts <- list(...)

  if (length(opts) == 0) {
    opts <- list()
    opts$outputType <- "summary"
  } else {
    opts <- list(...)
  }
  # expression matrix +
  em <- Biobase::exprs(es)
  pd <- Biobase::pData(es)

  res <- cbind(.checkEM(em, opts, self), .checkPD(pd, self), .checkBiosample(em,pd))
  return(res)
})



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

# expression matrix
.checkEM <- function(em, opts, self){
  # get feature set
  anno <- labkey.selectRows(
    baseUrl = self$config$labkey.url.base,
    folderPath = self$config$labkey.url.path,
    schemaName = "Microarray",
    queryName = "FeatureAnnotation",
    colSelect = c("FeatureId", "GeneSymbol"),
    maxRows = 20,
    showHidden = TRUE)

  # compare to gem rows
  anno <- ifelse(opts$outputType == "summary", anno$`Gene Symbol`, anno$`Feature Id`)
  anno_match <- all(anno %in% row.names(em))
  # check range (log2)
  expr_within_range <- all(0 < range(em) & range(em) < 30)
  # check num of genes
  min_genes <- ifelse(opts$outputType == 'summary', 10000, 20000)
  gene_num <- length(row.names(em)) >= min_genes
  res <- data.frame(outputType = opts$outputType, anno_match, expr_within_range, gene_num)
  return(res)
}

# Check pdata
.checkPD <- function(pd, self){
  cohort_type_col <- "cohort_type" %in% colnames(pd)
  ct_split <- do.call(rbind, strsplit(pd$cohort_type, "_", fixed = TRUE))
  # does cohort type cohort match pd$cohort
  cohort_match <- all(ct_split[,1] == pd$cohort)
  # does type == labkey lookup for cell type
  lk_smpl_type <- labkey.selectRows(baseUrl = self$config$labkey.url.base,
                                    folderPath = self$config$labkey.url.path,
                                    schemaName = "immport",
                                    queryName = "lk_sample_type",
                                    showHidden = TRUE,
                                    colSelect = "Name")
  type_match <- all(ct_split[,2] %in% lk_smpl_type$Name)
  res <- data.frame(cohort_type_col, cohort_match, type_match)
  return(res)
}

# Check biosamples from pdata and expression matrix
.checkBiosample <- function(em, pd){
  # change to all.equal fxn
  biosample_match<- all.equal(row.names(pd), colnames(em))
  if (all(biosample_match != TRUE)) {
    biosample_match <- FALSE
  }
  res <- data.frame(biosample_match)
  return(res)
}


