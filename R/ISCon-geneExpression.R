#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------

# List the available gene expression matrices
ISCon$set(
  which = "public",
  name = "listGEMatrices",
  value = function(verbose = FALSE, reload = FALSE, participantIds = NULL) {
    ## HELPERS
    ..getData <- function() {
      try(
        .getLKtbl(
          con = self,
          schema = "assay.ExpressionMatrix.matrix",
          query = "SelectedRuns",
          colNameOpt = "fieldname",
          viewName = "expression_matrices"
        ),
        silent = TRUE
      )
    }


    ## MAIN
    if (is.null(self$cache[[private$.constants$matrices]]) | reload) {
      if (verbose) {
        ge <- ..getData()
      } else {
        ge <- suppressWarnings(..getData())
      }

      if (inherits(ge, "try-error") || nrow(ge) == 0) {
        # No assay or no runs
        message("No gene expression data...")
        self$cache[[private$.constants$matrices]] <- NULL
      } else {
        # adding cols to allow for getGEMatrix() to update
        ge[, cacheinfo := ""][]
        setnames(ge, private$.munge(colnames(ge)))

        # adding cohort_type for use with getGEMatrix(cohort)
        samples <- labkey.executeSql(
          baseUrl = self$config$labkey.url.base,
          folderPath = self$config$labkey.url.path,
          schemaName = "study",
          sql = "
            SELECT DISTINCT expression_matrix_accession, cohort_type
            FROM HM_inputSamplesQuery
            GROUP BY expression_matrix_accession, cohort_type
          ",
          containerFilter = "CurrentAndSubfolders",
          colNameOpt = "fieldname"
        )
        ge$cohort_type <- samples$cohort_type[match(ge$name, samples$expression_matrix_accession)]

        # caching
        self$cache[[private$.constants$matrices]] <- ge
      }
    }
    if (is.null(participantIds)) {
      return(self$cache[[private$.constants$matrices]])
    } else {
      # Get matrices from participantIDs
      # Use assay.ExpressionMatrix.inputSamples in executeSQL to get distinct matrices

      sql <- paste0("SELECT DISTINCT Run.Name
                     FROM InputSamples
                     WHERE Biosample.participantId IN ('", paste0(participantIds, collapse = "','"), "')")

      matrixNames <- Rlabkey::labkey.executeSql(
        baseUrl = self$config$labkey.url.base,
        folderPath = self$config$labkey.url.path,
        schemaName = "assay.ExpressionMatrix.matrix",
        sql = sql,
        containerFilter = "CurrentAndSubfolders",
        colNameOpt = "fieldname"
      )

      return(self$cache[[private$.constants$matrices]][name %in% matrixNames$Name])
    }
  }
)


# List the available gene expression analyses
ISCon$set(
  which = "public",
  name = "listGEAnalysis",
  value = function() {
    GEA <- tryCatch(
      .getLKtbl(
        con = self,
        schema = "gene_expression",
        query = "gene_expression_analysis",
        showHidden = FALSE,
        colNameOpt = "rname"
      ),
      error = function(e) return(e)
    )

    if (length(GEA$message) > 0) {
      stop("Study does not have Gene Expression Analyses.")
    }

    GEA
  }
)


# Download a normalized gene expression matrix from ImmuneSpace
ISCon$set(
  which = "public",
  name = "getGEMatrix",
  value = function(matrixName = NULL,
                     cohortType = NULL,
                     outputType = "summary",
                     annotation = "latest",
                     reload = FALSE,
                     verbose = FALSE) {

    # Handle potential incorrect use of "ImmSig" annotation
    if (outputType == "summary" & annotation == "ImmSig") {
      stop("Not able to provide summary eSets for ImmSig annotated studies. Please use
           'raw' as outputType with ImmSig studies.")
    } else if (annotation == "ImmSig" & !grepl("IS1", self$config$labkey.url.path)) {
      stop("ImmSig annotation only allowable with IS1, no other studies")
    }

    # Handle use of cohortType instead of matrixName
    if (!is.null(cohortType)) {
      ct_name <- cohortType # can't use cohort = cohort in d.t
      if (all(ct_name %in% self$cache$GE_matrices$cohort_type)) {
        matrixName <- self$cache$GE_matrices[cohort_type %in% ct_name, name]
      } else {
        validCohorts <- self$cache$GE_matrices[, cohort_type]
        stop(paste(
          "No expression matrix for the given cohort_type.",
          "Valid cohort_types:", paste(validCohorts, collapse = ", ")
        ))
      }
    }

    # Return a combined eSet of all matrices by default
    if (is.null(matrixName)) {
      matrixName <- self$cache$GE_matrices$name
    }

    # Get matrix or matrices
    esetNames <- vapply(matrixName, function(name) {
      esetName <- .getEsetName(name, outputType, annotation)

      if (esetName %in% names(self$cache) & !reload) {
        message(paste0("Returning ", esetName, " from cache"))
      } else {
        self$cache[[esetName]] <- NULL
        private$.downloadMatrix(name, outputType, annotation, reload)
        private$.getGEFeatures(name, outputType, annotation, reload)
        private$.constructExpressionSet(name, outputType, annotation, verbose)

        # Add to cacheinfo
        cacheinfo_status <- self$cache$GE_matrices$cacheinfo[self$cache$GE_matrices$name == name]
        cacheinfo <- .getcacheinfo(outputType, annotation)
        if (!grepl(cacheinfo, cacheinfo_status)) {
          self$cache$GE_matrices$cacheinfo[self$cache$GE_matrices$name == name] <-
            paste0(
              cacheinfo_status,
              cacheinfo, ";"
            )
        }
      }
      return(esetName)
    },
    FUN.VALUE = "esetName"
    )

    # Combine if needed
    if (length(esetNames) > 1) {
      eset <- .combineEMs(self$cache[esetNames])
      # Handle cases where combineEMs() results in no return object
      if (dim(eset)[[1]] == 0) {
        warn <- "Returned ExpressionSet has 0 rows. No feature is shared across the selected runs or cohorts."
        if (outputType != "summary") {
          warn <- paste(warn, "Try outputType = 'summary' to merge matrices by gene symbol.")
        }
        warning(warn)
      }
    } else {
      eset <- self$cache[[esetNames]]
    }

    if (verbose == TRUE) {
      info <- Biobase::experimentData(eset)
      message("\nNotes:")
      dmp <- lapply(names(info@other), function(nm) {
        message(paste0(nm, ": ", info@other[[nm]]))
      })
      message("\n")
    }

    return(eset)
  }
)


# Retrieve a gene expression analysis
ISCon$set(
  which = "public",
  name = "getGEAnalysis",
  value = function(...) {
    GEAR <- tryCatch(
      .getLKtbl(
        con = self,
        schema = "gene_expression",
        query = "DGEA_filteredGEAR",
        viewName = "DGEAR",
        colNameOpt = "caption",
        ...
      ),
      error = function(e) return(e)
    )

    if (length(GEAR$message) > 0) {
      stop("Gene Expression Analysis not found for study.")
    }

    setnames(GEAR, private$.munge(colnames(GEAR)))

    GEAR
  }
)


# Retrieve gene expression inputs
ISCon$set(
  which = "public",
  name = "getGEInputs",
  value = function() {
    if (!is.null(self$cache[[private$.constants$matrix_inputs]])) {
      self$cache[[private$.constants$matrix_inputs]]
    } else {
      ge <- tryCatch(
        .getLKtbl(
          con = self,
          schema = "assay.Expressionmatrix.matrix",
          query = "InputSamples",
          viewName = "gene_expression_matrices",
          colNameOpt = "fieldname"
        ),
        error = function(e) return(e)
      )

      if (length(ge$message) > 0) {
        stop("Gene Expression Inputs not found for study.")
      }

      setnames(ge, private$.munge(colnames(ge)))
      self$cache[[private$.constants$matrix_inputs]] <- ge
    }
  }
)


# DEPRECATED
ISCon$set(
  which = "public",
  name = "getGEFiles",
  value = function(files, destdir = ".", quiet = FALSE) {
    .Deprecated("downloadGEFiles", old = "getGEFiles")
    self$downloadGEFiles(files, destdir)
  }
)


# Downloads the raw gene expression files to the local machine
#' @importFrom httr HEAD
ISCon$set(
  which = "public",
  name = "downloadGEFiles",
  value = function(files, destdir = ".") {
    stopifnot(file.exists(destdir))

    gef <- copy(con$getDataset(
      "gene_expression_files",
      colSelect = c("participantid", "file_info_name")
    ))
    gef[, study_accession := paste0("SDY", gsub("SUB\\d+.", "", participant_id))]
    gef[, participant_id := NULL]
    gef <- unique(gef[!is.na(file_info_name) & file_info_name %in% files])

    vapply(
      files,
      function(file) {
        if (!file %in% gef$file_info_name) {
          warning(
            file, " is not a valid file name. Skipping downloading this file..",
            call. = FALSE, immediate. = TRUE
          )
          return(FALSE)
        }

        link <- paste0(
          self$config$labkey.url.base,
          "/_webdav/Studies/",
          gef[file == file_info_name]$study_accession[1],
          "/%40files/rawdata/gene_expression/",
          gef[file == file_info_name]$file_info_name[1]
        )

        # check if file exists
        head <- HEAD(link, Rlabkey:::labkey.getRequestOptions())
        if (grepl("json", head$headers$`content-type`)) {
          warning(
            file, "does not exist. ",
            "Please contanct the ImmuneSpace team at immunespace@gmail.com",
            call. = FALSE, immediate. = TRUE
          )
          return(FALSE)
        }

        message("Downloading ", file, "..")
        GET(
          url = link,
          Rlabkey:::labkey.getRequestOptions(),
          write_disk(file.path(destdir, file), overwrite = TRUE)
        )

        TRUE
      },
      FUN.VALUE = logical(1)
    )
  }
)


# Add treatment information to the phenoData of an ExpressionSet
ISCon$set(
  which = "public",
  name = "addTreatment",
  value = function(expressionSet) {
    stopifnot(is(expressionSet, "ExpressionSet"))

    bsFilter <- makeFilter(
      c(
        "biosample_accession",
        "IN",
        paste(pData(expressionSet)$biosample_accession, collapse = ";")
      )
    )

    bs2es <- .getLKtbl(
      con = self,
      schema = "immport",
      query = "expsample_2_biosample",
      colFilter = bsFilter,
      colNameOpt = "rname"
    )

    esFilter <- makeFilter(
      c(
        "expsample_accession",
        "IN",
        paste(bs2es$expsample_accession, collapse = ";")
      )
    )

    es2trt <- .getLKtbl(
      con = self,
      schema = "immport",
      query = "expsample_2_treatment",
      colFilter = esFilter,
      colNameOpt = "rname"
    )

    trtFilter <- makeFilter(
      c(
        "treatment_accession",
        "IN",
        paste(es2trt$treatment_accession, collapse = ";")
      )
    )

    trt <- .getLKtbl(
      con = self,
      schema = "immport",
      query = "treatment",
      colFilter = trtFilter,
      colNameOpt = "rname"
    )

    bs2trt <- merge(bs2es, es2trt, by = "expsample_accession")
    bs2trt <- merge(bs2trt, trt, by = "treatment_accession")

    pData(expressionSet)$treatment <- bs2trt[
      match(
        pData(expressionSet)$biosample_accession,
        biosample_accession
      ),
      name
    ]

    expressionSet
  }
)


# Map the biosample ids in expressionSet object to Immunespace subject IDs
# concatenated with study time collected or experiment sample IDs
ISCon$set(
  which = "public",
  name = "mapSampleNames",
  value = function(EM = NULL, colType = "participant_id") {
    if (is.null(EM) || !is(EM, "ExpressionSet")) {
      stop("EM should be a valid ExpressionSet, as returned by getGEMatrix")
    }

    if (!all(grepl("^BS", sampleNames(EM)))) {
      stop("All sampleNames should be biosample_accession, as returned by getGEMatrix")
    }

    pd <- data.table(pData(EM))
    colType <- gsub("_.*$", "", tolower(colType))

    if (colType == "expsample") {
      bsFilter <- makeFilter(
        c(
          "biosample_accession",
          "IN",
          paste(pd$biosample_accession, collapse = ";")
        )
      )

      bs2es <- .getLKtbl(
        con = self,
        schema = "immport",
        query = "expsample_2_biosample",
        colFilter = bsFilter,
        colNameOpt = "rname"
      )

      pd <- merge(
        pd,
        bs2es[, list(biosample_accession, expsample_accession)],
        by = "biosample_accession"
      )

      sampleNames(EM) <-
        pData(EM)$expsample_accession <-
        pd[
          match(sampleNames(EM), pd$biosample_accession),
          expsample_accession
        ]
    } else if (colType == "participant") {
      pd[, nID := paste0(
        participant_id,
        "_",
        tolower(substr(study_time_collected_unit, 1, 1)),
        study_time_collected)
        ]
      sampleNames(EM) <- pd[match(sampleNames(EM), pd$biosample_accession), nID]
    } else if (colType == "biosample") {
      warning("Nothing done, the column names are already be biosample_accession numbers.")
    } else {
      stop("colType should be one of 'expsample_accession', 'biosample_accession', 'participant_id'.")
    }

    EM
  }
)



# PRIVATE ----------------------------------------------------------------------

# Download the gene expression matrix
#' @importFrom httr GET write_disk
#' @importFrom preprocessCore normalize.quantiles
ISCon$set(
  which = "private",
  name = ".downloadMatrix",
  value = function(matrixName,
                     outputType = "summary",
                     annotation = "latest",
                     reload = FALSE) {

    cache_name <- .getMatrixCacheName(matrixName, outputType, annotation)
    cacheinfo <- .getcacheinfo(outputType, annotation)

    # check if study has matrices
    if (nrow(subset(
      self$cache[[private$.constants$matrices]],
      name %in% matrixName
    )) == 0) {
      stop(sprintf("No matrix %s in study\n", matrixName))
    }

    # check if data in cache corresponds to current request
    # if it does, then no download needed.
    # Only use matrix from cache when
    #   a. outputType and annotation match cache OR
    #   b. outputType matches cache and is not summary
    # Otherwise, load a new matrix
    if (!reload) {
      if (grepl(cacheinfo, self$cache$GE_matrices$cacheinfo[self$cache$GE_matrices$name == matrixName])) {
        message(paste0("Returning ", outputType, " matrix from cache"))
        return()
      }
      if (outputType != "summary") {
        if (grepl(outputType, self$cache$GE_matrices$cacheinfo[self$cache$GE_matrices$name == matrixName])) {
          message(paste0("Returning ", outputType, " matrix from cache"))
          return()
        }
      }
    }

    if (annotation == "ImmSig") {
      fileSuffix <- ".immsig"
    } else {
      if (outputType == "summary") {
        fileSuffix <- switch(
          annotation,
          "latest" = ".summary",
          "default" = ".summary.orig"
        )
      } else {
        fileSuffix <- switch(
          outputType,
          "normalized" = "",
          "raw" = ".raw"
        )
      }
    }
    mxName <- paste0(matrixName, ".tsv", fileSuffix)

    # For HIPC studies, the matrix Import script generates subdirectories
    # based on the original runs table in /Studies/ with the format "Run123"
    # with the suffix being the RowId from the runs table. However, since
    # some of the original runs may have been deleted to fix issues found
    # later on, more complex logic must be used to find the correct flat file.
    if (grepl("HIPC", self$config$labkey.url.path)) {

      # get list of run sub-directories from webdav on /HIPC/ISx
      sdy <- regmatches(
        self$config$labkey.url.path,
        regexpr("IS\\d{1}", self$config$labkey.url.path)
      )
      folder_link <- paste0(
        self$config$labkey.url.base,
        "/_webdav/HIPC/",
        sdy,
        "/%40files/analysis/exprs_matrices?method=JSON"
      )
      runDirs <- unlist(lapply(folder_link, private$.listISFiles))
      runDirs <- grep("Run", runDirs, value = TRUE)

      # Map run sub-directories to the matrixNames passed to downloadMatrix
      id2MxNm <- vapply(runDirs, function(x) {
        run_folder <- gsub("exprs_matrices", paste0("exprs_matrices/", x), folder_link)
        fls <- unlist(lapply(run_folder, private$.listISFiles))
        newNm <- gsub("\\.tsv*", "", grep("tsv", fls, value = TRUE)[[1]])
      },
      FUN.VALUE = "newNm"
      )

      # Generate correct filepath in /HIPC/IS1/@files/analysis/exprs_matrices/
      runId <- names(id2MxNm)[ match(matrixName, id2MxNm) ]
      mxName <- paste0(runId, "/", mxName)
    }

    path <- ifelse(self$config$labkey.url.path == "/Studies/",
                   paste0("/Studies/", self$cache$GE_matrices[name == matrixName, folder], "/"),
                   gsub("^/", "", self$config$labkey.url.path))

    link <- URLdecode(
      file.path(
        gsub("/$", "", self$config$labkey.url.base),
        "_webdav",
        path,
        "@files/analysis/exprs_matrices",
        mxName
      )
    )

    localpath <- private$.localStudyPath(urlPath = link)
    runningLocally <- private$.isRunningLocally(localpath)
    if (runningLocally) {
      message("Reading local matrix")
      fl <- localpath
    } else {
      message("Downloading matrix..")
      opts <- self$config$curlOptions
      opts$options$netrc <- 1L
      fl <- tempfile()
      GET(url = link, config = opts, write_disk(fl))
    }

    EM <- data.table::fread(fl, sep = "\t", header = TRUE)

    if (nrow(EM) == 0) {
      stop("The downloaded matrix has 0 rows. Something went wrong.")
    }

    self$cache[[cache_name]] <- EM

    if (!runningLocally) {
      file.remove(fl)
    }
  }
)


# Get the gene expression features by matrix
ISCon$set(
  which = "private",
  name = ".getGEFeatures",
  value = function(matrixName,
                     outputType = "summary",
                     annotation = "latest",
                     reload = FALSE) {
    cacheinfo <- .getcacheinfo(outputType, annotation)
    cache_name <- paste0(matrixName, cacheinfo)

    if (!(matrixName %in% self$cache[[private$.constants$matrices]]$name)) {
      stop("Invalid gene expression matrix name")
    }

    cacheinfo_status <- self$cache$GE_matrices$cacheinfo[self$cache$GE_matrices$name == matrixName]

    # For raw or normalized, can reuse cached annotation
    if (!reload) {
      if (grepl(cacheinfo, cacheinfo_status)) {
        message(paste0("Returning ", annotation, " annotation from cache"))
        return()
      }
      if (outputType != "summary") {
        if (grepl(paste0("(raw_", annotation, ")|(normalized_", annotation, ")"), cacheinfo_status)) {
          message(paste0("Returning ", annotation, " annotation from cache"))
          return()
        }
      }
    }

    # ---- queries ------
    runs <- labkey.selectRows(
      baseUrl = self$config$labkey.url.base,
      folderPath = self$config$labkey.url.path,
      schemaName = "Assay.ExpressionMatrix.Matrix",
      queryName = "Runs",
      showHidden = TRUE
    )

    faSets <- labkey.selectRows(
      baseUrl = self$config$labkey.url.base,
      folderPath = self$config$labkey.url.path,
      schemaName = "Microarray",
      queryName = "FeatureAnnotationSet",
      showHidden = TRUE
    )

    fasMap <- labkey.selectRows(
      baseUrl = self$config$labkey.url.base,
      folderPath = self$config$labkey.url.path,
      schemaName = "Microarray",
      queryName = "FasMap",
      showHidden = TRUE
    )

    #--------------------

    # Map to correct annotation regardless of name of FAS at time of creation.
    # This is important because for legacy matrices, FAS name may not have '_orig'
    # even though it is the original annotation. 'ImmSig' anno only applies to IS1
    # as other ISx studies will use 'latest' from that study's container.

    if (annotation == "ImmSig") {
      sdy <- regmatches(matrixName, regexpr("SDY\\d{2,3}", matrixName))
      annoSetId <- faSets$`Row Id`[faSets$Name == paste0("ImmSig_", tolower(sdy))]
    } else {
      fasIdAtCreation <- runs$`Feature Annotation Set`[ runs$Name == matrixName ]
      idCol <- ifelse(annotation == "default", "Orig Id", "Curr Id")
      annoAlias <- gsub("_orig", "", faSets$Name[ faSets$`Row Id` == fasIdAtCreation ])
      annoSetId <- fasMap[ fasMap$Name == annoAlias, idCol ]
    }

    if (outputType != "summary") {
      if (paste0("featureset_", annoSetId) %in% names(self$cache)) {
        message(paste0("Returning ", annotation, " annotation from cache"))
      }

      message("Downloading Features..")
      featureAnnotationSetQuery <- sprintf(
        "SELECT * from FeatureAnnotation where FeatureAnnotationSetId='%s';",
        annoSetId
      )
      features <- labkey.executeSql(
        baseUrl = self$config$labkey.url.base,
        folderPath = self$config$labkey.url.path,
        schemaName = "Microarray",
        sql = featureAnnotationSetQuery,
        colNameOpt = "fieldname"
      )
      setnames(features, "GeneSymbol", "gene_symbol")

      # update cache$gematrices with correct fasId
      self$cache$GE_matrices$featureset[self$cache$GE_matrices$name == matrixName] <- annoSetId

      # push features to cache
      self$cache[[paste0("featureset_", annoSetId)]] <- features
    } else {
      # Get annotation from flat file b/c otherwise don't know order
      # NOTE: For ImmSig studies, this means that summaries use the latest
      # annotation even though that was not used in the manuscript for summarizing.
      features <- data.frame(
        FeatureId = self$cache[[cache_name]]$gene_symbol,
        gene_symbol = self$cache[[cache_name]]$gene_symbol
      )
    }
  }
)


# Constructs an expressionSet object with expression matrix (exprs),
# feature annotation data (fData), and subject metadata (pData)
ISCon$set(
  which = "private",
  name = ".constructExpressionSet",
  value = function(matrixName,
                   outputType,
                   annotation,
                   verbose) {

    cache_name <- .getMatrixCacheName(matrixName, outputType, annotation)
    esetName <- .getEsetName(matrixName, outputType, annotation)

    # expression matrix
    message("Constructing ExpressionSet")
    matrix <- self$cache[[cache_name]]

    # pheno
    runID <- self$cache$GE_matrices[name == matrixName, rowid]
    bs <- colnames(matrix)[ grep("^BS\\d{6}$", colnames(matrix))]
    pheno_filter <- makeFilter(
      c(
        "Run",
        "EQUAL",
        runID
      ),
      c(
        "biosample_accession",
        "IN",
        paste(bs, collapse = ";")
      )
    )

    pheno <- unique(
      .getLKtbl(
        con = self,
        schema = "study",
        query = "HM_inputSmplsPlusImmEx",
        containerFilter = "CurrentAndSubfolders",
        colNameOpt = "caption",
        colFilter = pheno_filter,
        showHidden = FALSE
      )
    )

    setnames(pheno, private$.munge(colnames(pheno)))
    pheno <- data.frame(pheno, stringsAsFactors = FALSE)

    # Need cohort for updateGEAR() mapping to arm_accession
    # Need cohort_type for modules
    pheno <- pheno[, colnames(pheno) %in% c(
      "biosample_accession",
      "participant_id",
      "cohort_type",
      "cohort",
      "study_time_collected",
      "study_time_collected_unit",
      "exposure_material_reported",
      "exposure_process_preferred"
    )]

    # ensure same order as GEM rownames
    rownames(pheno) <- pheno$biosample_accession
    order <- names(self$cache[[cache_name]])
    order <- order[-grep("feature_id|gene_symbol|X|V1", order)] # X|V1 kept for IS1 raw
    order <- order[ order != "BS694717.1" ] # rm SDY212 dup for the moment
    pheno <- pheno[match(order, row.names(pheno)), ]

    # handling multiple experiment samples per biosample (e.g. technical replicates)
    dups <- colnames(matrix)[duplicated(colnames(matrix))]
    if (length(dups) > 0) {
      for (dup in dups) {
        dupIdx <- grep(dup, colnames(matrix))
        matrix[ , dupIdx[[1]] ] <- rowMeans(matrix[, dupIdx, with = FALSE])
        matrix[ , (dupIdx[2:length(dupIdx)]) := NULL ]
      }
      if (verbose) {
        warning(
          "The matrix contains subjects with multiple measures per timepoint. ",
          "Averaging the expression values ..."
        )
      }
    }

    # gene features
    if (outputType == "summary") {
      fdata <- data.frame(
        FeatureId = matrix$gene_symbol,
        gene_symbol = matrix$gene_symbol
      )

      rownames(fdata) <- rownames(matrix) <- matrix$gene_symbol # exprs and fData must match
    } else {
      annoSetId <- self$cache$GE_matrices$featureset[self$cache$GE_matrices$name == matrixName]
      features <- self$cache[[paste0("featureset_", annoSetId)]][, c("FeatureId", "gene_symbol")]
      # IS1 matrices have not been standardized, otherwise all others should be 'feature_id'
      colnames(matrix)[[grep("feature_id|X|V1", colnames(matrix))]] <- "FeatureId"

      # Only known case is SDY300 for "2-Mar" and "1-Mar" which are
      # likely not actual probe_ids but mistransformed strings
      if (any(duplicated(matrix$FeatureId))) {
        matrix <- matrix[ !duplicated(matrix$FeatureId), ]
      }

      fdata <- data.frame(
        FeatureId = as.character(matrix$FeatureId),
        stringsAsFactors = FALSE
      )

      fdata <- merge(fdata, features, by = "FeatureId", all.x = TRUE)

      # exprs and fData must match
      rownames(fdata) <- fdata$FeatureId
      matrix <- matrix[ order(match(matrix$FeatureId, fdata$FeatureId)), ]
      rownames(matrix) <- matrix$FeatureId
    }

    # SDY212 has dbl biosample that is need for IS1 normalization, but later removed
    # in the IS1 report.
    if ("BS694717" %in% pheno$biosample_accession) {
      pheno["BS694717.1", ] <- pheno[pheno$biosample_accession == "BS694717", ]
      pheno$biosample_accession[rownames(pheno) == "BS694717.1"] <- "BS694717.1"
    }

    # Prep Eset and push
    # NOTES: At project level, InputSamples may be filtered
    matrix <- data.frame(matrix, stringsAsFactors = FALSE) # for when on rsT / rsP
    exprs <- matrix[, colnames(matrix) %in% pheno$biosample_accession] # rms gene_symbol!
    pheno <- pheno[colnames(exprs), ]

    # add processing information for user
    fasInfo <- .getLKtbl(
      con = self,
      schema = "Microarray",
      query = "FeatureAnnotationSet"
    )
    gemx <- self$cache$GE_matrices
    fasId <- gemx$featureset[ gemx$name == matrixName ]
    fasInfo <- fasInfo[ match(fasId, fasInfo$`Row Id`)]
    isRNA <- (fasInfo$Vendor == "NA" & !grepl("ImmSig", fasInfo$Name)) | grepl("SDY67", fasInfo$Name)
    if (fasInfo$Comment == "Do not update" | is.na(fasInfo$Comment)) {
      annoVer <- annotation
    } else {
      annoVer <- gsub(" ", "", strsplit(fasInfo$Comment, ":")[[1]][2])
    }

    processInfo <- list(
      normalization = ifelse(isRNA, "DESeq", "normalize.quantiles"),
      summarizeBy = ifelse(outputType == "summary", "mean", "none"),
      org.Hs.eg.db_version = annoVer,
      featureAnnotationSet = fasInfo$Name
    )

    self$cache[[esetName]] <- ExpressionSet(
      assayData = as.matrix(exprs),
      phenoData = AnnotatedDataFrame(pheno),
      featureData = AnnotatedDataFrame(fdata),
      experimentData = new("MIAME", other = processInfo)
    )
  }
)


# Get feature ID by matrix
ISCon$set(
  which = "private",
  name = ".getFeatureId",
  value = function(matrixName) {
    subset(self$cache[[private$.constants$matrices]], name %in% matrixName)[, featureset]
  }
)


# Rename the feature ID
ISCon$set(
  which = "private",
  name = ".mungeFeatureId",
  value = function(annotation_set_id) {
    sprintf("featureset_%s", annotation_set_id)
  }
)



# HELPER -----------------------------------------------------------------------

# Get the cache name of expression matrix by output type and annotation
.getMatrixCacheName <- function(matrixName, outputType, annotation) {
  outputSuffix <- switch(
    outputType,
    "summary" = "_sum",
    "normalized" = "_norm",
    "raw" = "_raw"
  )
  annotationSuffix <- switch(
    annotation,
    "latest" = "_latest",
    "default" = "_default",
    "ImmSig" = "_immsig"
  )

  matrixName <- paste0(matrixName, outputSuffix)
  if (annotation == "ImmSig" || outputType == "summary") {
    matrixName <- paste0(matrixName, annotationSuffix)
  }
  return(matrixName)
}

# Get the cache name for eset by output type and annotation
.getEsetName <- function(matrixName, outputType, annotation) {
  esetName <- paste0(matrixName, "_", outputType, "_", annotation, "_eset")
  return(esetName)
}

# Get cacheinfo string from output type and annotation
.getcacheinfo <- function(outputType, annotation) {
  cacheinfo <- paste0(outputType, "_", annotation)
  return(cacheinfo)
}

# Combine EMs and output only genes available in all EMs.
.combineEMs <- function(EMlist) {
  message("Combining ExpressionSets")

  fds <- lapply(EMlist, function(x) {
    droplevels(data.table(fData(x)))
  })

  fd <- Reduce(f = function(x, y) {
    merge(x, y, by = c("FeatureId", "gene_symbol"))
  }, fds)

  EMlist <- lapply(EMlist, "[", as.character(fd$FeatureId))

  for (i in seq_len(length(EMlist))) {
    fData(EMlist[[i]]) <- fd
  }

  Reduce(f = combine, EMlist)
}
