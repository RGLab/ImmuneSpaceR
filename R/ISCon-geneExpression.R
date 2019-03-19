#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------

# List the available gene expression matrices
ISCon$set(
  which = "public",
  name = "listGEMatrices",
  value = function(verbose = FALSE, reload = FALSE) {
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
        ge[, annotation := ""][, outputType := ""][] # see data.table #869
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

    self$cache[[private$.constants$matrices]]
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
        # SDY67 is special case. "Batch2" matrix is only day 0 and has overlapping
        # biosamples with full matrix "SDY67_HealthyAdults".  This causes
        # full matrix for cohort.
        if (grepl("SDY67", matrixName[[1]])) {
          matrixName <- matrixName[ grep("Batch2", matrixName, invert = T) ]
        }
      } else {
        validCohorts <- self$cache$GE_matrices[, cohort_type]
        stop(paste("No expression matrix for the given cohort_type.",
                   "Valid cohort_types:", paste(validCohorts, collapse = ", ")))
      }
    }

    cache_name <- .setCacheName(matrixName, outputType)
    esetName <- paste0(cache_name, "_eset")

    # Multiple matrices
    if (length(matrixName) > 1) {
      lapply(matrixName, private$.downloadMatrix, outputType, annotation, reload)
      lapply(matrixName, private$.getGEFeatures, outputType, annotation, reload)
      lapply(matrixName, private$.constructExpressionSet, outputType, annotation)
      ret <- .combineEMs(self$cache[esetName])

      # Handle cases where combineEMs() results in no return object
      if (dim(ret)[[1]] == 0) {
        warn <- "Returned ExpressionSet has 0 rows. No feature is shared across the selected runs or cohorts."
        if (outputType != "summary") {
          warn <- paste(warn, "Try outputType = 'summary' to merge matrices by gene symbol.")
        }
        warning(warn)
      }

      if (verbose == TRUE) {
        info <- Biobase::experimentData(ret)
        message("\nNotes:")
        dmp <- lapply(names(info@other), function(nm){
          message(paste0(nm, ": ", info@other[[nm]]))
        })
      }

      return(ret)

    # Single matrix
    } else {
      if (esetName %in% names(self$cache) & !reload) {
        message(paste0("returning ", esetName, " from cache"))
      } else {
        self$cache[[esetName]] <- NULL
        private$.downloadMatrix(matrixName, outputType, annotation, reload)
        private$.getGEFeatures(matrixName, outputType, annotation, reload)
        private$.constructExpressionSet(matrixName, outputType, annotation)
      }

      if (verbose == TRUE) {
        info <- Biobase::experimentData(self$cache[[esetName]])
        message("\nNotes:")
        dmp <- lapply(names(info@other), function(nm){
                 message(paste0(nm, ": ", info@other[[nm]]))
               })
      }

      return(self$cache[[esetName]])
    }
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


# Downloads the raw gene expression files to the local machine
ISCon$set(
  which = "public",
  name = "getGEFiles",
  value = function(files, destdir = ".", quiet = FALSE) {
    links <- paste0(
      self$config$labkey.url.base,
      "/_webdav/",
      self$config$labkey.url.path,
      "/%40files/rawdata/gene_expression/",
      files
    )

    sapply(
      links,
      function(x) {
        download.file(
          url = links[1],
          destfile = file.path(destdir, basename(x)),
          method = "curl",
          extra = "-n",
          quiet = quiet)
      }
    )
  }
)


# Add treatment information to the phenoData of an expression matrix available
# in the connection object.
ISCon$set(
  which = "public",
  name = "addTreatment",
  value = function(matrixName = NULL) {
    if (is.null(matrixName) || !matrixName %in% names(self$cache)) {
      stop(paste(matrixName, "is not a valid expression matrix."))
    }

    bsFilter <- makeFilter(
      c("biosample_accession",
        "IN",
        paste(pData(self$cache[[x]])$biosample_accession, collapse = ";")
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

    pData(self$cache[[x]])$treatment <- bs2trt[match(pData(self$cache[[x]])$biosample_accession,
                                                     biosample_accession), name]

    self$cache[[x]]
  }
)


# Map the expression set by the sample names
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
        pd[match(sampleNames(EM), pd$biosample_accession),
           expsample_accession]
    } else if (colType %in% c("participant", "subject")) {
      pd[, nID := paste0(participant_id,
                         "_",
                         tolower(substr(study_time_collected_unit, 1, 1)),
                         study_time_collected)
         ]
      sampleNames(EM) <- pd[match(sampleNames(EM), pd$biosample_accession), nID]
    } else if (colType == "biosample") {
      warning("Nothing done, the column names should already be biosample_accession numbers.")
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
    cache_name <- .setCacheName(matrixName, outputType)

    # check if study has matrices
    if (nrow(subset(self$cache[[private$.constants$matrices]],
                    name %in% matrixName)) == 0) {
      stop(sprintf("No matrix %s in study\n", matrixName))
    }

    # check if data in cache corresponds to current request
    # if it does, then no download needed.
    status <- self$cache$GE_matrices$outputtype[self$cache$GE_matrices$name == matrixName]
    if (status == outputType & reload != TRUE) {
      message(paste0("returning ", outputType, " matrix from cache"))
      return()
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
      sdy <- regmatches(self$config$labkey.url.path,
                        regexpr("IS\\d{1}", self$config$labkey.url.path))
      folder_link <- paste0(
        self$config$labkey.url.base,
        "/_webdav/HIPC/",
        sdy,
        "/%40files/analysis/exprs_matrices?method=JSON"
      )
      runDirs <- unlist(lapply(folder_link, private$.listISFiles))
      runDirs <- grep("Run", runDirs, value = TRUE)

      # Map run sub-directories to the matrixNames passed to downloadMatrix
      id2MxNm <- sapply(runDirs, function(x){
        run_folder <- gsub("exprs_matrices", paste0("exprs_matrices/", x), folder_link)
        fls <- unlist(lapply(run_folder, private$.listISFiles))
        newNm <- gsub("\\.tsv*", "", grep("tsv", fls, value = TRUE)[[1]])
      })

      # Generate correct filepath in /HIPC/IS1/@files/analysis/exprs_matrices/
      runId <- names(id2MxNm)[ match(matrixName, id2MxNm) ]
      mxName <- paste0(runId, "/", mxName)
    }

    if (self$config$labkey.url.path == "/Studies/") {
      path <- paste0("/Studies/", self$cache$GE_matrices[name == matrixName, folder], "/")
    } else {
      path <- gsub("^/", "", self$config$labkey.url.path)
    }

    link <- URLdecode(
      file.path(
        gsub("http:",
             "https:",
             gsub("/$", "", self$config$labkey.url.base)
        ),
        "_webdav",
        path,
        "@files/analysis/exprs_matrices",
        mxName
      )
    )

    localpath <- private$.localStudyPath(link)
    if (private$.isRunningLocally(localpath)) {
      message("Reading local matrix")
      self$cache[[cache_name]] <- read.table(
        localpath,
        header = TRUE,
        sep = "\t",
        quote = "\"",
        stringsAsFactors = FALSE
      )
    } else {
      opts <- self$config$curlOptions
      opts$options$netrc <- 1L

      message("Downloading matrix..")
      fl <- tempfile()
      GET(url = link, config = opts, write_disk(fl))

      # fread does not read correctly
      # SDY1289 has gene symbols in original version with single quote, therefore need
      # to change 'quote' so that only looks for double-quote
      EM <- read.table(
        fl,
        header = TRUE,
        sep = "\t",
        quote = "\"",
        stringsAsFactors = FALSE
      )

      if (nrow(EM) == 0) {
        stop("The downloaded matrix has 0 rows. Something went wrong.")
      }

      self$cache[[cache_name]] <- EM
      file.remove(fl)
    }

    # Be sure to note which output is already in cache. Colnames are "munged"
    self$cache$GE_matrices$outputtype[self$cache$GE_matrices$name == matrixName] <- outputType
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
    cache_name <- .setCacheName(matrixName, outputType)

    if (!(matrixName %in% self$cache[[private$.constants$matrices]]$name)) {
      stop("Invalid gene expression matrix name")
    }

    status <- self$cache$GE_matrices$annotation[self$cache$GE_matrices$name == matrixName]
    currOut <- self$cache$GE_matrices$outputtype[self$cache$GE_matrices$name == matrixName]
    if (status == annotation & reload != TRUE & currOut == outputType) {
      message(paste0("returning ", annotation, " annotation from cache"))
      return()
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
      idCol <- ifelse( annotation == "default", "Orig Id", "Curr Id")
      annoAlias <- gsub("_orig", "", faSets$Name[ faSets$`Row Id` == fasIdAtCreation ])
      annoSetId <- fasMap[ fasMap$Name == annoAlias, idCol ]
    }

    if (outputType != "summary") {
      message("Downloading Features..")
      featureAnnotationSetQuery <- sprintf("SELECT * from FeatureAnnotation
                                          where FeatureAnnotationSetId='%s';",
                                          annoSetId)
      features <- labkey.executeSql(
        baseUrl = self$config$labkey.url.base,
        folderPath = self$config$labkey.url.path,
        schemaName = "Microarray",
        sql = featureAnnotationSetQuery,
        colNameOpt = "fieldname"
      )
      setnames(features, "GeneSymbol", "gene_symbol")
    } else {
      # Get annotation from flat file b/c otherwise don't know order
      # NOTE: For ImmSig studies, this means that summaries use the latest
      # annotation even though that was not used in the manuscript for summarizing.
      features <- data.frame(
        FeatureId = self$cache[[cache_name]]$gene_symbol,
        gene_symbol = self$cache[[cache_name]]$gene_symbol
      )
    }

    # update cache$gematrices with correct fasId
    self$cache$GE_matrices$featureset[self$cache$GE_matrices$name == matrixName] <- annoSetId

    # Change ge_matrices$annotation
    self$cache$GE_matrices$annotation[self$cache$GE_matrices$name == matrixName] <- annotation

    # push features to cache
    self$cache[[paste0("featureset_", annoSetId)]] <- features
  }
)


# Constructs a expression set by matrix
ISCon$set(
  which = "private",
  name = ".constructExpressionSet",
  value = function(matrixName, outputType, annotation) {

    cache_name <- .setCacheName(matrixName, outputType)
    esetName <- paste0(cache_name, "_eset")

    # expression matrix
    message("Constructing ExpressionSet")
    matrix <- self$cache[[cache_name]]

    # pheno
    runID <- self$cache$GE_matrices[name == matrixName, rowid]
    bs <- colnames(matrix)[ grep("^BS\\d{6}$", colnames(matrix))]
    pheno_filter <- makeFilter(c("Run",
                                 "EQUAL",
                                 runID),
                               c("biosample_accession",
                                 "IN",
                                 paste(bs, collapse = ";")))

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
    # Need cohortType for modules
    pheno <- pheno[, colnames(pheno) %in% c("biosample_accession",
                                            "participant_id",
                                            "cohort_type",
                                            "cohort",
                                            "study_time_collected",
                                            "study_time_collected_unit",
                                            "exposure_material_reported",
                                            "exposure_process_preferred")]

    # ensure same order as GEM rownames
    rownames(pheno) <- pheno$biosample_accession
    order <- names(self$cache[[cache_name]])
    order <- order[-grep("feature_id|gene_symbol|X|V1", order)]
    order <- order[ order != "BS694717.1" ] # rm SDY212 dup for the moment
    pheno <- pheno[match(order, row.names(pheno)),]

    # handling multiple timepoints per subject
    dups <- colnames(matrix)[duplicated(colnames(matrix))]
    if (length(dups) > 0) {
      matrix <- data.table(matrix)
      for (dup in dups) {
        dupIdx <- grep(dup, colnames(matrix))
        newNames <- paste0(dup, 1:length(dupIdx))
        setnames(matrix, dupIdx, newNames)
        eval(substitute(matrix[, `:=`(dup,
                                      rowMeans(matrix[, dupIdx, with = FALSE]))],
                        list(dup = dup)))
        eval(substitute(matrix[, `:=`(newNames, NULL)], list(newNames = newNames)))
      }
      if (config$verbose) {
        warning("The matrix contains subjects with multiple measures per timepoint. Averaging expression values.")
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
      features <- self$cache[[paste0("featureset_", annoSetId)]][, c("FeatureId","gene_symbol")]

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
    fasInfo <- .getLKtbl(con = self,
                         schema = "Microarray",
                         query = "FeatureAnnotationSet")
    gemx <- self$cache$GE_matrices
    fasId <- gemx$featureset[ gemx$name == matrixName & gemx$outputtype == outputType ]
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
      featureAnnotationSet = fasInfo$Name)

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

# Set the cache name of expression matrix by output type
.setCacheName <- function(matrixName, outputType) {
  outputSuffix <- switch(
    outputType,
    "summary" = "_sum",
    "normalized" = "_norm",
    "raw" = "_raw"
  )

  paste0(matrixName, outputSuffix)
}


# Combine EMs and output only genes available in all EMs.
.combineEMs <- function(EMlist) {
  message("Combining ExpressionSets")

  fds <- lapply(EMlist, function(x) {
    droplevels(data.table(fData(x)))
  })

  fd <- Reduce(f = function(x, y) {merge(x, y, by = c("FeatureId", "gene_symbol"))}, fds)

  EMlist <- lapply(EMlist, "[", as.character(fd$FeatureId))

  for (i in 1:length(EMlist)) {
    fData(EMlist[[i]]) <- fd
  }

  Reduce(f = combine, EMlist)
}
