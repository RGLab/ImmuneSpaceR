#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------

# List the available gene expression matrices
ISCon$set(
  which = "public",
  name = "listGEMatrices",
  value = function(verbose = FALSE) {
    ## HELPERS
    ..getData <- function() {
      try(
        .getLKtbl(
          con = self,
          schema = "assay.ExpressionMatrix.matrix",
          query = "Runs",
          colNameOpt = "fieldname",
          viewName = "expression_matrices"
        ),
        silent = TRUE
      )
    }


    ## MAIN
    if (!is.null(self$cache[[private$.constants$matrices]])) {
      self$cache[[private$.constants$matrices]]
    } else {
      if (verbose) {
        ge <- ..getData()
      } else {
        ge <- suppressWarnings(..getData())
      }

      if (inherits(ge, "try-error") || nrow(ge) == 0) {
        # No assay or no runs
        message("No gene expression data")
        self$cache[[private$.constants$matrices]] <- NULL
      } else {
        # adding cols to allow for getGEMatrix() to update
        ge[, annotation := ""]
        ge[, outputType := ""]
        setnames(ge, private$.munge(colnames(ge)))
        self$cache[[private$.constants$matrices]] <- ge
      }
    }

    return(self$cache[[private$.constants$matrices]])
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
                   cohort = NULL,
                   outputType = "summary",
                   annotation = "latest",
                   reload = FALSE,
                   verbose = FALSE) {
    if (outputType == "summary" & annotation == "ImmSig") {
      stop("Not able to provide summary eSets for ImmSig annotated studies. Please use
           'raw' as outputType with ImmSig studies.")
    }

    cohort_name <- cohort # can't use cohort = cohort in d.t
    if (!is.null(cohort_name)) {
      if (all(cohort_name %in% self$cache$GE_matrices$cohort)) {
        matrixName <- self$cache$GE_matrices[cohort %in% cohort_name, name]
      } else {
        validCohorts <- self$cache$GE_matrices[, cohort]
        stop(paste("No expression matrix for the given cohort.",
                   "Valid cohorts:", paste(validCohorts, collapse = ", ")))
      }
    }

    cache_name <- .setCacheName(matrixName, outputType)
    esetName <- paste0(cache_name, "_eset")

    # length(x) > 1 means multiple cohorts
    if (length(matrixName) > 1) {
      lapply(matrixName, private$.downloadMatrix, outputType, annotation, reload)
      lapply(matrixName, private$.getGEFeatures, outputType, annotation, reload)
      lapply(matrixName, private$.constructExpressionSet, outputType, annotation)
      ret <- .combineEMs(self$cache[esetName])
      if (dim(ret)[[1]] == 0) {
        # No features shared
        warn <- "Returned ExpressionSet has 0 rows. No feature is shared across the selected runs or cohorts."
        if (outputType != "summary") {
          warn <- paste(warn,
                        "Try outputType = 'summary' to merge matrices by gene symbol.")
        }
        warning(warn)
      }

      if (verbose == TRUE) {
        print(Biobase::experimentData(ret))
      }

      return(ret)

    } else {
      if (esetName %in% names(self$cache) & !reload) {
        message(paste0("returning ", esetName, " from cache"))
      } else {
        self$cache[[esetName]] <- NULL
        private$.downloadMatrix(matrixName, outputType, annotation)
        private$.getGEFeatures(matrixName, outputType, annotation)
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

    if (annotation == "ImmSig"){
      fileSuffix <- ".immsig"
    }else{
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
        paste0(matrixName, ".tsv", fileSuffix)
      )
    )

    localpath <- private$.localStudyPath(link)
    if (private$.isRunningLocally(localpath)) {
      message("Reading local matrix")
      self$cache[[cache_name]] <- read.table(
        localpath,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    } else {
      opts <- self$config$curlOptions
      opts$options$netrc <- 1L

      message("Downloading matrix..")
      fl <- tempfile()
      GET(url = link, config = opts, write_disk(fl))

      # fread does not read correctly
      EM <- read.table(
        fl,
        header = TRUE,
        sep = "\t",
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
    if (status == annotation & reload != TRUE) {
      message(paste0("returning ", annotation, " annotation from cache"))
      return()
    }

    runs <- labkey.selectRows(
      baseUrl = self$config$labkey.url.base,
      folderPath = self$config$labkey.url.path,
      schemaName = "Assay.ExpressionMatrix.Matrix",
      queryName = "Runs",
      showHidden = TRUE
    )

    getOrigFasId <- function(config, matrixName) {
      # Get annoSet based on name of FeatureAnnotationSet + "_orig" tag
      faSets <- labkey.selectRows(
        baseUrl = config$labkey.url.base,
        folderPath = config$labkey.url.path,
        schemaName = "Microarray",
        queryName = "FeatureAnnotationSet",
        showHidden = TRUE
      )

      fasId <- runs$`Feature Annotation Set`[runs$Name == matrixName]
      fasNm <- faSets$Name[faSets$`Row Id` == fasId]

      # '_orig' annotation is the feature annotation copy generated by
      # updateAnno when run on the server.  In the case of SDY67, the copy is
      # not generated because it is has 'ImmSig'.  Therefore, no change to the
      # fasNm is needed.
      if (annotation == "default" & !grepl("_orig|ImmSig", fasNm)) {
        fasNm <- paste0(fasNm, "_orig")

        # This situation occurred with SDY400 where annotation was added after original
        # was generated.
      } else if (annotation == "latest" & grepl("_orig", fasNm, fixed = TRUE)) {
        # This situation occurred with SDY400 where annotation was added after
        # original was generated.
        fasNm <- gsub("_orig", "", fasNm, fixed = TRUE)
      }
      annoSetId <- faSets$`Row Id`[faSets$Name == fasNm]
    }

    # ImmuneSignatures data needs mapping from when microarray was read, not
    # 'original' when IS matrices were created.
    if (annotation == "ImmSig") {
      faSets <- labkey.selectRows(
        baseUrl = self$config$labkey.url.base,
        folderPath = self$config$labkey.url.path,
        schemaName = "Microarray",
        queryName = "FeatureAnnotationSet",
        showHidden = TRUE
      )

      sdy <- tolower(gsub("/Studies/", "", self$config$labkey.url.path))
      annoSetId <- faSets$`Row Id`[faSets$Name == paste0("ImmSig_", sdy)]
    } else if (annotation == "default") {
      annoSetId <- getOrigFasId(self$config, matrixName)
    } else if (annotation == "latest") {
      annoSetId <- runs$`Feature Annotation Set`[ runs$Name == matrixName]
    }

    if (outputType != "summary") {
      message("Downloading Features..")
      featureAnnotationSetQuery = sprintf("SELECT * from FeatureAnnotation
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
    pheno_filter <- makeFilter(c("Run",
                                 "EQUAL",
                                 runID),
                               c("Biosample/biosample_accession",
                                 "IN",
                                 paste(colnames(matrix), collapse = ";")))

    pheno <- unique(
      .getLKtbl(
        con = self,
        schema = "study",
        query = "HM_InputSamplesQuery",
        containerFilter = "CurrentAndSubfolders",
        colNameOpt = "caption",
        colFilter = pheno_filter,
        showHidden = FALSE
      )
    )

    setnames(pheno, private$.munge(colnames(pheno)))
    pheno <- data.frame(pheno, stringsAsFactors = FALSE)
    pheno <- pheno[, colnames(pheno) %in% c("biosample_accession",
                                            "participant_id",
                                            "cohort",
                                            "study_time_collected",
                                            "study_time_collected_unit")]
    rownames(pheno) <- pheno$biosample_accession

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

      # exprs and fData must match
      rownames(fdata) <- rownames(matrix) <- matrix$gene_symbol
    } else {
      annoSetId <- self$cache$GE_matrices$featureset[self$cache$GE_matrices$name == matrixName]

      features <- self$cache[[paste0("featureset_", annoSetId)]][, c("FeatureId","gene_symbol")]

      colnames(matrix)[[which(colnames(matrix) %in% c(" ", "V1", "X", "feature_id")) ]] <- "FeatureId"
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

    # SDY212 has dbl biosample that is removed for ImmSig, but needs to be
    # present for normalization, so needs to be included in eSet!
    if (runID == 469) {
      pheno["BS694717.1", ] <- pheno[pheno$biosample_accession == "BS694717", ]
      pheno$biosample_accession[rownames(pheno) == "BS694717.1"] <- "BS694717.1"
    }

    # Prep Eset and push
    # NOTES: At project level, InputSamples may be filtered
    matrix <- data.frame(matrix, stringsAsFactors = FALSE) # for when on rsT / rsP
    exprs <- matrix[, colnames(matrix) %in% pheno$biosample_accession] # rms gene_symbol!
    pheno <- pheno[colnames(exprs), ]

    # add processing information for user
    isRNA <- self$study %in% c("SDY888","SDY224","SDY67")
    fasInfo <- .getLKtbl(con = self,
                         schema = "Microarray",
                         query = "FeatureAnnotationSet")
    gemx <- self$cache$GE_matrices
    fasId <- gemx$featureset[ gemx$name == matrixName & gemx$outputtype == outputType ]
    fasInfo <- fasInfo[ match(fasId, fasInfo$`Row Id`)]
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
