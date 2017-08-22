#' @include ImmuneSpace.R
NULL

# helper
setCacheName <- function(matrixName, outputType){
  outputSuffix <- switch(outputType,
                         "summary" = "_sum",
                         "normalized" = "_norm",
                         "raw" = "_raw")
  return(paste0(matrixName, outputSuffix))
}

#' @importFrom RCurl getCurlHandle curlPerform basicTextGatherer
#' @importFrom preprocessCore normalize.quantiles
.ISCon$methods(
  downloadMatrix = function(matrixName, outputType = "summary", annotation = "latest"){
  
    cache_name <- setCacheName(matrixName, outputType)

    # check if study has matrices
    if( nrow( subset( data_cache[[ constants$matrices ]], 
                      name %in% matrixName) ) == 0 ){
      stop(sprintf("No matrix %s in study\n", matrixName))
    }

    # check if data in data_cache corresponds to current request
    # if it does, then no download needed
    if(outputType == "summary"){
      status <- data_cache$GE_matrices$comments[ data_cache$GE_matrices$name == matrixName ]
      if( !is.na(status) ){
        if( status == "Using default Annotation" & annotation == "default" ){
          message("default GE Matrix already in data_cache")
          return()
        }else if( status == "Using Updated Annotation" & annotation == "latest" ){
          message("Updated GE matrix already in data_cache")
          return()
        }
      }
    }

    fileSuffix <- if(outputType == "summary"){
      switch(annotation, "latest" = ".summary", "default" = ".summary.orig")
    }else{
      switch(outputType, "normalized" = "", "raw" = ".raw")
    }

    path <- if( config$labkey.url.path == "/Studies/" ){
      paste0("/Studies/", data_cache$GE_matrices[name == matrixName, folder], "/")
    }else{
      gsub("^/", "", config$labkey.url.path)
    }

    link <- URLdecode(file.path(gsub("http:",
                                     "https:",
                                     gsub("/$","",config$labkey.url.base)),
                                "_webdav",
                                path,
                                "@files/analysis/exprs_matrices",
                                paste0(matrixName, ".tsv", fileSuffix)))

    localpath <- .self$.localStudyPath(link)
    if( .self$.isRunningLocally(localpath) ){
      message("Reading local matrix")
      data_cache[[cache_name]] <<- read.table(localpath,
                                                    header = TRUE,
                                                    sep = "\t",
                                                    stringsAsFactors = FALSE)
    }else{
      opts <- config$curlOptions
      opts$netrc <- 1L
      handle <- getCurlHandle(.opts = opts)
      h <- basicTextGatherer()
      message("Downloading matrix..")
      curlPerform(url = link, curl = handle, writefunction = h$update)
      fl <- tempfile()
      write(h$value(), file = fl)
      EM <- read.table(fl,
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE) # fread does not read correctly
      if(nrow(EM) == 0){
        stop("The downloaded matrix has 0 rows. Something went wrong.")
      }
      data_cache[[cache_name]] <<- EM
      file.remove(fl)
    }
  }
)

.ISCon$methods(
  GeneExpressionFeatures = function(matrixName, 
                                    outputType = "summary", 
                                    annotation = "latest"){
    
    cache_name <- setCacheName(matrixName, outputType)

    if( !matrixName %in% data_cache[[constants$matrices]]$name ){
      stop("Invalid gene expression matrix name");
    }

    runs <- labkey.selectRows(baseUrl = config$labkey.url.base,
                              folderPath = config$labkey.url.path,
                              schemaName = "Assay.ExpressionMatrix.Matrix",
                              queryName = "Runs",
                              showHidden = TRUE)

    getOrigFasId <- function(config, matrixName){
      # Get annoSet based on name of FeatureAnnotationSet + "_orig" tag
      faSets <- labkey.selectRows(baseUrl = config$labkey.url.base,
                                  folderPath = config$labkey.url.path,
                                  schemaName = "Microarray",
                                  queryName = "FeatureAnnotationSet",
                                  showHidden = TRUE)

      fasId <- runs$`Feature Annotation Set`[ runs$Name == matrixName]
      fasNm <- faSets$Name[ faSets$`Row Id` == fasId]

      # Impt to use current fasNm if currAnno == T and also in case where
      # currAnno == F, but the matrix was made using the original annotation
      # AFTER a newer version of annotation was generated. E.g. SDY400 where
      # EH had to make matrices for ImmSig project using orig annotation at probe lvl.
      # Second logic makes it feasible to use the current annotation with such a study
      # if currAnno == T
      if(annotation == "default" & !grepl("_orig", fasNm, fixed = T)){
        fasNm <- paste0(fasNm, "_orig")
      }else if(annotation == "latest" & grepl("_orig", fasNm, fixed = T)){
        fasNm <- gsub("_orig", "", fasNm, fixed = T)
      }

      annoSetId <- faSets$`Row Id`[ faSets$Name == fasNm ]
    }
    
    # ImmuneSignatures data needs mapping from when microarray was was read, not
    # 'original' when IS matrices were created.
    if( annotation == "ImmSig" ){
      sdy <- gsub("/Studies/", "", config$labkey.url.path)
      annoSetId <- switch(sdy,
                          "SDY212" = 36,
                          "SDY63" = 37,
                          "SDY404" = 38,
                          "SDY400" = 39,
                          "SDY67" = 40)
    }else if( annotation == "default"){
      annoSetId <- getOrigFasId(config, matrixName)
    }else if( annotation == "latest"){
      annoSetId <- runs$`Feature Annotation Set`[ runs$Name == matrixName]
    }
    
    if(outputType != "summary"){
      message("Downloading Features..")
      featureAnnotationSetQuery = sprintf("SELECT * from FeatureAnnotation
                                        where FeatureAnnotationSetId='%s';",
                                          annoSetId);
      features <- labkey.executeSql(baseUrl = config$labkey.url.base,
                                    folderPath = config$labkey.url.path,
                                    schemaName = "Microarray",
                                    sql = featureAnnotationSetQuery,
                                    colNameOpt = "fieldname")
      setnames(features, "GeneSymbol", "gene_symbol")
    }else{
      # Get annotation from flat file b/c otherwise don't know order
      # NOTE: For ImmSig studies, this means that summaries use the latest
      # annotation even though that was not used in the manuscript for summarizing.
      features <- data.frame(FeatureId = data_cache[[cache_name]]$gene_symbol,
                             gene_symbol = data_cache[[cache_name]]$gene_symbol)
    }

    # update the data_cache$gematrices with correct fasId
    data_cache$GE_matrices$featureset[ data_cache$GE_matrices$name == matrixName ] <<- annoSetId

    # make comment to clarify
    data_cache$GE_matrices$comments[ data_cache$GE_matrices$name == matrixName ] <<- ifelse(
      annotation == "latest",
     "Using Updated Annotation",
     "Using Default Annotation")

    # push features to cache
    data_cache[[ paste0("featureset_", annoSetId) ]] <<- features
  }
)

.ISCon$methods(
  ConstructExpressionSet = function(matrixName, outputType){
    cache_name <- setCacheName(matrixName, outputType)

    # expression matrix
    message("Constructing ExpressionSet")
    matrix <- data_cache[[ cache_name ]]

    # Average expr vals for subs with multiple time points
    dups <- colnames(matrix)[ duplicated(colnames(matrix)) ]
    if(length(dups) > 0){
      matrix <- data.table(matrix)
      for(dup in dups){
        dupIdx <- grep(dup, colnames(matrix))
        newNames <- paste0(dup, 1:length(dupIdx))
        setnames(matrix, dupIdx, newNames)
        eval(substitute(matrix[, `:=`(dup,
                                      rowMeans(matrix[, dupIdx, with = FALSE]))],
                        list(dup = dup)))
        eval(substitute(matrix[, `:=`(newNames, NULL)], list(newNames = newNames)))
      }
      if(config$verbose){
        warning("The matrix contains subjects with multiple measures per timepoint. Averaging expression values.")
      }
    }

    # gene features
    if(outputType == "summary"){
      fdata <- data.frame(FeatureId = matrix$gene_symbol,
                          gene_symbol = matrix$gene_symbol)
      rownames(fdata) <- rownames(matrix) <- matrix$gene_symbol # exprs and fData must match
    }else{
      annoSetId <- data_cache$GE_matrices$featureset[ data_cache$GE_matrices$name == matrixName]
      
      features <- data_cache[[ paste0("featureset_", annoSetId)]][,c("FeatureId","gene_symbol")]
      
      colnames(matrix)[[ which(colnames(matrix) %in% c(" ", "V1", "X", "feature_id")) ]] <- "FeatureId"
      fdata <- data.frame(FeatureId = as.character(matrix$FeatureId), 
                          stringsAsFactors = FALSE)
      
      fdata <- merge(fdata, features, by = "FeatureId", all.x = TRUE)
      
      # exprs and fData must match
      rownames(fdata) <- fdata$FeatureId 
      matrix <- matrix[ order(match(matrix$FeatureId, fdata$FeatureId)), ]
      rownames(matrix) <- matrix$FeatureId
    }

    # pheno
    runID <- data_cache$GE_matrices[name == matrixName, rowid]
    pheno_filter <- makeFilter(c("Run",
                                 "EQUAL",
                                 runID),
                               c("Biosample/biosample_accession",
                                 "IN",
                                 paste(colnames(matrix), collapse = ";")))

    pheno <- unique(labkey.selectRows(baseUrl = config$labkey.url.base,
                                      folderPath = config$labkey.url.path,
                                      schemaName = "study",
                                      queryName = "HM_InputSamplesQuerySnapshot",
                                      containerFilter = "CurrentAndSubfolders",
                                      colNameOpt = "caption",
                                      colFilter = pheno_filter))
    

    colnames(pheno) <- sapply(colnames(pheno), .munge)
    keep <- c("biosample_accession",
              "participant_id",
              "cohort",
              "study_time_collected",
              "study_time_collected_unit")
    pheno <- pheno[ , colnames(pheno) %in% keep ]
    rownames(pheno) <- pheno$biosample_accession
    
    # SDY212 has dbl biosample that is removed for ImmSig, but needs to be 
    # present for normalization, so needs to be included in eSet!
    if(runID == 469){
      pheno[ "BS694717.1", ] <- pheno[ pheno$biosample_accession == "BS694717", ]
      pheno$biosample_accession[ rownames(pheno) == "BS694717.1"] <- "BS694717.1"
    }
    
    # Prep Eset and push
    # NOTES: At project level, InputSamples may be filtered
    matrix <- data.frame(matrix) # for when on rsT / rsP
    posNames <- pheno$biosample_accession
    
    if(runID == 469){ posNames <- c(posNames, "BS694717.1") } # for SDY212
    exprs <- matrix[, colnames(matrix) %in% posNames] # rms gene_symbol!
    pheno <- pheno[ colnames(exprs), ]

    data_cache[[cache_name]] <<- ExpressionSet(assayData = as.matrix(exprs),
                                               phenoData = AnnotatedDataFrame(pheno),
                                               featureData = AnnotatedDataFrame(fdata))
  }
)

# Downloads a normalized gene expression matrix from ImmuneSpace.
.ISCon$methods(
  getGEMatrix = function(matrixName = NULL,
                         cohort = NULL,
                         outputType = "summary",
                         annotation = "latest",
                         reload = FALSE){

    "Downloads a normalized gene expression matrix from ImmuneSpace.\n
    `x': A `character'. The name of the gene expression matrix to download.\n
    `cohort': A `character'. The name of a cohort that has an associated gene
    expression matrix. Note that if `cohort' isn't NULL, then `x' is ignored.
    `summary': A `logical'. If set to TRUE, Downloads a matrix with expression
    averaged by gene symbol.
    `reload': A `logical'. If set to TRUE, the matrix will be downloaded again,
    even if a cached cop exist in the ImmuneSpaceConnection object."

    cohort_name <- cohort #can't use cohort = cohort in d.t
    if( !is.null(cohort_name) ){
      if( all(cohort_name %in% data_cache$GE_matrices$cohort) ){
        matrixName <- data_cache$GE_matrices[cohort %in% cohort_name, name]
      } else{
        validCohorts <- data_cache$GE_matrices[, cohort]
        stop(paste("No expression matrix for the given cohort.",
                   "Valid cohorts:", paste(validCohorts, collapse = ", ")))
      }
    }

    cache_name <- setCacheName(matrixName, outputType)

    # length(x) > 1 means multiple cohorts
    if( length(matrixName) > 1 ){
      data_cache[cache_name] <<- NULL
      data_cache$GE_matrices$comments <<- NA
      lapply(matrixName, downloadMatrix, outputType, annotation)
      lapply(matrixName, GeneExpressionFeatures, outputType, annotation)
      lapply(matrixName, ConstructExpressionSet, outputType)
      ret <- .combineEMs(data_cache[cache_name])
      if(dim(ret)[[1]] == 0){
        # No features shared
        warn <- "Returned ExpressionSet has 0 rows. No feature is shared across the selected runs or cohorts."
        if(outputType != "summary"){
          warn <- paste(warn, 
                        "Try outputType = 'summary' to merge matrices by gene symbol.")
        }
        warning(warn)
      }

      return(ret)

    }else{

      .getNewEset = function(matrixName, outputType, annotation){
        data_cache[[cache_name]] <<- NULL
        downloadMatrix(matrixName, outputType, annotation)
        GeneExpressionFeatures(matrixName, outputType, annotation)
        ConstructExpressionSet(matrixName, outputType)
      }

      if( cache_name %in% names(data_cache) & !reload ){
        # check if data in data_cache corresponds to current request
        # if it does, then no download needed
        if(outputType == "summary"){
          status <- data_cache$GE_matrices$comments[ data_cache$GE_matrices$name == matrixName ]
          if( !is.na(status) ){
            if( status == "Using Default Annotation" & annotation == "default" ){
              message("Returning Default expressionSet from data_cache")
            }else if( status == "Using Updated Annotation" & annotation == "latest" ){
              message("Returning Updated expressionSet from data_cache")
            }else{
              .getNewEset(matrixName, outputType, annotation)
            }
          }
        }

      }else{
        .getNewEset(matrixName, outputType, annotation)
      }
      
      return(data_cache[[cache_name]])
    }
  }
)

# Combine EMs and output only genes available in all EMs.
#' @importFrom Biobase fData
.combineEMs <- function(EMlist){
  message("Combining ExpressionSets")
  fds <- lapply(EMlist, function(x){ droplevels(data.table(fData(x))) })
  fd <- Reduce(f = function(x, y){ merge(x, y, by = c("FeatureId", "gene_symbol"))}, fds)
  EMlist <- lapply(EMlist, "[", as.character(fd$FeatureId))
  for(i in 1:length(EMlist)){ fData(EMlist[[i]]) <- fd}
  res <- Reduce(f = combine, EMlist)
}

# Add treatment information to the phenoData of an expression matrix available in the connection object.
.ISCon$methods(
  addTreatment=function(matrixName = NULL){
    "Add treatment information to the phenoData of an expression matrix
    available in the connection object.\n
    x: A character. The name of a expression matrix that has been downloaded 
    from the connection."

    if(is.null(matrixName) | !matrixName %in% names(data_cache)){
      stop(paste(matrixName, "is not a valid expression matrix."))
    } 

    bsFilter <- makeFilter(c("biosample_accession",
                             "IN",
                             paste(pData(data_cache[[matrixName]])$biosample_accession, collapse = ";")))
    bs2es <- data.table(labkey.selectRows(baseUrl = config$labkey.url.base,
                                          folderPath = config$labkey.url.path,
                                          schemaName = "immport",
                                          queryName = "expsample_2_biosample",
                                          colFilter = bsFilter,
                                          colNameOpt = "rname"))

    esFilter <- makeFilter(c("expsample_accession",
                             "IN",
                             paste(bs2es$expsample_accession, collapse = ";")))
    es2trt <- data.table(labkey.selectRows(baseUrl = config$labkey.url.base,
                                           folderPath = config$labkey.url.path,
                                           schemaName = "immport",
                                           queryName = "expsample_2_treatment",
                                           colFilter = esFilter,
                                           colNameOpt = "rname"))

    trtFilter <- makeFilter(c("treatment_accession",
                              "IN",
                              paste(es2trt$treatment_accession, collapse = ";")))
    trt <- data.table(labkey.selectRows(baseUrl = config$labkey.url.base,
                                        folderPath = config$labkey.url.path,
                                        schemaName = "immport",
                                        queryName = "treatment",
                                        colFilter = trtFilter,
                                        colNameOpt = "rname"))

    bs2trt <- merge(bs2es, es2trt, by = "expsample_accession")
    bs2trt <- merge(bs2trt, trt, by = "treatment_accession")

    pData(data_cache[[matrixName]])$treatment <<- bs2trt[match(pData(data_cache[[matrixName]])$biosample_accession, biosample_accession), name]

    return(data_cache[[matrixName]])
  }
)

.ISCon$methods(
  .getFeatureId=function(matrixName){
    subset(data_cache[[constants$matrices]],name%in%matrixName)[, featureset]
  }
)

.ISCon$methods(
  .mungeFeatureId=function(annotation_set_id){
    return(sprintf("featureset_%s",annotation_set_id))
  }
)


#' @importFrom Biobase pData sampleNames
.ISCon$methods(
  EMNames=function(EM = NULL, colType = "participant_id"){
    "Change the sampleNames of an ExpressionSet fetched by getGEMatrix using the
    information in the phenodData slot.\n
    x: An ExpressionSet, as returned by getGEMatrix.\n
    colType: A character. The type of column names. Valid options are 'expsample_accession'
    and 'participant_id'."

    if(is.null(EM) | !is(EM, "ExpressionSet")){
      stop("EM should be a valid ExpressionSet, as returned by getGEMatrix")
    }

    if(!all(grepl("^BS", sampleNames(EM)))){
      stop("All sampleNames should be biosample_accession, as returned by getGEMatrix")
    }

    pd <- data.table(pData(EM))
    colType <- gsub("_.*$", "", tolower(colType))

    if(colType == "expsample"){
      bsFilter <- makeFilter(c("biosample_accession",
                               "IN",
                               paste(pd$biosample_accession, collapse = ";")))
      bs2es <- data.table(labkey.selectRows(baseUrl = config$labkey.url.base,
                                            folderPath = config$labkey.url.path,
                                            schemaName = "immport",
                                            queryName = "expsample_2_biosample",
                                            colFilter = bsFilter,
                                            colNameOpt = "rname"))
      pd <- merge(pd,
                  bs2es[, list(biosample_accession, expsample_accession)],
                  by = "biosample_accession")

      es <- pd[match(sampleNames(EM), pd$biosample_accession), expsample_accession]
      
      sampleNames(EM) <- pData(EM)$expsample_accession <- es

    } else if(colType %in% c("participant", "subject")){
      pd <- pd[, nID := paste0(participant_id,
                               "_",
                               tolower(substr(study_time_collected_unit, 1, 1)),
                               study_time_collected)]
      sampleNames(EM) <- pd[ match(sampleNames(EM), pd$biosample_accession), nID]

    } else if(colType == "biosample"){
      warning("Nothing done, the column names should already be biosample_accession numbers.")

    } else {
      stop("colType should be one of 'expsample_accession', 'biosample_accession', 'participant_id'.")
    }

    return(EM)
  }
)