#' @include ImmuneSpace.R
NULL


#' @importFrom RCurl getCurlHandle curlPerform basicTextGatherer
.ISCon$methods(
  downloadMatrix=function(x, summary = FALSE){
    cache_name <- paste0(x, ifelse(summary, "_sum", ""))
#     if(is.null(data_cache[[x]])){
    if(is.null(data_cache[[cache_name]])){
      if(nrow(subset(data_cache[[constants$matrices]], name%in%x)) == 0){
        stop(sprintf("No matrix %s in study\n", x))
      }
      summary <- ifelse(summary, ".summary", "")
      link <- URLdecode(file.path(gsub("http:","https:", gsub("/$","",config$labkey.url.base)),
                                  "_webdav", gsub("^/","",config$labkey.url.path),
                                  "@files/analysis/exprs_matrices",
                                  paste0(x, ".tsv", summary)))
      localpath <- .self$.localStudyPath(link)
      if(.self$.isRunningLocally(localpath)){
        fl <- localpath
        message("Reading local matrix")
#         data_cache[[x]]<<-fread(fl,header=TRUE)
        data_cache[[cache_name]] <<- fread(fl, header = TRUE)
      }else{
        opts <- config$curlOptions
        opts$netrc <- 1L
        #opts$httpauth <- 1L
        handle <- getCurlHandle(.opts=opts)
        h <- basicTextGatherer()
        message("Downloading matrix..")
        curlPerform(url = link, curl = handle, writefunction = h$update)
        fl <- tempfile()
        write(h$value(), file = fl)
        EM <- fread(fl, header = TRUE)
        if(nrow(EM) == 0){
          stop("The downloaded matrix has 0 rows. Something went wrong")
        }
        data_cache[[cache_name]] <<-EM
        file.remove(fl)
      }
      
    }else{
      data_cache[[cache_name]]
    }
  }
)

.ISCon$methods(
  GeneExpressionFeatures=function(matrix_name,summary=FALSE){
    cache_name <- paste0(matrix_name, ifelse(summary, "_sum", ""))
    if(!matrix_name %in% data_cache[[constants$matrices]][, name]){
      stop("Invalid gene expression matrix name");
    }
    annotation_set_id <- .self$.getFeatureId(matrix_name)
    if(is.null(data_cache[[.self$.mungeFeatureId(annotation_set_id)]])){
      if(!summary){
        message("Downloading Features..")
        featureAnnotationSetQuery=sprintf("SELECT * from FeatureAnnotation where FeatureAnnotationSetId='%s';",annotation_set_id);
        features<-labkey.executeSql(config$labkey.url.base,config$labkey.url.path,schemaName = "Microarray",sql = featureAnnotationSetQuery ,colNameOpt = "fieldname")
        setnames(features, "GeneSymbol", "gene_symbol")
      }else{
        features<-data.frame(FeatureId=data_cache[[cache_name]][,gene_symbol],
                             gene_symbol=data_cache[[cache_name]][,gene_symbol])
      }
      data_cache[[.self$.mungeFeatureId(annotation_set_id)]]<<-features
    }
  }
)

.ISCon$methods(
  ConstructExpressionSet=function(matrix_name, summary){
    cache_name <- paste0(matrix_name, ifelse(summary, "_sum", ""))
    #matrix
    message("Constructing ExpressionSet")
    matrix <- data_cache[[cache_name]]
    #features
    features <- data_cache[[.self$.mungeFeatureId(.self$.getFeatureId(matrix_name))]][,c("FeatureId","gene_symbol")]
    #inputs
    pheno_filter <- makeFilter(c("Run/DataOutputs/Name", "EQUAL", paste0(matrix_name, ".tsv")),
                               c("Biosample/biosample_accession", "IN", paste(colnames(matrix), collapse = ";")))
    pheno <- unique(data.table(labkey.selectRows(
      config$labkey.url.base, config$labkey.url.path,
      "assay.ExpressionMatrix.matrix", "InputSamples", "gene_expression_matrices",
      colNameOpt = "rname", colFilter = pheno_filter)))
    setnames(pheno, colnames(pheno), gsub("^biosample_", "", .self$.munge(colnames(pheno))))
    pheno <- pheno[, list(biosample_accession, subject_accession, arm_name,
                          study_time_collected, study_time_collected_unit)]
    
    if(summary){
      fdata <- data.frame(FeatureId = matrix$gene_symbol, gene_symbol = matrix$gene_symbol, row.names = matrix$gene_symbol)
      fdata <- AnnotatedDataFrame(fdata)
    } else{
      try(setnames(matrix," ","FeatureId"),silent=TRUE)
      fdata <- data.table(FeatureId = as.character(matrix$FeatureId))
      fdata <- merge(fdata, features, by = "FeatureId", all.x = TRUE)
      fdata <- as.data.frame(fdata)
      rownames(fdata) <- fdata$FeatureId
      fdata <- AnnotatedDataFrame(fdata)
    }
    pheno <- data.frame(pheno)
    rownames(pheno) <- pheno$biosample_accession
    pheno<-pheno[colnames(matrix)[-1L],]
    ad_pheno<-AnnotatedDataFrame(data=pheno)
    es<-ExpressionSet(assayData=as.matrix(matrix[,-1L,with=FALSE]),phenoData=ad_pheno,featureData=fdata)
    data_cache[[cache_name]]<<-es
  }
)

#' @title get Gene Expression Matrix
#'
#' @description 
#' Downloads a normalized gene expression matrix from ImmuneSpace.
#'
#' con$getGEMatrix*x = NULL, cohort = NULL, summary = FALSE, reload = FALSE)
#'
#' @param x  A \code{character}. The name of the Gene Expression Matrix to download.
#' @param cohort A \code{character}. The name of a cohort that has an associated
#'  gene expression matrix. Note that if cohort is not NULL, then x will be ignored.
#' @param summary A \code{logical}. If set to TRUE. Downloads a matrix with
#'  expression averaged by gene symbol. By default the feature level matrix is
#'  downloaded.
#' @param reload A \code{logical}. If set to TRUE, the matrix will be downloaded
#'  again, even if a cached copy exist in the \code{ImmuneSpaceConnection} object.
#'  
#' @return An \code{ExpressionSet}.
#' @aliases getGEMatrix
#' @name ImmuneSpaceConnection_getGEMatrix
#' @examples
#' labkey.url.base = "https://www.immunespace.org"
#' labkey.url.path = "/Studies/SDY269"
#' sdy269 <- CreateConnection("SDY269")
#' sdy269$getGEMatrix("TIV_2008")
.ISCon$methods(
  getGEMatrix=function(x = NULL, cohort = NULL, summary = FALSE, reload=FALSE){
    "Downloads a s normalized gene expression matrix from ImmuneSpace.\n
    `x': A `character'. The name of the gene expression matrix to download.\n
    `cohort': A `character'. The name of a cohort that has an associated gene
    expression matrix. Note that if `cohort' isn't NULL, then `x' is ignored."
    cohort_name <- cohort #can't use cohort = cohort in d.t
    if(!is.null(cohort_name)){
      if(all(cohort_name %in% data_cache$GE_matrices$cohort)){
        x <- data_cache$GE_matrices[cohort == cohort_name, name]
      } else{
        validCohorts <- data_cache$GE_matrices[, cohort]
        stop(paste("No expression matrix for the given cohort.",
                   "Valid cohorts:", paste(validCohorts, collapse = ", ")))
      }
    }
    cache_name <- paste0(x, ifelse(summary, "_sum", ""))
    if(length(x) > 1){
      data_cache[cache_name] <<- NULL
      lapply(x, downloadMatrix, summary)
      lapply(x, GeneExpressionFeatures,summary)
      lapply(x, ConstructExpressionSet, summary)
      return(Reduce(f=combine, data_cache[cache_name]))
    } else{
      if (cache_name %in% names(data_cache) && !reload) {
        data_cache[[cache_name]]
      }
      else {
        data_cache[[cache_name]] <<- NULL
        downloadMatrix(x, summary)
        GeneExpressionFeatures(x, summary)
        ConstructExpressionSet(x, summary)
        data_cache[[cache_name]]
      }
    }
  }
)

.ISCon$methods(
  .getFeatureId=function(matrix_name){
    subset(data_cache[[constants$matrices]],name%in%matrix_name)[, featureset]
  }
)

.ISCon$methods(
  .mungeFeatureId=function(annotation_set_id){
    return(sprintf("featureset_%s",annotation_set_id))
  }
)