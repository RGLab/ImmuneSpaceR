#'@docType package
#'@title A Thin Wrapper Around ImmuneSpace.
#'@description ImmuneSpaceR provides a convenient API for accessing data sets
#'within the ImmuneSpace database.
#'
#'@details Uses the Rlabkey package to connect to ImmuneSpace. Implements
#'caching, and convenient methods for accessing data sets.
#'
#'@name ImmuneSpaceR-package
#'@aliases ImmuneSpaceR
#'@author Greg Finak
#'@import data.table Rlabkey methods Biobase
NULL

.CreateConnection = function(study = NULL
                             , labkey.url.base
                             , labkey.url.path
                             , labkey.user.email
                             , curlOptions
                             , verbose
                             , ...){
  labkey.url.path<-try(get("labkey.url.path",.GlobalEnv),silent=TRUE)
  if(inherits(labkey.url.path,"try-error")){
    if(is.null(study)){
      stop("study cannot be NULL")
    }
    labkey.url.path <- paste0("/Studies/",study)
  }else if(!is.null(study)){
    labkey.url.path <- file.path(dirname(labkey.url.path),study)
  }
  config <- list(labkey.url.base = labkey.url.base,
                  labkey.url.path = labkey.url.path,
                  labkey.user.email = labkey.user.email,
                  curlOptions = curlOptions,
                  verbose = verbose)
  
  .ISCon(config = config)
}

#'@name ImmuneSpaceConnection
#'@aliases ImmuneSpaceConnection-class
#'@aliases ImmuneSpace
#'@rdname ImmuneSpaceConnection-class
#'@docType class
#'@title The ImmuneSpaceConnection class
#'@description Instantiate this class to access a study
#'@details
#' Uses global variables \code{labkey.url.base}, and \code{labkey.url.path}, to
#' access a study. \code{labkey.url.base} should be
#' \code{https://www.immunespace.org/}. \code{labkey.url.path} should be
#' \code{/Studies/studyname}, where 'studyname' is the accession number of the
#' study.
#' The ImmunespaceConnection will initialize itself, and look for a
#' \code{.netrc} file in \code{"~/"} the user's home directory. The
#' \code{.netrc} file should contain a \code{machine}, \code{login}, and
#' \code{password} entry to allow access to ImmuneSpace, where \code{machine} is
#' the host name like "www.immunespace.org".
#' 
#'@seealso \code{\link{ImmuneSpaceR-package}} \code{\link{ImmuneSpaceConnection_getGEMatrix}}  \code{\link{ImmuneSpaceConnection_getDataset}}  \code{\link{ImmuneSpaceConnection_listDatasets}}
#'@exportClass ImmuneSpaceConnection
#'@examples
#'labkey.url.base <- "https://www.immunespace.org"
#'labkey.url.path <- "/Studies/SDY269"
#'labkey.user.email <- 'gfinak at fhcrc.org'
#'sdy269 <- CreateConnection("SDY269")
#'sdy269
#'@return An instance of an ImmuneSpaceConnection for a study in `labkey.url.path`
.ISCon <- setRefClass(Class = "ImmuneSpaceConnection",
            fields = list(study = "character", config="list",
                          available_datasets = "data.table",
                          data_cache="list",constants="list")
)

#' @importFrom gtools mixedsort
.ISCon$methods(
  checkStudy=function(verbose = FALSE){
    validStudies <- mixedsort(grep("^SDY", basename(lsFolders(getSession(config$labkey.url.base, "Studies"))), value = TRUE))
    req_study <- basename(config$labkey.url.path)
    if(!req_study %in% validStudies){
      if(!verbose){
        stop(paste0(req_study, " is not a valid study"))
      } else{
        stop(paste0(req_study, " is not a valid study\nValid studies: ",
                    paste(validStudies, collapse=", ")))
      }
    }
  }
)

.ISCon$methods(
  getAvailableDataSets=function(){
    if(length(available_datasets)==0){
      dataset_filter <- makeFilter(c("showbydefault", "EQUAL", TRUE))
      df <- labkey.selectRows(baseUrl = config$labkey.url.base
                        , config$labkey.url.path
                        , schemaName = "study"
                        , queryName = "DataSets"
                        , colFilter = dataset_filter)
      
      available_datasets <<- data.table(df)[,list(Label,Name,Description,`Key Property Name`)]
    }
  }
)

.ISCon$methods(
  .munge=function(x){
    new <- tolower(gsub(" ","_",basename(x)))
    idx <- which(duplicated(new) | duplicated(new, fromLast = TRUE))
    if(length(idx)>0)
      new[idx] <- .munge(gsub("(.*)/.*$", "\\1", x[idx]))
    return(new)
  }
)
.ISCon$methods(
  GeneExpressionInputs=function(){
    if(!is.null(data_cache[[constants$matrix_inputs]])){
      data_cache[[constants$matrix_inputs]]
    }else{
      ge<-data.table(labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "assay.ExpressionMatrix.matrix",queryName = "InputSamples",colNameOpt = "fieldname",viewName = "gene_expression_matrices",showHidden=TRUE))
      setnames(ge,.self$.munge(colnames(ge)))
      data_cache[[constants$matrix_inputs]]<<-ge
    }
  }
)

.ISCon$methods(
  GeneExpressionFeatures=function(matrix_name,summary=FALSE){
    cache_name <- paste0(matrix_name, ifelse(summary, "_sum", ""))
#     if(!any((data_cache[[constants$matrices]][,"name"] %in% cache_name))){
    if(!matrix_name %in% data_cache[[constants$matrices]][, name]){
      stop("Invalid gene expression matrix name");
    }
    annotation_set_id<-.self$.getFeatureId(matrix_name)
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
  GeneExpressionMatrices=function(verbose = FALSE){
    if(!is.null(data_cache[[constants$matrices]])){
      data_cache[[constants$matrices]]
    }else{
      if(verbose){
        ge <- try(data.table(
          labkey.selectRows(baseUrl = config$labkey.url.base,
                            config$labkey.url.path,
                            schemaName = "assay.ExpressionMatrix.matrix",
                            queryName = "Runs",
                            colNameOpt = "fieldname",
                            showHidden = TRUE,
                            viewName = "expression_matrices")),
        )
      } else {
        suppressWarnings(
          ge <- try(data.table(
            labkey.selectRows(baseUrl = config$labkey.url.base,
                              config$labkey.url.path,
                              schemaName = "assay.ExpressionMatrix.matrix",
                              queryName = "Runs",
                              colNameOpt = "fieldname",
                              showHidden = TRUE,
                              viewName = "expression_matrices")),
          )
        )
      }
      if(inherits(ge, "try-error") || nrow(ge) == 0){
        #No assay or no runs
        message("No gene expression data")
        data_cache[[constants$matrices]] <<- NULL
      } else{
        setnames(ge,.self$.munge(colnames(ge)))
        data_cache[[constants$matrices]]<<-ge
      }
    }
    return(data_cache[[constants$matrices]])
  }
)

#' @importFrom RCurl getCurlHandle curlPerform basicTextGatherer
.ISCon$methods(
  downloadMatrix=function(x, summary = FALSE){
    cache_name <- paste0(x, ifelse(summary, "_sum", ""))
#     if(is.null(data_cache[[x]])){
    if(is.null(data_cache[[cache_name]])){
      if(nrow(subset(data_cache[[constants$matrices]],name%in%x))==0){
        stop(sprintf("No matrix %s in study\n",x))
      }
      summary <- ifelse(summary, ".summary", "")
      link <- URLdecode(file.path(gsub("http:","https:", gsub("/$","",config$labkey.url.base)),
                                  "_webdav", gsub("^/","",config$labkey.url.path),
                                  "@files/analysis/exprs_matrices",
                                  paste0(x, ".tsv", summary)))
      localpath <- .self$.localStudyPath(link)
      if(.self$.isRunningLocally(localpath)){
        fl<-localpath
        message("Reading local matrix")
#         data_cache[[x]]<<-fread(fl,header=TRUE)
        data_cache[[cache_name]]<<-fread(fl,header=TRUE)
      }else{
        opts <- config$curlOptions
        opts$netrc <- 1L
        #opts$httpauth <- 1L
        handle<-getCurlHandle(.opts=opts)
        h<-basicTextGatherer()
        message("Downloading matrix..")
        curlPerform(url=link,curl=handle,writefunction=h$update)
        fl<-tempfile()
        write(h$value(),file=fl)
        EM <- fread(fl,header=TRUE)
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
  ConstructExpressionSet=function(matrix_name, summary){
    cache_name <- paste0(matrix_name, ifelse(summary, "_sum", ""))
    #matrix
    message("Constructing ExpressionSet")
    matrix<-data_cache[[cache_name]]
    #features
    features<-data_cache[[.self$.mungeFeatureId(.self$.getFeatureId(matrix_name))]][,c("FeatureId","gene_symbol")]
    #inputs
    #pheno<-unique(data_cache[[constants$matrix_inputs]][biosample_accession %in% colnames(matrix)][, c("biosample_accession","subject_accession","arm_name","study_time_collected", "study_time_collected_unit")])
    pheno_filter <- makeFilter(c("Run/DataOutputs/Name", "EQUAL", paste0(matrix_name, ".tsv")),
                               c("Biosample/biosample_accession", "IN", paste(colnames(matrix), collapse = ";")))
    pheno <- unique(data.table(labkey.selectRows(
      config$labkey.url.base, config$labkey.url.path,
      "assay.ExpressionMatrix.matrix", "InputSamples", "gene_expression_matrices",
      colNameOpt = "rname", colFilter = pheno_filter)))
    setnames(pheno, colnames(pheno), gsub("^biosample_", "", .self$.munge(colnames(pheno))))
    pheno <- pheno[, list(biosample_accession, subject_accession, arm_name,
                          study_time_collected, study_time_collected_unit)]
    #colnames(pheno) <- gsub("^biosample_", "", .self$.munge(colnames(pheno)))
    
    
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
      #setkey(matrix,FeatureId)
      #rownames(features)<-features$FeatureId
      #features<-features[matrix$FeatureId,]#order feature info
      #fdata <- AnnotatedDataFrame(features)
    }
    pheno <- data.frame(pheno)
    rownames(pheno) <- pheno$biosample_accession
    pheno<-pheno[colnames(matrix)[-1L],]
    ad_pheno<-AnnotatedDataFrame(data=pheno)
    es<-ExpressionSet(assayData=as.matrix(matrix[,-1L,with=FALSE]),phenoData=ad_pheno,featureData=fdata)
    data_cache[[cache_name]]<<-es
  }
)

#'@title get Gene Expression Matrix
#'@aliases getGEMatrix
#'@param x \code{"character"} name of the Gene Expression Matrix
#'@details Returns an `ExpressionSet` from the matrix named 'x', downloads it if it is not already cached.
#'@return an \code{ExpressionSet}
#'@name ImmuneSpaceConnection_getGEMatrix
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269$getGEMatrix("TIV_2008")
.ISCon$methods(
  getGEMatrix=function(x = NULL, cohort = NULL, summary = FALSE, reload=FALSE){
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
.ISCon$methods(
  .isRunningLocally=function(path){
    file.exists(path)
  }
)
.ISCon$methods(
  .localStudyPath=function(urlpath){
    #LOCALPATH<-"/shared/silo_researcher/Gottardo_R/immunespace"
    LOCALPATH<-"/share/files/"
    PRODUCTION_HOST<-"www.immunespace.org"
    TEST_HOST<-"test.immunespace.org"
    PRODUCTION_PATH<-""
    #if(grepl(PRODUCTION_HOST, urlpath)){
    #  PROCESS<-PRODUCTION_PATH
    #}else if(grepl(TEST_HOST, urlpath)){
    #  PROCESS <- ""
    #}else{
    #  stop("Can't determine if we are running on immunespace (production) or posey (staging)")
    #}
    #gsub(file.path(gsub("/$","",config$labkey.url.base), "_webdav"), file.path(LOCALPATH,PROCESS), urlpath)
    gsub(file.path(gsub("/$","",config$labkey.url.base), "_webdav"), file.path(LOCALPATH), urlpath)
  }
)

#'@title get a dataset
#'@aliases getDataset
#'@param x A \code{character}. The name of the dataset
#'@param original_view A \code{logical}. If set to TRUE, download the ImmPort
#' view. Else, download the default grid view. Note: Once data is cached,
#' changing value of this argument won't have effect on the subsequent calls
#' unless \code{reload} is set to 'TRUE'.
#'@param ... Arguments to be passed to the underlying \code{labkey.selectRows}.
#'@param reload A \code{logical}. Clear the cache. If set to TRUE, download the 
#'  dataset, whether a cached version exist or not.
#'@details
#'Returns the dataset named 'x', downloads it if it is not already cached. Note
#'that if additional arguments (...) are passed, the dataset will not be
#'reloaded.
#'@return a \code{data.table}
#'@name ImmuneSpaceConnection_getDataset
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269$getDataset("hai")
.ISCon$methods(
    getDataset = function(x, original_view = FALSE, reload=FALSE, ...){
      if(nrow(available_datasets[Name%in%x])==0){
        wstring <- paste0(study, " has invalid data set: ",x)
        if(config$verbose){
          wstring <- paste0(wstring, "\n",
                            "Vali datasets for ", study, ": ",
                            paste(available_datasets$Name, collapse = ", "), ".")
        }
        warning(wstring)
        NULL
      }else{
        cache_name <- paste0(x, ifelse(original_view, "_full", ""))
        if(!is.null(data_cache[[cache_name]]) & !reload & length(list(...)) == 0){
          data_cache[[cache_name]]
        }else{
          viewName <- NULL
          if(original_view){
            viewName <- "full"
          }
          
          data_cache[[cache_name]] <<- data.table(labkey.selectRows(baseUrl = config$labkey.url.base
                                                           ,config$labkey.url.path
                                                           ,schemaName = "study"
                                                           , queryName = x
                                                           , viewName = viewName
                                                           , colNameOpt = "caption"
                                                           , ...)
          )
          setnames(data_cache[[cache_name]],
                   .self$.munge(colnames(data_cache[[cache_name]])))
          data_cache[[cache_name]]
        }
      }
    })


#'@title list available datasets
#'@aliases listDatasets
#'@details Prints the names of the available datasets
#'@return Doesn't return anything, just prints to console.
#'@name ImmuneSpaceConnection_listDatasets
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269$listDatasets()
.ISCon$methods(
    listDatasets=function(){
      cat("datasets\n")
      
      for(i in 1:nrow(available_datasets)){
        cat(sprintf("\t%s\n",available_datasets[i,Name]))
      }
      if(!is.null(data_cache[[constants$matrices]])){
        cat("Expression Matrices\n")
        for(i in 1:nrow(data_cache[[constants$matrices]])){
          cat(sprintf("\t%s\n",data_cache[[constants$matrices]][i, name]))
        }
      }
    })


#'@title list available gene expression analysis
#'@aliases listGEAnalysis
#'@details Prints the table of differential expression analysis
#'@return A \code{data.frame}. The list of gene expression analysis.
#'@name ImmuneSpaceConnection_listGEAnalysis
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269$listGEAnalysis()
.ISCon$methods(
    listGEAnalysis = function(){
      GEA <- data.table(labkey.selectRows(config$labkey.url.base,
                                          config$labkey.url.path,
                                          "gene_expression",
                                          "gene_expression_analysis",
                                          colNameOpt = "rname"))
      return(GEA)
    })

#' Get gene expression analysis
#' 
#' Download the result of a Gene epxression analysis experiment from the
#' gene_expression schema.
#' 
#' @param ... A \code{list} of arguments to be passed to \code{labkey.selectRows}.
#'  
#' @return A \code{data.table} containing the requested gene expression analysis
#'  results.
#' 
#' @import data.table
#' @aliases getGEAnalysis
#' @name ImmuneSpaceConnection_getGEAnalysis
.ISCon$methods(
  getGEAnalysis = function(...){
    GEAR <- data.table(labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
        "gene_expression", "gene_expression_analysis_results",  colNameOpt = "fieldname", ...))
    setnames(GEAR, .self$.munge(colnames(GEAR)))
    return(GEAR)
  }
)

#' Get HAI response at peak immunogenicity
#' 
#' Peak immunogenicity is defined as the timepoint with the maximum average fold
#' change to baseline. It is calculated per cohort.
#' 
#' @return A \code{data.table} with columns subject_accession, response and arm name
#' 
#' @aliases getHAIResponse
#' @name ImmuneSpaceConnection_getHAIResponse
.ISCon$methods(
  getHAIResponse = function(reload){
    hai <- .self$getDataset("hai", reload = TRUE)
    hai <- hai[, list(name, study_time_collected, study_time_collected_unit,
                        response = value_reported/mean(value_reported[study_time_collected<=0])),
                 by = "virus_strain,subject_accession"]
    hai <- hai[, mr := mean(response), by="study_time_collected"]
    hai <- hai[, ma := max(mr), by = "name"]
    peak <- unique(hai[mr ==ma, list(study_time_collected, name)])
    hai <- merge(hai, peak, by=c("study_time_collected", "name"))
    hai <- hai[, list(response=log2(max(response)), name), by="subject_accession"]   
    return(hai)    
  }
)

.ISCon$methods(
  clear_cache = function(){
    data_cache[grep("^GE", names(data_cache), invert = TRUE)] <<- NULL
  }
)
.ISCon$methods(
    show=function(){
      cat(sprintf("Immunespace Connection to study %s\n",study))
      cat(sprintf("URL: %s\n",file.path(gsub("/$","",config$labkey.url.base),gsub("^/","",config$labkey.url.path))))
      cat(sprintf("User: %s\n",config$labkey.user.email))
      cat("Available datasets\n")
      for(i in 1:nrow(available_datasets)){
        cat(sprintf("\t%s\n",available_datasets[i,Name]))
      }
      if(!is.null(data_cache[[constants$matrices]])){
        cat("Expression Matrices\n")
        for(i in 1:nrow(data_cache[[constants$matrices]])){
          cat(sprintf("\t%s\n",data_cache[[constants$matrices]][i, name]))
        }
      }
    }
)

#'@title get Gene Expression Files
#'@param files A \code{character}. The name of the files to download
#'@param destdir A \code{character}. The destination directory
#'
#'@details getGEFiles makes calls to base function \code{download.file} which
#' in turn makes system calls to curl.
#'@return An \code{integer} vector of the same length as the \code{files}
#' parameter. See the Value section of \code{?download.file} for more
#' information.
#'
#'@name ImmuneSpaceConnection_getGEMatrix
#'@aliases getGEFiles
#'
#'@examples
#'#Downloads CEL files in the current directory
#'sdy269<-CreateConnection("SDY269")
#'sdy269$getGEFiles(c("GSM733843.CEL", "GSM733844.CEL"))
.ISCon$methods(
  getGEFiles=function(files, destdir = "."){
    links <- paste0(config$labkey.url.base, "/_webdav/",
                    config$labkey.url.path,
                    "/%40files/rawdata/gene_expression/", files)
    sapply(links, function(x){
      download.file(url = links[1], destfile = file.path(destdir, basename(x)),
                    method = "curl", extra = "-n")
    })
  }
)

.ISCon$methods(
  initialize=function(..., config = NULL){
    
    #invoke the default init routine in case it needs to be invoked 
    #(e.g. when using $new(object) to construct the new object based on the exiting object)
    callSuper(...)
    
    constants <<- list(matrices="GE_matrices",matrix_inputs="GE_inputs")
    
    if(!is.null(config))
      config <<- config

    study <<- basename(config$labkey.url.path)
    if(config$verbose){
      checkStudy(config$verbose)
    }
    
    getAvailableDataSets()

    gematrices_success <- GeneExpressionMatrices(verbose = FALSE)
    
  }
)
