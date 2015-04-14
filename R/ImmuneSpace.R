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
#'@seealso \code{\link{ImmuneSpaceR-package}} 
#'\code{\link{ImmuneSpaceConnection_getGEMatrix}}
#'\code{\link{ImmuneSpaceConnection_getDataset}}
#'\code{\link{ImmuneSpaceConnection_listDatasets}}
#'\code{\link{ImmuneSpaceConnection_getGEAnalysis}}
#'\code{\link{ImmuneSpaceConnection_listGEAnalysis}}
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

.ISCon$methods(
  .isRunningLocally=function(path){
    file.exists(path)
  }
)
.ISCon$methods(
  .localStudyPath=function(urlpath){
    LOCALPATH <- "/share/files/"
    PRODUCTION_HOST <- "www.immunespace.org"
    TEST_HOST <- "test.immunespace.org"
    gsub(file.path(gsub("/$","",config$labkey.url.base), "_webdav"), file.path(LOCALPATH), urlpath)
  }
)

#' @title Get a dataset
#' 
#' @description 
#' Downloads a dataset and cache the result in the connection object.
#' 
#' con$getDataset(x, original_view = FALSE, reload = FALSE, ...)
#' 
#' @param x A \code{character}. The name of the dataset.
#' @param original_view A \code{logical}. If set to TRUE, download the ImmPort
#'  view. Else, download the default grid view. Note: Once data is cached,
#'  changing value of this argument won't have effect on the subsequent calls
#'  unless \code{reload} is set to 'TRUE'.
#' @param reload A \code{logical}. Clear the cache. If set to TRUE, download the 
#'  dataset, whether a cached version exist or not.
#' @param ... Arguments to be passed to the underlying \code{labkey.selectRows}.
#' @details
#' Returns the dataset named 'x', downloads it if it is not already cached. Note
#' that if additional arguments (...) are passed, the dataset will be reloaded,
#' even if a cached copy exist.
#' 
#' @return a \code{data.table}
#' @name ImmuneSpaceConnection_getDataset
#' @examples
#' labkey.url.base = "https://www.immunespace.org"
#' labkey.url.path = "/Studies/SDY269"
#' sdy269 <- CreateConnection("SDY269")
#' sdy269$getDataset("hai")
.ISCon$methods(
    getDataset = function(x, original_view = FALSE, reload=FALSE, ...){
      "Get a dataset form the connection\n
      original_view: A logical. If set tot TRUE, download the ImmPort view.
      Else, download the default grid view.\n
      reload: A logical. Clear the cache. If set to TRUE, download the dataset,
      whether a cached version exist or not."
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


#' List available datasets
#' 
#' @description 
#' List the datasets available in the study or studies of the connection.
#'  
#' con$listDatasets()
#'  
#' @details Prints the names of the available datasets
#' @return Doesn't return anything, just prints to console.
#' @name ImmuneSpaceConnection_listDatasets
#' @examples
#' labkey.url.base = "https://www.immunespace.org"
#' labkey.url.path = "/Studies/SDY269"
#' sdy269 <- CreateConnection("SDY269")
#' sdy269$listDatasets()
.ISCon$methods(
    listDatasets=function(){
      "List the datasets available in the study or studies of the connection."
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


#' @title List available gene expression analysis
#'
#' @description 
#' List the available gene expression analysis runs.
#' 
#' con$listGEAnalysis() 
#'
#' @details Prints the table of differential expression analysis
#' @return A \code{data.frame}. The list of gene expression analysis.
#' @name ImmuneSpaceConnection_listGEAnalysis
#' @examples
#' labkey.url.base="https://www.immunespace.org"
#' labkey.url.path="/Studies/SDY269"
#' sdy269<-CreateConnection("SDY269")
#' sdy269$listGEAnalysis()
.ISCon$methods(
    listGEAnalysis = function(){
      "List available gene expression analysis for the connection."
      GEA <- data.table(labkey.selectRows(config$labkey.url.base,
                                          config$labkey.url.path,
                                          "gene_expression",
                                          "gene_expression_analysis",
                                          colNameOpt = "rname"))
      return(GEA)
    })

#' @title Get gene expression analysis
#' 
#' @description 
#' Download the result of a Gene epxression analysis experiment from the
#' gene_expression schema.
#' 
#' con$getGEAnalysis(...)
#' 
#' @param ... A \code{list} of arguments to be passed to \code{labkey.selectRows}.
#'  
#' @return A \code{data.table} containing the requested gene expression analysis
#'  results.
#' 
#' @import data.table
#' @seealso labkey.selectRows
#' @name ImmuneSpaceConnection_getGEAnalysis
#' @examples
#' 
#' labkey.url.base = "https://www.immunespace.org"
#' labkey.url.path = "/Studies/SDY269"
#' sdy269 <- CreateConnection("SDY269")
#' sdy269$listGEAnalysis()
.ISCon$methods(
  getGEAnalysis = function(...){
    "Downloads data from the gene expression analysis results table"
    GEAR <- data.table(labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
        "gene_expression", "gene_expression_analysis_results",  colNameOpt = "fieldname", ...))
    setnames(GEAR, .self$.munge(colnames(GEAR)))
    return(GEAR)
  }
)

# Get HAI response at peak immunogenicity
# 
# Peak immunogenicity is defined as the timepoint with the maximum average fold
# change to baseline. It is calculated per cohort.
# 
# @return A \code{data.table} with columns subject_accession, response and arm name
# 
# @aliases getHAIResponse
# @name ImmuneSpaceConnection_getHAIResponse
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

#' @title get Gene Expression Files
#' @param files A \code{character}. The name of the files to download
#' @param destdir A \code{character}. The destination directory
#' 
#' @details getGEFiles makes calls to base function \code{download.file} which
#'  in turn makes system calls to curl.
#' @return An \code{integer} vector of the same length as the \code{files}
#'  parameter. See the Value section of \code{?download.file} for more
#'  information.
#' 
#' @name ImmuneSpaceConnection_getGEFiles
#' 
#' @examples
#' #Downloads CEL files in the current directory
#' sdy269<-CreateConnection("SDY269")
#' sdy269$getGEFiles(c("GSM733843.CEL", "GSM733844.CEL"))
.ISCon$methods(
  getGEFiles=function(files, destdir = "."){
    "Download gene expression raw data files.\n
    files: A character. Filenames as shown on the gene_expression_files dataset.\n
    destdir: A character. The loacal path to store the downloaded files."
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
