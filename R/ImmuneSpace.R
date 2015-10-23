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
        "gene_expression", "gene_expression_analysis_results",  colNameOpt = "caption", ...))
    setnames(GEAR, .self$.munge(colnames(GEAR)))
    return(GEAR)
  }
)

# Get HAI response at peak immunogenicity
# 
# Peak immunogenicity is defined as the timepoint with the maximum average fold
# change to baseline. It is calculated per cohort.
# 
# @return A \code{data.table} with columns participant_id, response and arm name
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
    hai <- hai[, list(response=log2(max(response)), name), by="participant_id"]   
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

