.ISCon$methods(
  .munge=function(x){
    new <- tolower(gsub(" ","_",basename(x)))
    idx <- which(duplicated(new) | duplicated(new, fromLast = TRUE))
    if(length(idx)>0)
      new[idx] <- .munge(gsub("(.*)/.*$", "\\1", x[idx]))
    return(new)
  }
)

# Returns TRUE if the connection is at project level ("/Studies")
.ISCon$methods(
  .isProject=function()
    if(config$labkey.url.path == "/Studies/"){
      TRUE
    } else{
      FALSE
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

.ISCon$methods(
    listDatasets=function(which = c("datasets", "expression")){
      "List the datasets available in the study or studies of the connection."
      
      if("datasets" %in% which){
        cat("datasets\n")
        for(i in 1:nrow(available_datasets)){
          cat(sprintf("\t%s\n",available_datasets[i,Name]))
        }
      }
      if("expression" %in% which){
        if(!is.null(data_cache[[constants$matrices]])){
          cat("Expression Matrices\n")
          for(i in 1:nrow(data_cache[[constants$matrices]])){
            cat(sprintf("\t%s\n",data_cache[[constants$matrices]][i, name]))
          }
        }
      }
    })


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

.ISCon$methods(
  getGEAnalysis = function(...){
    "Downloads data from the gene expression analysis results table.\n
    '...': A list of arguments to be passed to labkey.selectRows."
    GEAR <- data.table(labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
        "gene_expression", "DGEA_filteredGEAR",  "DGEAR", colNameOpt = "caption", ...))
    setnames(GEAR, .self$.munge(colnames(GEAR)))
    return(GEAR)
  }
)

.ISCon$methods(
  clear_cache = function(){
  "Clear the data_cache. Remove downloaded datasets and expression matrices."
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

# Returns a named logical where TRUE marks files that are accessible.
# This function is used for administrative purposes to check that the flat files
# are properly loaded and accessible to the users.
#' @importFrom RCurl url.exists
#' @import parallel
.ISCon$methods(
  .test_files=function(what = c("gene_expression_files", "fcs_sample_files", "protocol")){
    
    labkey.url.base <- "https://www.immunespace.org"
    ret <- list()
    what <- tolower(what)
    
    
    
    for(i in what){
      df <- .self$getDataset(i, original_view = TRUE)
      df <- df[!is.na(file_info_name)]
      
      # handle gene expr / fcs separately bc similar link paths
      if(i == "gene_expression_files" | i == "fcs_sample_files"){ 
        df <- unique(df[, list(study_accession, file_info_name)])
        
        numrow <- nrow(df)
        link_text <- ""
        
        if(i == "gene_expression_files"){
          link_text <- "gene_expression"
        }else if( i == "fcs_sample_files"){
          link_text <- "flow_cytometry"
        }
        
        links <- paste0(labkey.url.base, "/_webdav/", "/Studies/", 
                        gef$study_accession, "/%40files/rawdata/",link_text,
                        sapply(gef$file_info_name, URLencode))
        
        all_res <- unlist(mclapply(links, link_test, mc.cores = detectCores()))
        bound_res <- cbind(gef,all_res)

        print(paste0(length(all_res), "/", numrow, " ", i, " with valid links."))
        
        ret[[i]] <- bound_res
        
        #handle protocols alone 
      }else{
        #TODO: Get list of all studies
        # make list of links
        # mcapply to test links like above
        
        study_string <- strsplit(.self$config$labkey.url.path, "/")
        study <- study_string[[1]][3]
        res <- link_test(paste0(labkey.url.base, "/_webdav/Studies/",
                                 study, "/%40files/protocols/", folder,
                                 "_protocol.zip"))
        ret[[i]] <- res
        }
      }
      
    return(ret)
  }
)

link_test <- function(link){
  
  info <- GET(link)
  status <- info$status_code
  return(status)
}

#    if("protocol" %in% what){
#if(.self$.isProject()){
#  message("Project level, all protocols will be checked.")
#  # Get study list and check em all. Bonus points for returning named logical.
#} else{
# https://www.immunespace.org/_webdav/Studies/SDY112/%40files/protocols/SDY112_protocol.zip
#}