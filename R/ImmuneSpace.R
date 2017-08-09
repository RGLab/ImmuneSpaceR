.ISCon$methods(
  .munge=function(x){
    new <- tolower(gsub(" ","_",basename(x)))
    idx <- which(duplicated(new) | duplicated(new, fromLast = TRUE))
    if( length(idx) > 0 ){ new[idx] <- .munge(gsub("(.*)/.*$", "\\1", x[idx])) }
    return(new)
  }
)

# Returns TRUE if the connection is at project level ("/Studies")
.ISCon$methods(
  .isProject=function()
    res <- config$labkey.url.path == "/Studies/"
)

.ISCon$methods(
  GeneExpressionInputs=function(){
    if(!is.null(data_cache[[constants$matrix_inputs]])){
      data_cache[[constants$matrix_inputs]]
    }else{
      ge <- .getLKtbl(con = .self, 
                      schema = "assay.Expressionmatrix.matrix",
                      query = "InputSamples",
                      viewName = "gene_expression_matrices",
                      colNameOpt = "fieldname")
      setnames(ge,.self$.munge(colnames(ge)))
      data_cache[[constants$matrix_inputs]] <<- ge
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
    }
)

.ISCon$methods(
    listGEAnalysis = function(){
      "List available gene expression analysis for the connection."
      GEA <- .getLKtbl(con = .self, 
                       schema = "gene_expression",
                       query = "gene_expression_analysis",
                       showHidden = FALSE,
                       colNameOpt = "rname")
    }
)

.ISCon$methods(
  getGEAnalysis = function(...){
    "Downloads data from the gene expression analysis results table.\n
    '...': A list of arguments to be passed to labkey.selectRows."
    GEAR <- .getLKtbl(con = .self, 
                      schema = "gene_expression",
                      query = "DGEA_filteredGEAR",
                      viewName = "DGEAR",
                      colNameOpt = "caption",
                      ...)
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
    "Display information about the object."
      cat(sprintf("Immunespace Connection to study %s\n",study))
      
      cat(sprintf("URL: %s\n",
                  file.path(gsub("/$","",config$labkey.url.base),
                            gsub("^/","",config$labkey.url.path)))
          )
      
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
      download.file(url = links[1], 
                    destfile = file.path(destdir, basename(x)),
                    method = "curl", 
                    extra = "-n")
    })
  }
)

# Returns a list of data frames where TRUE in file_exists column marks files that are accessible.
# This function is used for administrative purposes to check that the flat files
# are properly loaded and accessible to the users.
#' @importFrom RCurl getURL url.exists
#' @importFrom rjson fromJSON
#' @importFrom parallel mclapply detectCores
.ISCon$methods(
  .test_files=function(what = c("gene_expression_files", "fcs_sample_files", "protocols", "ge_matrices")){
    list_files <- function(link){
      response <- NULL
      if (url.exists(url = link, netrc = TRUE)){
        response_json <- fromJSON(getURL(url = link, netrc = TRUE))
        response <- unlist(lapply(response_json$files, function(x){ x$text } ))
      }
      response
    }
    
    check_links <- function (dataset, folder){
      res <- data.frame(file_info_name = NULL, 
                        study_accession = NULL, 
                        file_link = NULL, 
                        file_exists = NULL, 
                        stringsAsFactors = FALSE)
      
      if (dataset %in% .self$available_datasets$Name){
        temp <- .self$getDataset(dataset, original_view = TRUE)
        temp <- temp[!is.na(file_info_name)]
        temp <- unique(temp[, list(study_accession, file_info_name)])
        
        file_link <- paste0(config$labkey.url.base, "/_webdav/Studies/", 
                            temp$study_accession, "/%40files/rawdata/", folder, "/", 
                            sapply(temp$file_info_name, URLencode))
        
        folder_link <- paste0(config$labkey.url.base, "/_webdav/Studies/", 
                              unique(temp$study_accession), "/%40files/rawdata/", folder, "?method=JSON")
        
        file_list <- unlist(mclapply(folder_link, list_files, mc.cores = detectCores()))
        
        file_exists <- temp$file_info_name %in% file_list
        
        res <- data.frame(study = temp$study_accession, 
                          file_link = file_link, 
                          file_exists = file_exists, 
                          stringsAsFactors = FALSE) 
        
        print(paste0(sum(res$file_exists), "/", nrow(res), " ", dataset, " with valid links."))
      }
      res
    }
    
    ret <- list()
    what <- tolower(what)
    
    if("gene_expression_files" %in% what){
      ret$gene_expression_files <- check_links("gene_expression_files", "gene_expression")
    }
    if("fcs_sample_files" %in% what){
      ret$fcs_sample_files <- check_links("fcs_sample_files", "flow_cytometry")
    }
    if("protocols" %in% what){
      if(.self$.isProject()){
        folders_list <- labkey.getFolders(baseUrl = config$labkey.url.base, 
                                          folderPath = "/Studies/")
        folders <- folders_list[, 1]
        folders <- folders[!folders %in% c("SDY_template","Studies")]
      } else{
        folders <- basename(config$labkey.url.path)
      }
      
      file_link <- paste0(config$labkey.url.base, "/_webdav/Studies/", folders, 
                              "/%40files/protocols/", folders, "_protocol.zip")
      file_exists <- unlist(mclapply(file_link, url.exists, netrc = TRUE, mc.cores = detectCores()))
      
      res <- data.frame(study = folders, 
                        file_link = file_link, 
                        file_exists = file_exists, 
                        stringsAsFactors = FALSE)      
      
      print(paste0(sum(res$file_exists), "/", nrow(res), " protocols with valid links."))
      ret$protocols <- res
    }
    if ("ge_matrices" %in% what){
      matrix_queries <- labkey.getQueries(baseUrl = config$labkey.url.base, 
                                          folderPath = config$labkey.url.path, 
                                          schemaName = "assay.ExpressionMatrix.matrix")
      
      if ("OutputDatas" %in% matrix_queries$queryName) {

        ge <-.getLKtbl(con = .self, 
                       schema = "assay.ExpressionMatrix.matrix", 
                       query = "OutputDatas", 
                       colNameOpt = "rname", 
                       viewName = "links")
        
        output <- lapply(ge[4], function(x){
                          gsub("@", 
                               "%40", 
                               gsub("file:/share/files", 
                                    paste0(config$labkey.url.base, "/_webdav"), 
                                    x))} 
                        )
        
        file_exists <- unlist(mclapply(output$data_datafileurl, 
                                       url.exists, 
                                       netrc = TRUE, 
                                       mc.cores = detectCores()))
        
        res <- data.frame(file_link = output$data_datafileurl, 
                          file_exists = file_exists, 
                          stringsAsFactors = FALSE)
        
        print(paste0(sum(res$file_exists), "/", nrow(res), " ge_matrices with valid links."))
        
      } else {
        res <- data.frame(file_link = NULL, 
                          file_exists = NULL, 
                          stringsAsFactors = FALSE)
      }
      
      ret$ge_matrices <- res
    }
    return(ret)
  }
)

###############################################################################
# -------------- PARTICIPANT FILTERING METHODS -------------------------------#
###############################################################################

# Ensure only one valid user email is returned or possible errors handled 
.validateUser <- function(con){
  # First Check for apiKey and netrc file
  api <- Rlabkey:::ifApiKey()
  sink(tempfile())
  validNetrc <- tryCatch({
    check_netrc()
  }, error = function(e){
    return(NULL)
  })
  sink()
  
  # Case 1: if apiKey, but no Netrc (possible in Integrated RStudio session UI)
  # get user email from global environment, which is set by labkey.init and
  # handle unexpected case of labkey.user.email not being found.
  if( !is.null(api) & is.null(validNetrc) ){
    user <- try( get("labkey.user.email"), silent = TRUE )
    if( inherits(labkey.user.email, "try-error") ){
      stop("labkey.user.email not found, please set")
    }
    
    # Case 2: valid netrc file (with or without ApiKey)
    # To mimic LK.get() method - pull first login and use that one.
    # Need to differentiate between test / prod too.
  }else if( !is.null(validNetrc) ){
    netrc <- readLines(validNetrc)
    machine <- gsub("https://", "machine ", con$config$labkey.url.base)
    loc <- which(netrc == machine)[1]
    login <- netrc[ loc + 1 ]
    user <- gsub("login ", "", login)
  }
  
  return(user)
}

# less verbose wrapper for LK.selectRows calls
.getLKtbl <- function(con, schema, query, showHidden = TRUE, ...){
  data.table(labkey.selectRows(baseUrl = con$config$labkey.url.base,
                               folderPath = con$config$labkey.url.path,
                               schemaName = schema,
                               queryName = query,
                               showHidden = showHidden,
                               ...),
             stringsAsFactors = FALSE)
}

.ISCon$methods(
  listParticipantGroups = function(){
    "returns a dataframe with all saved Participant Groups on ImmuneSpace.\n"
    
    if(config$labkey.url.path != "/Studies/"){
      stop("labkey.url.path must be /Studies/. Use CreateConnection with all studies.")
    }
    
    user <- .validateUser(con = .self)
    
    # get userId from email, need separate call b/c of different schema
    userSql <- sprintf("SELECT UserId FROM Users WHERE Email='%s'", user)
    userId <- labkey.executeSql(config$labkey.url.base,
                                config$labkey.url.path,
                                schemaName = "core",
                                sql = userSql,
                                colNameOpt = "fieldname")
    
    # get rowId (for ParticipantCategory) with userId
    pcatSql <- sprintf("SELECT rowId FROM ParticipantCategory WHERE OwnerId='%s'", userId)
    
    # get GroupId (rowId in ParticipantGroup), label from CategoryId (rowId in PartCat)
    pgrpSql <- paste0("SELECT RowId, Label FROM ParticipantGroup WHERE CategoryId IN (", pcatSql, ")")
    pgrpIds <- labkey.executeSql(config$labkey.url.base,
                                config$labkey.url.path,
                                schemaName = "study",
                                sql = pgrpSql,
                                colNameOpt = "fieldname",
                                showHidden = TRUE)
    
    # get number of participants in group
    grpStr <- paste( pgrpIds$RowId, collapse = ",")
    sumSql <- paste0("SELECT GroupId FROM ParticipantGroupMap WHERE GroupId IN (", grpStr, ")")
    sumGrps <- data.table(labkey.executeSql(config$labkey.url.base,
                                 config$labkey.url.path,
                                 schemaName = "study",
                                 sql = sumSql,
                                 colNameOpt = "fieldname"))

    # Create results table with label, groupId, and number of participants
    sumGrps <- sumGrps[, Subjects := .N, by = "GroupId" ]
    sumGrps <- data.frame(sumGrps[ !duplicated(sumGrps), ])
    result <- merge(sumGrps, pgrpIds, by.x = "GroupId", by.y = "RowId")
    colnames(result) <- c("groupId", "subjects", "groupName")
    
    if( nrow(result) == 0 ){
      warning(paste0("No participant groups found for user email: ", user))
    }

    return(result)
  }
)

.ISCon$methods(
  getParticipantData = function(group, dataType){
    "returns a dataframe with ImmuneSpace data subset by groupId.\n
    group: Use listParticipantGroups() to find Participant groupId or groupName.\n
    dataType: Use con$available_datasets to see possible dataType inputs.\n"
    
    if(config$labkey.url.path != "/Studies/"){
      stop("labkey.url.path must be /Studies/. Use CreateConnection with all studies.")
    }
    
    if(!(dataType %in% .self$available_datasets$Name)){
      stop("DataType must be in ", paste(.self$available_datasets$Name, collapse = ", "))
    }
    
    # Checking to ensure user is owner of group on correct machine
    # Note that groupId is used for actually sql statement later regardless of user input
    user <- .validateUser(con = .self)
    validGrps <- .self$listParticipantGroups()
    if(typeof(group) == "double"){
      col <- "groupId"
      groupId <- group
    }else if( typeof(group) == "character"){
      col <- "groupName"
      groupId <- validGrps$groupId[ validGrps$groupName == group ]
    }else{
      stop("Group ID or Name not interpretable as double or character. Please reset")
    }
    
    if( !( group %in% validGrps[[col]] ) ){
      stop(paste0("group ", group,
                  " is not in the set of ", col,
                  " created by ", user,
                  " on ", .self$config$labkey.url.base))
    }

    # Get data
    dt <- tolower(dataType)

    # assay data + participantGroup + demographic data
    sqlAssay <- paste0("SELECT ",
                  dt, ".*,",
                  "pmap.GroupId AS groupId, ",
                  "demo.age_reported AS age_reported, ",
                  "demo.race AS race, ",
                  "demo.gender AS gender, ",
                  "FROM ",
                  dt,
                  " JOIN ParticipantGroupMap AS pmap ",
                  "ON ", dt, ".ParticipantId = pmap.ParticipantId",
                  " JOIN Demographics AS demo ",
                  "ON ", dt, ".ParticipantId = demo.ParticipantId",
                  " WHERE groupId = ", as.character(groupId))
    
    assayData <- labkey.executeSql(baseUrl = .self$config$labkey.url.base,
                                 folderPath = .self$config$labkey.url.path,
                                 schemaName = "study",
                                 sql = sqlAssay,
                                 colNameOpt = "fieldname")
    
    # Arm Data needs to come from different schema so need second call
    sqlArm <- "SELECT name, arm_accession FROM arm_or_cohort"
    armData <- labkey.executeSql(baseUrl = .self$config$labkey.url.base,
                                 folderPath = .self$config$labkey.url.path,
                                 schemaName = "immport",
                                 sql = sqlArm,
                                 colNameOpt = "fieldname")

    # Add Arm data to assayData
    assayData$arm_name <- armData$name[ match(assayData$arm_accession, armData$arm_accession) ]

    # SelectRows = easiest way to grab original columns since DataType could be one of
    # a number of possibilities. Leave one row otherwise have to handle empty df warning msg
    defaultCols <- labkey.selectRows(baseUrl = .self$config$labkey.url.base,
                                     folderPath = .self$config$labkey.url.path,
                                     schemaName = "study",
                                     queryName = dataType,
                                     colNameOpt = "fieldname",
                                     maxRows = 1)

    # b/c of lookups defaultCols has versions of demo / arm fieldnames with path
    possibleNames <- c( colnames(defaultCols), "race", "gender", "age_reported", "arm_name" )

    filtData <- assayData[ , (colnames(assayData) %in% possibleNames) ]
    
    return(filtData)
  }
)
