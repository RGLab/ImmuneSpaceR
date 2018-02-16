########################################################################
###                        IS CON HELPER FN                          ###
########################################################################

.ISCon$methods(
  .munge = function(x) {
    new <- tolower(gsub(" ", "_", basename(x)))
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
      ge <- tryCatch(.getLKtbl(con = .self, 
                      schema = "assay.Expressionmatrix.matrix",
                      query = "InputSamples",
                      viewName = "gene_expression_matrices",
                      colNameOpt = "fieldname"),
                     error = function(e) return(e))
      
      if (length(ge$message) > 0) {
        stop("Gene Expression Inputs not found for study.")
      }
      
      setnames(ge,.self$.munge(colnames(ge)))
      data_cache[[constants$matrix_inputs]] <<- ge
    }
  }
)

.ISCon$methods(
  .isRunningLocally = function(path) {
    file.exists(path)
  }
)

.ISCon$methods(
  .localStudyPath = function(urlpath) {
    LOCALPATH <- "/share/files/"
    PRODUCTION_HOST <- "www.immunespace.org"
    TEST_HOST <- "test.immunespace.org"

    gsub(file.path(gsub("/$", "", config$labkey.url.base), "_webdav"),
         file.path(LOCALPATH),
         urlpath)
  }
)

.ISCon$methods(
  listDatasets = function(output = c("datasets", "expression")) {
    "List the datasets available in the study or studies of the connection."

    if (!all(output %in% c("datasets", "expression"))) {
      stop("output other than datasets and expressions not allowed")
    }

    if ("datasets" %in% output) {
      cat("datasets\n")
      for (i in 1:nrow(available_datasets)) {
        cat(sprintf("\t%s\n", available_datasets[i, Name]))
      }
    }

    if ("expression" %in% output) {
      if (!is.null(data_cache[[constants$matrices]])) {
        cat("Expression Matrices\n")
        for (i in 1:nrow(data_cache[[constants$matrices]])) {
          cat(sprintf("\t%s\n", data_cache[[constants$matrices]][i, name]))
        }
      } else {
        cat("No Expression Matrices Available")
      }
    }
  }
)

.ISCon$methods(
    listGEAnalysis = function(){
      "List available gene expression analysis for the connection."
      GEA <- tryCatch(.getLKtbl(con = .self, 
                       schema = "gene_expression",
                       query = "gene_expression_analysis",
                       showHidden = FALSE,
                       colNameOpt = "rname"),
                      error = function(e) return(e) )
      
      if( length(GEA$message) > 0 ){
        stop("Study does not have Gene Expression Analyses.")
      }
      
      return(GEA)
    }
)

.ISCon$methods(
  getGEAnalysis = function(...) {
    "Downloads data from the gene expression analysis results table.\n
    '...': A list of arguments to be passed to labkey.selectRows."

    GEAR <- tryCatch(.getLKtbl(con = .self, 
                      schema = "gene_expression",
                      query = "DGEA_filteredGEAR",
                      viewName = "DGEAR",
                      colNameOpt = "caption",
                      ...),
                     error = function(e) return(e))
    
    if(length(GEAR$message) > 0) {
      stop("Gene Expression Analysis not found for study.")
    }

    setnames(GEAR, .self$.munge(colnames(GEAR)))

    return(GEAR)
  }
)

.ISCon$methods(
  clear_cache = function() {
    "Clear the data_cache. Remove downloaded datasets and expression matrices."
    data_cache[grep("^GE", names(data_cache), invert = TRUE)] <<- NULL
  }
)

.ISCon$methods(
  show = function() {
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
  getGEFiles = function(files, destdir = ".", quiet = FALSE) {
    "Download gene expression raw data files.\n
    files: A character. Filenames as shown on the gene_expression_files dataset.\n
    destdir: A character. The local path to store the downloaded files."
    links <- paste0(config$labkey.url.base, "/_webdav/",
                    config$labkey.url.path,
                    "/%40files/rawdata/gene_expression/", files)

    sapply(links, function(x) {
      download.file(url = links[1],
                    destfile = file.path(destdir, basename(x)),
                    method = "curl",
                    extra = "-n",
                    quiet = quiet)
    })
  }
)

#########################################################################
###                 NON-CON HELPER FN                                 ###
#########################################################################

# sort studies by number
.spSort <- function(vec){
  if( length(vec) > 0 ){
    vec <- sort(as.numeric(gsub("SDY","",vec)))
    vec <- paste0("SDY", as.character(vec))
  }else{
    vec <- NULL
  }
}

# get vector of study folders
.getSdyVec <- function(.self){
  studies <- labkey.getFolders(baseUrl = .self$config$labkey.url.base, 
                               folderPath = "/Studies/")[, 1]
  studies <- studies[ !studies %in% c("SDY_template","Studies") ]
}

# helper to get sdy ids from subids
.subidsToSdy <- function(subids){
  sdys <- unique(gsub("^SUB.+", NA, unlist(strsplit(subids, split = "\\."))))
  sdys <- sdys[ !is.na(sdys) ]
  sdys <- paste0("SDY", sdys)
}

# get names of files in a single folder from webdav link
.listISFiles <- function(link, .self) {
  opts <- .self$config$curlOptions
  opts$options$netrc <- 1L
  
  response <- NULL
  
  res <- GET(url = link, config = opts)
  if (!http_error(res)) {
    response_json <- httr::content(res)
    response <- unlist(lapply(response_json$files, function(x) x$text))
  }
  
  response
}

.urlExists <- function(url, .self) {
  opts <- .self$config$curlOptions
  opts$options$netrc <- 1L
  
  res <- HEAD(url, config = opts)
  
  if (http_error(res)) {
    ret <- FALSE
  } else {
    if (http_type(res) == "application/json") {
      res <- GET(url, config = opts)
      cont <- httr::content(res)
      ret <- is.null(cont$exception)
    } else {
      ret <- TRUE
    }
  }
  
  ret
}

# Generate named list of files in either rawdata or analysis/exprs_matrices folders
.getGEFileNms <- function(.self, rawdata) {
  studies <- .getSdyVec(.self)
  
  # check webdav folder for presence of rawdata
  file_list <- lapply(studies, FUN = function(sdy) {
    suffix <- ifelse(rawdata == TRUE,
                     "/%40files/rawdata/gene_expression?method=JSON",
                     "/%40files/analysis/exprs_matrices?method=JSON")
    dirLink <-  paste0(.self$config$labkey.url.base,
                       "/_webdav/Studies/",
                       sdy,
                       suffix)
    files <- .listISFiles(dirLink, .self)
    if( rawdata == TRUE ){
      if( !is.null(files) ){
        files <- length(files) > 0
      } 
    }
    return(files)
  })
  
  names(file_list) <- studies
  
  return(file_list)
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
    if( inherits(user, "try-error") ){
      stop("labkey.user.email not found, please set")
    }

    # Case 2: valid netrc file (with or without ApiKey)
    # To mimic LK.get() method - use first login for correct machine
  }else if( !is.null(validNetrc) ){
    machine <- gsub("https://", "", con$config$labkey.url.base)
    netrc <- unlist(strsplit(readLines(validNetrc), split = " "))
    user <- netrc[ grep(machine, netrc) + 2 ][[1]]
  }

  siteUsers <- labkey.selectRows(baseUrl = con$config$labkey.url.base,
                                 folderPath = con$config$labkey.url.path,
                                 schemaName = "core",
                                 queryName = "siteUsers")

  # Admin user pulls whole table
  if( dim(siteUsers)[[1]] > 1 ){
    userId <- siteUsers$`User Id`[ siteUsers$Email == user ]

  # Non-Admin sees only self and no email
  }else{
    userId <- siteUsers$`User Id`[[1]]
  }

  return(userId)
}


#################################################################################
###                     MAINTAINENANCE METHODS                                ###
#################################################################################

# Returns a list of data frames where TRUE in file_exists column marks files that are accessible.
# This function is used for administrative purposes to check that the flat files
# are properly loaded and accessible to the users.
#' @importFrom rjson fromJSON
#' @importFrom parallel mclapply detectCores
.ISCon$methods(
  .test_files = function(what = c("gene_expression_files",
                                  "fcs_sample_files",
                                  "fcs_control_files",
                                  "protocols",
                                  "ge_matrices")) {
    
    # HELPER fn ----------------------------------
    check_links <- function (dataset, folder) {
      res <- data.frame(file_info_name = NULL,
                        study_accession = NULL,
                        file_link = NULL,
                        file_exists = NULL,
                        stringsAsFactors = FALSE)
      
      if( dataset %in% .self$available_datasets$Name ){
        temp <- .self$getDataset(dataset, original_view = TRUE)

        if( dataset == "fcs_control_files" ){
          temp <- temp[, file_info_name := control_file]
          temp <- temp[, c("pid", "sid") := data.table::tstrsplit(participant_id, "\\.")]
          temp <- temp[, study_accession := paste0("SDY", sid)]
        }

        temp <- temp[ !is.na(file_info_name) ]
        temp <- unique(temp[, list(study_accession, file_info_name)])
        
        file_link <- paste0(.self$config$labkey.url.base,
                            "/_webdav/Studies/",
                            temp$study_accession,
                            "/%40files/rawdata/",
                            folder,
                            "/",
                            sapply(temp$file_info_name, URLencode))
        
        studies <- unique(temp$study_accession)
        folder_link <- paste0(.self$config$labkey.url.base,
                              "/_webdav/Studies/",
                              studies,
                              "/%40files/rawdata/",
                              folder,
                              "?method=JSON")

        file_list <- unlist(
          mclapply(
            folder_link,
            .listISFiles,
            .self,
            mc.cores = detectCores()
          )
        )

        file_exists <- temp$file_info_name %in% file_list

        res <- data.frame(study = temp$study_accession,
                          file_link = file_link,
                          file_exists = file_exists,
                          stringsAsFactors = FALSE)

        print(paste0(sum(res$file_exists),
                     "/",
                     nrow(res),
                     " ",
                     dataset,
                     " with valid links."))
      }
      res
    }
    
    # MAIN fn ------------------------------
    ret <- list()
    what <- tolower(what)
    
    if ("gene_expression_files" %in% what) {
      ret$gene_expression_files <- check_links("gene_expression_files",
                                               "gene_expression")
    }

    if ("fcs_sample_files" %in% what) {
      ret$fcs_sample_files <- check_links("fcs_sample_files", "flow_cytometry")
    }

    if ("fcs_control_files" %in% what) {
      ret$fcs_control_files <- check_links("fcs_control_files", "flow_cytometry")
    }

    if ("protocols" %in% what) {
      if (.self$.isProject()) {
        folders_list <- labkey.getFolders(baseUrl = config$labkey.url.base,
                                          folderPath = "/Studies/")
        folders <- folders_list[, 1]
        folders <- folders[!folders %in% c("SDY_template","Studies")]
      } else {
        folders <- basename(config$labkey.url.path)
      }

      file_link <- paste0(config$labkey.url.base,
                          "/_webdav/Studies/",
                          folders,
                          "/%40files/protocols/",
                          folders,
                          "_protocol.zip")

      file_exists <- unlist(
        mclapply(
          file_link,
          .urlExists,
          .self,
          mc.cores = detectCores()
        )
      )

      print(paste0(sum(file_exists),
                   "/",
                   length(file_exists),
                   " protocols with valid links."))

      ret$protocols <- data.frame(study = folders,
                                  file_link = file_link,
                                  file_exists = file_exists,
                                  stringsAsFactors = FALSE)
    }

    if ("ge_matrices" %in% what) {
        mx <- .getLKtbl(con = .self,
                        schema = "assay.ExpressionMatrix.matrix",
                        query = "Runs",
                        colNameOpt = "rname")

        mxLinks <- paste0(.self$config$labkey.url.base,
                          "/_webdav/Studies/",
                          mx$folder_name,
                          "/@files/analysis/exprs_matrices/",
                          mx$name,
                          ".tsv")

        file_exists <- unlist(
          mclapply(
            mxLinks,
            .urlExists,
            .self,
            mc.cores = detectCores()
          )
        )

        print(paste0(sum(file_exists),
                     "/",
                     length(file_exists),
                     " ge_matrices with valid links."))

        ret$ge_matrices <- data.frame(file_link = mxLinks,
                                      file_exists = file_exists,
                                      stringsAsFactors = FALSE)

    } else {
      ret$ge_matrices <- data.frame(file_link = NULL,
                          file_exists = NULL, 
                          stringsAsFactors = FALSE)
    }
    
    return(ret)
  }
)

#-------------------GEM CLEANUP-----------------------------------
# This function outputs a list of lists.  There are six possible sub-lists >
# 1. $gemAndRaw = study has gene expression matrix flat file and raw data
# 2. $gemNoRaw = study has gene expression matrix flat file but no raw data (unlikely)
# 3. $rawNoGem = study has raw data but no gene expression flat file, which is likely
#                in the case that no annotation pkg is availble in bioconductor or
#                it is RNAseq and had trouble being processed.
# 4. $gefNoGem = con$getDataset("gene_expression_files") reports raw data available, but
#                no gene expression matrix flat file has been generated. Similar to rawNoGem.
# 5. $gefNoRaw = con$getDataset() reports files being available, but no rawdata found.  This
#                may be due to files being in GEO, but not having been downloaded to ImmPort
#                and ImmuneSpace.
# 6. $rawNoGef = This would be unexpected. Rawdata present on server, but no gene expression
#                files found by con$getDataset().
.ISCon$methods(
  .sdysWithoutGems = function() {

    res <- list()

    # get list of matrices and determine which sdys they represent
    gems <- .self$data_cache$GE_matrices
    withGems <- unique(gems$folder)

    file_list <- .getGEFileNms(.self = .self, rawdata = TRUE)
    file_list <- file_list[file_list != "NULL"]
    emptyFolders <- names(file_list)[file_list == FALSE]
    withRawData <- names(file_list)[file_list == TRUE]

    # Compare lists
    res$gemAndRaw <- .spSort(intersect(withRawData, withGems))
    res$gemNoRaw <- .spSort(setdiff(withGems, withRawData))
    res$rawNoGem <- ..spSort(setdiff(withRawData, withGems))

    # Check which studies without gems have gef in IS
    ge <- con$getDataset("gene_expression_files")
    geNms <- unique(ge$participant_id)
    gefSdys <- unique(sapply(geNms, FUN = function(x) {
      res <- strsplit(x, ".", fixed = TRUE)
      return(res[[1]][2])
    }))
    gefSdys <- paste0("SDY", gefSdys)

    res$gefNoGem <- .spSort(gefSdys[!(gefSdys %in% withGems)])
    res$gefNoRaw <- .spSort(setdiff(res$gefNoGem, res$rawNoGem))
    res$rawNoGef <- .spSort(setdiff(res$rawNoGem, res$gefNoGem))

    return(res)
  }
)

# Remove gene expression matrices that do not correspond to a run currently on prod or test
# in the query assay.ExpressionMatrix.matrix.Runs. NOTE: Important to change the labkey.url.base
# variable depending on prod / test to ensure you are not deleting any incorrectly.
#' @importFrom httr DELETE HEAD http_type http_error
.ISCon$methods(
  .rmOrphanGems = function() {

    .getNoRunPres <- function(.self, runs) {
      # get flat files list with appropriate names (sdy)
      emFls <- .getGEFileNms(.self = .self, rawdata = FALSE)
      emFls <- emFls[emFls != "NULL"]
      tmpNms <- rep(x = names(emFls), times = lengths(emFls))
      emFls <- unlist(emFls)
      names(emFls) <- tmpNms

      # get names from emFls for comparison
      emNms <- sapply(emFls, FUN = function(x) {
        return(strsplit(x, "\\.tsv")[[1]][1])
      })

      # check runs for file names - assume if not in runs$Name then ok to delete!
      noRunPres <- emNms[!(emNms %in% runs$Name)]
      noRunPres <- noRunPres[!duplicated(noRunPres)]
    }

    # helper curl FN
    curlDelete <- function(baseNm, sdy, .self) {
      opts <- .self$config$curlOptions
      opts$options$netrc <- 1L
      
      tsv <-  paste0(.self$config$labkey.url.base,
                     "/_webdav/Studies/",
                     sdy,
                     "/%40files/analysis/exprs_matrices/",
                     baseNm,
                     ".tsv")
      smry <- paste0(tsv, ".summary")
      
      tsvRes <- DELETE(url = tsv, config = opts)
      smryRes <- DELETE(url = smry, config = opts)
      
      list(tsv = tsvRes, summary = smryRes)
    }

    # Double check we are working at project level and on correct server!
    if (!.self$.isProject()) stop("Can only be run at project level")
    chkBase <- readline(prompt = paste0("You are working on ",
                                        .self$config$labkey.url.base,
                                        ". Continue? [T / f] "))
    if (!(chkBase %in% c("T", "t", ""))) return("Operation Aborted.")

    # get runs listed in the proper table
    runs <- data.table(labkey.selectRows(baseUrl = .self$config$labkey.url.base,
                                         folderPath = .self$config$labkey.url.path,
                                         schemaName = "assay.ExpressionMatrix.matrix",
                                         queryName = "Runs",
                                         showHidden = TRUE))

    noRunPres <- .getNoRunPres(.self = .self, runs = runs)

    # if files to-be-rm, confirm rm ok, attempt delete, check results and report
    if (length(noRunPres) == 0) {
      return("No orphans found. Hurray!")
    } else {
      print(noRunPres)
    }
    ok2rm <- readline(prompt = "Ok to remove all files listed above? [Y / n] ")
    if (toupper(ok2rm) == "Y" | ok2rm == "") {
      for (i in 1:length(noRunPres)) {
        curlDelete(baseNm = noRunPres[i],
                   sdy = names(noRunPres)[i],
                   .self = .self)
      }
      noRunPresPost <- .getNoRunPres(.self = .self, runs = runs)
      if (length(noRunpresPost) == 0) {
        return("No orphans found after removal. Success!")
      } else {
        print("Problems Occurred. Remaining Files with No Runs Present")
        print(noRunPresPost)
        print("************")
      }
    } else {
      print("Operation aborted.")
    }
  }
)

#----------------STUDY-CHECK-------------------------------------------------
# This method allows admin to check which studies are compliant with the 
# following modules/data files (GEF, RAW, GEO, GEM, DE, GEE, IRP, GSEA)
.ISCon$methods(
  .studiesComplianceCheck = function(filterNonGE = TRUE, showAllCols = FALSE, onlyShowNonCompliant = TRUE){
    
    # For labkey.executeSql calls
    baseUrl <- .self$config$labkey.url.base
    
    # res table
    if( .self$study != "Studies" ){
      Sdys <- c(.self$study) # single study
    } else {
      Sdys <- .spSort(.getSdyVec(.self)) # all studies
    }
    
    mods <- c("GEF", "RAW", "GEO", "GEM", "DE_implied", "GEE_implied", "GSEA_implied", "IRP_implied")
    compDF <- data.frame(matrix(nrow = length(Sdys),
                                ncol = length(mods)),
                         row.names = Sdys)
    colnames(compDF) <- mods
    
    # GEM
    withGems <- unique(.self$data_cache$GE_matrices$folder)
    compDF$GEM <- rownames(compDF) %in% withGems
    
    # RAW
    file_list <- .getGEFileNms(.self = .self, rawdata = TRUE)
    file_list <- file_list[ file_list != "NULL" ]
    compDF$RAW <- rownames(compDF) %in% names(file_list)[ file_list == TRUE ]
    
    # GEF
    gef <- .self$getDataset("gene_expression_files")
    compDF$GEF <- rownames(compDF) %in% .subidsToSdy(gef$participant_id)
    
    # GEO
    geoGef <- gef[ !is.na(gef$geo_accession), ]
    compDF$GEO <- rownames(compDF) %in% .subidsToSdy(geoGef$participant_id)
    
    # GEE - Studies with subjects having both GEM and response data for any timepoints
    # NOTE: when GEE is changed to allow NAb, can uncomment nab / respSubs lines
    hai <- .self$getDataset("hai")
    # nab <- .self$getDataset("neut_ab_titer")
    # resp <- union(resp$participant_id, nab$participant_id)
    inputSmpls <- labkey.selectRows(baseUrl = baseUrl,
                                    folderPath = "/Studies",
                                    schemaName = "study",
                                    queryName = "HM_InputSamplesQuery",
                                    containerFilter = "CurrentAndSubfolders")
    exprResp <- merge(inputSmpls, hai,
                      by.x = c("Participant Id", "Study Time Collected"),
                      by.y = c("participant_id", "study_time_collected"))
    compDF$GEE_implied <- rownames(compDF) %in% .subidsToSdy(unique(exprResp$`Participant Id`))
    
    # GSEA - studies with subjects having GEAR results (i.e. compared multiple GEM timepoints)
    gear <- sapply(withGems, FUN = function(sdy){
      res <- tryCatch(labkey.executeSql(baseUrl = baseUrl,
                                        folderPath = paste0("/Studies/", sdy),
                                        schemaName = "gene_expression",
                                        sql = "SELECT COUNT (*) FROM gene_expression_analysis_results"),
                      error = function(e){ return( NA ) })
      output <- res[[1]] > 0 & !is.na(res)
    })
    compDF$GSEA_implied <- rownames(compDF) %in% names(gear)[gear == TRUE]

    # IRP - Studies with subjects from multiple cohorts with GEM data at both target
    # timepoint and baseline + response data
    nab <- .self$getDataset("neut_ab_titer")
    resp <- rbind(nab, hai, fill = TRUE) # nab has a col that hai does not
    inputSmpls$study <- gsub("^SUB\\d{6}\\.", "SDY", inputSmpls$`Participant Id`)
    library(dplyr)
    # NOTE: At least SDY180 has overlapping study_time_collected for both hours and days
    # so it is important to group by study_time_collected_unit as well. This is reflected
    # in IRP_timepoints_hai/nab.sql
    geCohortSubs <- inputSmpls %>%
                group_by(study, `Study Time Collected`, `Study Time Collected Unit`) %>%
                  # need multiple cohorts per timePoint (from IRP.js)
                filter( (length(unique(Cohort)) > 1 ) == TRUE ) %>% 
                ungroup() %>%
                group_by(study, Cohort, `Study Time Collected Unit`) %>%
                  # need baseline + later timePoint (from IRP.Rmd)
                filter( (length(unique(`Study Time Collected`)) > 1 ) == TRUE )
    geRespSubs <- geCohortSubs[ geCohortSubs$`Participant Id` %in% unique(resp$participant_id), ]
    compDF$IRP_implied <- rownames(compDF) %in% .subidsToSdy(geRespSubs$`Participant Id`)
    studyTimepoints <- geRespSubs %>%
                        group_by(study) %>%
                        summarize(timepoints = paste(unique(`Study Time Collected`), collapse = ","))
    compDF$IrpTimepoints <- studyTimepoints$timepoints[ match(rownames(compDF), studyTimepoints$study) ]

    # DE
    deSets <- c("Neutralizing antibody titer",
                "Enzyme-linked immunosorbent assay (ELISA)",
                "Enzyme-Linked ImmunoSpot (ELISPOT)",
                "Hemagglutination inhibition (HAI)",
                "Polymerisation chain reaction (PCR)",
                "Flow cytometry analyzed results",
                "Multiplex bead array asssay")
    compDF$DE_implied <- sapply(rownames(compDF), FUN = function(sdy){
      res <- suppressWarnings(tryCatch(labkey.executeSql(baseUrl = baseUrl,
                                                         folderPath = paste0("/Studies/", sdy),
                                                         schemaName = "study",
                                                         sql = "SELECT Label FROM ISC_datasets"),
                                       error = function(e){ return( NA ) })
      )
      ret <- any(res[[1]] %in% deSets)
    })

    # Validation based on modules being turned on
    getModSdys <- function(name){
      url <- paste0("https://www.immunespace.org/immport/studies/containersformodule.api?name=", name)
      res <- unlist( lapply(fromJSON(Rlabkey:::labkey.get(url))[[1]], function(x){ x[["name"]]} ) )
      res <- .spSort( res[ res != "Studies" & res != "SDY_template"] )
    }

    compDF$DE_actual <- rownames(compDF) %in% getModSdys("DataExplorer")
    compDF$GEE_actual <- rownames(compDF) %in% getModSdys("GeneExpressionExplorer")
    compDF$GSEA_actual <- rownames(compDF) %in% getModSdys("GeneSetEnrichmentAnalysis")
    compDF$IRP_actual <- rownames(compDF) %in% getModSdys("ImmuneResponsePredictor")

    # Differences between implied and actual
    compDF$DE_diff <- compDF$DE_implied == compDF$DE_actual
    compDF$GEE_diff <- compDF$GEE_implied == compDF$GEE_actual
    compDF$GSEA_diff <- compDF$GSEA_implied == compDF$GSEA_actual
    compDF$IRP_diff <- compDF$IRP_implied == compDF$IRP_actual

    colOrder <- c("RAW",
                  "GEF",
                  "GEO",
                  "GEM",
                  "DE_implied",
                  "DE_actual",
                  "DE_diff",
                  "GEE_implied",
                  "GEE_actual",
                  "GEE_diff",
                  "GSEA_implied",
                  "GSEA_actual",
                  "GSEA_diff",
                  "IRP_implied",
                  "IRP_actual",
                  "IRP_diff")

    compDF <- compDF[ , order(match(colnames(compDF), colOrder)) ]

    # Filter out studies that don't have GE since this is basis for everything
    if( filterNonGE == TRUE ){
      compDF <- compDF[ compDF$GEO == TRUE | compDF$GEF == TRUE, ]
    }

    # Subset to only show problematic studies
    if( onlyShowNonCompliant == TRUE ){
      redux <- compDF[ , grep("diff", colnames(compDF)) ]
      idx <- which(apply(redux, 1, all))
      compDF <- compDF[-(idx), ]
    }

    # Defaults to showing only the actual module status and the difference with the implied
    if( showAllCols == FALSE ){
      compDF <- compDF[ , grep("diff|act|time", colnames(compDF)) ]
      message("NOTE: \n Return object has columns marked with 'diff' that indicate a difference between the actual column that was generated by a call to the module url and the believed column that \n was created by looking at the filesystem.  The 'act' column is included as a reference.")
    }

    return(compDF)

  }
)
