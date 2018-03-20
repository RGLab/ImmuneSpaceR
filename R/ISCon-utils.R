#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------

# Print the ISCon object
ISCon$set(
  which = "public",
  name = "print",
  value = function() {
    cat(sprintf("Immunespace Connection to study %s\n", self$study))

    cat(sprintf("URL: %s\n",
                file.path(gsub("/$","", self$config$labkey.url.base),
                          gsub("^/","", self$config$labkey.url.path)))
    )

    cat(sprintf("User: %s\n", self$config$labkey.user.email))

    cat("Available datasets\n")

    for (i in 1:nrow(self$availableDatasets)) {
      cat(sprintf("\t%s\n", self$availableDatasets[i, Name]))
    }

    if (!is.null(self$cache[[private$.constants$matrices]])) {
      cat("Expression Matrices\n")
      for (i in 1:nrow(self$cache[[private$.constants$matrices]])) {
        cat(sprintf("\t%s\n", self$cache[[private$.constants$matrices]][i, name]))
      }
    }
  }
)


# Clear the cache field
ISCon$set(
  which = "public",
  name = "clearCache",
  value = function() {
    self$cache[grep("^GE", names(self$cache), invert = TRUE)] <- NULL
  }
)



# PRIVATE ----------------------------------------------------------------------

# Check if study is valid
#' @importFrom gtools mixedsort
ISCon$set(
  which = "private",
  name = ".checkStudy",
  value = function(verbose = FALSE) {
    sdyNm <- basename(self$config$labkey.url.path)
    dirNm <- dirname(self$config$labkey.url.path)
    gTerm <- ifelse(dirNm == "/HIPC", "^IS\\d{1,3}$", "^SDY\\d{2,4}$")

    # adjust for "" connection
    if (sdyNm == "Studies") {
      sdyNm <- ""
      dirNm <- "/Studies"
    }

    folders <- labkey.getFolders(self$config$labkey.url.base, dirNm)
    subdirs <- gsub(paste0(dirNm, "/"), "", folders$folderPath)
    validSdys <- mixedsort(subdirs[grep(gTerm, subdirs)])

    if (!(sdyNm %in% c("", validSdys))) {
      if (verbose == FALSE) {
        stop(paste0(sdyNm, " is not a valid study. \n Use `verbose = TRUE` to see list of valid studies."))
      } else {
        stop(paste0(sdyNm, " is not a valid study\nValid studies: ",
                    paste(validSdys, collapse=", ")))
      }
    }
  }
)


# Clean up the column names
ISCon$set(
  which = "private",
  name = ".munge",
  value = function(x) {
    new <- tolower(gsub(" ", "_", basename(x)))
    idx <- which(duplicated(new) | duplicated(new, fromLast = TRUE))

    if (length(idx) > 0) {
      new[idx] <- private$.munge(gsub("(.*)/.*$", "\\1", x[idx]))
    }

    new
  }
)


# Check if the connection is at project level ("/Studies")
ISCon$set(
  which = "private",
  name = ".isProject",
  value = function() {
    self$config$labkey.url.path == "/Studies/"
  }
)


# Check if the coonection is running locally (prod/test)
ISCon$set(
  which = "private",
  name = ".isRunningLocally",
  value = function(path) {
    file.exists(path)
  }
)


# Replace url path to a local path
ISCon$set(
  which = "private",
  name = ".localStudyPath",
  value = function(urlpath) {
    LOCALPATH <- "/share/files/"
    PRODUCTION_HOST <- "www.immunespace.org"
    TEST_HOST <- "test.immunespace.org"

    gsub(
      file.path(gsub("/$", "", self$config$labkey.url.base), "_webdav"),
      file.path(LOCALPATH),
      urlpath
    )
  }
)



# HELPER -----------------------------------------------------------------------

# less verbose wrapper for labkey.selectRows function
.getLKtbl <- function(con, schema, query, showHidden = TRUE, ...) {
  data.table(
    labkey.selectRows(
      baseUrl = con$config$labkey.url.base,
      folderPath = con$config$labkey.url.path,
      schemaName = schema,
      queryName = query,
      showHidden = showHidden,
      ...
    ),
    stringsAsFactors = FALSE
  )
}


# Ensure only one valid user email is returned or possible errors handled
.validateUser <- function(con) {
  # First Check for apiKey and netrc file
  api <- Rlabkey:::ifApiKey()
  sink(tempfile())
  validNetrc <- tryCatch({
    check_netrc()
  }, error = function(e) {
    return(NULL)
  })
  sink()

  # Case 1: if apiKey, but no Netrc (possible in Integrated RStudio session UI)
  # get user email from global environment, which is set by labkey.init and
  # handle unexpected case of labkey.user.email not being found.
  if (!is.null(api) & is.null(validNetrc)) {
    user <- try(get("labkey.user.email"), silent = TRUE)
    if (inherits(user, "try-error")) {
      stop("labkey.user.email not found, please set")
    }

    # Case 2: valid netrc file (with or without ApiKey)
    # To mimic LK.get() method - use first login for correct machine
  } else if (!is.null(validNetrc)) {
    machine <- gsub("https://", "", con$config$labkey.url.base)
    netrc <- unlist(strsplit(readLines(validNetrc), split = " "))
    user <- netrc[grep(machine, netrc) + 2][[1]]
  }

  siteUsers <- labkey.selectRows(
    baseUrl = con$config$labkey.url.base,
    folderPath = con$config$labkey.url.path,
    schemaName = "core",
    queryName = "siteUsers"
  )

  # Admin user pulls whole table
  if (dim(siteUsers)[[1]] > 1) {
    userId <- siteUsers$`User Id`[ siteUsers$Email == user ]

    # Non-Admin sees only self and no email
  } else {
    userId <- siteUsers$`User Id`[[1]]
  }

  return(userId)
}
