#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------

# Print the ISCon object
ISCon$set(
  which = "public",
  name = "print",
  value = function() {
    cat("<ImmuneSpaceConnection>\n")

    cat(sprintf("  Study: %s\n", self$study))

    cat(
      sprintf(
        "  URL: %s\n",
        file.path(
          gsub("/$", "", self$config$labkey.url.base),
          gsub("^/", "", self$config$labkey.url.path)
        )
      )
    )

    cat(sprintf("  User: %s\n", self$config$labkey.user.email))

    cat(" ", nrow(self$availableDatasets), "Available Datasets\n")

    for (i in seq_len(nrow(self$availableDatasets))) {
      cat(sprintf("    - %s\n", self$availableDatasets[i, Name]))
    }

    GEMs <- self$cache[[private$.constants$matrices]]
    if (!is.null(GEMs)) {
      cat(" ", nrow(GEMs), "Available Expression Matrices\n")
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
        stop(paste0(
          sdyNm, " is not a valid study\nValid studies: ",
          paste(validSdys, collapse = ", ")
        ))
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


# Get names of files in a single folder from webdav link
ISCon$set(
  which = "private",
  name = ".listISFiles",
  value = function(link) {
    response <- NULL
    res <- tryCatch(
      Rlabkey:::labkey.get(link),
      warning = function(w) return(w),
      error = function(e) return(NULL)
    )
    if (!is.null(res)) {
      tmp <- rjson::fromJSON(res)
      response <- sapply(tmp$files, function(x) {
        return(x$text)
      }) # basename only
    }
    response
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
