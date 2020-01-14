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
ISCon$set(
  which = "private",
  name = ".checkStudy",
  value = function(verbose = FALSE) {
    folders <- labkey.getFolders(
      baseUrl = self$config$labkey.url.base,
      folderPath = "",
      includeSubfolders = TRUE,
      includeEffectivePermissions = TRUE
    )
    folders <- folders[ grepl("IS\\d{1}|SDY\\d{2,4}", folders$name), ]
    study <- basename(self$config$labkey.url.path)

    if (!(study %in% c("Studies", folders$name))) {
      msg <- ifelse(verbose == FALSE,
        " is not a valid study. \n Use `verbose = TRUE` to see list of valid studies.",
        paste0(" is not a valid study\nValid studies: ", paste(folders$name, collapse = ", "))
      )
      stop(paste0(study, msg))
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


# Replace url with webdav to a local path on server or dev machine
ISCon$set(
  which = "private",
  name = ".localStudyPath",
  value = function(link) {
    # If running on individual dev machine, must have symlink set for '/share' to
    # ~/release<ver>/build/deploy
    gsub(
      file.path(gsub("/$", "", self$config$labkey.url.base), "(|/)_webdav"),
      file.path("/share/files"),
      link
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
      warning = function(w) {
        return(w)
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (!is.null(res)) {
      tmp <- fromJSON(res, simplifyDataFrame = FALSE)
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


.mixedsort <- function(x) {
  x[order(as.integer(gsub("[A-z]+", "", x)))]
}
