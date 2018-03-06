# Decode filters in case the user used column labels instead of column names
#' @importFrom utils URLdecode URLencode
.check_filter <- function(con, schema, query, view = "", colFilter) {
  # helper functions
  extractNames <- function(colFilter) {
    tolower(unlist(lapply(gsub("~.*$", "", colFilter), URLdecode)))
  }

  getColnames <- function(colNameOpt) {
    colnames(
      .getLKtbl(
        con = con,
        schema = schema,
        query = query,
        viewName = view,
        maxRows = 0,
        colNameOpt = colNameOpt,
        showHidden = FALSE
      )
    )
  }

  fixColFilter <- function(colFilter, labels, names) {
    if (any(old %in% labels)) {
      idx <- which(old %in% labels)
      new <- unlist(lapply(names[match(old[idx], labels)], URLencode))
      colFilter[idx] <- paste0(paste0(new, "~"), gsub("^.*~", "", colFilter[idx]))
    }

    colFilter
  }

  # step 1:
  # Extract the colnames used in the filter
  old <- extractNames(colFilter)
  old[old == "participant_id"] <- "participant id"

  # get colnames by fieldname and caption
  suppressWarnings({
    fieldname <- getColnames(colNameOpt = "fieldname")
    caption <- tolower(getColnames(colNameOpt = "caption"))
  })

  # Fix the colFilter with caption
  colFilter <- fixColFilter(colFilter, caption, fieldname)

  # step 2:
  # Extract the colnames in the modified filter
  old <- extractNames(colFilter)

  # get lookups columns and labels to fix them by
  lookups <- fieldname[grep("/", fieldname)]
  labels <- gsub("\\w+/", "", lookups)

  # Fix the colFilter with lookups columns
  # (e.g., "DataSets/demographics/age_reported")
  colFilter <- fixColFilter(colFilter, labels, lookups)

  colFilter
}

## TODO: Need a way to translate Curl operators into math
## (so that we can filter the saved tables using colFilter).
filterDT <- function(table, col, value) {
  table[eval(as.name(col)) == value]
}

filter_cached_copy <- function(filters, data) {
  decoded <- unlist(lapply(filters, URLdecode))

  cols <- gsub(".*/", "", gsub("~.*", "", decoded))
  values <- gsub(".*=", "", decoded)

  for (i in 1:length(filters)) {
    data <- filterDT(table, cols[i], values[i])
  }

  data
}

# Downloads a dataset and cache the result in the connection object.
ISCon$set(
  which = "public",
  name = "getDataset",
  value = function(x, original_view = FALSE, reload = FALSE, colFilter = NULL, ...) {
    if (nrow(self$available_datasets[Name %in% x]) == 0) {
      wstring <- paste0(study, " has invalid data set: ", x)
      if (config$verbose) {
        wstring <- paste0(wstring, "\n",
                          "Valid datasets for ", study, ": ",
                          paste(self$available_datasets$Name, collapse = ", "), ".")
      }
      stop(wstring)
    } else {
      cache_name <- paste0(x, ifelse(original_view, "_full", ""))
      nOpts <- length(list(...))

      if (!is.null(self$cache[[cache_name]]) &&
          !reload &&
          is.null(colFilter) &&
          nOpts == 0) { # Serve cache
        data <- self$cache[[cache_name]]
        #if(!is.null(colFilter)){
        #  data <- filter_cached_copy(colFilter, data)
        #  return(data)
        #} else{
        #}

      } else { # Download the data
        viewName <- NULL
        if (original_view) {
          viewName <- "full"
        }

        if (!is.null(colFilter)) {
          colFilter <- .check_filter(
            con = self,
            schema = "study",
            query = x,
            view = viewName,
            colFilter = colFilter
          )
          cache <- FALSE
        } else if (length(nOpts) > 0) {
          cache <- FALSE
        } else{
          cache <- TRUE
        }

        data <- .getLKtbl(
          con = self,
          schema = "study",
          query = x,
          viewName = viewName,
          colNameOpt = "caption",
          colFilter = colFilter,
          showHidden = FALSE,
          ...
        )
        setnames(data, private$.munge(colnames(data)))

        if (cache) {
          self$cache[[cache_name]] <- data
        }
      }

      if (!is.null(self$config$use.data.frame) && self$config$use.data.frame) {
        data <- data.frame(data)
      }

      data
    }
  }
)
