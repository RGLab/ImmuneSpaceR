#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------

# Lists available datasets and expressions
ISCon$set(
  which = "public",
  name = "listDatasets",
  value = function(output = c("datasets", "expression")) {
    if (!all(output %in% c("datasets", "expression"))) {
      stop("output other than datasets and expressions not allowed")
    }

    if ("datasets" %in% output) {
      cat("datasets\n")
      for (i in seq_len(nrow(self$availableDatasets))) {
        cat(sprintf("\t%s\n", self$availableDatasets[i, Name]))
      }
    }

    if ("expression" %in% output) {
      if (!is.null(self$cache[[private$.constants$matrices]])) {
        cat("Expression Matrices\n")
        for (i in seq_len(nrow(self$cache[[private$.constants$matrices]]))) {
          cat(sprintf("\t%s\n", self$cache[[private$.constants$matrices]][i, name]))
        }
      } else {
        cat("No Expression Matrices Available")
      }
    }
  }
)


# Downloads a dataset and cache the result in the connection object.
#' @importFrom digest digest
ISCon$set(
  which = "public",
  name = "getDataset",
  value = function(x, original_view = FALSE, reload = FALSE, colFilter = NULL, transformMethod = "none", ...) {
    if (nrow(self$availableDatasets[Name %in% x]) == 0) {
      wstring <- paste0(
        "Empty data.table was returned.",
        " `", x, "` is not a valid dataset for ", self$study
      )
      if (self$config$verbose) {
        wstring <- paste0(
          wstring, "\n",
          "Valid datasets for ", self$study, ":\n  - ",
          paste(self$availableDatasets$Name, collapse = "\n  - ")
        )
      }
      warning(wstring, immediate. = TRUE)
      return(data.table())
    }

    # build a list of arguments to digest and compare
    args <- list(
      x = x,
      original_view = original_view,
      reload = reload,
      colFilter = colFilter,
      transformMethod = transformMethod,
      ...
    )

    # retrieve dataset from cache if arguments match
    digestedArgs <- digest(args)
    if (digestedArgs %in% names(self$cache) && !reload) {
      return(self$cache[[digestedArgs]]$data)
    }

    viewName <- NULL
    if (original_view) {
      viewName <- "full"
    }

    if (!is.null(colFilter)) {
      colFilter <- private$.checkFilter(
        colFilter = colFilter,
        schema = "study",
        query = x,
        view = viewName
      )
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

    # transform data if needed
    data <- .transformData(data, x, transformMethod)

    # caching
    self$cache[[digestedArgs]] <- list(args = args, data = data)

    data
  }
)



# PRIVATE ----------------------------------------------------------------------

# Set available datasets
ISCon$set(
  which = "private",
  name = ".setAvailableDatasets",
  value = function() {
    if (length(self$availableDatasets) == 0) {
      self$availableDatasets <- .getLKtbl(
        con = self,
        schema = "study",
        query = "ISC_study_datasets"
      )
    }
  }
)


# Retrieve column names
ISCon$set(
  which = "private",
  name = ".getColnames",
  value = function(colNameOpt, schema, query, view = "") {
    colnames(
      .getLKtbl(
        con = self,
        schema = schema,
        query = query,
        viewName = view,
        maxRows = 0,
        colNameOpt = colNameOpt,
        showHidden = FALSE
      )
    )
  }
)


# Decode filters in case the user used column labels instead of column names
ISCon$set(
  which = "private",
  name = ".checkFilter",
  value = function(colFilter, schema, query, view = "") {
    stopifnot(is.matrix(colFilter))

    holdsGroupFilter <- grepl("ParticipantId/(\\S+)~eq=\\1$", colFilter[1, 1])
    if (nrow(colFilter) == 1 && holdsGroupFilter) {
      return(colFilter)
    }

    # step 1:
    # Extract the column names used in the column filter
    old <- .extractNames(colFilter)
    old[old == "participant_id"] <- "participant id"

    # Retrieve the column names by `fieldname` and `caption`
    suppressWarnings({
      fieldname <- private$.getColnames("fieldname", schema, query, view)
      caption <- tolower(private$.getColnames("caption", schema, query, view))
    })

    # Fix the column filter with `caption`
    colFilter <- .fixColFilter(colFilter, old, caption, fieldname)

    # step 2:
    # Extract the column names in the modified column filter from step 1
    old <- .extractNames(colFilter)

    # Extract the lookups columns and the labels to fix the column names by
    lookups <- fieldname[grep("/", fieldname)]
    labels <- gsub("\\w+/", "", lookups)

    # Fix the colulm filter with lookups columns
    # (e.g., "DataSets/demographics/age_reported")
    colFilter <- .fixColFilter(colFilter, old, labels, lookups)

    colFilter
  }
)



# HELPER -----------------------------------------------------------------------

# Extract the column names from the column filter
#' @importFrom utils URLdecode
.extractNames <- function(colFilter) {
  tolower(unlist(lapply(gsub("~.*$", "", colFilter), URLdecode)))
}


# Fix the column filter
#' @importFrom utils URLencode
.fixColFilter <- function(colFilter, old, labels, names) {
  if (any(old %in% labels)) {
    idx <- which(old %in% labels)
    new <- unlist(lapply(names[match(old[idx], labels)], URLencode))
    colFilter[idx] <- paste0(paste0(new, "~"), gsub("^.*~", "", colFilter[idx]))
  }

  colFilter
}


# Transform data
.transformData <- function(data, dataType, transformMethod = "auto") {
  if (transformMethod == "none") {
    return(data)
  }

  # List of datasets with transformable column names (current as of March 2019)
  # and default transformation methods (selected by Raphael)
  datasets <- list(
    "elisa" = list(column = "value_reported", method = "log"),
    "elispot" = list(column = "spot_number_reported", method = "log1p"),
    "fcs_analyzed_result" = list(column = "population_cell_number", method = "sqrt"),
    "hai" = list(column = "value_preferred", method = "log"),
    "neut_ab_titer" = list(column = "value_preferred", method = "log")
  )

  if (!dataType %in% names(datasets)) {
    warning(
      "'", dataType, "' is not a dataset that can be transformed. ",
      "`transformMethod` argument is ignored."
    )
    return(data)
  }

  if (!transformMethod %in% c("auto", "log", "log1p", "sqrt")) {
    warning(
      "'", transformMethod, "' is not a valid transformation method. ",
      "Please use 'log', 'log1p', or 'sqrt'. ",
      "`transformMethod` argument is ignored."
    )
    return(data)
  }

  # retrieve transformation function
  if (transformMethod == "auto") {
    transformMethod <- datasets[[dataType]][["method"]]
  }
  transformFunction <- get(transformMethod)

  # retrieve column name
  column <- datasets[[dataType]][["column"]]

  message(
    "Transformed '", column, "' column of '", dataType, "' dataset with ",
    "'", transformMethod, "' transformation function."
  )
  if (dataType == "fcs_analyzed_result") {
    data[, population_cell_number := as.numeric(population_cell_number)]
  }
  data[, (column) := lapply(.SD, function(x) transformFunction(x)),
    .SDcols = grep(column, colnames(data))
  ][]
}
