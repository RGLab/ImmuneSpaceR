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
      for (i in 1:nrow(self$availableDatasets)) {
        cat(sprintf("\t%s\n", self$availableDatasets[i, Name]))
      }
    }

    if ("expression" %in% output) {
      if (!is.null(self$cache[[private$.constants$matrices]])) {
        cat("Expression Matrices\n")
        for (i in 1:nrow(self$cache[[private$.constants$matrices]])) {
          cat(sprintf("\t%s\n", self$cache[[private$.constants$matrices]][i, name]))
        }
      } else {
        cat("No Expression Matrices Available")
      }
    }
  }
)


# Downloads a dataset and cache the result in the connection object.
ISCon$set(
  which = "public",
  name = "getDataset",
  value = function(x, original_view = FALSE, reload = FALSE, colFilter = NULL, transformMethod = "none", ...) {
    if (nrow(self$availableDatasets[Name%in%x]) == 0) {
      wstring <- paste0(
        "Empty data frame was returned.",
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
      return(data.frame())
    }

    cache_name <- paste0(x, ifelse(original_view, "_full", ""))
    nOpts <- length(list(...))

    if (!is.null(self$cache[[cache_name]]) &&
        !reload &&
        is.null(colFilter) &&
        nOpts == 0) { # Serve cache
      data <- self$cache[[cache_name]]
      # if(!is.null(colFilter)) {
      #   data <- .filterCachedCopy(colFilter, data)
      #   return(data)
      # } else {
      # }

    } else { # Download the data
      viewName <- NULL
      if (original_view) {
        viewName <- "full"
      }

      if (!is.null(colFilter)) {
        colFilter <- private$.checkFilter(
          schema = "study",
          query = x,
          colFilter = colFilter,
          view = viewName
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

      noTrx <- c("pcr", "mbaa", "hla_typing", "kir_typing", "gene_expression_files")

      if (transformMethod != "none" && !(x %in% noTrx)) {

        if (x == "fcs_analyzed_result"){
          data[, population_cell_number := as.numeric(population_cell_number)]
        }

        # colNames current as of 8/2018
        if (!(x %in% c("elispot", "fcs_analyzed_result"))) {
          cNm <- "value_preferred"
        } else if (x == "elispot") {
          cNm <- "spot_number_reported"
        } else if (x == "fcs_analyzed_result"){
          cNm <- "population_cell_number"
        }

        if (transformMethod == "auto") {
          # Transformation options selected by RG
          if (x %in% c("hai", "neut_ab_titer", "elisa")) {
            tFun <- log
          } else if (x %in% c("elispot")) {
            tFun <- log1p
          } else if (x %in% c("fcs_analyzed_result")){
            tFun <- sqrt
          }

        }else{
          if (transformMethod %in% c("log", "log1p", "sqrt")) {
            tFun <- get(transformMethod)
            data[, (cNm) := lapply(.SD, function(x) tFun(x)), .SDcols=grep(cNm, colnames(data))]
          } else {
            warning(paste0("'transformMethod' ", transformMethod, " not recognized. Please use 'log', 'log1p', or 'sqrt'. 'transformMethod' ignored." ))
          }
        }

      } else if (transformMethod != "none" && x %in% noTrx) {
        message(paste0(x, " is not a dataset that can be transformed. \n 'transformMethod' ignored."))
      }

      if (cache) {
        self$cache[[cache_name]] <- data
      }
    }

    if (!is.null(self$config$use.data.frame) && self$config$use.data.frame) {
      data <- data.frame(data)
    }

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
      .getLKtbl(
        con = self,
        schema = "study",
        query = "ISC_study_datasets"
      )
    }
  }
)


# Decode filters in case the user used column labels instead of column names
#' @importFrom utils URLdecode URLencode
ISCon$set(
  which = "private",
  name = ".checkFilter",
  value = function(schema, query, colFilter, view = "") {
    ## HELPERS
    extractNames <- function(colFilter) {
      tolower(unlist(lapply(gsub("~.*$", "", colFilter), URLdecode)))
    }

    getColnames <- function(colNameOpt) {
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

    fixColFilter <- function(colFilter, labels, names) {
      if (any(old %in% labels)) {
        idx <- which(old %in% labels)
        new <- unlist(lapply(names[match(old[idx], labels)], URLencode))
        colFilter[idx] <- paste0(paste0(new, "~"), gsub("^.*~", "", colFilter[idx]))
      }

      colFilter
    }


    ## MAIN
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
)



# HELPER -----------------------------------------------------------------------

# Filter the data table object by column
# TODO: Need a way to translate Curl operators into math
# (so that we can filter the saved tables using colFilter).
.filterDT <- function(table, col, value) {
  table[eval(as.name(col)) == value]
}


# Filter data by the prodivded filters
.filterCachedCopy <- function(filters, data) {
  decoded <- unlist(lapply(filters, URLdecode))

  cols <- gsub(".*/", "", gsub("~.*", "", decoded))
  values <- gsub(".*=", "", decoded)

  for (i in 1:length(filters)) {
    data <- .filterDT(table, cols[i], values[i])
  }

  data
}
