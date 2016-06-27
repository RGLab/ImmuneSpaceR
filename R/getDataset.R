# Decode filters in case the user used column labels instead of column names
#' @importFrom RCurl curlUnescape curlEscape
.check_filter <- function(lub, lup, schema, query, view = "", colFilter){
  # Get the names used in the filter
  old <- tolower(curlUnescape(gsub("~.*$", "", colFilter)))
  suppressWarnings({
    labels <- tolower(colnames(labkey.selectRows(lub, lup, schema, query, view, maxRows = 0, colNameOpt = "caption")))
    names <- colnames(labkey.selectRows(lub, lup, schema, query, view, maxRows = 0, colNameOpt = "fieldname"))
  })
  # Get the new names
  idx <- which(old %in% labels)
  new <- curlEscape(names[match(old[idx], labels)])
  colFilter[idx] <- paste0(paste0(new, "~"), gsub("^.*~", "", colFilter[idx]))
  return(colFilter)
}

## TODO: Need a way to translate Curl operators into math (so that we can filter the saved tables using colFilter).
filterDT <- function(table, col, value){
  table <- table[eval(as.name(col)) == value]
  return(table)
}
filter_cached_copy <- function(filters, data){
  decoded <- curlUnescape(filters)
  cols <- gsub(".*/", "", gsub("~.*", "", decoded))
  values <- gsub(".*=", "", decoded)
  for(i in 1:length(filters)){
    data <- filterDT(table, cols[i], values[i])
  }
  return(data)
}

# Downloads a dataset and cache the result in the connection object.
.ISCon$methods(
  getDataset = function(x, original_view = FALSE, reload = FALSE, colFilter = NULL, ...){
    "Get a dataset form the connection\n
    original_view: A logical. If set tot TRUE, download the ImmPort view.
    Else, download the default grid view.\n
    reload: A logical. Clear the cache. If set to TRUE, download the dataset,
    whether a cached version exist or not.\n
    colFilter: A character. A filter as returned by Rlabkey's makeFilter function.\n
    '...': Extra arguments to be passed to labkey.selectRows."
    if(nrow(available_datasets[Name%in%x])==0){
      wstring <- paste0(study, " has invalid data set: ",x)
      if(config$verbose){
        wstring <- paste0(wstring, "\n",
                          "Valid datasets for ", study, ": ",
                          paste(available_datasets$Name, collapse = ", "), ".")
      }
      warning(wstring)
      NULL
    } else{
      cache_name <- paste0(x, ifelse(original_view, "_full", ""))
      if(!is.null(data_cache[[cache_name]]) & !reload & is.null(colFilter)){ # Serve cache
        data <- data_cache[[cache_name]]
        #if(!is.null(colFilter)){
        #  data <- filter_cached_copy(colFilter, data)
        #  return(data)
        #} else{
        #}
      } else{ # Download the data
        viewName <- NULL
        if(original_view){
          viewName <- "full"
        }
        if(!is.null(colFilter)){
          colFilter <- .check_filter(config$labkey.url.base, 
                                     config$labkey.url.path, 
                                     "study", x, viewName, colFilter)
          cache <- FALSE
        } else{
          cache <- TRUE
        }
        data <- data.table(
          labkey.selectRows(baseUrl = config$labkey.url.base,
                            config$labkey.url.path,
                            schemaName = "study",
                            queryName = x,
                            viewName = viewName,
                            colNameOpt = "caption",
                            colFilter = colFilter,
                            ...))
        setnames(data, .self$.munge(colnames(data)))
        if(cache){
          data_cache[[cache_name]] <<- data
        }
      }
      if(!is.null(config$use.data.frame) & config$use.data.frame){
        data <- data.frame(data)
      }
      return(data)
    }
  })