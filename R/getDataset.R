#' @title Get a dataset
#' 
#' @description 
#' Downloads a dataset and cache the result in the connection object.
#' 
#' con$getDataset(x, original_view = FALSE, reload = FALSE, ...)
#' 
#' @param x A \code{character}. The name of the dataset.
#' @param original_view A \code{logical}. If set to TRUE, download the ImmPort
#'  view. Else, download the default grid view. Note: Once data is cached,
#'  changing value of this argument won't have effect on the subsequent calls
#'  unless \code{reload} is set to 'TRUE'.
#' @param reload A \code{logical}. Clear the cache. If set to TRUE, download the 
#'  dataset, whether a cached version exist or not.
#' @param ... Arguments to be passed to the underlying \code{labkey.selectRows}.
#' @details
#' Returns the dataset named 'x', downloads it if it is not already cached. Note
#' that if additional arguments (...) are passed, the dataset will be reloaded,
#' even if a cached copy exist.
#' 
#' @return a \code{data.table}
#' @name ImmuneSpaceConnection_getDataset
#' @examples
#' labkey.url.base = "https://www.immunespace.org"
#' labkey.url.path = "/Studies/SDY269"
#' sdy269 <- CreateConnection("SDY269")
#' sdy269$getDataset("hai")
.ISCon$methods(
  getDataset = function(x, original_view = FALSE, reload = FALSE, colFilter = NULL, ...){
    "Get a dataset form the connection\n
    original_view: A logical. If set tot TRUE, download the ImmPort view.
    Else, download the default grid view.\n
    reload: A logical. Clear the cache. If set to TRUE, download the dataset,
    whether a cached version exist or not."
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
      filter_state <- paste0("filter_state_", cache_name)
      if(!is.null(data_cache[[cache_name]]) &
         !reload &
         is.null(colFilter) &
         (is.null(data_cache[[filter_state]]) || !data_cache[[filter_state]])){
        data_cache[[cache_name]]
      } else{
        viewName <- NULL
        if(original_view){
          viewName <- "full"
        }
        if(!is.null(colFilter)){
          colFilter <- .check_filter(config$labkey.url.base, 
                                     config$labkey.url.path, 
                                     "study", x, viewName, colFilter)
          data_cache[[filter_state]] <<- TRUE 
        } else{
          data_cache[[filter_state]] <<- FALSE
        }
        data_cache[[cache_name]] <<- data.table(
          labkey.selectRows(baseUrl = config$labkey.url.base,
                            config$labkey.url.path,
                            schemaName = "study",
                            queryName = x,
                            viewName = viewName,
                            colNameOpt = "caption",
                            colFilter = colFilter,
                            ...))
        setnames(data_cache[[cache_name]],
                 .self$.munge(colnames(data_cache[[cache_name]])))
        data_cache[[cache_name]]
      }
    }
  })