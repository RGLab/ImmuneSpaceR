#' Load an ImmuneSpaceConnection/ImmuneSpaceConnectionList object from disk
#' 
#' @rdname loadConnection
#' 
#' @param file the file name to be saved to or loaded from
#' @return ImmuneSpaceConnection or ImmuneSpaceConnectionList object
#'@export
loadConnection <- function(file){
  con <- readRDS(file = file)
  conType <- class(con)
  if(conType == 'ImmuneSpaceConnection') 
    labkey.url.base <- con$config$labkey.url.base
  else if(conType == 'ImmuneSpaceConnectionList')
    labkey.url.base <- con$connections[[1]]$config$labkey.url.base
  else
    stop("invalid ImmuneSpaceConnection object!")
  
  #init labkey.setCurlOptions
  labkey.setCurlOptions(ssl.verifyhost = 2, sslversion=1)
  con
}

#' Save an ImmuneSpaceConnection/ImmuneSpaceConnectionList object to disk
#' 
#' @param con An \code{ImmuneSpaceConnection}. The connection to save to file. 
#'  To be loaded later using \code{loadConnection}.
#' 
#' @rdname loadConnection
#' @export
saveConnection <- function(con, file){
  saveRDS(con, file = file)
}

#' ImmuneSpace palette
#' 
#' Create a color gradient of the selected length that matches the ImmuneSpace
#' theme.
#' 
#' @param n A \code{numeric}. The length of the desired palette.
#' @return A \code{character} vector colors in hexadecimal code of length
#'  \code{n}.
#' 
#' @importFrom gplots colorpanel
#' @export
ISpalette <- function(n){
  colorpanel(n, low = "#268bd2", mid = "#fdf6e3", high = "#dc322f")
}

#' @importFrom RCurl curlUnescape curlEscape
.check_filter <- function(lub, lup, schema, query, view = "", colFilter){
  # Get the names used in the filter
  old <- tolower(curlUnescape(gsub("~.*$", "", colFilter)))
  # 
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
