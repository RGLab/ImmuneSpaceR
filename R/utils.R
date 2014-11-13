CreateConnection = function(study = NULL, verbose = FALSE){
  if(length(study) <= 1)
    .CreateConnection(study = study, verbose = verbose)
  else
  {
    conList <- sapply(study, .CreateConnection, verbose = verbose)
    .ISConList(connections = conList)
  }
    
  
}
#' Load an ImmuneSpaceConnection/ImmuneSpaceConnectionList object from disk
#' 
#' @rdname loadConnection
#' 
#' @param file the file name to be saved to or loaded from
#' @return ImmuneSpaceConnection or ImmuneSpaceConnectionList object
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
  if(gsub("https://", "", labkey.url.base) == "www.immunespace.org"){
    labkey.setCurlOptions(ssl.verifyhost = 2, ssl.cipher.list="ALL")
  } else{
    labkey.setCurlOptions(ssl.verifyhost = 2, sslversion=1)
  }
  con
}

#' Save an ImmuneSpaceConnection/ImmuneSpaceConnectionList object to disk
#' 
#' @rdname loadConnection
saveConnection <- function(con, file){
  saveRDS(con, file = file)
}

#' @importFrom gplots colorpanel
#' @export
ISpalette <- function(n){
  colorpanel(n, low = "#268bd2", mid = "#fdf6e3", high = "#dc322f")
}
