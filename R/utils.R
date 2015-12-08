#' Load an ImmuneSpaceConnection object from disk
#' 
#' @rdname loadConnection
#' 
#' @param file the file name to be saved to or loaded from
#' @return An ImmuneSpaceConnection object
#' @export
loadConnection <- function(file){
  con <- readRDS(file = file)
  conType <- class(con)
  if(conType == 'ImmuneSpaceConnection') 
    labkey.url.base <- con$config$labkey.url.base
  else
    stop("invalid ImmuneSpaceConnection object!")
  
  #init labkey.setCurlOptions
  labkey.setCurlOptions(ssl.verifyhost = 2, sslversion=1)
  con
}

#' Save an ImmuneSpaceConnection object to disk
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