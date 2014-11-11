CreateConnection = function(study = NULL, verbose = FALSE){
  if(length(study) <= 1)
    .CreateConnection(study = study, verbose = verbose)
  else
  {
    conList <- sapply(study, .CreateConnection, verbose = verbose)
    .ISConList(connections = conList)
  }
    
  
}


#' @importFrom gplots colorpanel
#' @export
ISpalette <- function(n){
  colorpanel(n, low = "#268bd2", mid = "#fdf6e3", high = "#dc322f")
}
