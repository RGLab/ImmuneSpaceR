#'@title CreateConnection
#'@name CreateConnection
#'@param study \code{"character"} vector naming the study.
#' @param verbose \code{"logical"} wehther to print the extra details for troubleshooting. 
#'@description Constructor for \code{ImmuneSpaceConnection} class
#'@details Instantiates and \code{ImmuneSpaceConnection} for \code{study}
#'The constructor will try to take the values of the various `labkey.*` parameters from the global environment.
#'If they don't exist, it will use default values. These are assigned to `options`, which are then used by the \code{ImmuneSpaceConnection} class.
#'@export CreateConnection
#'@return an instance of an \code{ImmuneSpaceConnection} or \code{ImmuneSpaceConnectionList}
CreateConnection = function(study = NULL, verbose = FALSE){
  # try to parse labkey options from global environment 
  # which really should have been done through option()/getOption() mechanism
  # Here we do this to be compatible to labkey online report system 
  # that automatically assigns these variables in global environment
  labkey.url.base<-try(get("labkey.url.base",.GlobalEnv),silent=TRUE)
  if(inherits(labkey.url.base,"try-error"))
    labkey.url.base<-"https://www.immunespace.org"
  labkey.url.base<-gsub("http:","https:",labkey.url.base)
  if(length(grep("^https://", labkey.url.base)) == 0)
    labkey.url.base <- paste0("https://", labkey.url.base)
  labkey.user.email<-try(get("labkey.user.email",.GlobalEnv),silent=TRUE)
  if(inherits(labkey.user.email,"try-error"))
    labkey.user.email="unknown_user at not_a_domain.com"
  
  # set options to avoid the dealing with .GlobalEnv
  
  
  
  
  # set curoption for Rlabkey package
  #
  # Rlabkey stores the Curl options in its package environment through labkey.setCurlOptions call.
  # So in theory we need to reset it prior to each Rlabkey query 
  # because  multiple connections created by user indiviudally (not as ImmuneSystemConnectionList)
  # may have different different urls and ssl settings. 
  # (Ideally labkey.selectRows should optionally parse the options from its argument besides package environment)
  # 
  # for now we assume they all share the same setting and init it only once here
  curlOptions <- labkey.setCurlOptions(ssl.verifyhost = 2, sslversion=1)
  
  if(length(study) <= 1)
    .CreateConnection(study = study
                      , labkey.url.base=labkey.url.base
                      , labkey.user.email=labkey.user.email
                      , verbose = verbose
                      , curlOptions = curlOptions
                      )
  else
  {
    conList <- sapply(study
                      , .CreateConnection
                      , labkey.url.base=labkey.url.base
                      , labkey.user.email=labkey.user.email
                      , verbose = verbose
                      , curlOptions = curlOptions
                      )
    
    .ISConList(connections = conList)
  }
    
  
}
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
#' @rdname loadConnection
#' @export
saveConnection <- function(con, file){
  saveRDS(con, file = file)
}

#' @importFrom gplots colorpanel
#' @export
ISpalette <- function(n){
  colorpanel(n, low = "#268bd2", mid = "#fdf6e3", high = "#dc322f")
}
