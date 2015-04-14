.onLoad <- function(libname, pkgname){
  if(.Platform$OS.type == "windows")#set ca bundle file path for windows
    options(RCurlOptions = list(cainfo = system.file("ssl_certs/ca-bundle.crt", package = pkgname)))  
}
