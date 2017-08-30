.onLoad <- function(libname, pkgname) {
  def_netrc <- ifelse(.Platform$OS.type == "windows", "~/_netrc", "~/.netrc")
  if (!file.exists(def_netrc)) {
    packageStartupMessage("A .netrc file is required to connect to ImmuneSpace. For more information on how to create one, refer to the Configuration section of the introduction vignette.")
  }

  if (.Platform$OS.type == "windows") # set ca bundle file path for windows
    options(RCurlOptions = list(cainfo = system.file("ssl_certs/ca-bundle.crt",
                                                     package = pkgname)))
}
