.onAttach <- function(libname, pkgname) {
  netrc <- ifelse(.Platform$OS.type == "windows", "~/_netrc", "~/.netrc")

  if (!file.exists(netrc) &&
      !exists("labkey.sessionCookieName") &&
      !exists("apiKey", where = Rlabkey:::.lkdefaults) &&
      Sys.getenv("IS_login") == "") {
    packageStartupMessage("A .netrc file is required to connect to ImmuneSpace. For more information on how to create one, refer to the Configuration section of the introduction vignette.")
  }
}
