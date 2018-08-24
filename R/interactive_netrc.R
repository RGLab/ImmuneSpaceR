

## Get OS
## Get login
## Get Password
## Print .netrc or _netrc based on operating system

fileConn <- file(filepath)
writeLines(c(login, password))
close(fileConn)

## TODO
## Figure out file paths for windows
## write .netrc to home directory
## put together in a full function
# 
# C:/Users/Tyler/Documents"

# machine www.immunespace.org
# login myUser@mySite.com
# password superSecretPassword

.get_path <- function() {
  os <- .Platform[[1]]
  home <- Sys.getenv("HOME")
  netrc <- ifelse (os == "unix", ".netrc", "_netrc")
  filepath <- paste0(home, "/", netrc)
  return(filepath)
}

## ADD CHECK NETRC FOR EXISTING FILE??

write_netrc_lw <- function() {
  # generate netrc path
  filepath <- .get_path()
  # check if netrc exists
  if(file.exists(filepath)) {
    check_netrc()[1]
    overwrite <- menu(c("yes", "no"), title = cat("A netrc file already exists! \nOverwrite?"))
  }
  if(overwrite != 2) {
    login <- readline("What is your ImmuneSpace login email?  ")
    password <- readline("What is your ImmuneSpace password   ")
    cat("machine www.immunespace.org\nlogin ", login, "\npassword ", password, "\n", file = filepath)
  }
  
}
