#' Write a netrc file
#'
#' Write a netrc file that is valid for accessing ImmuneSpace
#'
#' @param login A \code{character}. The email address used for loging in on
#'  ImmuneSpace.
#' @param password A \code{character}. The password associated with the login.
#' @param file A \code{character}. The credentials will be written into that
#'  file. If left NULL, the netrc will be written into a temporary file.
#' @export
#' @return A character vector containing the file paths for netrc
#' @examples
#' write_netrc("immunespaceuser@gmail.com", "mypassword")
#'
write_netrc <- function(login,
                        password,
                        machine = "www.immunespace.org",
                        file = NULL) {
  string <- paste("machine", machine,
                  "login", login,
                  "password", password)
  if (is.null(file)) {
    file <- tempfile()
  } else if(file.exists(file)) {
    stop("The file you are trying to write to already exists. Remove manually if you wish to overwrite.")
  }
  write(string, file)
  file
}

#' Check netrc file
#'
#' Check that there is a netrc file with a valid entry for ImmuneSpace.
#'
#' @return The name of the netrc file
#'
#' @details
#' In order to connect to ImmuneSpace, you will need a `.netrc` file in your
#' contains a `machine` name (hostname of ImmuneSpace), and `login` and
#' `password`. See [here](https://www.labkey.org/wiki/home/Documentation/page.view?name=netrc)
#' for more information. By default \code{RCurl} will look for the file in your
#' home directoty.
#'
#' If no netrc is available or it is not formatted properly, \code{write_netrc}
#' can be used to write one. Otherwise, when specifying login and password in
#' \code{CreateConnection}, a temporary file will be created for that connection.
#'
#' @seealso CreateConnection write_netrc
#' @examples
#' try(check_netrc())
#'
#' @export
check_netrc <- function(){
  if(exists("labkey.netrc.file", .GlobalEnv)){
    netrc_file <- labkey.netrc.file
  } else{
    netrc_file <- "~/.netrc"
  }
  if(!file.exists(netrc_file)){
    stop("There is no netrc file. Use `write_netrc`")
  } else{
    print(paste("netrc file found at", netrc_file))
  }
  lines <- readLines(netrc_file)
  lines <- gsub("http.*//", "", lines)
  if(length(grep("machine\\swww.immunespace.org", lines)) == 0){
    stop("No entry found for www.immunespace.org in the netrc file.")
  }
  print("The netrc looks valid.")
  return(netrc_file)
}

# Get (and create) temporary netrc file from environment variables
get_env_netrc <- function() {
  ISR_login <- Sys.getenv("ISR_login")
  ISR_pwd <- Sys.getenv("ISR_pwd")
  ISR_machine <- ifelse(Sys.getenv("ISR_machine") == "",
                        "www.immunespace.org",
                        Sys.getenv("ISR_machine"))
  if (ISR_login != ""  &  ISR_pwd != "") {
    write_netrc(login = ISR_login, password = ISR_pwd, machine = ISR_machine)
  }
}

# Get labkey.url.base from environment variable
get_env_url <- function() {
  ifelse(Sys.getenv("ISR_machine") == "",
         "https://www.immunespace.org",
         paste0("https://", Sys.getenv("ISR_machine")))
}
