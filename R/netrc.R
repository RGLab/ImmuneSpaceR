#' Interactively write a netrc file
#'
#' Write a netrc file that is valid for accessing ImmuneSpace
#'
#' @return A netrc file that is verified to connect to ImmuneSpace
#' @examples
#' \dontrun{
#' interactive_netrc()
#' }
#'
#' @export
interactive_netrc <- function() {
  # generate netrc path
  filepath <- .get_path()
  overwrite <- ""

  # check if netrc exists
  if (file.exists(filepath)) {
    message("A netrc file already exists!")
    message("***Printing existing netrc to console***")
    cat(readChar(filepath, nchars = 10000))
    cat("\n\n")
    overwrite <- readline(prompt = "Overwrite current netrc? [Y / n] ")
  }

  # write netrc
  if (toupper(overwrite) == "Y" | overwrite == "") {
    login <- readline("What is your ImmuneSpace login email?  ")
    password <- readline("What is your ImmuneSpace password?  ")
    chk <- .check_con(login, password)

    if (chk == TRUE) {
      message("writing netrc to ", filepath)
      cat("machine www.immunespace.org\nlogin ", login, "\npassword ", password, "\n", file = filepath)
    }
  } else {
    # don't overwrite - validate available netrc
    chk <- .check_con()
  }

  filepath
}


#' Write a netrc file
#'
#' Write a netrc file that is valid for accessing ImmuneSpace
#'
#' @param login A \code{character}. The email address used for loging in on
#'  ImmuneSpace.
#' @param password A \code{character}. The password associated with the login.
#' @param machine A \code{character}. The server to connect.
#' @param file A \code{character}. The credentials will be written into that
#'  file. If left NULL, the netrc will be written into a temporary file.
#' @return A character vector containing the file paths for netrc
#' @examples
#' write_netrc("immunespaceuser@gmail.com", "mypassword")
#' @export
write_netrc <- function(login,
                        password,
                        machine = "www.immunespace.org",
                        file = NULL) {
  string <- paste0(
    "machine ", machine, "\n",
    "login ", login, "\n",
    "password ", password, "\n"
  )
  if (is.null(file)) {
    file <- tempfile()
  } else if (file.exists(file)) {
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
#' @export
check_netrc <- function() {
  if (exists("labkey.netrc.file", .GlobalEnv)) {
    netrc_file <- labkey.netrc.file
  } else {
    netrc_file <- "~/.netrc"
  }
  if (!file.exists(netrc_file)) {
    stop("There is no netrc file. Use `write_netrc`")
  } else {
    print(paste("netrc file found at", netrc_file))
  }
  lines <- readLines(netrc_file)
  lines <- gsub("http.*//", "", lines)
  if (length(grep("machine\\swww.immunespace.org", lines)) == 0) {
    stop("No entry found for www.immunespace.org in the netrc file.")
  }
  print("The netrc looks valid.")
  return(netrc_file)
}


# Get (and create) temporary netrc file from environment variables
.get_env_netrc <- function() {
  ISR_login <- Sys.getenv("ISR_login")
  ISR_pwd <- Sys.getenv("ISR_pwd")
  ISR_machine <- ifelse(Sys.getenv("ISR_machine") == "",
    "www.immunespace.org",
    Sys.getenv("ISR_machine")
  )
  if (ISR_login != "" & ISR_pwd != "") {
    write_netrc(login = ISR_login, password = ISR_pwd, machine = ISR_machine)
  }
}


# Get labkey.url.base from environment variable
# Ensure secure connection for server
# Allow ISCon.R `.get_url.base` method to handle local
.get_env_url <- function() {
  machine <- Sys.getenv("ISR_machine")
  # if blank, then use production
  if (machine == ""){
    return("https://www.immunespace.org")

  # if local assume http, no ssl
  } else if (grepl("^10\\.", machine)){
    return(machine)

  # assume format for test / prod as "(test|www).immunespace.org"
  } else {
    return(paste0("https://", machine))
  }
}

# get the path to where a netrc file should be
.get_path <- function() {
  home <- Sys.getenv("HOME")
  netrc <- ifelse(.Platform[[1]] == "unix", ".netrc", "_netrc")
  filepath <- paste0(home, "/", netrc)
  return(filepath)
}

# check that user can connect to IS with netrc file
.check_con <- function(login = NULL, passwd = NULL) {
  # try to connect to IS -- if no connection return NA
  if (!is.null(login)) {
    message("Validating netrc ...")
    con <- tryCatch(CreateConnection(study = "SDY269", login = login, password = passwd),
      error = function(e) {
        return(NULL)
      }
    )
  } else {
    message("Validating netrc ...")
    con <- tryCatch(CreateConnection(study = "SDY269"),
      error = function(e) {
        return(NULL)
      }
    )
  }
  if (is.null(con)) {
    message("Cannot connect to ImmuneSpace with current netrc information -- check login and password for errors")
    return(FALSE)
  } else {
    message("Ability to connect to ImmuneSpace with netrc confirmed")
    return(TRUE)
  }
}
