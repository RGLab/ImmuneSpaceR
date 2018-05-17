#' @title CreateConnection
#'
#' @name CreateConnection
#'
#' @param study A \code{"character"} vector naming the study.
#' @param login A \code{"character"}. Optional argument. If there is no netrc
#'  file a temporary one can be written by passing login and password of an
#'  active ImmuneSpace account.
#' @param password A \code{"character"}. Optional. The password for the selected
#'  login.
#' @param use.data.frame A \code{"logical"}. If set to TRUE, the functions will
#'  return \code{data.frame} objects instead of \code{data.table}.
#' @param verbose A \code{"logical"} whether to print the extra details for
#' troubleshooting.
#' @param onTest A \code{"logical"} whether to connect to the test server
#' (https://test.immunespace.org/) instead of the production server
#' (https://www.immunespace.org/).
#'
#' @description Constructor for \code{ImmuneSpaceConnection} class.
#'
#' @details Instantiates an \code{ImmuneSpaceConnection} for \code{study}
#' The constructor will try to take the values of the various `labkey.*`
#' parameters from the global environment. If they don't exist, it will use
#' default values. These are assigned to `options`, which are then used by the
#' \code{ImmuneSpaceConnection} class.
#'
#' @return an instance of an \code{ImmuneSpaceConnection}
#'
#' @seealso \code{\link{ImmuneSpaceConnection}}
#'
#' @examples
#' \dontrun{
#' # Single study
#' con <- CreateConnection("SDY269")
#' # Cross study
#' con <- CreateConnection("")
#' }
#'
#' sdy <- try(CreateConnection("SDY269"))
#' if (inherits(sdy, "try-error")) {
#'   print("Read the Introduction vignette for more information on how to set up
#'   a .netrc file.")
#' }
#'
#' @export
#' @importFrom utils packageVersion
CreateConnection <- function(study = NULL,
                             login = NULL,
                             password = NULL,
                             use.data.frame = FALSE,
                             verbose = FALSE,
                             onTest = FALSE) {
  # check internet connection
  if (!curl::has_internet()) {
    stop("No internet connection. Please connect to internet and try again.")
  }

  # check if the portal is up
  url <- ifelse(onTest, "test.immunespace.org", "www.immunespace.org")
  if (is.null(curl::nslookup(url, error = FALSE))) {
    stop("The portal is currently down. Try again later.")
  }


  # Try to parse labkey options from global environment
  # which really should have been done through option()/getOption() mechanism
  # Here we do this to be compatible to labkey online report system
  # that automatically assigns these variables in global environment
  labkey.url.base <- try(get("labkey.url.base", .GlobalEnv), silent = TRUE)
  if (inherits(labkey.url.base, "try-error")) {
    labkey.url.base <- paste0("https://", url)
  }

  labkey.url.base <- gsub("http:", "https:", labkey.url.base)
  if (length(grep("^https://", labkey.url.base)) == 0) {
    labkey.url.base <- paste0("https://", labkey.url.base)
  }

  labkey.user.email <- try(get("labkey.user.email", .GlobalEnv), silent = TRUE)
  if (inherits(labkey.user.email, "try-error")) {
    labkey.user.email <- "unknown_user at not_a_domain.com"
  }

  # set curoption for Rlabkey package
  #
  # Rlabkey stores the Curl options in its package environment through labkey.setCurlOptions call.
  # So in theory we need to reset it prior to each Rlabkey query
  # because  multiple connections created by user indiviudally (not as ImmuneSystemConnectionList)
  # may have different urls and ssl settings.
  # (Ideally labkey.selectRows should optionally parse the options from its argument besides package environment)
  #
  # for now we assume they all share the same setting and init it only once here
  if (!is.null(login) & is.null(password)) {
    stop(
      "login = ",
      login,
      " given without password. Please try again with password"
    )
  } else if (!is.null(login) & !is.null(password)) {
    nf <- write_netrc(login, password)
  } else {
    nf <- try(get("labkey.netrc.file", .GlobalEnv), silent = TRUE)
  }

  useragent <- paste0(
    "R/", R.version$major, ".", R.version$minor,
    " (", Sys.info()["sysname"], " ", Sys.info()["machine"], ")",
    " Rlabkey/", packageVersion("Rlabkey"),
    " ImmuneSpaceR/", packageVersion("ImmuneSpaceR")
  )

  if (!inherits(nf, "try-error") && !is.null(nf)) {
    curlOptions <- labkey.setCurlOptions(
      ssl_verifyhost = 2,
      sslversion = 1,
      netrc_file = nf,
      useragent = useragent
    )
  } else {
    curlOptions <- labkey.setCurlOptions(
      ssl_verifyhost = 2,
      sslversion = 1,
      useragent = useragent
    )
  }

  if (length(study) <= 1) {
    .CreateConnection(
      study = study,
      labkey.url.base = labkey.url.base,
      labkey.user.email = labkey.user.email,
      use.data.frame = use.data.frame,
      verbose = verbose,
      curlOptions = curlOptions
    )
  } else {
    stop("For multiple studies, use an empty string and filter the connection.")
  }
}


.CreateConnection <- function(study = NULL,
                              labkey.url.base,
                              labkey.user.email,
                              use.data.frame,
                              curlOptions,
                              verbose,
                              ...) {
  labkey.url.path <- try(
    get("labkey.url.path", .GlobalEnv),
    silent = TRUE
  )

  if (inherits(labkey.url.path,"try-error")) {
    if (is.null(study)) {
      stop("study cannot be NULL")
    }
    pathStr <- ifelse(
      grepl("^IS\\d{1,3}$", study),
      "/HIPC/",
      "/Studies/"
    )
    labkey.url.path <- paste0(pathStr, study)
  } else if (!is.null(study)) {
    labkey.url.path <- file.path(dirname(labkey.url.path), study)
  }

  if (class(use.data.frame) != "logical") {
    warning("use.data.frame should be of class `logical`. Setting it to FALSE.")
    use.data.frame <- FALSE
  }

  # check credential
  if (verbose) message("Checking credential...")
  res <- GET(
    url = paste0(labkey.url.base, "/login-whoami.view"),
    config = Rlabkey:::labkey.getRequestOptions()
  )
  if (res$status_code == 200) {
    if (grepl("json", res$headers$`content-type`)) {
      parsed <- httr::content(res)

      if (parsed$displayName == "guest") {
        stop("Invalid credential or deactivated account. Check your account in the portal.")
      }
    } else {
      stop("Something went wrong. Check the portal and try again.")
    }
  } else if (res$status_code == 401) {
    stop("Invalid credential or deactivated account. Check your account in the portal.")
  } else if (res$status_code == 403) {
    stop("The portal is in admin-only mode. Please try again later.")
  } else {
    stop("Something went wrong. Check the portal and try again.")
  }

  config <- list(
    labkey.url.base = labkey.url.base,
    labkey.url.path = labkey.url.path,
    labkey.user.email = labkey.user.email,
    use.data.frame = use.data.frame,
    curlOptions = curlOptions,
    verbose = verbose
  )

  ISCon$new(config = config)
}
