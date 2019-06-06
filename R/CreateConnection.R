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
#'   warning("Read the Introduction vignette for more information on how to set
#' up a .netrc file.")
#' }
#' @export
#' @importFrom utils packageVersion
#' @importFrom curl has_internet nslookup
CreateConnection <- function(study = NULL,
                             login = NULL,
                             password = NULL,
                             verbose = FALSE,
                             onTest = FALSE) {
  ISCon$new(study, login, password, verbose, onTest)
}
