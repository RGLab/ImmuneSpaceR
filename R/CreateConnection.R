################################################################################
# Main class and instantiation
################################################################################

#' @title CreateConnection
#' @name CreateConnection
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
#' @description Constructor for \code{ImmuneSpaceConnection} class
#' @details Instantiates an \code{ImmuneSpaceConnection} for \code{study}
#' The constructor will try to take the values of the various `labkey.*`
#' parameters from the global environment. If they don't exist, it will use
#' default values. These are assigned to `options`, which are then used by the
#' \code{ImmuneSpaceConnection} class.
#' @export CreateConnection
#' @return an instance of an \code{ImmuneSpaceConnection}
#' @seealso ImmuneSpaceConnection
#' @examples
#' \dontrun{
#'   # Single study
#'   con <- CreateConnection("SDY269")
#'   # Cross study
#'   con <- CreateConnection("")
#' }
#'
#' sdy <- try(CreateConnection("SDY269"))
#' if(inherits(sdy, "try-error")){
#'   print("Read the Introduction vignette for more information on how to set up
#'   a .netrc file.")
#' }
#' @importFrom utils packageVersion
CreateConnection <- function(study = NULL,
                             login = NULL,
                             password = NULL,
                             use.data.frame = FALSE,
                             verbose = FALSE,
                             onTest = FALSE) {
  # Try to parse labkey options from global environment
  # which really should have been done through option()/getOption() mechanism
  # Here we do this to be compatible to labkey online report system
  # that automatically assigns these variables in global environment
  labkey.url.base <- try(get("labkey.url.base", .GlobalEnv), silent = TRUE)
  if (inherits(labkey.url.base, "try-error")){
    labkey.url.base <- ifelse(onTest,
                              "https://test.immunespace.org",
                              "https://www.immunespace.org")
  }
  labkey.url.base <- gsub("http:", "https:", labkey.url.base)
  if (length(grep("^https://", labkey.url.base)) == 0){
    labkey.url.base <- paste0("https://", labkey.url.base)
  }
  labkey.user.email <- try(get("labkey.user.email", .GlobalEnv), silent = TRUE)
  if (inherits(labkey.user.email, "try-error")){
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
    stop(paste0("login = ",
                login,
                " given without password. Please try again with password"))
  } else if (!is.null(login) & !is.null(password)) {
    nf <- write_netrc(login, password)
  } else {
    nf <- try(get("labkey.netrc.file", .GlobalEnv), silent = TRUE)
  }

  useragent <- paste("ImmuneSpaceR", packageVersion("ImmuneSpaceR"))
  if (!inherits(nf, "try-error") && !is.null(nf)) {
    curlOptions <- labkey.setCurlOptions(ssl_verifyhost = 2,
                                         sslversion = 1,
                                         netrc_file = nf,
                                         useragent = useragent)
  } else {
    curlOptions <- labkey.setCurlOptions(ssl_verifyhost = 2,
                                         sslversion = 1,
                                         useragent = useragent)
  }

  if (length(study) <= 1) {
    .CreateConnection(study = study,
                      labkey.url.base = labkey.url.base,
                      labkey.user.email = labkey.user.email,
                      use.data.frame = use.data.frame,
                      verbose = verbose,
                      curlOptions = curlOptions)
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
  labkey.url.path <- try(get("labkey.url.path", .GlobalEnv), silent = TRUE)
  if (inherits(labkey.url.path,"try-error")) {
    if (is.null(study)) { stop("study cannot be NULL") }
    pathStr <- ifelse( grepl("^IS\\d{1,3}$", study),
                       "/HIPC/",
                       "/Studies/")
    labkey.url.path <- paste0(pathStr, study)
  } else if (!is.null(study)) {
    labkey.url.path <- file.path(dirname(labkey.url.path),study)
  }

  if (class(use.data.frame) != "logical") {
    warning("use.data.frame should be of class `logical`. Setting it to FALSE.")
    use.data.frame <- FALSE
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

#'@rdname ImmuneSpaceConnection-class
#'@aliases
#'ImmuneSpace
#'ImmuneSpaceConnection
#'ImmuneSpaceConnection-class
#'getGEMatrix
#'getDataset
#'listDatasets
#'getGEAnalysis
#'listGEAnalysis
#'addTrt
#'EMNames
#'quick_plot
#'@title The ImmuneSpaceConnection class
#'
#'@description
#'A connection respresents a study or a set of studies available on ImmuneSpace.
#'It provides function to download and display the data within these studies.
#'
#'@field study A \code{character}. The study accession number. Use an empty
#'string ("") to create a connection at the project level.
#'@field config A \code{list}. Stores configuration of the connection object
#'such as URL, path and username.
#'@field available_datasets A \code{data.table}. The table of datasets available
#'in the connection object.
#'@field data_cache A \code{list}. Stores the data to avoid downloading the same
#'tables multiple times.
#'@field constants A \code{list}. Used to store information regarding
#'gene-expression data.
#'
#'@details
#' Uses global variables \code{labkey.url.base}, and \code{labkey.url.path}, to
#' access a study. \code{labkey.url.base} should be
#' \code{https://www.immunespace.org/}. \code{labkey.url.path} should be
#' \code{/Studies/studyname}, where 'studyname' is the accession number of the
#' study.
#' The ImmunespaceConnection will initialize itself, and look for a
#' \code{.netrc} file in \code{"~/"} the user's home directory. The
#' \code{.netrc} file should contain a \code{machine}, \code{login}, and
#' \code{password} entry to allow access to ImmuneSpace, where \code{machine} is
#' the host name like "www.immunespace.org".
#'
#' @seealso
#'  \code{\link{CreateConnection}}
#'  \code{\link{ImmuneSpaceR-package}}
#' @exportClass ImmuneSpaceConnection
#' @examples
#' \dontrun{
#'   sdy269 <- CreateConnection("SDY269")
#'   sdy269
#' }
#'@return An instance of an ImmuneSpaceConnection for a study in `labkey.url.path`
ISCon <- R6Class(
  classname = "ImmuneSpaceConnection",
  public = list(
    study = character(),
    config = list(),
    available_datasets = data.table(),
    data_cache = list(),
    constants = list()
  )
)

# Functions used in initialize need to be declared ahead of it
#' @importFrom gtools mixedsort
ISCon$set(
  which = "public",
  name = "checkStudy",
  value = function(verbose = FALSE) {
    sdyNm <- basename(self$config$labkey.url.path)
    dirNm <- dirname(self$config$labkey.url.path)
    gTerm <- ifelse(dirNm == "/HIPC", "^IS\\d{1,3}$", "^SDY\\d{2,4}$")

    # adjust for "" connection
    if (sdyNm == "Studies") {
      sdyNm <- ""
      dirNm <- "/Studies"
    }

    folders <- labkey.getFolders(self$config$labkey.url.base, dirNm)
    subdirs <- gsub(paste0(dirNm, "/"), "", folders$folderPath)
    validSdys <- mixedsort(subdirs[grep(gTerm, subdirs)])

    if (!(sdyNm %in% c("", validSdys))) {
      if (verbose == FALSE) {
        stop(paste0(sdyNm, " is not a valid study. \n Use `verbose = TRUE` to see list of valid studies."))
      } else {
        stop(paste0(sdyNm, " is not a valid study\nValid studies: ",
                    paste(validSdys, collapse=", ")))
      }
    }
  }
)

ISCon$set(
  which = "public",
  name = "setAvailableDatasets",
  value = function() {
    if (length(self$available_datasets) == 0) {
      .getLKtbl(
        con = self,
        schema = "study",
        query = "ISC_study_datasets"
      )
    }
  }
)

ISCon$set(
  which = "public",
  name = "GeneExpressionMatrices",
  value = function(verbose = FALSE) {
    getData <- function() {
      try(
        .getLKtbl(
          con = self,
          schema = "assay.ExpressionMatrix.matrix",
          query = "Runs",
          colNameOpt = "fieldname",
          viewName = "expression_matrices"
        ),
        silent = TRUE
      )
    }

    if (!is.null(self$data_cache[[self$constants$matrices]])) {
      self$data_cache[[self$constants$matrices]]
    } else {
      if (verbose) {
        ge <- getData()
      } else {
        ge <- suppressWarnings(getData())
      }

      if (inherits(ge, "try-error") || nrow(ge) == 0) {
        # No assay or no runs
        message("No gene expression data")
        self$data_cache[[self$constants$matrices]] <- NULL
      } else {
        # adding cols to allow for getGEMatrix() to update
        ge[, annotation := ""]
        ge[, outputType := ""]
        setnames(ge, self$.munge(colnames(ge)))
        self$data_cache[[self$constants$matrices]] <- ge
      }
    }

    return(self$data_cache[[self$constants$matrices]])
  }
)

ISCon$set(
  which = "public",
  name = "initialize",
  value = function(..., config = NULL) {

    # invoke the default init routine in case it needs to be invoked
    # (e.g. when using $new(object) to construct the new object based on the existing object)
    # callSuper(...) # THIS MIIGHT NOT APPLICABLE IN R6

    self$constants <- list(
      matrices = "GE_matrices",
      matrix_inputs = "GE_inputs"
    )

    if (!is.null(config)) {
      self$config <- config
    }

    self$study <- basename(config$labkey.url.path)

    self$checkStudy(self$config$verbose)

    self$available_datasets <- self$setAvailableDatasets()

    gematrices_success <- self$GeneExpressionMatrices()
  }
)
