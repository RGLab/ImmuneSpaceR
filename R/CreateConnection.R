################################################################################
# Main class and instantiation
################################################################################

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

#' @title The ImmuneSpaceConnection class
#' @name ImmuneSpaceConnection
#' @aliases
#' getGEMatrix
#' getDataset
#' listDatasets
#' getGEAnalysis
#' listGEAnalysis
#' addTreatmentt
#' EMNames
#' quick_plot
#'
#' @description
#' A connection respresents a study or a set of studies available on ImmuneSpace.
#' It provides function to download and display the data within these studies.
#'
#' @details
#' Uses global variables \code{labkey.url.base}, and \code{labkey.url.path}, to
#' access a study. \code{labkey.url.base} should be
#' \code{https://www.immunespace.org/}. \code{labkey.url.path} should be
#' \code{/Studies/studyname}, where 'studyname' is the accession number of the
#' study. The ImmuneSpaceConnection will initialize itself, and look for a
#' \code{.netrc} file in \code{"~/"} the user's home directory. The
#' \code{.netrc} file should contain a \code{machine}, \code{login}, and
#' \code{password} entry to allow access to ImmuneSpace, where \code{machine} is
#' the host name like "www.immunespace.org".
#'
#' @return An instance of an ImmuneSpaceConnection for a study in
#' \code{labkey.url.path}.
#'
#' @section Constructor:
#' \code{\link{CreateConnection}}
#'
#' @section Fields:
#' \describe{
#'   \item{\code{study}}{
#'     A \code{character}. The study accession number. Use an empty string ("")
#'     to create a connection at the project level.
#'   }
#'   \item{\code{config}}{
#'     A \code{list}. Stores configuration of the connection object such as
#'     URL, path and username.
#'   }
#'   \item{\code{available_datasets}}{
#'     A \code{data.table}. The table of datasets available in the connection
#'     object.
#'   }
#'   \item{\code{data_cache}}{
#'     A \code{list}. Stores the data to avoid downloading the same tables
#'     multiple times.
#'   }
#'   \item{\code{constants}}{
#'     A \code{list}. Used to store information regarding gene-expression data.
#'   }
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{initialize(..., config = NULL)}}{
#'     Initialize \code{ImmuneSpaceConnection} class.
#'     See \code{\link{CreateConnection}}.
#'   }
#'   \item{\code{print()}}{
#'     Print \code{ImmuneSpaceConnection} class.
#'   }
#'   \item{\code{getDataset(x, original_view = FALSE, reload = FALSE,
#'   colFilter = NULL, ...)}}{
#'     Get a dataset form the connection
#'
#'     \code{original_view}: A logical. If set tot TRUE, download the ImmPort
#'     view. Else, download the default grid view.
#'
#'     \code{reload}: A logical. Clear the cache. If set to TRUE, download the
#'     dataset, whether a cached version exist or not.
#'
#'     \code{colFilter}: A character. A filter as returned by Rlabkey's
#'     \code{makeFilter} function.
#'
#'     \code{...}: Extra arguments to be passed to \code{labkey.selectRows}.
#'   }
#'   \item{\code{addTreatment(matrixName = NULL)}}{
#'     Adds treatment information to the phenoData of an expression matrix
#'     available in the connection object.
#'
#'     \code{x}: A character. The name of a expression matrix that has been
#'      downloaded from the connection.
#'   }
#'   \item{\code{getGEMatrix(matrixName = NULL, cohort = NULL,
#'   outputType = "summary", annotation = "latest", reload = FALSE)}}{
#'     Downloads a normalized gene expression matrix from ImmuneSpace.
#'
#'     \code{x}: A character. The name of the gene expression matrix to download.
#'
#'     \code{cohort}: A character. The name of a cohort that has an associated
#'     gene expression matrix. Note that if `cohort` isn't NULL, then `x` is
#'     ignored.
#'
#'     \code{outputType}: one of 'raw', 'normalized' or 'summary'. If 'raw'
#'     then returns an expression matrix of non-normalized values by probe.
#'     'normalized' returns normalized values by probe.  'summary' returns
#'     normalized values averaged by gene symbol.
#'
#'     \code{annotation}: one of 'default', 'latest', or 'ImmSig'. Determines
#'     which feature annotation set is used.  'default' uses the fas from when
#'     the matrix was generated. latest' uses a recently updated fas based on
#'     the original.'ImmSig' is specific to studies involved in the
#'     ImmuneSignatures project and uses the annotation from when the
#'     meta-study's manuscript was created.
#'
#'     \code{reload}: A logical. If set to TRUE, the matrix will be downloaded
#'     again, even if a cached cop exist in the ImmuneSpaceConnection object.
#'   }
#'   \item{\code{EMNames(EM = NULL, colType = "participant_id")}}{
#'     Change the sampleNames of an ExpressionSet fetched by \code{getGEMatrix}
#'     using the information in the phenodData slot.
#'
#'     \code{x}: An ExpressionSet, as returned by \code{getGEMatrix}.
#'
#'     \code{colType}: A character. The type of column names. Valid options are
#'     'expsample_accession' and 'participant_id'.
#'   }
#'   \item{\code{listDatasets(output = c("datasets", "expression"))}}{
#'     Lists the datasets available in the study or studies of the connection.
#'   }
#'   \item{\code{listGEAnalysis()}}{
#'     Lists available gene expression analysis for the connection.
#'   }
#'   \item{\code{getGEAnalysis(...)}}{
#'     Downloads data from the gene expression analysis results table.
#'
#'     \code{...}: A list of arguments to be passed to \code{labkey.selectRows}.
#'   }
#'   \item{\code{clear_cache()}}{
#'     Clears the data_cache. Removes downloaded datasets and expression
#'     matrices.
#'   }
#'   \item{\code{getGEFiles(files, destdir = ".", quiet = FALSE)}}{
#'     Download gene expression raw data files.
#'
#'     \code{files}: A character. Filenames as shown on the
#'     gene_expression_files dataset.
#'
#'     \code{destdir}: A character. The local path to store the downloaded
#'     files.
#'   }
#'   \item{\code{listParticipantGroups()}}{
#'     Returns a dataframe with all saved Participant Groups on ImmuneSpace.
#'   }
#'   \item{\code{getParticipantData(group, dataType, original_view = FALSE,
#'   ...)}}{
#'     Returns a dataframe with ImmuneSpace data subset by groupId.
#'
#'     \code{group}: Use con$listParticipantGroups() to find Participant groupId
#'     or groupName.
#'
#'     \code{dataType}: Use \code{con$listDatasets('datasets')} to see possible
#'     dataType inputs.
#'   }
#'   \item{\code{quick_plot(...)}}{
#'     "Plots a selected dataset. This is the function used by the DataExplorer
#'     module on ImmuneSpace.
#'
#'     \code{dataset}: A character. The name of the dataset to plot, as
#'     displayed by the listDataset method.
#'
#'     \code{normalize_to_baseline}: A logical. If set to TRUE, the values are
#'     plotted as log2 fold-change from baseline.
#'
#'     \code{type}: A character. The type of plot. Valid choices are 'auto',
#'     'heatmap', 'boxplot', 'lineplot', 'violinplot'. If set to 'auto', the
#'     function will select an appropriate plot type for the selected data.
#'
#'     \code{filter}: A filter as created by the makeFilter function from
#'     Rlabkey.
#'
#'     \code{facet}: The facetting for ggplot2 based plots. Valid choices are
#'     'grid' and 'wrap'.
#'
#'     \code{text_size}: The size of all text elements in the plot.
#'
#'     \code{legend}: A character. Columns of the dataset or demographics to be
#'     added as legend on the heatmap. This argument is ignored if the plot type
#'     isn't heatmap.
#'
#'     \code{show_virus_strain}: A logical. Should all the virus strains be
#'     shown or should the values be averaged. Only used when dataset = 'hai'.
#'
#'     \code{interactive}: A logical. If set to TRUE, an interactive plot will
#'     be created. The default is FALSE.
#'
#'     \code{...}: Extra argument to be passed to ggplot. e.g: shape = 'Age',
#'     color = 'Race'.
#'   }
#' }
#'
#' @seealso \code{\link{CreateConnection}} \code{\link{ImmuneSpaceR-package}}
#'
#' @examples
#' \dontrun{
#' # Create a connection (Initiate a ImmuneSpaceConnection object)
#' sdy269 <- CreateConnection("SDY269")
#' sdy269
#' }
#'
#' @docType class
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
  which = "private",
  name = ".checkStudy",
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
  which = "private",
  name = ".setAvailableDatasets",
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
        setnames(ge, private$.munge(colnames(ge)))
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

    private$.checkStudy(self$config$verbose)

    self$available_datasets <- private$.setAvailableDatasets()

    gematrices_success <- self$GeneExpressionMatrices()
  }
)
