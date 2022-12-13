# PUBLIC -----------------------------------------------------------------------

#' @title The ImmuneSpaceConnection class
#' @name ImmuneSpaceConnection
#' @aliases
#' ImmuneSpaceConnection
#' listDatasets
#' listGEMatrices
#' listGEAnalysis
#' listParticipantGroups
#' listParticipantMatrices
#' listWorkspaces
#' listGatingSets
#' summarizeCyto
#' summarizeGatingSet
#' loadGatingSet
#' getDataset
#' getGEMatrix
#' getGEAnalysis
#' getGEFiles
#' getGEInputs
#' getParticipantData
#' getParticipantGEMatrix
#' addTreatmentt
#' mapSampleNames
#'
#' @description
#' A connection respresents a study or a set of studies available on ImmuneSpace.
#' It provides function to download and display the data within these studies.
#'
#' @details
#' The ImmuneSpaceConnection will initialize itself, and look for a
#' \code{.netrc} file in \code{"~/"} the user's home directory. The
#' \code{.netrc} file should contain a \code{machine}, \code{login}, and
#' \code{password} entry to allow access to ImmuneSpace, where \code{machine} is
#' the host name like "datatools.immunespace.org".
#'
#' It can also use global variables \code{labkey.url.base}, and
#' \code{labkey.url.path}, to access a study. \code{labkey.url.base} should be
#' \code{https://datatools.immunespace.org/}. \code{labkey.url.path} should be
#' \code{/Studies/studyname}, where 'studyname' is the accession number of the
#' study.
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
#'   \item{\code{availableDatasets}}{
#'     A \code{data.table}. The table of datasets available in the connection
#'     object.
#'   }
#'   \item{\code{cache}}{
#'     A \code{list}. Stores the data to avoid downloading the same tables
#'     multiple times.
#'   }
#'   \item{\code{config}}{
#'     A \code{list}. Stores configuration of the connection object such as
#'     URL, path and username.
#'   }
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{initialize()}}{
#'     Initialize \code{ImmuneSpaceConnection} class.
#'     See \code{\link{CreateConnection}}.
#'   }
#'   \item{\code{print()}}{
#'     Print \code{ImmuneSpaceConnection} class.
#'   }
#'   \item{\code{listDatasets(output = c("datasets", "expression"))}}{
#'     Lists the datasets available in the study or studies of the connection.
#'   }
#'   \item{\code{listGEMatrices(verbose = FALSE, reload = FALSE, participantIds = NULL)}}{
#'     Lists available gene expression matrices for the connection.
#'
#'     \code{verbose}: A logical. If TRUE, whether to print the extra details
#'     for troubleshooting.
#'
#'     \code{reload}: A logical. If TRUE, retrieve the table of available gene
#'     expression matrices whether a cached version exist or not.
#'
#'     \code{participantIds}: A character vector of participant ids to filter
#'     by. Only matrices with data from \code{participantIds} will be returned.
#'     If \code{NULL}, all matrices are returned.
#'   }
#'   \item{\code{listGEAnalysis()}}{
#'     Lists available gene expression analysis for the connection.
#'   }
#'   \item{\code{listParticipantGroups()}}{
#'     Lists available participant groups on the ImmuneSpace portal.
#'   }
#'   \item{\code{listParticipantGEMatrices(group, verbose = FALSE)}}{
#'     Lists available gene expression matrices for participants in \code{group}.
#'
#'     \code{group}: A character or integer.
#'     Call \code{con$listParticipantGroups()} to see available participants
#'     groups. Use group_id or group_name as input.
#'
#'     \code{verbose}: A logical. If TRUE, whether to print the extra details
#'     for troubleshooting.
#'
#'   }
#'   \item{\code{listWorkspaces(reload = FALSE)}}{
#'     Lists available workspaces for the connection.
#'
#'     \code{reload}: A logical. If TRUE, download the table whether a cached
#'     version exist or not.
#'   }
#'   \item{\code{listGatingSets(reload = FALSE)}}{
#'     Lists available gating sets for the connection.
#'
#'     \code{reload}: A logical. If TRUE, download the table whether a cached
#'     version exist or not.
#'   }
#'   \item{\code{summarizeCyto()}}{
#'     Prints a summary of cytometry data for the connection.
#'   }
#'   \item{\code{summarizeGatingSet(gatingSet)}}{
#'     Prints a summary of a gating set. Note that this method currently works
#'     only in the ImmuneSpace RStudio session.
#'
#'     \code{gatingSet}: A character. The name of the gating set to summarize.
#'   }
#'   \item{\code{loadGatingSet(gatingSet)}}{
#'     Loads a gating set via \code{flowWorkspace::load_gs} to the
#'     current environment. Note that this method currently works only in the
#'     ImmuneSpace RStudio Docker session.
#'
#'     \code{gatingSet}: A character. The name of the gating set to load.
#'   }
#'   \item{\code{getDataset(x, original_view = FALSE, reload = FALSE,
#'   colFilter = NULL, ...)}}{
#'     Get a dataset form the connection.
#'
#'     \code{x}: A character. The name of the dataset to download.
#'
#'     \code{original_view}: A logical. If TRUE, download the original ImmPort
#'     view; else, download the default grid view.
#'
#'     \code{reload}: A logical. If TRUE, download the dataset whether a cached
#'     version exist or not.
#'
#'     \code{colFilter}: A character. A filter as returned by Rlabkey's
#'     \code{makeFilter} function.
#'
#'     \code{...}: Extra arguments to be passed to \code{labkey.selectRows}.
#'   }
#'   \item{\code{getGEMatrix(matrixName = NULL, cohortType = NULL,
#'   outputType = "summary", annotation = "latest", reload = FALSE,
#'   verbose = FALSE)}}{
#'     Downloads a probe-level or gene-symbol summarized expression matrix from
#'     ImmuneSpace and constructs an ExpressionSet. Use \code{experimentData()}
#'     on the resulting ExpressionSet object to see version info for annotation.
#'
#'     \code{matrixName}: A character. The name of the gene expression matrix
#'     to download.
#'
#'     \code{cohortType}: A character. The name of a cohortType that has an
#'     associated gene expression matrix. Note that if this argument is not
#'     NULL, then \code{matrixName} is ignored. CohortType is a concatenation of
#'     "cohort" and "cell type" that allows the user to specify a matrix for the
#'     cell type subset of a cohort.
#'
#'     \code{outputType}: A character. one of 'raw', 'normalized' or 'summary'.
#'     If 'raw', returns an expression matrix of non-normalized values by probe.
#'     'normalized' returns normalized values by probe. 'summary' returns
#'     normalized values averaged by gene symbol.
#'
#'     \code{annotation}: A character. one of 'default', 'latest', or 'ImmSig'.
#'     Determines which feature annotation set (FAS) is used. 'default' uses the
#'     FAS from when the matrix was generated. latest' uses a recently updated
#'     FAS based on the original. 'ImmSig' is specific to studies involved in
#'     the ImmuneSignatures project and uses the annotation from when the
#'     meta-study's manuscript was created.
#'
#'     \code{reload}: A logical. If set to TRUE, the matrix will be downloaded
#'     again, even if a cached copy exists in the ImmuneSpaceConnection object.
#'
#'     \code{verbose}: A logical. If set to TRUE, notes on how the expressionSet
#'     object was created will be printed, including normalization, summarization,
#'     feature_annotation_set, and alias2symbol mapping version of org.Hs.eg.db.
#'   }
#'   \item{\code{getGEAnalysis(...)}}{
#'     Downloads data from the gene expression analysis results table.
#'
#'     \code{...}: A list of arguments to be passed to \code{labkey.selectRows}.
#'   }
#'   \item{\code{getGEInputs()}}{
#'     Downloads data from the gene expression input samples table.
#'   }
#'   \item{\code{getParticipantData(group, dataType, original_view = FALSE,
#'   ...)}}{
#'     Returns a data.table with data subset by participant group.
#'
#'     \code{group}: A character or integer.
#'     Call \code{con$listParticipantGroups()} to see available participants
#'     groups. Use group_id or group_name as input.
#'
#'     \code{dataType}: A character. Use \code{con$availableDatasets} to see
#'     available dataset names.
#'   }
#'   \item{\code{getParticipantGEMatrix(group, outputType = "summary",
#'   annotation = "latest", reload = FALSE)}}{
#'     Downloads probe-level or gene-symbol summarized expression matrices for all
#'     participants within \code{group} from ImmuneSpace and constructs an
#'     ExpressionSet containing observations for each participant in \code{group}
#'     where gene expression data is available.
#'
#'     \code{group}: A character or integer.
#'     Call \code{con$listParticipantGroups()} to see available participants
#'     groups. Use group_id or group_name as input.
#'
#'     \code{outputType}: A character. one of 'raw', 'normalized' or 'summary'.
#'     If 'raw', returns an expression matrix of non-normalized values by probe.
#'     'normalized' returns normalized values by probe. 'summary' returns
#'     normalized values averaged by gene symbol.
#'
#'     \code{annotation}: A character. one of 'default', 'latest', or 'ImmSig'.
#'     Determines which feature annotation set (FAS) is used. 'default' uses the
#'     FAS from when the matrix was generated. latest' uses a recently updated
#'     FAS based on the original. 'ImmSig' is specific to studies involved in
#'     the ImmuneSignatures project and uses the annotation from when the
#'     meta-study's manuscript was created.
#'
#'     \code{reload}: A logical. If set to TRUE, matrices will be downloaded
#'     again, even if a cached copy exists in the ImmuneSpaceConnection object.
#'
#'   }
#'   \item{\code{downloadGEFiles(files, destdir = ".")}}{
#'     Downloads gene expression raw data files.
#'
#'     \code{files}: A character. Filenames as shown on the
#'     gene_expression_files dataset.
#'
#'     \code{destdir}: A character. The local path to store the downloaded
#'     files.
#'   }
#'   \item{\code{addTreatment(expressionSet)}}{
#'     Adds treatment information to the phenoData of an ExpressionSet.
#'
#'     \code{expressionSet}: An ExpressionSet. The ExpressionSet object that has
#'     been downloaded from the connection.
#'   }
#'   \item{\code{mapSampleNames(EM = NULL, colType = "participant_id")}}{
#'     Changes the sampleNames of an ExpressionSet fetched by \code{getGEMatrix}
#'     using the information in the phenodData slot.
#'
#'     \code{EM}: An ExpressionSet, as returned by \code{getGEMatrix}.
#'
#'     \code{colType}: A character. The type of column names. Valid options are
#'     'expsample_accession' and 'participant_id'.
#'   }
#'   \item{\code{plot(...)}}{
#'     Visualizes a selected dataset. This method is used by the DataExplorer
#'     module on the ImmuneSpace portal.
#'
#'     \code{dataset}: A character. The name of the dataset to plot, as
#'     displayed by the listDataset method.
#'
#'     \code{normalize_to_baseline}: A logical. If TRUE, the values are plotted
#'     as log2 fold-change from baseline.
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
#'     \code{interactive}: A logical. If TRUE, an interactive plot will be
#'     created. The default is FALSE.
#'
#'     \code{...}: Extra argument to be passed to ggplot. e.g: shape = 'Age',
#'     color = 'Race'.
#'   }
#'   \item{\code{clearCache()}}{
#'     Clears the cache. Removes downloaded datasets and expression matrices.
#'   }
#' }
#'
#' @seealso \code{\link{CreateConnection}} \code{\link{ImmuneSpaceR-package}}
#'
#' @examples
#' \dontrun{
#' # Create a connection (Initiate a ImmuneSpaceConnection object)
#' sdy269 <- CreateConnection("SDY269")
#'
#' # Print the connection object
#' sdy269
#'
#' # Retrieve the HAI dataset
#' HAI <- sdy269$getDataset("hai")
#'
#' # Fetch a summarized gene expresssion matrix with latest annotation
#' LAIV <- sdy269$getGEMatrix("LAIV_2008")
#'
#' # Visualize the ELISA dataset
#' sdy269$plot("elisa")
#' }
#'
#' sdy <- try(CreateConnection("SDY269"))
#' if (inherits(sdy, "try-error")) {
#'   warning("Read the Introduction vignette for more information on how to set
#' up a .netrc file.")
#' }
#' @docType class
ISCon <- R6Class(
  classname = "ImmuneSpaceConnection",
  public = list(
    study = character(),
    config = list(),
    availableDatasets = data.table(),
    cache = list()
  ),
  private = list(
    .constants = list()
  )
)


# Initialize the ISCon object
ISCon$set(
  which = "public",
  name = "initialize",
  value = function(study = NULL,
                   login = NULL,
                   password = NULL,
                   verbose = FALSE,
                   onTest = FALSE) {
    if (length(study) > 1) {
      stop("For multiple studies, use an empty string and filter the connection.")
    }
    .check_internet()
    .check_portal(onTest)

    # fetch config variables
    labkey.url.base <- .get_url_base(onTest, verbose)
    labkey.user.email <- .get_user_email()
    labkey.url.path <- .get_url_path(study)
    curlOptions <- .set_curl_options(login, password)

    # set fields
    self$config <- list(
      labkey.url.base = labkey.url.base,
      labkey.url.path = labkey.url.path,
      labkey.user.email = labkey.user.email,
      curlOptions = curlOptions,
      verbose = verbose
    )
    self$study <- basename(self$config$labkey.url.path)
    private$.constants <- list(
      matrices = "GE_matrices",
      matrix_inputs = "GE_inputs"
    )

    # validate connection
    .check_credential(labkey.url.base, verbose)
    private$.checkStudy(self$config$verbose)

    # fetch available datasets and expression matrices
    private$.setAvailableDatasets()
    self$listGEMatrices()

    self
  }
)



# PRIVATE ----------------------------------------------------------------------



# HELPER -----------------------------------------------------------------------

# check internet connection
.check_internet <- function() {
  if (!has_internet()) {
    stop("No internet connection. Please connect to internet and try again.")
  }
}


.get_host <- function(onTest = FALSE) {
  ifelse(onTest, "datatools-dev.immunespace.org", "datatools.immunespace.org")
}


# check if the portal is up
.check_portal <- function(onTest = FALSE) {
  host <- .get_host(onTest)
  if (is.null(nslookup(host, error = FALSE))) {
    stop("The portal is currently down. Try again later.")
  }
}


# Try to parse labkey options from global environment
# which really should have been done through option()/getOption() mechanism
# Here we do this to be compatible to labkey online report system
# that automatically assigns these variables in global environment
.get_url_base <- function(onTest = FALSE, verbose = FALSE) {
  labkey.url.base <- try(get("labkey.url.base", .GlobalEnv), silent = TRUE)
  if (inherits(labkey.url.base, "try-error")) {
    labkey.url.base <- paste0("https://", .get_host(onTest))
  }

  if (!is.null(getOption("labkey.baseUrl"))) {
    labkey.url.base <- getOption("labkey.baseUrl")
  }

  # Allow http (no SSL) for local
  if (grepl("(immunespace|lji)", labkey.url.base)) {
    labkey.url.base <- gsub("http:", "https:", labkey.url.base)
    if (length(grep("^https://", labkey.url.base)) == 0) {
      labkey.url.base <- paste0("https://", labkey.url.base)
    }
  } else {
    if (!grepl(":8080/*$", labkey.url.base)) {
      labkey.url.base <- paste0(labkey.url.base, ":8080")
    }
  }

  if (isTRUE(verbose)) message(sprintf("URL base: %s", labkey.url.base))

  labkey.url.base
}


.get_user_email <- function() {
  labkey.user.email <- try(get("labkey.user.email", .GlobalEnv), silent = TRUE)
  if (inherits(labkey.user.email, "try-error")) {
    labkey.user.email <- "unknown_user at not_a_domain.com"
  }

  if (!is.null(getOption("labkey.user.email"))) {
    labkey.user.email <- getOption("labkey.user.email")
  }

  labkey.user.email
}


.get_url_path <- function(study) {
  if (is.null(study)) {
    if (exists("labkey.url.path", .GlobalEnv)) {
      labkey.url.path <- get("labkey.url.path", .GlobalEnv)
    } else {
      stop("study cannot be NULL")
    }
  } else {
    base_path <- ifelse(
      grepl("^IS\\d{1,3}$", study),
      "/HIPC/",
      "/Studies/"
    )
    labkey.url.path <- paste0(base_path, study)
  }

  labkey.url.path
}


# set curoption for Rlabkey package
# Rlabkey stores the Curl options in its package environment through
# labkey.setCurlOptions call, so in theory we need to reset it prior to each
# Rlabkey query because  multiple connections created by user individually may
# have different urls and ssl settings. Ideally labkey.selectRows should
# optionally parse the options from its argument besides package environment.
# For now, we assume they all share the same setting and init it only once here
.set_curl_options <- function(login = NULL, password = NULL) {
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
      netrc_file = nf,
      useragent = useragent
    )
  } else {
    curlOptions <- labkey.setCurlOptions(
      useragent = useragent
    )
  }

  curlOptions
}


# check credential
.check_credential <- function(labkey.url.base, verbose = FALSE) {
  if (verbose) message("Checking credential...")

  res <- GET(
    url = paste0(labkey.url.base, "/login-whoami.view"),
    config = Rlabkey:::labkey.getRequestOptions()
  )

  if (res$status_code == 200) {
    if (grepl("json", res$headers$`content-type`)) {
      parsed <- httr::content(res)

      if (parsed$displayName == "guest" & !grepl("8080", labkey.url.base)) {
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
}
