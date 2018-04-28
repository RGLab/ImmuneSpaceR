# PUBLIC -----------------------------------------------------------------------

#' @title The ImmuneSpaceConnection class
#' @name ImmuneSpaceConnection
#' @aliases
#' ImmuneSpaceConnection
#' getGEMatrix
#' getDataset
#' listDatasets
#' getGEAnalysis
#' listGEAnalysis
#' addTreatmentt
#' mapSampleNames
#' plot
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
#' the host name like "www.immunespace.org".
#'
#' It can also use global variables \code{labkey.url.base}, and
#' \code{labkey.url.path}, to access a study. \code{labkey.url.base} should be
#' \code{https://www.immunespace.org/}. \code{labkey.url.path} should be
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
#'   \item{\code{initialize(..., config = NULL)}}{
#'     Initialize \code{ImmuneSpaceConnection} class.
#'     See \code{\link{CreateConnection}}.
#'   }
#'   \item{\code{print()}}{
#'     Print \code{ImmuneSpaceConnection} class.
#'   }
#'   \item{\code{listDatasets(output = c("datasets", "expression"))}}{
#'     Lists the datasets available in the study or studies of the connection.
#'   }
#'   \item{\code{listGEMatrices()}}{
#'     Lists available gene expression matrices for the connection.
#'   }
#'   \item{\code{listGEAnalysis()}}{
#'     Lists available gene expression analysis for the connection.
#'   }
#'   \item{\code{listParticipantGroups()}}{
#'     Lists available participant groups on the ImmuneSpace portal.
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
#'   \item{\code{getGEMatrix(matrixName = NULL, cohort = NULL,
#'   outputType = "summary", annotation = "latest", reload = FALSE, verbose = FALSE)}}{
#'     Downloads a normalized gene expression matrix from ImmuneSpace.
#'
#'     \code{matrixName}: A character. The name of the gene expression matrix
#'     to download.
#'
#'     \code{cohort}: A character. The name of a cohort that has an associated
#'     gene expression matrix. Note that if this argument is not NULL, then
#'     \code{matrixName} is ignored.
#'
#'     \code{outputType}: one of 'raw', 'normalized' or 'summary'. If 'raw'
#'     then returns an expression matrix of non-normalized values by probe.
#'     'normalized' returns normalized values by probe. 'summary' returns
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
#'   \item{\code{getGEFiles(files, destdir = ".", quiet = FALSE)}}{
#'     Downloads gene expression raw data files.
#'
#'     \code{files}: A character. Filenames as shown on the
#'     gene_expression_files dataset.
#'
#'     \code{destdir}: A character. The local path to store the downloaded
#'     files.
#'   }
#'   \item{\code{getGEInputs()}}{
#'     Downloads data from the gene expression input samples table.
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
#'   \item{\code{addTreatment(matrixName = NULL)}}{
#'     Adds treatment information to the phenoData of an expression matrix
#'     available in the connection object.
#'
#'     \code{matrixName}: A character. The name of a expression matrix that has
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
#'     Clears the cache. Removes downloaded datasets and expression
#'     matrices.
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
  value = function(..., config = NULL) {
    # invoke the default init routine in case it needs to be invoked
    # (e.g. when using $new(object) to construct the new object based on the existing object)
    # callSuper(...) # THIS MIIGHT NOT APPLICABLE IN R6

    private$.constants <- list(
      matrices = "GE_matrices",
      matrix_inputs = "GE_inputs"
    )

    if (!is.null(config)) {
      self$config <- config
    }

    self$study <- basename(config$labkey.url.path)

    private$.checkStudy(self$config$verbose)

    self$availableDatasets <- private$.setAvailableDatasets()

    gematrices_success <- self$listGEMatrices()
  }
)



# PRIVATE ----------------------------------------------------------------------



# HELPER -----------------------------------------------------------------------
