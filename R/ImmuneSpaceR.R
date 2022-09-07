#' @title A Thin Wrapper Around ImmuneSpace
#'
#' @description ImmuneSpaceR provides a convenient API for accessing data sets
#' within the ImmuneSpace database.
#'
#' @details Uses the Rlabkey package to connect to ImmuneSpace. Implements
#' caching, and convenient methods for accessing data sets.
#'
#' @name ImmuneSpaceR-package
#' @aliases ImmuneSpaceR
#' @seealso \code{\link{CreateConnection}}
#' @import data.table Rlabkey Biobase
#' @importFrom R6 R6Class
NULL


# globalVariables to remove RCC NOTES due to data.table scoping
globalVariables(
  c(
    "biosample_accession",
    "study_time_collected",
    "study_time_collected_unit",
    "arm_name",
    "response",
    "analyte",
    "analyte_name",
    "value_preferred",
    "spot_number_reported",
    "cell_number_reported",
    "threshold_cycles",
    "entrez_gene_id",
    "mfi",
    "concentration_value",
    "population_cell_number",
    "population_name_reported",
    "x", "y", "err", # ggplot2 aes
    "ID", # qpHeatmap
    "stc", "stcu", "time_str", # standardize_time
    "virus_strain", "virus", "cohort",
    "participant_id",
    "PROBE_ID",
    "geo_accession",
    "expsample_accession",
    "feature_id",
    "labkey.netrc.file",
    "cell_number_preferred"
  )
)
