# Load depenendcies ------------------------------------------------------------
suppressPackageStartupMessages(library(ImmuneSpaceR))
suppressPackageStartupMessages(library(data.table))


# Declare test-wide variables --------------------------------------------------
connections <- list(
  SDY269 = try(CreateConnection("SDY269", verbose = FALSE)),
  SDY180 = try(CreateConnection("SDY180", verbose = FALSE)),
  SDY28 = try(CreateConnection("SDY28", verbose = FALSE))
  # IS1 = try(CreateConnection("IS1", verbose = FALSE))
)
commonColumns <- c("age_reported", "gender", "race", "participant_id")
specificColumnsSet <- list(
  hai = data.table(
    name = c("value_preferred", "virus"),
    type = c("numeric", "character")
  ),
  neut_ab_titer = data.table(
    name = c("value_preferred", "virus"),
    type = c("numeric", "character")
  ),
  hla_typing = data.table(
    name = c("allele_1", "allele_2", "locus_name"),
    type = c("character", "character", "character")
  ),
  elisa = data.table(
    name = c("value_reported", "analyte"),
    type = c("numeric", "character")
  ),
  elispot = data.table(
    name = c("spot_number_reported", "analyte"),
    type = c("numeric", "character")
  ),
  pcr = data.table(
    name = c("value_reported", "entrez_gene_id"),
    type = c("numeric", "character")
  ),
  gene_expression_files = data.table(
    name = c("file_info_name", "cohort"),
    type = c("character", "character")
  ),
  mbaa = data.table(
    name = c("analyte_name", "concentration_value")
  ),
  fcs_analyzed_result = data.table(
    name = c("population_cell_number", "population_definition_reported"),
    type = c("character", "character")
  ),
  fcs_sample_files = data.table(
    name = c("file_info_name"),
    type = c("character")
    ),
  fcs_control_files = data.table(
    name = c("sample_file", "control_file"),
    type = c("character", "character")
  )
)


# Define helper test functions -------------------------------------------------
test_getDataset <- function(study, dataset) {
  test_that(paste(study, dataset), {
    con <- connections[[study]]
    specificColumns <- specificColumnsSet[[dataset]]

    data <- try(con$getDataset(dataset, reload = TRUE))

    # Dataset is of the right class and not empty
    expect_is(data, "data.table")
    expect_gt(nrow(data), 0)

    # All required columns are here
    expect_true(all(commonColumns %in% colnames(data)))
    expect_true(all(specificColumns$name %in% colnames(data)))

    # Important columns have no NAs
    expect_false(all(is.na(data[, specificColumns$name, with = FALSE])))

    # Important columns have the right class
    columns <- data[, specificColumns$name, with = FALSE]
    columnTypes <- vapply(columns, class, character(1), USE.NAMES = FALSE)
    expect_equal(columnTypes, specificColumns$type)
  })
}
