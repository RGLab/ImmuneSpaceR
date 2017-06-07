# 
#  Global variables
#

# Environment ------------------------------------------------
TEST <- "https://test.immunespace.org"
PROD <- "https://www.immunespace.org"
labkey.url.base <- PROD


#------DATASET-COLUMN-NAMES------------------------------------
demo <- c("age_reported", "gender", "race")
idCols <- c("participant_id")
common_cols <- c(demo, idCols)


haiCols     <- data.frame(name = c("value_reported", "virus"),
                          type = c("numeric", "character"))
nabCols     <- data.frame(name = c("value_reported", "virus"),
                          type = c("numeric", "character"))
hlaCols     <- data.frame(name = c("allele_1", "allele_2", "locus_name"),
                          type = c("character", "character", "character"))
elisaCols   <- data.frame(name = c("value_reported", "analyte"),
                          type = c("numeric", "character"))
elispotCols <- data.frame(name = c("spot_number_reported", "analyte"),
                          type = c("numeric", "character"))
pcrCols     <- data.frame(name = c("value_reported", "entrez_gene_id"),
                          type = c("numeric", "character"))
gefCols     <- data.frame(name = c("file_info_name", "cohort"),
                          type = c("character", "character"))
mbaaCols    <- data.frame(name = c("analyte_name", "concentration_value"))
farCols     <- data.frame(name = c(#"population_cell_number", # declared as VARCHAR(500)
  "population_definition_reported"),
                          type = c(#"numeric",
                            "character"))
fcsCols     <- data.frame(name = c("file_info_name"),
                          type = c("character"))
fccCols     <- data.frame(name = c("sample_file", "control_file"),
                          type = c("character", "character"))

#------VALID-ESETS-----------------------------------------------
# these are studies with valid expression matrices
valid_esets <- c("SDY63", "SDY212", "SDY400", "SDY404")
