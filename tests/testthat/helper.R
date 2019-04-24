# Load depenendcies ------------------------------------------------------------
suppressPackageStartupMessages(library(ImmuneSpaceR))
suppressPackageStartupMessages(library(data.table))


# Declare global test variables ------------------------------------------------
if (!any(file.exists("~/.netrc", "~/_netrc"))) {
  assign("labkey.netrc.file", ImmuneSpaceR:::.get_env_netrc(), .GlobalEnv)
  assign("labkey.url.base", ImmuneSpaceR:::.get_env_url(), .GlobalEnv)
}

CONNECTIONS <- list(
  ALL = suppressWarnings(try(CreateConnection(""), silent = TRUE)),
  SDY28 = suppressWarnings(try(CreateConnection("SDY28"), silent = TRUE)),
  SDY180 = suppressWarnings(try(CreateConnection("SDY180"), silent = TRUE)),
  SDY269 = suppressWarnings(try(CreateConnection("SDY269"), silent = TRUE))
  # IS1 = try(CreateConnection("IS1"), silent = TRUE)
)
