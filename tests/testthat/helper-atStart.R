# Load dependencies ------------------------------------------------------------
suppressPackageStartupMessages(library(ImmuneSpaceR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rlabkey))


# Declare global test variables ------------------------------------------------

# Use variables put in global environment through .Renviron config file
if (!any(file.exists("~/.netrc", "~/_netrc"))) {
  assign("labkey.netrc.file", ImmuneSpaceR:::.get_env_netrc(), .GlobalEnv)
  assign("labkey.url.base", ImmuneSpaceR:::.get_env_url(), .GlobalEnv)
}

# Initialize connections for studies used throughout testing
ALL    <- suppressMessages(CreateConnection(""))
SDY269 <- suppressMessages(CreateConnection("SDY269"))
SDY67  <- suppressMessages(CreateConnection("SDY67"))
SDY180 <- suppressMessages(CreateConnection("SDY180"))
SDY28  <- suppressMessages(CreateConnection("SDY28"))
SDY87  <- suppressMessages(CreateConnection("SDY87"))
IS1    <- suppressMessages(CreateConnection("IS1"))
SDY404 <- suppressMessages(CreateConnection("SDY404"))

