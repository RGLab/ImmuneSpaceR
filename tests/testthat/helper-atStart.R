# Load dependencies ------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rlabkey))
suppressPackageStartupMessages(library(ImmuneSpaceR))

# Declare global test variables ------------------------------------------------

# Move netrc file out the way to test against login and password in .Renviron
if (!any(file.exists("~/.netrc", "~/_netrc"))) {
  assign("labkey.netrc.file", ImmuneSpaceR:::.get_env_netrc(), .GlobalEnv)
}

# Regardless of netrc presence, test against `ISR_machine` in .Renviron
assign("labkey.url.base", ImmuneSpaceR:::.get_env_url(), .GlobalEnv)

# Initialize connections for studies used throughout testing
ALL <- suppressMessages(CreateConnection(""))
IS1 <- suppressMessages(CreateConnection("IS1"))
SDY28 <- suppressMessages(CreateConnection("SDY28"))
SDY67 <- suppressMessages(CreateConnection("SDY67"))
SDY87 <- suppressMessages(CreateConnection("SDY87"))
SDY180 <- suppressMessages(CreateConnection("SDY180"))
SDY269 <- suppressMessages(CreateConnection("SDY269"))
SDY404 <- suppressMessages(CreateConnection("SDY404"))
