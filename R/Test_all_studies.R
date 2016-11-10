# Evan Henrich
# ImmuneSpace Test Script
# Last Updated: 11/8/2016

# Purpose: Ensure all studies have necessary files


require(Rlabkey)
require(ImmuneSpaceR)
require(parallel)

# Create Output file and switch printing to file from stdout
ts <-paste(format(Sys.time(), "%Y_%m_%d"), "pdf", sep = ".")
filename <- paste0("Test_IS_Links_", ts)
sink(filename)

# Hit LK DB to get list of available studies
labkey.url.base <- "https://www.immunespace.org"
#folders <- labkey.getFolders(baseUrl = labkey.url.base, folderPath = "/Studies")
#SDY_list <- folders[1]
#SDY_list <- SDY_list[SDY_list != "Studies"]
#SDY_list <- SDY_list[SDY_list != "SDY_template"]

# Start Time
starttime <- strftime(Sys.time(),format = "%T")
cat(starttime)
cat("\n")

# Open connection
con <- CreateConnection("")
tryCatch({
  con$.test_files(mc=TRUE)
  },
  error = function(e){
    message(e)
  },
  warning = function(w){
    message(w)
  })


# Loop through all studies using RSauteraud's script(See ln133>ln179 in ImmuneSpace.R)
#for(i in SDY_list){
#  cat("Working on Study: ", i)
#  cat("\n")
#  con <- CreateConnection(i)
#  tryCatch({
#    con$.test_files(mc=TRUE)
#  },
#  error = function(e){
#    message(e)
#  },
#  warning = function(w){
#    message(w)
#  })
#  cat("\n")
#  cat("\n")
#}

# End Time
endtime <- strftime(Sys.time(),format = "%T")
cat(endtime)
