# Evan Henrich
# ImmuneSpace Test Script
# Last Updated: 11/8/2016

# Purpose: Ensure all studies have necessary files
require(Rlabkey)
require(ImmuneSpaceR)

# Create Output file
ts <-paste(format(Sys.time(), "%Y_%m_%d"), "pdf", sep = ".")
filename <- paste0("Test_IS_Links_", ts)
sink(filename)

# Hit LK DB to get list of available studies
labkey.url.base = "https://www.immunespace.org"
folders <- labkey.getFolders(baseUrl = labkey.url.base, folderPath = "/Studies")
SDY_list <- folders[1]
SDY_list <- SDY_list[SDY_list != "Studies"]
SDY_list <- SDY_list[SDY_list != "SDY_template"]

# Loop through all studies using RSauteraud's script(See ln133>ln179 in ImmuneSpace.R)
for(i in SDY_list){
  cat("Working on Study: ", i)
  cat("\n")
  con <- CreateConnection(i)
  con$.test_files()
  cat("\n")
  cat("\n")
  
}
