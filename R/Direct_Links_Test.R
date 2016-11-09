# Evan Henrich
# ImmuneSpace Test Script
# Last Updated: 11/8/2016

# Purpose: Ensure all links work correctly

# Pull data from labkey
get_data <- function(filetype, filename_col,link_text){
  df <- labkey.selectRows(baseUrl = "https://www.immunespace.org",
                         "/Studies",
                         schemaName = "study",
                         queryName = filetype,
                         colNameOpt = "caption")
  
  #iterate through the df to check link
  for(i in df){
  
  # pull out sdy ID as by splitting participant ID
  SDYid <- i$participant_id.split(".")[1]
  filename <- i$file_column
  
  # build link
  link <- paste0("https://www.immunespace.org/_webdav/Studies/SDY",SDYid,"/@files/rawdata/",link_text,"/",filename)
  
  #test_link
  url.exists(link)
  
  # if the link fails, keep track of it in error log
  }
}
  

get_data("fcs_sample_files","File Info Name","fcs")
get_data("gene_expression_files","File Info Name","gene_expression")

#TODO: figure out protocols ... different method of data extraction?

