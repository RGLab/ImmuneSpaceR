# Evan Henrich
# ImmuneSpace Test Script
# Last Updated: 11/8/2016

# Purpose: Ensure all links for gene expression, fcs, and protocols
# work correctly for ImmuneSpace, both test / production

#--------Requirements--------------------------------------
require(Rlabkey)
require(plyr)
require(httr)
require(parallel)

#---------Helper-Functions-----------------------------------
# Pull data from labkey
get_data <- function(filetype){
  df <- labkey.selectRows(baseUrl = "https://www.immunespace.org",
                          "/Studies",
                          schemaName = "study",
                          queryName = filetype,
                          colNameOpt = "caption")
  df <- df[!is.na(df$`File Info Name`),]
  df <- df[!duplicated(df$`File Info Name`),]
  
  numrow <- nrow(df)
  blurb <- paste0("There are ", numrow, " files to check in ", filetype)
  print(blurb)
  
  return(df)
}

# splitting function for SDYID
get_sdyid <- function(subjectid){
  id_table <- read.table(text = subjectid, sep = ".", as.is = T, fill = T)
  sdyid <- id_table[1,2]
  return(sdyid)
}

# parse SDYids into a list so they can be merged easily later
parse_ids <- function(raw_data_file){
  Sdyids <- apply(MARGIN=1, raw_data_file, FUN = function(x) get_sdyid(x))
  return(Sdyids)
}

# Test a link with GET b/c url.exists throws FALSE if certain options aren't specified
link_test <- function(filename,study_id,link_text){
 
    filename <- URLencode(filename)
    link <- paste0("https://www.immunespace.org/_webdav/Studies/SDY",study_id,"/@files/rawdata/",link_text,"/",filename) 
    info <- GET(link)
    status <- info$status_code
    return(status)
}

# does analysis and outputs table of just bad links with study id, filename, and error code from GET (impt. for knowing if actual bad link)
analyzer <- function(files,link_text){
  
  raw_data <- get_data(files)
  sdyids <- parse_ids(raw_data)
  filenames <- raw_data[ , "File Info Name"]
  link_test_results <- mcmapply(FUN=link_test, filename=filenames, study_id=sdyids, link_text=link_text,mc.cores = detectCores())
  link_status <- unname(link_test_results)
  final_df <- cbind(sdyids,filenames,link_status)
  bad_links_table <- subset(final_df, link_status != "200")
  
  numrow <- nrow(bad_links_table)
  blurb <- paste0("There are ", numrow, " bad links in ", files)
  print(blurb)
  
  return(bad_links_table)
}

#-------Execution------------------------------------------

result <- list()

result$fcs <- analyzer(files = "fcs_sample_files",link_text = "flow_cytometry")
result$ge <- analyzer(files = "gene_expression_files", link_text = "gene_expression")

return(result)

#TODO: figure out protocols ... different method of data extraction?


