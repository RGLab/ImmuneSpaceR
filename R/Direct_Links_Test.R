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
require(R2HTML)

#---------Helper-Functions-----------------------------------
# Pull data from labkey
get_data <- function(filetype){
  df <- labkey.selectRows(baseUrl = "https://www.immunespace.org",
                          "/Studies",
                          schemaName = "study",
                          queryName = filetype,
                          colNameOpt = "caption")
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

# Test a link and pass errors / warnings (b/c timeout is a problem for some reason)
link_test <- function(filename,study_id,link_text){
  result <- "init"
  
  if(is.na(filename)){
    result <- "NA"
  }else if(filename == "<NA>"){
    result <- "NA"
  }else{
    link <- paste0("https://www.immunespace.org/_webdav/Studies/SDY",study_id,"/@files/rawdata/",link_text,"/",filename) 
    result <- tryCatch({
      info <- GET(link)
      status <- info$status_code
      return(status)
    },
    error = function(e){
      return(e)
    },
    warning = function(w){
      return(w)
    })
  }
  return(result)
}

# does analysis and outputs table of just bad links with SDYID / filename
analyzer <- function(info_set){
  files <- info_set[1]
  filename_col <- info_set[2]
  link_text <- info_set[3]
  
  raw_data <- get_data(files)
  Sdyids <- parse_ids(raw_data)
  Filenames <- raw_data[ , filename_col]
  link_test_results <- mcmapply(FUN=link_test, filename=Filenames, study_id=Sdyids, link_text=link_text,mc.cores = detectCores())
  Link_Status <- unname(link_test_results)
  final_df <- cbind(Sdyids,Filenames,Link_Status)
  bad_links_table <- subset(final_df, Link_Status != "200")
  
  return(bad_links_table)
}

# Gathers performance info and packages with bad links table for pdf output
bad_links_to_pdf <- function(name, info_set){
  start <- strftime(Sys.time(),format = "%T")
  bad_links <- analyzer(info_set)
  end <- strftime(Sys.time(),format = "%T")
  
  start_string <- paste0(name," run started at: ", start)
  end_string <- paste0(name," run ended at: ", end)
  
  #if badlinks are found ....
  if(nrow(bad_links) > 1){
    ts <-paste(format(Sys.time(), "%Y_%m_%d %T"), "html", sep = ".")
    htmlname <- paste0(name," ", ts)
    HTML(bad_links, file=htmlname)
    print(htmlname)
  }
  
  print("Run information")
  print(start_string)
  print(end_string)
  
}

#-------Execution------------------------------------------
fcs <- c("fcs_sample_files","File Info Name","fcs")
ge <- c("gene_expression_files","File Info Name","gene_expression")

bad_links_to_pdf("Gene Expression",ge)
bad_links_to_pdf("FCS",fcs)

#TODO: figure out protocols ... different method of data extraction?


