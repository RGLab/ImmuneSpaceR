# Evan Henrich
# ImmuneSpace Test Script
# Last Updated: 11/8/2016

# Purpose: Ensure all links for gene expression, fcs, and protocols
# work correctly for ImmuneSpace, both test / production

#--------Requirements--------------------------------------
require(Rlabkey)
require(plyr)
require(grid)
require(gridExtra)
require(httr)

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
}

# parse SDYids into a list so they can be merged easily later
parse_ids <- function(raw_data_file){
  Sdyids <- apply(MARGIN=1, raw_data_file, FUN = function(x) get_sdyid(x))
  return(Sdyids)
}

# Test a link and pass errors / warnings (b/c timeout is a problem for some reason)
link_test <- function(study_id,filename,link_text){
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

# Parse link test results into list for merging
parse_results <- function(Link_Test_Results){
  final_Res_List = list()
  for(i in Link_Test_Results){
    final_Res_List[i] <- Link__Test_Results[[i]]
  }
  return(final_Res_List)
} 

# does analysis and outputs table of just bad links with SDYID / filename
analyzer <- function(info_set){
  files <- info_set[1]
  filename_col <- info_set[2]
  link_text <- info_set[3]
  
  raw_data <- get_data(files)
  Sdyids <- parse_ids(raw_data)
  Filenames <- raw_data[ , filename_col]
  Link_Test_Results <- mapply(link_test, Filenames, Sdyids, link_text)
  final_Res_List <- parse_results(Link_Test_Results)
  final_df <- cbind(Sdyids,Filenames,final_Res_List)
  
  num_rows <- nrow(final_df)
  bad_links_table <- final_df[which (final_df$Results != 200), ]
  
  return(bad_links_table)
}

# Gathers performance info and packages with bad links table for pdf output
bad_links_to_pdf <- function(name, info_set){
  start <- strftime(Sys.time(),format = "%T")
  bad_links <- analyzer(info_set)
  end <- strftime(Sys.time(),format = "%T")
  
  ts <-paste(format(Sys.time(), "%Y_%m_%d %T"), "pdf", sep = ".")
  pdfname <- paste0(name," ", ts)
  pdf(pdfname)
  grid.table(bad_links)
  dev.off()
}

#-------Execution------------------------------------------
fcs <- c("fcs_sample_files","File Info Name","fcs")
ge <- c("gene_expression_files","File Info name","gene_expression")

bad_links_to_pdf("Gene Expression",ge)
bad_links_to_pdf("FCS",fcs)

#TODO: figure out protocols ... different method of data extraction?


