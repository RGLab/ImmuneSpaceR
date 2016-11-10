# Evan Henrich
# ImmuneSpace Test Script
# Last Updated: 11/8/2016

# Purpose: Ensure all links work correctly
require(plyr)
require(grid)
require(gridExtra)
require(httr)


#splitting function for SDYID
get_sdyid <- function(subjectid){
  id_table <- read.table(text = subjectid, sep = ".", as.is = T, fill = T)
  sdyid <- id_table[1,2]
}

# test link function
test_link <- function(study_id,filename,link_text){
  
  if(is.na(filename)){
    result <- "NA"
  }else if(filename == "<NA>"){
    result <- "NA"
  }else{
    link <- paste0("https://www.immunespace.org/_webdav/Studies/SDY",study_id,"/@files/rawdata/",link_text,"/",filename) 
    tryCatch({
      status <- GET(link)
      result <- status$status_code
    },
    error = function(e){
      result <- "Error Occurred in GET(link)"
    },
    warning = function(w){
      message(w)
    })
  }
  return(result)
}

# Pull data from labkey
get_data <- function(filetype, filename_col,link_text){
  df <- labkey.selectRows(baseUrl = "https://www.immunespace.org",
                         "/Studies",
                         schemaName = "study",
                         queryName = filetype,
                         colNameOpt = "caption")
  
  # get SDY ids into a table
  Sdyids <- apply(MARGIN=1, df, FUN = function(x) get_sdyid(x))
  
  # get filenames into separate table / list
  Filenames <- df[filename_col]
  
  # test link and output to temp table
  Results <- mapply(test_link, Sdyids, Filenames, link_text)
  
  # merge temp table with comb df as third column
  #final_df <- cbind(Sdyids,Filenames,Results)
  
  # return table
  return(Results)
}

  

#fcs_table <- get_data("fcs_sample_files","File Info Name","fcs")
gene_expr_table <- get_data("gene_expression_files","File Info Name","gene_expression")

#sorted_ge <- gene_expr_table[order(Results),]
#sort by result

#TODO: figure out protocols ... different method of data extraction?

# Print tables to pdf
ts <-paste(format(Sys.time(), "%Y_%m_%d %T"), "pdf", sep = ".")
pdfname <- paste0("Test_IS_Links_", ts)
pdf(pdfname)
grid.table(gene_expr_table)
dev.off()



