# Returns a named logical where TRUE marks files that are accessible.
# This function is used for administrative purposes to check that the flat files
# are properly loaded and accessible to the users.
#' @import httr
#' @import parallel
#' @import Rlabkey
.ISCon$methods(
  .test_files=function(what = c("gene_expression_files", "fcs_sample_files", "protocol")){
    
    labkey.url.base <- .self$config$labkey.url.path
    ret <- list()
    what <- tolower(what)
    
    for(i in what){
      
      # handle gene expr / fcs separately from protocols bc similar link construction
      if(i == "gene_expression_files" | i == "fcs_sample_files"){ 
        df <- .self$getDataset(i, original_view = TRUE)
        df <- df[!is.na(file_info_name)]
        df <- unique(df[, list(study_accession, file_info_name)])
        
        link_text <- ""
        if(i == "gene_expression_files"){
          link_text <- "gene_expression"
        }else if( i == "fcs_sample_files"){
          link_text <- "flow_cytometry"
        }
        
        links <- paste0(labkey.url.base, "/_webdav/", "/Studies/", 
                        df$study_accession, "/%40files/rawdata/",link_text,
                        sapply(df$file_info_name, URLencode))
        
        numrow <- nrow(df)
        
        bound_res <- res_table_maker(links_to_test = links,info_table = df, numrow = numrow, filetype = i)
        
        ret[[i]] <- bound_res
        
        #handle protocols alone 
      }else{
        
        # if all studies, then pull links from folders list.  Assumption is that each SDY folder should have a protocol.
        if(.self$.isProject()){
          folders_list <- labkey.getFolders(baseUrl = labkey.url.base, folderPath = "/Studies/")
          studies <- unlist(folders_list[1])
          studies <- studies [! studies %in% c("SDY_template","Studies")]
          links <- lapply(studies, make_link)
          numrow <- length(studies)
          
          bound_res <- res_table_maker(links_to_test = links,info_table = studies, numrow = numrow, filetype = i)
          
          ret[[i]] <- bound_res
          
          #if single study, then id study from the path value from CreateConnection
        }else{
          study_string <- strsplit(.self$config$labkey.url.path, "/")
          study <- study_string[[1]][3]
          res <- link_test(paste0(labkey.url.base, "/_webdav/Studies/",
                                  study, "/%40files/protocols/", folder,
                                  "_protocol.zip"))
          ret[[i]] <- res
        }
      }
    }
    return(ret)
  }
)

#-----helper functions for test_files()-----------
link_test <- function(link){
  info <- GET(link)
  status <- info$status_code
  return(status)
}

make_link <- function(study){
  link <- paste0(labkey.url.base, "/_webdav/Studies/",
                 study, "/%40files/protocols/", study,
                 "_protocol.zip")
  return(link)
}

res_table_maker <- function(filetype,numrow,links_to_test,info_table){
  http_status <- unlist(mclapply(links_to_test, link_test, mc.cores = detectCores()))
  bound_res <- cbind(info_table,http_status)
  num_good_links <- length(which(http_status == 200))
  print(paste0(num_good_links, "/", numrow, " ", filetype, " with valid links."))
  return(bound_res)
}