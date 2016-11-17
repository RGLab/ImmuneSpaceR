test_files=function(what = c("gene_expression_files", "fcs_sample_files", "protocol")){
  
    labkey.url.base <- "https://www.immunespace.org"
    ret <- list()
    what <- tolower(what)
       
    for(i in what){
         
    # handle gene expr / fcs separately from protocols bc similar link construction
      if(i == "gene_expression_files" | i == "fcs_sample_files"){ 
        df <- .self$getDataset(i, original_view = TRUE)
        df <- df[!is.na(file_info_name)]
        df <- unique(df[, list(study_accession, file_info_name)])
                   
        numrow <- nrow(df)
        link_text <- ""
                       
        if(i == "gene_expression_files"){
          link_text <- "gene_expression"
        }else if( i == "fcs_sample_files"){
          link_text <- "flow_cytometry"
        }
                       
        links <- paste0(labkey.url.base, "/_webdav/", "/Studies/", 
                      df$study_accession, "/%40files/rawdata/",link_text,
                      sapply(df$file_info_name, URLencode))
                    
        http_status <- unlist(mclapply(links, link_test, mc.cores = detectCores()))
        num_good_links <- length(which(http_status == 200))
        bound_res <- cbind(df,http_status)
                   
        print(paste0(num_good_links, "/", numrow, " ", i, " with valid links."))
                           
        ret[[i]] <- bound_res
                           
        #handle protocols alone 
      }else{

        # if all studies, then pull links from folders list.  Assumption is that each SDY folder should have a protocol.
        if(.self$.isProject()){
          folders_list <- labkey.getFolders(baseUrl = labkey.url.base, folderPath = "/Studies/")
          studies <- unlist(folders_list[1])
          studies <- studies [! studies %in% c("SDY_template","Studies")]
          links <- lapply(studies, make_link)
          http_status <- unlist(mclapply(links, link_test, mc.cores = detectCores()))
          num_good_links <- length(which(http_status == 200))
          bound_res <- cbind(studies,http_status)
                                      
          print(paste0(num_good_links, "/", length(studies), " ", i, " with valid links."))
                                       
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
