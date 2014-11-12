#'@name ImmuneSpaceConnectionList
#'@aliases ImmuneSpaceConnectionList-class
#'@aliases ImmuneSpaceConnectionList
#'@rdname ImmuneSpaceConnectionList-class
#'@docType class
#'@title The ImmuneSpaceConnectionList class
#'@description Instantiate this class to access a list of studies
#'@seealso \code{\link{ImmuneSpaceConnection-class}}
#'@exportClass ImmuneSpaceConnectionList
#'@examples
#'labkey.url.base <- "https://www.immunespace.org"
#'labkey.url.path <- "/Studies/SDY269"
#'labkey.user.email <- 'gfinak at fhcrc.org'
#'cons <- CreateConnection(c("SDY269", "SDY180"))
#'cons
#'@return An instance of an ImmuneSpaceConnection for a study in `labkey.url.path`
.ISConList <- setRefClass(Class = "ImmuneSpaceConnectionList", fields = list(connections="list"))

#' print names of the studies associated with the connection list.
#' @examples
#' \dontrun{
#' cons <- CreateConnection(c("SDY269", "SDY180"))
#' cons$studies
#' }
.ISConList$methods(
  studyNames = function(){
    names(connections)
  })



.ISConList$methods(
  show=function(){
    cat(length(connections), "Immunespace Connections:\n\t")
    cat(studyNames(),sep = "\n\t")
    cat("use study('xxx') method to access each individual ImmuneSpaceConnection object.")
    }
)
.ISConList$methods(
  study = function(name){
    connections[[name]]
  })

.ISConList$methods(
  listDatasets=function(){
#       browser()
      cat("Datasets\n")
      # get ds from each con
      dsList <- lapply(connections, function(con)con$available_datasets[, Name])
      # get union set from all cons
      all_ds <- unique(unlist(dsList))
      # get 0,1 marks for each con
      df <- ldply(dsList, function(ds)as.integer(all_ds%in%ds))
      # transpose table to display
      rownames(df) <- df[, ".id"]
      df[, ".id"] <- NULL
      df <- t(df)
      rownames(df) <- all_ds
      print(df)
      
    
      cat("\nExpression Matrices\n")
      for(con in connections)
      {
          mat_name <- con$constants$matrices
          gem_names <- con$data_cache[[mat_name]][,"name"]
          cat("\tStudy:", con$study, "\n\t\t\t")
          cat(gem_names, sep = "\n\t\t\t")
      }
    
      
    
  })

.ISConList$methods(
  listGEAnalysis = function(merge = TRUE){
    GEA_list <- lapply(connections, function(con)con$listGEAnalysis())
    if(merge)
      ldply(GEA_list, function(gea)gea,.id = "study")
    else
      GEA_list
      })

.ISConList$methods(
    getDataset = function(x, merge = TRUE, ...){
    
    dsList <- lapply(connections, function(con){
                  ds <- con$getDataset(x = x, ...)
                  if(!is.null(ds))
                    ds[, study:= con$study]
                  ds
      })
    if(merge)
      rbindlist(dsList)
    else
      dsList
  })

.ISConList$methods(
  clear_cache = function(){
    for(con in connections)
      con$clear_cache()
  }
)

.ISConList$methods(
  getGEMatrix=function(x, merge = TRUE, ...){
#     browser()
    eSetList <- unlist(lapply(connections, function(con){
      
                      unlist(lapply(x, function(thisName){
                        
                            eSet <- try(con$getGEMatrix(x = thisName, ...), silent = T)  
                            if(class(eSet) == "try-error")
                              NULL
                            else
                              eSet
                          }))
                    }))
    
    
    
  }
)                                                    