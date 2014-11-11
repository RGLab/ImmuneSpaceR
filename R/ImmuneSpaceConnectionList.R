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
  size = function(){
    length(connections)
  })

.ISConList$methods(
  show=function(){
    cat(size(), "Immunespace Connections:\n")
    cat(studyNames(),sep = "\n")
    cat("use study('xxx') method to access each connection object.")
    }
)
.ISConList$methods(
  study = function(name){
    connections[[name]]
  })

.ISConList$methods(
  listDatasets=function(){
#       browser()
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
      
    
#       cat("Expression Matrices\n")
#       connections[[2]]$data_cache[[con$constants$matrices]][,"name"]))
      
    
  })

.ISConList$methods(
  listGEAnalysis = function(){
    ldply(connections, function(con)con$listGEAnalysis(), .id = "study")
      })