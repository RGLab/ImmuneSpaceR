#'@docType package
#'@title A Thin Wrapper Around Immunspace.
#'@description RImmuneSpace provides a convenient API for accessing data sets within the ImmuneSpace database.
#'
#'@details Uses the Rlabkey package to connect to ImmuneSpace. Implements caching, and convenient methods for accessing data sets.
#'
#'@name RImmuneSpace-package
#'@aliases RImmuneSpace
#'@author Greg Finak
#'@import data.table Rlabkey
NULL

#mycon<-new("ImmuneSpaceConnection")
#mycon
#mycon$listDatasets()
# mycon$getDataset("hai")
# mycon$getDataset("elisa_mbaa")
# mycon$getDataset("demographics")
# mycon$getDataset("gene_expression_files")
# mycon$getDataset("cohort_membership")
#mycon$getGEMatrix("LAIV_2008")

#'@name ImmuneSpaceConnection
#'@aliases ImmuneSpace
#'@docType class
#'@title The ImmuneSpaceConnection class
#'@description Instantiate this class to access a study
#'@details Uses global variables \code{labkey.url.base}, and \code{labkey.url.path}, to access a study.
#'\code{labkey.url.base} should be \code{https://www.immunespace.org/}.
#'\code{labkey.url.path} should be \code{/Studies/studyname}, where 'studyname' is the name of the study.
#'The ImmunspaceConnection will initialize itself, and look for a \code{.netrc} file in \code{"~/"} the user's home directory.
#'The \code{.netrc} file should contain a \code{machine}, \code{login}, and \code{password} entry to allow access to ImmuneSpace,
#'where \code{machine} is the host name like "www.immunespace.org".
#'@usage new("ImmuneSpaceConnection")
#'@seealso \code{\link{RImmuneSpace-package}} \code{\link{ImmuneSpaceConnection_getGEMatrix}}  \code{\link{ImmuneSpaceConnection_getDataset}}  \code{\link{ImmuneSpaceConnection_listDatasets}}
#'@exportClass ImmuneSpaceConnection
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-new("ImmuneSpaceConnection")
#'sdy269
#'@return An instance of an ImmuneSpaceConnection for a study in `labkey.url.path`
setRefClass(Class = "ImmuneSpaceConnection",fields = list(study="character",config="list",available_datasets="data.table",data_cache="list"),
            methods=list(
              initialize=function(){
                .AutoConfig()
                .GeneExpressionMatrices()
                .GeneExpressionInputs()
              },
              .AutoConfig=function(){
                labkey.url.base<-get("labkey.url.base",.GlobalEnv)
                labkey.url.path<-get("labkey.url.path",.GlobalEnv)
                labkey.user.email<-get("labkey.user.email",.GlobalEnv)
                study<<-basename(labkey.url.path)
                config<<-list(labkey.url.base=labkey.url.base,labkey.url.path=labkey.url.path,labkey.user.email=labkey.user.email)
                .getDataSets();
              },
            show=function(){
              cat(sprintf("Immunspace Connection to study %s\n",study))
              cat(sprintf("URL: %s\n",file.path(gsub("/$","",config$labkey.url.base),gsub("^/","",config$labkey.url.path))))
              cat(sprintf("User: %s\n",config$labkey.user.email))
              cat("Available datasets\n")
              for(i in 1:nrow(available_datasets)){
                cat(sprintf("\t%s\n",available_datasets[i,Name]))
              }
              if(!is.null(data_cache[["GE_matrices"]])){
                cat("Expression Matrices\n")
                for(i in 1:nrow(data_cache[["GE_matrices"]])){
                  cat(sprintf("%s\n",data_cache[["GE_matrices"]][i,name]))
                }
              }
            },
            .getDataSets=function(){
              if(length(available_datasets)==0){
                available_datasets<<-data.table(labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "study",queryName = "DataSets"))[,list(Label,Name,Description,`Key Property Name`)]
              }
            },
            getDataset=function(x,reload=FALSE){
              if(nrow(available_datasets[Name%in%x])==0){
                stop(sprintf("Invalid data set: %s",x))
              }else{
                if(!is.null(data_cache[[x]])&!reload){
                  data_cache[[x]]
                }else{
                  data_cache[[x]]<<-data.table(labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "study",queryName = x,colNameOpt = "fieldname")) 
                  setnames(data_cache[[x]],.munge(colnames(data_cache[[x]])))
                  data_cache[[x]]
                }                
              }
            },
            listDatasets=function(){
              for(i in 1:nrow(available_datasets)){
                cat(sprintf("\t%s\n",available_datasets[i,Name]))
              }
              if(!is.null(data_cache[["GE_matrices"]])){
                cat("Expression Matrices\n")
                for(i in 1:nrow(data_cache[["GE_matrices"]])){
                  cat(sprintf("%s\n",data_cache[["GE_matrices"]][i,name]))
                }
              }
            },
            .munge=function(x){
              tolower(gsub(" ","_",basename(x)))
            },
            .GeneExpressionInputs=function(){
              if(!is.null(data_cache[["GE_inputs"]])){
                data_cache[["GE_inputs"]]
              }else{
                ge<-data.table(labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "assay.ExpressionMatrix.matrix",queryName = "InputSamples",colNameOpt = "fieldname",viewName = "gene_expression_matrices"))
                setnames(ge,.munge(colnames(ge)))
                data_cache[["GE_inputs"]]<<-ge
              }
            },
            .GeneExpressionMatrices=function(){
              if(!is.null(data_cache[["GE_matrices"]])){
                data_cache[["GE_matrices"]]
              }else{
                ge<-data.table(labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "assay.ExpressionMatrix.matrix",queryName = "Runs",colNameOpt = "fieldname",showHidden = TRUE, viewName = "Expression Matrices"))
                setnames(ge,.munge(colnames(ge)))
                data_cache[["GE_matrices"]]<<-ge
              }
            },
            .downloadMatrix=function(x){
              if(is.null(data_cache[[x]])){
                if(nrow(data_cache[["GE_matrices"]][name%in%x,])==0){
                  stop(sprintf("No matrix %s in study\n",x))
                }
                pb <- txtProgressBar(style=3)
                setTxtProgressBar(pb,0)
                link<-URLdecode(file.path(gsub("www.","",gsub("http:","https:",gsub("/$","",config$labkey.url.base))),gsub("^/","",data_cache[["GE_matrices"]][name%in%x,downloadlink])))
                opts<-curlOptions(.opts=list(netrc=TRUE,ssl.verifyhost=FALSE,httpauth=1L,ssl.verifypeer=FALSE,followlocation=TRUE,verbose=FALSE))
                handle<-getCurlHandle(.opts=opts)
                h<-basicTextGatherer()
                curlPerform(url=link,curl=handle,writefunction=h$update)
                setTxtProgressBar(pb,0.33)
                con<-textConnection(h$value())
                setTxtProgressBar(pb,0.66)
                data_cache[[x]]<<-data.table(read.csv(con,sep="\t"))
                setTxtProgressBar(pb,1)
                close(pb)
              }else{
                data_cache[[x]]
              }
            },
            getGEMatrix=function(x){
              if(x%in%names(data_cache)){
                data_cache[[x]]              
              }else{
                .downloadMatrix(x)
                data_cache[[x]]
              }
            }
))

#'@title get Gene Expression Matrix
#'@usage getGEMatrix(x)
#'@aliases getGEMatrix
#'@param x \code{"character"} name of the Gene Expression Matrix
#'@details Returns the gene expression matrix named 'x', downloads it if it is not already cached.
#'@return a \code{data.table}
#'@name ImmuneSpaceConnection_getGEMatrix
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-new("ImmuneSpaceConnection")
#'sdy269$getGEMatrix("TIV_2008")
NULL

#'@title get a dataset
#'@aliases getDataset
#'@usage getDataset(x)
#'@param x \code{"character"} name of the dataset
#'@details Returns the dataset named 'x', downloads it if it is not already cached.
#'@return a \code{data.table}
#'@name ImmuneSpaceConnection_getDataset
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-new("ImmuneSpaceConnection")
#'sdy269$getDataset("hai")
NULL

#'@title list available datasets
#'@aliases listDatasets
#'@usage listDatasets()
#'@details Prints the names of the available datasets
#'@return Doesn't return anything, just prints to console.
#'@name ImmuneSpaceConnection_listDatasets
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-new("ImmuneSpaceConnection")
#'sdy269$listDatasets()
NULL