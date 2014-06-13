#'@docType package
#'@title A Thin Wrapper Around Immunspace.
#'@description ImmuneSpaceR provides a convenient API for accessing data sets within the ImmuneSpace database.
#'
#'@details Uses the Rlabkey package to connect to ImmuneSpace. Implements caching, and convenient methods for accessing data sets.
#'
#'@name ImmuneSpaceR-package
#'@aliases ImmuneSpaceR
#'@author Greg Finak
#'@import data.table Rlabkey methods Biobase
NULL

#'@title CreateConnection
#'@name CreateConnection
#'@param study \code{"character"} vector naming the study.
#'@description Constructor for \code{ImmuneSpaceConnection} class
#'@details Instantiates and \code{ImmuneSpaceConnection} for \code{study}
#'The constructor will try to take the values of the various `labkey.*` parameters from the global environment.
#'If they don't exist, it will use default values. These are assigned to `options`, which are then used by the \code{ImmuneSpaceConnection} class.
#'@export CreateConnection
#'@return an instance of an \code{ImmuneSpaceConnection}
CreateConnection = function(study=NULL){
  if(is.null(study)){
    stop("study cannot be NULL");
  }
  
  labkey.url.path<-try(get("labkey.url.path",.GlobalEnv),silent=TRUE)
  if(inherits(labkey.url.path,"try-error")){
    labkey.url.path<-paste0("/Studies/",study)
  }else{
    labkey.url.path<-file.path(dirname(labkey.url.path),study)
  }
  labkey.url.base<-try(get("labkey.url.base",.GlobalEnv),silent=TRUE)
  if(inherits(labkey.url.base,"try-error"))
    labkey.url.base<-"https://www.immunespace.org"
  labkey.url.base<-gsub("http:","https:",labkey.url.base)
  labkey.user.email<-try(get("labkey.user.email",.GlobalEnv),silent=TRUE)
  if(inherits(labkey.user.email,"try-error"))
    labkey.user.email="unknown_user at not_a_domain.com"
  
  
  options(labkey.url.base=labkey.url.base)
  options(labkey.url.path=labkey.url.path)
  options(labkey.user.email=labkey.user.email)
  
  new("ImmuneSpaceConnection")
}

#'@name ImmuneSpaceConnection
#'@aliases ImmuneSpaceConnection-class
#'@aliases ImmuneSpace
#'@rdname ImmuneSpaceConnection-class
#'@docType class
#'@title The ImmuneSpaceConnection class
#'@description Instantiate this class to access a study
#'@details Uses global variables \code{labkey.url.base}, and \code{labkey.url.path}, to access a study.
#'\code{labkey.url.base} should be \code{https://www.immunespace.org/}.
#'\code{labkey.url.path} should be \code{/Studies/studyname}, where 'studyname' is the name of the study.
#'The ImmunspaceConnection will initialize itself, and look for a \code{.netrc} file in \code{"~/"} the user's home directory.
#'The \code{.netrc} file should contain a \code{machine}, \code{login}, and \code{password} entry to allow access to ImmuneSpace,
#'where \code{machine} is the host name like "www.immunespace.org".
#'@seealso \code{\link{ImmuneSpaceR-package}} \code{\link{ImmuneSpaceConnection_getGEMatrix}}  \code{\link{ImmuneSpaceConnection_getDataset}}  \code{\link{ImmuneSpaceConnection_listDatasets}}
#'@exportClass ImmuneSpaceConnection
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269
#'@return An instance of an ImmuneSpaceConnection for a study in `labkey.url.path`
setRefClass(Class = "ImmuneSpaceConnection",fields = list(study="character",config="list",available_datasets="data.table",data_cache="list",constants="list"),
            methods=list(
              initialize=function(){
                constants<<-list(matrices="GE_matrices",matrix_inputs="GE_inputs")
                .AutoConfig()
                gematrices_success<-try(.GeneExpressionMatrices(),silent=TRUE)
                geinputs_success<-try(.GeneExpressionInputs(),silent=TRUE)
                if(inherits(gematrices_success,"try-error")){
                  message("No gene expression data")
                }
              },
              
              .AutoConfig=function(){
                #should use options
                labkey.url.base<-getOption("labkey.url.base")
                labkey.url.path<-getOption("labkey.url.path")
                labkey.user.email<-getOption("labkey.user.email")
                
                study<<-basename(labkey.url.path)
                config<<-list(labkey.url.base=labkey.url.base,labkey.url.path=labkey.url.path,labkey.user.email=labkey.user.email)
                .getAvailableDataSets();
              },
              
            show=function(){
              cat(sprintf("Immunspace Connection to study %s\n",study))
              cat(sprintf("URL: %s\n",file.path(gsub("/$","",config$labkey.url.base),gsub("^/","",config$labkey.url.path))))
              cat(sprintf("User: %s\n",config$labkey.user.email))
              cat("Available datasets\n")
              for(i in 1:nrow(available_datasets)){
                cat(sprintf("\t%s\n",available_datasets[i,Name]))
              }
              if(!is.null(data_cache[[constants$matrices]])){
                cat("Expression Matrices\n")
                for(i in 1:nrow(data_cache[[constants$matrices]])){
                  cat(sprintf("%s\n",data_cache[[constants$matrices]][i,"name"]))
                }
              }
            },
            
            .getAvailableDataSets=function(){
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
              if(!is.null(data_cache[[constants$matrices]])){
                cat("Expression Matrices\n")
                for(i in 1:nrow(data_cache[[constants$matrices]])){
                  cat(sprintf("%s\n",data_cache[[constants$matrices]][i,"name"]))
                }
              }
            },
            
            .munge=function(x){
              tolower(gsub(" ","_",basename(x)))
            },
            
            .GeneExpressionInputs=function(){
              if(!is.null(data_cache[[constants$matrix_inputs]])){
                data_cache[[constants$matrix_inputs]]
              }else{
                ge<-labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "assay.ExpressionMatrix.matrix",queryName = "InputSamples",colNameOpt = "fieldname",viewName = "gene_expression_matrices",showHidden=TRUE)
                setnames(ge,.munge(colnames(ge)))
                data_cache[[constants$matrix_inputs]]<<-ge
              }
            },
            
            .GeneExpressionFeatures=function(matrix_name){
              if(!any((data_cache[[constants$matrices]][,"name"]%in%matrix_name))){
                stop("Invalid gene expression matrix name");
              }
              annotation_set_id<-.getFeatureId(matrix_name)
              if(is.null(data_cache[[.mungeFeatureId(annotation_set_id)]])){
                message("Downloading Features..")
                featureAnnotationSetQuery=sprintf("SELECT * from FeatureAnnotation where FeatureAnnotationSetId='%s';",annotation_set_id);
                features<-labkey.executeSql(config$labkey.url.base,config$labkey.url.path,schemaName = "Microarray",sql = featureAnnotationSetQuery ,colNameOpt = "fieldname")
                data_cache[[.mungeFeatureId(annotation_set_id)]]<<-features
              }
            },
            
            .GeneExpressionMatrices=function(){
              if(!is.null(data_cache[[constants$matrices]])){
                data_cache[[constants$matrices]]
              }else{
                ge<-labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "assay.ExpressionMatrix.matrix",queryName = "Runs",colNameOpt = "fieldname",showHidden = TRUE, viewName = "Expression Matrices")
                setnames(ge,.munge(colnames(ge)))
                data_cache[[constants$matrices]]<<-ge
              }
            },
            
            .downloadMatrix=function(x){
              if(is.null(data_cache[[x]])){
                if(nrow(subset(data_cache[[constants$matrices]],name%in%x))==0){
                  stop(sprintf("No matrix %s in study\n",x))
                }
                link<-URLdecode(file.path(gsub("www.","",gsub("http:","https:",gsub("/$","",config$labkey.url.base))),gsub("^/","",subset(data_cache[[constants$matrices]],name%in%x)[,"downloadlink"])))
                opts<-curlOptions(.opts=list(netrc=TRUE,ssl.verifyhost=FALSE,httpauth=1L,ssl.verifypeer=FALSE,followlocation=TRUE,verbose=FALSE))
                handle<-getCurlHandle(.opts=opts)
                h<-basicTextGatherer()
                message("Downloading matrix..")
                curlPerform(url=link,curl=handle,writefunction=h$update)
                fl<-tempfile()
                write(h$value(),file=fl)
                data_cache[[x]]<<-fread(fl,header=TRUE)
                file.remove(fl)
              }else{
                data_cache[[x]]
              }
            },
            
            getGEMatrix=function(x){
              if(x%in%names(data_cache)){
                data_cache[[x]]              
              }else{
                .downloadMatrix(x)
                .GeneExpressionFeatures(x)
                .ConstructExpressionSet(x)
                data_cache[[x]]
              }
            },
          .ConstructExpressionSet=function(matrix_name){
            #matrix
            message("Constructing ExpressionSet")
            matrix<-data_cache[[matrix_name]]
            #features
            features<-data_cache[[.mungeFeatureId(.getFeatureId(matrix_name))]][,c("FeatureId","GeneSymbol")]
            #inputs
            pheno<-unique(subset(data_cache[[constants$matrix_inputs]],biosample_accession%in%colnames(matrix))[,c("biosample_accession","subject_accession","arm_name","study_time_collected")])
            try(setnames(matrix," ","FeatureId"),silent=TRUE)
            setkey(matrix,FeatureId)
            rownames(features)<-features$FeatureId
            features<-features[matrix$FeatureId,]#order feature info
            rownames(pheno)<-pheno$biosample_accession
            pheno<-pheno[colnames(matrix)[-1L],]
            ad_pheno<-AnnotatedDataFrame(data=pheno)
            ad_features<-AnnotatedDataFrame(features)
            es<-ExpressionSet(assayData=as.matrix(matrix[,-1L,with=FALSE]),phenoData=ad_pheno,featureData=ad_features)
            data_cache[[matrix_name]]<<-es
          },
          .getFeatureId=function(matrix_name){
            subset(data_cache[[constants$matrices]],name%in%matrix_name)[,"featureset"]
          },
          .mungeFeatureId=function(annotation_set_id){
            return(sprintf("featureset_%s",annotation_set_id))
          }
))

#'@title get Gene Expression Matrix
#'@aliases getGEMatrix
#'@param x \code{"character"} name of the Gene Expression Matrix
#'@details Returns an `ExpressionSet` from the matrix named 'x', downloads it if it is not already cached.
#'@return an \code{ExpressionSet}
#'@name ImmuneSpaceConnection_getGEMatrix
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269$getGEMatrix("TIV_2008")
NULL

#'@title get a dataset
#'@aliases getDataset
#'@param x \code{"character"} name of the dataset
#'@details Returns the dataset named 'x', downloads it if it is not already cached.
#'@return a \code{data.table}
#'@name ImmuneSpaceConnection_getDataset
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269$getDataset("hai")
NULL

#'@title list available datasets
#'@aliases listDatasets
#'@details Prints the names of the available datasets
#'@return Doesn't return anything, just prints to console.
#'@name ImmuneSpaceConnection_listDatasets
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269$listDatasets()
NULL