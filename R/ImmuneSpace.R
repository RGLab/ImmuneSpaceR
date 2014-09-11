#'@docType package
#'@title A Thin Wrapper Around ImmuneSpace.
#'@description ImmuneSpaceR provides a convenient API for accessing data sets within the ImmuneSpace database.
#'
#'@details Uses the Rlabkey package to connect to ImmuneSpace. Implements caching, and convenient methods for accessing data sets.
#'
#'@name ImmuneSpaceR-package
#'@aliases ImmuneSpaceR
#'@author Greg Finak
#'@import data.table Rlabkey methods Biobase gtools
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
CreateConnection = function(study=NULL, verbose = FALSE){
  labkey.url.path<-try(get("labkey.url.path",.GlobalEnv),silent=TRUE)
  if(inherits(labkey.url.path,"try-error")){
    if(is.null(study)){
      stop("study cannot be NULL")
    }
    labkey.url.path<-paste0("/Studies/",study)
  }else if(!is.null(study)){
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
  options(ISverbose = verbose)
  
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
#'The ImmunespaceConnection will initialize itself, and look for a \code{.netrc} file in \code{"~/"} the user's home directory.
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
setRefClass(Class = "ImmuneSpaceConnection",
            fields = list(study = "character", config="list",
                          available_datasets = "data.table",
                          data_cache="list",constants="list"),
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
                verbose <- getOption("ISverbose")
                curlOptions <- curlOptions(netrc=TRUE,ssl.verifyhost=FALSE,
                                            httpauth=1L,ssl.verifypeer=FALSE,
                                            followlocation=TRUE,verbose=FALSE,
                                            ssl.cipher.list="ALL")
                study<<-basename(labkey.url.path)
                config<<-list(labkey.url.base=labkey.url.base,
                              labkey.url.path=labkey.url.path,
                              labkey.user.email=labkey.user.email,
                              curlOptions = curlOptions,
                              verbose = verbose)
                .checkStudy(config$verbose)
                .getAvailableDataSets();
              },
              
            show=function(){
              cat(sprintf("Immunespace Connection to study %s\n",study))
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
            
            .checkStudy = function(verbose = FALSE){
              if(length(available_datasets)==0){
                validStudies <- mixedsort(grep("^SDY", basename(lsFolders(getSession(config$labkey.url.base, "Studies"))), value = TRUE))
                req_study <- basename(config$labkey.url.path) 
                if(!req_study %in% validStudies){
                  if(!verbose){
                    stop(paste0(req_study, " is not a valid study"))
                  } else{
                    stop(paste0(req_study, " is not a valid study\nValid studies: ",
                                paste(validStudies, collapse=", ")))
                  }
                }
              }
            },
            .getAvailableDataSets=function(){
              if(length(available_datasets)==0){
                dataset_filter <- makeFilter(c("showbydefault", "EQUAL", TRUE))
                available_datasets<<-data.table(labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "study",queryName = "DataSets", colFilter = dataset_filter))[,list(Label,Name,Description,`Key Property Name`)]
              }
            },
            
            getDataset=function(x, original_view = FALSE, reload=FALSE, ...){
              if(nrow(available_datasets[Name%in%x])==0){
                stop(sprintf("Invalid data set: %s",x))
              }else{
                if(!is.null(data_cache[[x]])&!reload){
                  data_cache[[x]]
                }else{
                  viewName <- NULL
                  if(original_view){
                    viewName <- "full"
                  }
                  data_cache[[x]] <<- data.table(labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "study", queryName = x, viewName = viewName, colNameOpt = "fieldname", ...)) 
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
            
            .GeneExpressionFeatures=function(matrix_name,summary=FALSE){
              if(!any((data_cache[[constants$matrices]][,"name"]%in%matrix_name))){
                stop("Invalid gene expression matrix name");
              }
              annotation_set_id<-.getFeatureId(matrix_name)              
              if(is.null(data_cache[[.mungeFeatureId(annotation_set_id)]])){
                if(!summary){
                  message("Downloading Features..")
                  featureAnnotationSetQuery=sprintf("SELECT * from FeatureAnnotation where FeatureAnnotationSetId='%s';",annotation_set_id);
                  features<-labkey.executeSql(config$labkey.url.base,config$labkey.url.path,schemaName = "Microarray",sql = featureAnnotationSetQuery ,colNameOpt = "fieldname")
                }else{
                  features<-data.frame(FeatureId=con$data_cache[[matrix_name]][,gene_symbol],GeneSymbol=con$data_cache[[matrix_name]][,gene_symbol])
                }
                data_cache[[.mungeFeatureId(annotation_set_id)]]<<-features
              }
            },
            
            .GeneExpressionMatrices=function(){
              if(!is.null(data_cache[[constants$matrices]])){
                data_cache[[constants$matrices]]
              }else{
                ge<-labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "assay.ExpressionMatrix.matrix",queryName = "Runs",colNameOpt = "fieldname",showHidden = TRUE, viewName = "expression_matrices")
                setnames(ge,.munge(colnames(ge)))
                data_cache[[constants$matrices]]<<-ge
              }
            },
            
            .downloadMatrix=function(x, summary = FALSE){
              if(is.null(data_cache[[x]])){
                if(nrow(subset(data_cache[[constants$matrices]],name%in%x))==0){
                  stop(sprintf("No matrix %s in study\n",x))
                }
                summary <- ifelse(summary, ".summary", "")
                #link<-URLdecode(file.path(gsub("www.","",gsub("http:","https:",gsub("/$","",config$labkey.url.base))), paste0(gsub("^/","",subset(data_cache[[constants$matrices]],name%in%x)[,"downloadlink"]),summary)))
               
                #shouldn't be removing the www reported by labkey. Fix your netrc entry instead
                link<-URLdecode(file.path(gsub("http:","https:",gsub("/$","",config$labkey.url.base)), 
                                          "_webdav", gsub("^/","",config$labkey.url.path), "@files/analysis/exprs_matrices",
                                          paste0(x, ".tsv", summary)))
                localpath<-.localStudyPath(link)
                if(.isRunningRemotely(localpath)){
                  fl<-localpath
                  message("Reading local matrix")
                  data_cache[[x]]<<-fread(fl,header=TRUE)                
                }else{                
                  opts <- config$curlOptions
                  handle<-getCurlHandle(.opts=opts)
                  h<-basicTextGatherer()
                  message("Downloading matrix..")
                  curlPerform(url=link,curl=handle,writefunction=h$update)
                  fl<-tempfile()
                  write(h$value(),file=fl)
                  data_cache[[x]]<<-fread(fl,header=TRUE)
                  file.remove(fl)
                }

              }else{
                data_cache[[x]]
              }
            },
            
            getGEMatrix=function(x, summary = FALSE){
              if(x%in%names(data_cache)){
                data_cache[[x]]              
              }else{
                .downloadMatrix(x, summary)
                
                  .GeneExpressionFeatures(x,summary)
                
                .ConstructExpressionSet(x, summary)
                data_cache[[x]]
              }
            },
          .ConstructExpressionSet=function(matrix_name, summary){
            #matrix
            message("Constructing ExpressionSet")
            matrix<-data_cache[[matrix_name]]
            #features
              features<-data_cache[[.mungeFeatureId(.getFeatureId(matrix_name))]][,c("FeatureId","GeneSymbol")]
            #inputs
            pheno<-unique(subset(data_cache[[constants$matrix_inputs]],biosample_accession%in%colnames(matrix))[,c("biosample_accession","subject_accession","arm_name","study_time_collected")])
            
            if(summary){
              fdata <- data.frame(FeatureId = matrix$gene_symbol, gene_symbol = matrix$gene_symbol, row.names = matrix$gene_symbol)
              fdata <- AnnotatedDataFrame(fdata)
            } else{
              try(setnames(matrix," ","FeatureId"),silent=TRUE)
              setkey(matrix,FeatureId)
              rownames(features)<-features$FeatureId
              features<-features[matrix$FeatureId,]#order feature info
              fdata <- AnnotatedDataFrame(features)
            }
            rownames(pheno)<-pheno$biosample_accession
            pheno<-pheno[colnames(matrix)[-1L],]
            ad_pheno<-AnnotatedDataFrame(data=pheno)
            es<-ExpressionSet(assayData=as.matrix(matrix[,-1L,with=FALSE]),phenoData=ad_pheno,featureData=fdata)
            data_cache[[matrix_name]]<<-es
          },
          .getFeatureId=function(matrix_name){
            subset(data_cache[[constants$matrices]],name%in%matrix_name)[,"featureset"]
          },
          .mungeFeatureId=function(annotation_set_id){
            return(sprintf("featureset_%s",annotation_set_id))
          },
          .isRunningRemotely=function(path){
            file.exists(path)
          },
          .localStudyPath=function(urlpath){
            LOCALPATH<-"/shared/silo_researcher/Gottardo_R/immunespace"
            PRODUCTION_HOST<-"immunespace.org"
            STAGING_HOST<-"posey.fhcrc.org"
            PRODUCTION_PATH<-"production/files"
            STAGING_PATH<-"staging/files"
            if(grepl(PRODUCTION_HOST,urlpath)){
              PROCESS<-PRODUCTION_PATH
            }else if(grepl(STAGING_HOST,urlpath)){
              PROCESS<-STAGING_PATH
            }else{
              stop("Can't determine if we are running on immunespace (production) or posey (staging)")
            }
            gsub(file.path(gsub("/$","",config$labkey.url.base), "_webdav"), file.path(LOCALPATH,PROCESS), urlpath)
          },
          listGEAnalysis = function(){
            GEA <- labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
                                     "gene_expression", "gene_expression_analysis",
                                     colNameOpt = "rname")
            print(GEA)
          },
          getGEAnalysis = function(analysis_accession){
            "Get gene expression analysis resluts from a connection"
            if(missing(analysis_accession)){
              stop("Missing analysis_accession argument.
                   Use listGEAnalysis to get a list of available 
                   analysis_accession numbers")
            }
            AA_filter <- makeFilter(c("analysis_accession", "IN", analysis_accession))
            GEAR <- labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
                                     "gene_expression", "gene_expression_analysis_results",
                                     colFilter = AA_filter)
            colnames(GEAR) <- .munge(colnames(GEAR))
            return(GEAR)
          },
          quick_plot = function(dataset, normalize_to_baseline = TRUE,
                                type = "auto", filter = NULL, ...){
            palette <- rev(brewer.pal(n = 11, name = "RdYlBu"))
            addPar <- c("gender", "age_reported", "race")
            annoCols <- c("name", "subject_accession", "study_time_collected", addPar)
            toKeep <- c("response", "analyte", annoCols)
            
            e <- try({
              dt <- con$getDataset(dataset, reload = TRUE, colFilter = filter)
              if(length(grep("analyte",colnames(dt)))==0){
                if(type == "auto"){
                  type <- "boxplot"
                }
                dt <- dt[, analyte := ""]
              } else{
                if(type == "auto"){
                  type <- "heatmap"
                }
              }
              if(dataset == "elispot"){
                dt <- dt[, value_reported := (spot_number_reported +1) / cell_number_reported]
              } else if(dataset == "pcr"){
                if(all(is.na(dt[, threshold_cycles]))){
                  plot(1:10)
                  stop("PCR results cannot be displayed for studies that do not use threshold cycles")
                }
                dt <- dt[, analyte := entrez_gene_id]
              }
              dt <- dt[, response := mean(log2(value_reported), na.rm = TRUE),
                       by = "name,subject_accession,analyte,study_time_collected"]
              dt <- unique(dt[, toKeep, with = FALSE])
              
              if(normalize_to_baseline){
                dt <- dt[,response:=response-response[study_time_collected==0],
                         by="name,subject_accession,analyte"][study_time_collected!=0]
                ylab <- "Response normalized to baseline"
                scale <- "none"
              } else{
                ylab <- "Response (log2)"
                scale <- "row"
              }
            })
            
            if(inherits(e, "try-error")){
              type <- "error"
              error_string <- attr(e, "condition")$message
            }
            
            #dt: analyte, subject, timepoint, response, arm + addPar
            if(type == "heatmap"){
              mat <- acast(dt, analyte ~ name + study_time_collected + subject_accession, value.var = "response")
              anno <- data.frame(unique(dt[, annoCols, with = FALSE]))
              rownames(anno) <- paste(anno$name, anno$study_time_collected, anno$subject_accession, sep = "_")
              anno <- anno[, c("study_time_collected", "name")]
              max <- max(abs(mat))
              show_rnames <- TRUE
              if(scale == "none"){
                breaks <- seq(-max, max, length.out = length(palette) + 1)
              } else{
                breaks <- NA
              }
              cluster_rows <- ifelse(nrow(mat) > 2, TRUE, FALSE)
              p <- list(mat = mat, annotation = anno, show_colnames = FALSE,
                            show_rownames = show_rnames, cluster_cols = FALSE,
                            cluster_rows = cluster_rows, color = palette,
                            scale = scale, breaks = breaks)
              do.call("pheatmap", p)
            } else if(type == "boxplot"){
              p <- qplot(as.factor(study_time_collected), response, data = dt,
                         facets = analyte~name, geom = c("boxplot", "jitter"),
                         xlab = "time", ylab = ylab, 
                         group = study_time_collected, ...)
              print(p)
            } else if(type == "line"){
              p <- qplot(study_time_collected, response, data = dt,
                         facets = analyte~name, geom=c("line", "point"),
                         xlab = "time", ylab = ylab, group = subject_accession,
                         ...)
              print(p)
#               error_string <- "Lineplot has not been enabled yet"
#               data <- data.frame(x = 0, y = 0, err = error_string)
#               p <- ggplot(data = data) + geom_text(aes(x, y, label = err))
            } else{#} if(type == "error"){
              #grid.text(error_string)
              #grid.newpage()
              #grid.text(
              error_string <- "The datset you are trying to visualize does not follow HIPC standards.
                        Please use the built-in LabKey charts for this data."
              data <- data.frame(x = 0, y = 0, err = error_string)
              p <- ggplot(data = data) + geom_text(aes(x, y, label = err))
              print(p)
            }
            return(p)
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
#'@param x A \code{character}. The name of the dataset
#'@param original_view A \code{logical}. If set to TRUE, download the ImmPort view. Else,
#'  download the default grid view.
#'@param reload A \code{logical}. Clear the cache. If set to TRUE, download the 
#'  dataset, whether a cached versio exist or not.
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


#'@title list available gene expression analysis
#'@aliases listGEAnalysis
#'@details Prints the table of differential expression analysis
#'@return A \code{data.frame}. The list of gene expression analysis.
#'@name ImmuneSpaceConnection_listGEAnalysis
#'@examples
#'labkey.url.base="https://www.immunespace.org"
#'labkey.url.path="/Studies/SDY269"
#'labkey.user.email='gfinak at fhcrc.org'
#'sdy269<-CreateConnection("SDY269")
#'sdy269$listGEAnalysis()
NULL
