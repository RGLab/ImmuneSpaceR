#' @include ImmuneSpace.R
NULL

#' Detect whether the code is running on an immunespace server
#' 
#' @return A \code{logical}. TRUE if the code is executed on one of the
#' immunespace servers. FALSE otherwise.
#' 
#' @name isRunningOnServer
#' 
.ISCon$methods(
  isRunningOnServer = function(){
    onServer <- ifelse(Sys.info()["user"] == "immunespace", TRUE, FALSE)
    names(onServer) <- NULL
    return(onServer)
  }
)

#' makeMatrix
#' 
#' Create a standard expression matrix from a given set of files.
#'
#' @param gef A \code{data.table}. The gene_expression_files dataset, filtered
#'  for the selected cohort and relvant files.
#' 
#' @name makeMatrix 
#' @importFrom tools file_ext
#' 
.ISCon$methods(
  makeMatrix = function(gef){
    inputFiles <- unique(gef$file_info_name)
    cohort <- unique(gef$arm_name)
    ext <- unique(file_ext(inputFiles))
    if(length(cohort) > 1){
      message("There are more than one cohort selected in this HIPCMatrix run")
    }
    
    inputFiles <- file.path("/share/files/", config$labkey.url.path, "@files/rawdata/gene_expression", inputFiles)
    # Filetypes
    # After this step norm_exprs is a matrix with features as rownames and expsample as colnames
    if(length(ext) > 1){
      stop(paste("There is more than one file extension:", paste(ext, collapse = ",")))
    } else if(ext == "CEL"){
      norm_exprs <- .process_CEL(con, gef, inputFiles)
    } else if(ext == "txt" | study %in% c("SDY162", "SDY180", "SDY212")){
      norm_exprs <- .process_others(gef, inputFiles)
    } else if(ext == "tsv"){
      norm_exprs <- .process_TSV(gef, inputFiles)
    } else{
      stop("File extension not supported.")
    }
    if(!is(norm_exprs, "data.table")){
      norm_exprs <- data.table(norm_exprs, keep.rownames = TRUE)
      setnames(norm_exprs, "rn", "feature_id")
    }
    # This step should eventually be removed as we move from biosample to expsample
    norm_exprs <- .es2bs(.self, norm_exprs)
    return(norm_exprs)
  }
)

# Process CEL files
# @param gef A \code{data.table} the gene_expression_files table or a subset of
#  it.
# @param inputFiles A \code{character}. The filenames.
# @return A \code{matrix} with biosample_accession as cols and feature_id as rownames
#' @importFrom affy ReadAffy rma
.process_CEL <- function(con, gef, inputFiles){
  affybatch <- ReadAffy(filenames = inputFiles)                                 
  eset <- rma(affybatch)                                                        
  norm_exprs <- exprs(eset)                                                     
  if (all(file_ext(colnames(norm_exprs)) == "CEL")) {#If filenames used as samplenames
    colnames(norm_exprs) <- gef[match(colnames(norm_exprs), gef$file_info_name), biosample_accession]
  }
  return(norm_exprs)
}

# This will work for files that follow the standards from immport
# Eventually, all tsv files should be rewritten to follow this standard.
# @return A matrix with biosample_accession as cols and feature_id as rownames
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom reshape2 acast
.process_TSV <- function(gef, inputFiles){
  exprs <- fread(inputFiles, header = TRUE)
  exprs <- .clean_colnames(exprs)
  
  if(!all(c("target_id", "raw_signal") %in% colnames(exprs))){
    stop("The file does not follow HIPC standards.")
  }
  try(setnames(exprs, "experiment_sample_accession", "expsample_accession"))
  if("expsample_accession" %in% colnames(exprs)){
    EorB <- "expsample_accession"
  } else if("biosample_accession" %in% colnames(exprs)){
    EorB <- "biosample_accession"
  } else{
    stop("The input file must contain either biosample_accession or expsample_accession")
  }
  exprs <- acast(exprs, formula = paste("target_id ~", EorB), value.var = "raw_signal")
  cnames <- colnames(exprs)
  rnames <- rownames(exprs)
  exprs <- preprocessCore::normalize.quantiles(exprs)
  colnames(exprs) <- cnames 
  rownames(exprs) <- rnames 
  
  norm_exprs <- log2(pmax(exprs, 1))
  norm_exprs <- norm_exprs[, c(colnames(norm_exprs) %in% gef$expsample_accession)]
  return(norm_exprs)
}

#biosample_accession as colnames
#Works for SDY212 & 162
#' @importFrom preprocessCore normalize.quantiles
.process_others <- function(gef, inputFiles){
  exprs <- fread(inputFiles, header = TRUE)
  sigcols <- grep("Signal", colnames(exprs), value = TRUE)
  rnames <- exprs[, PROBE_ID]
  if(length(sigcols) > 0){
    exprs <- exprs[, sigcols, with = FALSE]
    #exprs <- exprs[, c("PROBE_ID", sigcols), with = FALSE]
    setnames(exprs, gsub(".AVG.*$", "", colnames(exprs)))
    #try(setnames(exprs, "PROBE_ID", "feature_id"))
    #setnames(exprs, c("PROBE_ID", colnames(exprs)),
    #         c("feature_id", gsub(".AVG.*$", "", colnames(exprs))))
  } else if(all(gef$biosample_accession %in% colnames(exprs))){
    exprs <- exprs[, gef$biosample_accession, with = FALSE]
  } else{
    stop("Unknown format: check data and add code if needed.")
  }
  
  cnames <- colnames(exprs)
  #rnames <- rownames(exprs)
  exprs <- as.matrix(exprs)
  exprs <- preprocessCore::normalize.quantiles(exprs)
  colnames(exprs) <- cnames
  rownames(exprs) <- rnames
  
  norm_exprs <- log2(pmax(exprs, 1))
  norm_exprs <- norm_exprs[, c(colnames(norm_exprs) %in% gef$biosample_accession)]
  return(norm_exprs)
}
  


.clean_colnames <- function(table){
  setnames(table, tolower(chartr(" ", "_", names(table))))
}

.es2bs <- function(con, table){
  ess <- grep("^ES", colnames(table), value = TRUE)
  if(length(ess) > 0){
    esFilter <- makeFilter(c("expsample_accession", "IN", paste0(ess, collapse = ";")))
    bs2es <- data.table(labkey.selectRows(con$config$labkey.url.base,
                                          con$config$labkey.url.path,
                                          "immport", "biosample_2_expsample",
                                          colFilter = esFilter,
                                          colNameOpt = "rname"))
    bss <- bs2es[match(ess, bs2es$expsample_accession), biosample_accession]
    setnames(table, ess, bss)
  }
  return(table)
}
