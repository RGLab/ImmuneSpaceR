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
    
    # Filetypes
    if(length(ext) > 1){
      stop(paste("There is more than one file extension:", paste(ext, collapse = ",")))
    } else if(ext == "CEL"){
      process_CEL(gef, inputFiles)
    } else if(ext == "tsv"){
      #process_TSV(gef, inputFiles)
    }
    return(cohort)
  }
)
    

# Process CEL files
# @param gef A \code{data.table} the gene_expression_files table or a subset of
#  it.
# @param inputFiles A \code{character}. The filenames.
# @return A \code{matrix} with biosample_accession as cols and feature_id as rownames
process_CEL <- function(gef, inputFiles){
  library(affy)
  affybatch <- ReadAffy(filenames = inputFiles)                                 
  eset <- rma(affybatch)                                                        
  norm_exprs <- exprs(eset)                                                     
  if (all(file_ext(colnames(norm_exprs)) == "CEL")) {#If filenames used as samplenames
    colnames(norm_exprs) <- gef[match(colnames(norm_exprs), gef$file_info_name), biosample_accession]
  }
}














# 
# input_csv <- "~/Dropbox (Gottardo Lab)/ImmPort/Microarray_values_noPGP.csv"
# sdy <- "SDY296"
# 
# library(Rlabkey)
# library(ImmuneSpaceR)
# con <- CreateConnection(sdy)
# data <- fread(input_csv)
# setnames(data, tolower(gsub(" ", "_", colnames(data))))
# 
# pdata <- con$getDataset("gene_expression_files", original_view = TRUE)
# 
# 
# es <- data.table(labkey.selectRows("test.immunespace.org", paste0("Studies/", sdy), "immport",
#                              "expsample", colNameOpt="rname"))
# es2f <- data.table(labkey.selectRows("test.immunespace.org", paste0("Studies/", sdy), "immport",
#                                    "expsample_2_file_info", colNameOpt="rname"))
# 
# #Merge expsample_accession w/ data
# #acast data and put it in exprs of an eSet
# #dcast.data.table(data, formula="target_id ~ experiment_sample_user-defined_id", value.var = "raw_signal")
