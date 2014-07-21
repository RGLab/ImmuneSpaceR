#' Quick plot of a data set from a study
#' 
#' Generates a quick plot of a data set automatically. Based on ggplot2's qplot.
#' @param dt \code{data.table} the data to be plotted.
#' @param normalize_to_baseline \code{logical} defaults to \code{TRUE}.
#' @param type \code{character} one of \code{c("auto","heatmap","boxplot")}. Heatmap chosen automatically if there are more than 5 analytes.
#' @param ... additional arguments passed to qplot.
#' @importFrom ggplot2 qplot
#' @return \code{ggplot2} structure.
#' @usage quick_plot(dt, normalize_to_baseline = TRUE, type="auto", ...)
#' @export
#' @examples
#' \dontrun{ 
#' study <- CreateConnection("SDY269")
#' dt <- study$getDataset("elisa_mbaa")
#' quick_plot(dt, normalize_to_baseline = F, type="auto")
#' quick_plot(dt, normalize_to_baseline = F, type="boxplot")
#' }
quick_plot <- function(dt, normalize_to_baseline=TRUE, type="auto", ...)
{
  # Add a dummy analyte for consistency
  if(length(grep("analyte",colnames(dt)))==0)
    dt <- dt[,analyte:=""]
  # Compute summaries over all repeated measures (e.g. multiple virus strains)
  dt_unique <- dt[,list(response=mean(log2(value_reported), na.rm=TRUE)),by="arm_accession,subject_accession,analyte,study_time_collected"]
  
  if(type=="auto" & length(unique(dt_unique$analyte))>5)
    type <- "heatmap"
  else if(type=="auto")
    type <- "boxplot"
  
  # Compute fold changes
  if(normalize_to_baseline==TRUE)
  {
    # Remove the time zero response
    dt_unique <- dt_unique[,response:=response-response[study_time_collected==0],by="arm_accession,subject_accession,analyte"][study_time_collected!=0]
  }
  
  # TODO: change arm_accession with cohort, once it's there
  if(type=="boxplot")
  {
    p <- qplot(as.factor(study_time_collected), response, data=dt_unique, facets=analyte~arm_accession,geom=c("boxplot","jitter"), ..., xlab = "time", ylab="response (log2)")
  }
  else if(type=="line")
  {
    p <- qplot(as.factor(study_time_collected), response, data=dt_unique, facets=analyte~arm_accession,geom=c("line","point"), ..., xlab = "time", ylab="response (log2)", group=subject_accession)  
  }
  else if(type=="heatmap")
  {
    p <- qplot(as.factor(study_time_collected), analyte, data=dt_unique, facets=~arm_accession,geom=c("raster"), ..., xlab = "time", fill=response) + scale_fill_gradient2(high = "#a50026", mid="#ffffbf", low="#313695")
  }
  p+theme_minimal()
}

