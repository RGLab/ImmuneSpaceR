#' @include ImmuneSpace.R
NULL

.ISCon$methods(
  quick_plot = function(...){
    .quick_plot(.self, ...)
  }
)

# @importFrom ggthemr ggthemr
#Currently we have to depend on ggthemr because it depends on ggplot2
#' @import ggthemr 
#' @importFrom ggplot2 facet_grid facet_wrap geom_text element_blank
#' @importFrom Biobase pData
#' @importFrom reshape2 melt
.quick_plot <- function(con, dataset, normalize_to_baseline = TRUE,
                      type = "auto", filter = NULL,
                      facet = "grid", text_size = 15,
                      legend = NULL, ...){
  logT <- TRUE #By default, log transform the value_reported
  message_out <- ""
  extras <- list(...)
  
  # legend
  if(!is.null(legend)){
    legend <- unique(paste0(toupper(substring(legend, 1, 1)),
                      substring(tolower(gsub("_.*$", "", legend)), 2)
                    )
    )
  }
  
  # Datasets
  e <- try({
    dt <- .getDataToPlot(con, dataset, filter = filter)
    dt <- standardize_time(dt, keep_units = TRUE)
    if(logT){
      dt <- dt[, response := mean(log2(response+1), na.rm = TRUE),
               by = "arm_name,subject_accession,analyte,time_str"]
    } else{
      dt <- dt[, response := mean(response, na.rm = TRUE),
               by = "arm_name,subject_accession,analyte,time_str"]
    }
    dt <- unique(dt)
    
    if(normalize_to_baseline){
      dt <- dt[,response:=response-response[study_time_collected==0],
               by="arm_name,subject_accession,analyte"][study_time_collected!=0]
      ylab <- "Response normalized to baseline"
    } else{
      ylab <- "Response (log2)"
    }
    if(type == "auto"){
      if(length(unique(dt$analyte)) < 10){
        type <- "boxplot"
      } else{
        type <- "heatmap"
      }
    }
  })
  
  if(inherits(e, "try-error")){
    type <- "error"
    error_string <- attr(e, "condition")$message
  }
  
  # Plot
  if(facet == "grid"){
    facet <- facet_grid(aes(analyte, arm_name), scales = "free")
  } else if(facet == "wrap"){
    facet <- facet_wrap(~arm_name + analyte, scales = "free")
  }
  if(type == "heatmap"){
    p <- .qpHeatmap2(dt, normalize_to_baseline, legend, text_size)
  } else if(type %in% c("boxplot", "violin")){
    .qpBoxplotViolin(dt, type, facet, ylab, text_size, extras, ...)
  } else if(type == "line"){
    .qpLineplot(dt, facet, ylab, text_size, extras, ...)
  } else{#} if(type == "error"){
    data <- data.frame(x = 0, y = 0, err = error_string)
    p <- ggplot(data = data) + geom_text(aes(x, y, label = err), size = 10) +
      theme(line = element_blank(), text = element_blank())
    print(p)
  }
}

# Select appropriate unit and merge equivalent timepoints
# Order timepoints by time (not numeric or alphabetic)
# @param data A \code{data.table} containing a study_time_collected and a 
#  study_time_collected_unit columns
# @value Returns a data.table with an additional time_str column of ordered
#  factors
standardize_time <- function(data, keep_units = TRUE){
  data <- data[, study_time_collected_unit := gsub("s$", "", tolower(data$study_time_collected_unit))]
  if(!all(unique(data$study_time_collected_unit) %in% c("hour", "day"))){
    stop("Time should be expressed in Days or Hours")
  }
  if(keep_units){
    data <- data[, stcu := gsub("s$", "", tolower(data$study_time_collected_unit))]
    data <- data[, stc := study_time_collected]
    data <- data[stcu == "day", stc := stc * 24]
    # Merge equivalent TP
    data <- data[, stcu := ifelse(abs(stc) < 24, "hour", "day")]
    # Get levels 
    ut <- sort(unique(data$stc))
    levs <- ifelse(abs(ut) < 24, paste(ut, "hour"), paste(ut/24, "day"))
    # Concatenate time and unit
    data <- data[, stc := ifelse(abs(stc) < 24, stc, stc/24)]
    data <- data[, time_str := factor(paste(stc, stcu), levels = levs)]
    # Cleanup
    data <- data[, c("study_time_collected", "study_time_collected_unit") := 
                   list(stc, stcu)]
    data <- data[, c("stc", "stcu") := NULL]
  } else{
    data <- data[study_time_collected_unit == "day", study_time_collected := study_time_collected * 24]
    if(max(abs(data$study_time_collected)) > 24){
      data[, study_time_collected_unit := "Days"]
      data[, study_time_collected := study_time_collected/24]
    } else{
      data[, study_time_collected_unit := "Hours"]
    }
  }
  return(data)
}


#' @importFrom pheatmap pheatmap
#' @importFrom reshape2 acast
.qpHeatmap2 <- function(dt, normalize_to_baseline, legend, text_size){
  palette <- ISpalette(20)
  
  dt <- dt[, ID := paste(arm_name, time_str, subject_accession, sep = "_")]
  mat <- acast(data = dt, formula = formula("analyte ~ ID"), value.var = "response")
  
  if(ncol(mat) > 2 & nrow(mat) > 1){
    mat <- mat[rowSums(apply(mat, 2, is.na)) < ncol(mat),, drop = FALSE]
  }
  
  annos <- .heatmapAnnotations(dt, legend)
  anno <- annos[[1]]
  anno_color <- annos[[2]]
  mat <- mat[, rownames(anno), drop = FALSE]
  
  # pheatmap parameters
  if(normalize_to_baseline){
    scale <- "none"                                                 
    max <- max(abs(mat), na.rm = TRUE)
    breaks <- seq(-max, max, length.out = length(palette))
  } else{                                                           
    scale <- "row"                                                  
    breaks <- NA
  }
  
  show_rnames <- ifelse(nrow(mat) < 50, TRUE, FALSE)
  cluster_rows <- ifelse(nrow(mat) > 2 & ncol(mat) > 2, TRUE, FALSE)
  
  e <- try({
    p <- pheatmap(mat = mat, annotation = anno, show_colnames = FALSE,
                  show_rownames = show_rnames, cluster_cols = FALSE,
                  cluster_rows = cluster_rows, color = palette,
                  scale = scale, breaks = breaks,
                  fontsize = text_size, annotation_colors = anno_color)
  })
  if(inherits(e, "try-error")){
    p <- pheatmap(mat = mat, annotation = anno, show_colnames = FALSE,
                  show_rownames = show_rnames, cluster_cols = FALSE,
                  cluster_rows = FALSE, color = palette,
                  scale = scale, breaks = breaks,
                  fontsize = text_size, annotation_colors = anno_color)
  }
  return(p)
}


#' @importFrom ggplot2 ggplot geom_violin geom_boxplot geom_jitter
#' @importFrom ggplot2 theme element_text aes_string aes xlab ylab
.qpBoxplotViolin <- function(dt, type, facet, ylab, text_size, extras, ...){
  ggthemr("solarized")
  if(type == "violin"){
    geom_type <- geom_violin() #+ stat_summary(fun.y="median", geom="point")
  } else{
    geom_type <- geom_boxplot(outlier.size = 0)
  }
  p <- ggplot(data = dt, aes(as.factor(study_time_collected), response)) +
  geom_type + xlab("Time") + ylab(ylab) + facet + 
  theme(text = element_text(size = text_size), axis.text.x = element_text(angle = 45))
  if(!is.null(extras[["size"]])){                                           
    p <- p + geom_jitter(aes_string(...))                                   
  } else{                                                                   
    p <- p + geom_jitter(size = 3, aes_string(...))                         
  }                                                                         
  print(p)
}

#' @importFrom ggplot2 ggplot geom_line geom_point
#' @importFrom ggplot2 theme element_text aes_string aes xlab ylab
.qpLineplot <- function(dt, facet, ylab, text_size, extras, ...){
  ggthemr("solarized")
  p <- ggplot(data = dt, aes(study_time_collected, response, group = subject_accession)) +
  geom_line(size = 1, aes_string(...)) +                                            
  xlab("Time") + ylab(ylab) + facet + 
  theme(text = element_text(size = text_size), axis.text.x = element_text(angle = 45))
  if(!is.null(extras[["size"]])){                                           
    p <- p + geom_point(aes_string(...))                                    
  } else{                                                                   
    p <- p + geom_point(size = 3, aes_string(...))                          
  }                                                                         
  print(p)                                                                  
}

# Get the data
# Add standard columns for analyte and response
# NOTE: No need for analyte_name check (Fixed in DR12)
.getDataToPlot <- function(con, dataset, filter = NULL){
  # All columns that can potentially be used
  demo_cols <- c("gender", "age_reported", "race")
  out_cols <- c("study_time_collected", "study_time_collected_unit", "arm_name", "subject_accession")
  out_cols <- c(c("response", "analyte"), demo_cols, out_cols)
  if(dataset != "gene_expression_analysis_results"){
    dt <- copy(con$getDataset(dataset, colFilter = filter, reload = TRUE))
    if(!"analyte" %in% colnames(dt)){
      dt <- dt[, analyte := ""]
    }
  }
  
  if(dataset == "elispot"){
    dt <- dt[, value_reported := (spot_number_reported) / cell_number_reported]
  } else if(dataset == "pcr"){
    if(all(is.na(dt[, threshold_cycles]))){
      stop("PCR results cannot be displayed for studies that do not use threshold cycles.
           Use LabKey Quick Chart interface to plot this dataset.")
    }
    dt <- dt[, value_reported := threshold_cycles]
    dt <- dt[, analyte := entrez_gene_id]
    logT <- FALSE #Threshold cycle is already log transformed
  } else if(dataset == "mbaa"){
    if(all(dt$concentration_value ==0) || all(is.na(dt$concentration_value))){
      if(any(!is.na(dt$mfi)) && any(dt$mfi != 0)){
        dt <- dt[, value_reported := as.numeric(mfi)]
      }else{
        stop("Plotting MBAA requires either concentration or MFI values")
      }
    } else{
      dt <- dt[, value_reported := as.numeric(concentration_value)]
    }
  } else if(dataset == "fcs_analyzed_result"){
    dt <- dt[, value_reported := as.numeric(population_cell_number)]
    dt <- dt[, analyte := population_name_reported]
  } else if(dataset == "gene_expression_analysis_results"){
    logT <- FALSE #Matrices are already log2 transformed
    dt <- copy(con$getGEAnalysis(colFilter = filter))
    uarm <- unique(dt$arm_name)
    ugenes <- unique(dt$gene_symbol)
    ugenes <- ugenes[ ugenes != "NA"]
    EM <- con$getGEMatrix(cohort = uarm, summary = TRUE)
    EM <- EM[ugenes,]
    pd <- data.table(pData(EM))
    demo <- con$getDataset("demographics")
    dt <- data.table(melt(exprs(EM)))
    setnames(dt, c("analyte", "biosample_accession", "value_reported"))
    dt <- merge(dt, pd, by = "biosample_accession", all.x = TRUE) # Add s_t_c, s_t_c_u, arm 
    dt <- merge(dt, demo, by = "subject_accession", all.x = TRUE) # Add race, gender, age
    setkey(dt, NULL)
  }
  dt <- dt[, response := ifelse(value_reported <0, 0, value_reported)]
  dt <- dt[, out_cols, with = FALSE]
  setnames(dt, demo_cols, c("Gender", "Age", "Race"))
  if(nrow(dt) != nrow(unique(dt))){
    print(paste("There are: ", nrow(dt) - nrow(unique(dt)), " duplicates."))
  }
  return(dt)
}

# dt has ID and all relevant columns
.heatmapAnnotations <- function(dt, legend){
  annoCols <- c("arm_name", "time_str", legend)
  # Annotations  
  anno <- data.table(unique(dt[, c("ID", annoCols), with = FALSE]))
  # Order: Arm > Time > Age > Gender > Race
  order_cols <- c("arm_name", "time_str", legend)
  setorderv(anno, order_cols)
  setcolorder(anno, c("ID", rev(legend), "time_str", "arm_name"))
  # Set colors
  setnames(anno, c("arm_name", "time_str"), c("Arm", "Time"))
  anno_color <- list(Time = colorpanel(n = length(levels(anno$Time)),
                                       low = "white", high = "black"))
  names(anno_color$Time) <- levels(anno$Time)
  if("Age" %in% legend){
    anno_color$Age <- c("yellow", "red")
  }
  # data.frame for pheatmap call
  anno <- data.frame(anno, row.names = anno$ID)
  anno$ID <- NULL
  return(list(anno, anno_color))
}