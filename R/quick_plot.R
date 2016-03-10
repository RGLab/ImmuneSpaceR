#' @include ImmuneSpace.R
NULL

.ISCon$methods(
  quick_plot = function(...){
    "Plots a selected dataset. This is the function used by the DataExplorer 
    module on ImmuneSpace.\n
    dataset: A character. The name of the dataset to plot, as displayed by the
    listDataset method.\n
    normalize_to_baseline: A logical. If set to TRUE, the values are plotted as
    log2 fold-change from baseline.\n
    type: A character. The type of plot. Valid choices are 'auto', 'heatmap',
    'boxplot', 'lineplot', 'violinplot'. If set to 'auto', the function will
    select an appropriate plot type for the selected data.\n
    filter: A filter as created by the makeFilter function from Rlabkey.\n
    facet: The facetting for ggplot2 based plots. Valid choices are 'grid' and 
    'wrap'.\n
    text_size: The size of all text elements in the plot.\n
    legend: A character. Columns of the dataset or demographics to be added as
    legend on the heatmap. This argument is ignored if the plot type isn't 
    heatmap.\n
    show_virus_strain: A logical. Should all the virus strains be shown or should
    the values be averaged. Only used when dataset = 'hai'.\n
    '...': Extra argument to be passed to ggplot. e.g: shape = 'Age', color =
    'Race'.\n"
    .quick_plot(.self, ...)
  }
)

#' @importFrom ggplot2 facet_grid facet_wrap geom_text element_blank
#' @importFrom Biobase pData
.quick_plot <- function(con, dataset, normalize_to_baseline = TRUE,
                      type = "auto", filter = NULL,
                      facet = "grid", text_size = 15,
                      legend = NULL, show_virus_strain = FALSE, ...){
  logT <- TRUE #By default, log transform the value_reported
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
    if(tolower(dataset) %in% c("ge", "dgea_filteredgear", "gene_expression", "gene_expression_analysis_results")){
      dataset <- "gene_expression"
    }
    dt <- .getDataToPlot(con, dataset, filter = filter, show_virus_strain)
    dt <- .standardize_time(dt)
    ylab <- .format_lab(dataset, normalize_to_baseline)
    if(logT){
      dt <- dt[, response := mean(log2(response+1), na.rm = TRUE),
               by = "cohort,participant_id,analyte,time_str"]
    } else{
      dt <- dt[, response := mean(response, na.rm = TRUE),
               by = "cohort,participant_id,analyte,time_str"]
    }
    dt <- unique(dt)
    if(normalize_to_baseline){
      dt <- dt[,response:=response-response[study_time_collected <= 0],
               by="cohort,participant_id,analyte"][study_time_collected > 0]
      if(nrow(dt) == 0){
        stop("All data points are <= 0. Cannot normalize to baseline.")
      }
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
    facet <- facet_grid(aes(analyte, cohort), scales = "free")
  } else if(facet == "wrap"){
    facet <- facet_wrap(~cohort + analyte, scales = "free")
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
.standardize_time <- function(data){
  data <- data[, stcu := gsub("s$", "", tolower(data$study_time_collected_unit))]
  if(length(unique(data$stcu)) > 1){
    if(!all(unique(data$stcu) %in% c("hour", "day"))){
      stop("Time should be expressed in Days or Hours")
    }
    data <- data[, stc := study_time_collected]
    data <- data[stcu == "day", stc := stc * 24]
    # Merge equivalent TP
    data <- data[, stcu := ifelse(abs(stc) < 24, "Hour", "Day")]
    # Get levels 
    ut <- sort(unique(data$stc))
    levs <- ifelse(abs(ut) < 24, paste("Hour", ut), paste("Day", ut/24))
    # Concatenate time and unit
    data <- data[, stc := ifelse(abs(stc) < 24, stc, stc/24)]
    data <- data[, time_str := factor(paste(stcu, stc), levels = levs)]
    # Cleanup
    data <- data[, c("study_time_collected", "study_time_collected_unit") := 
                   list(stc, stcu)]
    data <- data[, c("stc", "stcu") := NULL]
  } else{
    ut <- sort(unique(data$study_time_collected))
    levs <- paste(unique(data$study_time_collected_unit), ut)
    data <- data[, time_str := factor(paste(study_time_collected_unit, study_time_collected), levels = levs)]
  }
  return(data)
}


#' @importFrom pheatmap pheatmap
#' @importFrom reshape2 acast
.qpHeatmap2 <- function(dt, normalize_to_baseline, legend, text_size){
  palette <- ISpalette(20)
  
  dt <- dt[, ID := paste(cohort, time_str, participant_id, sep = "_")]
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
# @importFrom ggplot2 scale_color_manual scale_color_gradient
.qpBoxplotViolin <- function(dt, type, facet, ylab, text_size, extras, ...){
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
  p <- p + theme_IS()
  print(p)
}

#' @importFrom ggplot2 ggplot geom_line geom_point
#' @importFrom ggplot2 theme element_text aes_string aes xlab ylab
.qpLineplot <- function(dt, facet, ylab, text_size, extras, ...){
  p <- ggplot(data = dt, aes(study_time_collected, response, group = participant_id)) +
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
.getDataToPlot <- function(con, dataset, filter = NULL, show_virus_strain = FALSE){
  # All columns that can potentially be used
  demo_cols <- c("gender", "age_reported", "race")
  out_cols <- c("study_time_collected", "study_time_collected_unit", "cohort", "participant_id")
  out_cols <- c(c("response", "analyte"), demo_cols, out_cols)
  if(dataset != "gene_expression"){
    dt <- copy(con$getDataset(dataset, colFilter = filter, reload = TRUE))
    if(!"analyte" %in% colnames(dt)){
      dt <- dt[, analyte := ""]
    }
  }
  
  if(dataset == "elispot"){
    dt <- dt[, value_reported := (spot_number_reported) / cell_number_reported]
  } else if(dataset %in% c("hai", "neut_ab_titer")){
    if(isTRUE(show_virus_strain)){
      dt <- dt[, analyte := virus]
    }
  } else if(dataset == "pcr"){
    if(all(is.na(dt[, threshold_cycles]))){
      stop("PCR results cannot be displayed for studies that do not use threshold cycles.")
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
  } else if(dataset == "gene_expression"){
    logT <- FALSE #Matrices are already log2 transformed
    dt <- copy(con$getGEAnalysis(colFilter = filter))
    if(!is.null(filter) & any(sapply(filter, function(x){gsub("~.*$", "", x)}) == "cohort")){
      uarm <- unique(dt$cohort)
    } else{
    uarm <- labkey.selectRows(
      baseUrl=con$config$labkey.url.base, folderPath=con$config$labkey.url.path,
      schemaName="assay.ExpressionMatrix.matrix",
      queryName="SelectedRuns",
      viewName="expression_matrices",
      colFilter=NULL,
      containerFilter= "CurrentAndSubfolders")$Cohort
    }
    ugenes <- unique(dt$gene_symbol)
    ugenes <- ugenes[ ugenes != "NA"]
    EM <- con$getGEMatrix(cohort = uarm, summary = TRUE)
    ugenes <- ugenes[ugenes %in% featureNames(EM)]
    EM <- EM[ugenes,]
    pd <- data.table(pData(EM))
    demo <- con$getDataset("demographics")
    demo <- unique(demo[, c("participant_id", demo_cols), with = FALSE])
    dt <- data.table(melt(exprs(EM)))
    setnames(dt, c("analyte", "biosample_accession", "value_reported"))
    dt <- merge(dt, pd, by = "biosample_accession", all.x = TRUE) # Add s_t_c, s_t_c_u, arm 
    dt <- merge(dt, demo, by = "participant_id", all.x = TRUE) # Add race, gender, age
    setkey(dt, NULL)
  }
  dt <- dt[, response := ifelse(value_reported <0, 0, value_reported)]
  dt <- dt[, out_cols, with = FALSE]
  setnames(dt, demo_cols, c("Gender", "Age", "Race"))
  return(dt)
}

# dt has ID and all relevant columns
.heatmapAnnotations <- function(dt, legend){
  annoCols <- c("cohort", "time_str", legend)
  # Annotations  
  anno <- data.table(unique(dt[, c("ID", annoCols), with = FALSE]))
  # Order: Arm > Time > Age > Gender > Race
  order_cols <- c("cohort", "time_str", legend)
  setorderv(anno, order_cols)
  setcolorder(anno, c("ID", rev(legend), "time_str", "cohort"))
  # Set colors
  setnames(anno, c("cohort", "time_str"), c("Cohort", "Time"))
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

# Display name for the Y-axis
# TODO: Add units. But they are in dt.
.format_lab <- function(dataset, normalize_to_baseline){
  lab <- switch(dataset,
                "hai" = "HAI",
                "elisa" = "Concentration",
                "elispot" =  "Spot count",
                "mbaa" = "Concentration",
                "fcs_analyzed_result" = "Cell number",
                "gene_expression" = "")
  if(normalize_to_baseline){
    lab <- paste(lab, "normalized to baseline")
  }
  return(lab)       
}
