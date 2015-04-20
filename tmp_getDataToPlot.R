# Get the data
# Add standard columns for analyte and response
# NOTE: No need for analyte_name check (Fixed in DR12)
.getDataToPlot <- function(con, dataset, filter = NULL){
  # All columns that can potentially be used
  demo <- c("gender", "age_reported", "race")
  out_cols <- c("study_time_collected", "study_time_collected_unit", "arm_name", "subject_accession")
  out_cols <- c(c("response", "analyte"), demo, out_cols)
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
  if(nrow(dt) != nrow(unique(dt))){
    print(paste("There are: ", nrow(dt) - nrow(unique(dt)), " duplicates."))
  }
  return(dt)
}

test_dat <- function(dt){
  demo <- c("age_reported", "race", "gender")
  others <- c("study_time_collected", "study_time_collected_unit", "arm_name", "subject_accession")
  expected_cn <- c(c("value_reported", "analyte"), demo, others)
  if(nrow(dt) == 0){
    warning("Empty dataset!")
  }
  if(!all(expected_cn %in% colnames(dt))){
    print(paste0(
      "Missing column(s):", expected_cn[!expected_cn %in% colnames(dt)]
    ))
  }
}

# dt has ID and all relevant columns
.createAnnotations <- function(dt, legend){
  annoCols <- c("arm_name", "time_str", legend)
  # Annotations  
  anno <- data.table(unique(dt[, c("ID", annoCols), with = FALSE]))
  # Order: Arm > Time > Age > Gender > Race
  order_cols <- c("arm_name", "time_str", legend)
  setorderv(anno, order_cols)
  setcolorder(anno, c("ID", rev(legend), "time_str", "arm_name"))
  # Set colors
  setnames(anno, "time_str", "Time")
  anno_color <- list(Time = colorpanel(n = length(levels(anno$Time)),
                                       low = "white", high = "black"))
  names(anno_color$Time) <- levels(anno$Time)
  if("age_reported" %in% legend){
    setnames(anno, "age_reported", "Age")
    anno_color$Age <- c("yellow", "red")
  }
  # data.frame for pheatmap call
  anno <- data.frame(anno, row.names = anno$ID)
  anno$ID <- NULL
  return(list(anno, anno_color))
}