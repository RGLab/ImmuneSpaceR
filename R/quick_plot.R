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
    ggthemr("solarized")
    addPar <- c("Gender", "Age", "Race")
    annoCols <- c("arm_name", "subject_accession", "study_time_collected", addPar)
    toKeep <- c("response", "analyte", annoCols)
    logT <- TRUE #By default, log transform the value_reported
    message_out <- ""
    extras <- list(...)
    
    e <- try({
      if(dataset != "gene_expression_analysis_results"){
        dt <- copy(con$getDataset(dataset, colFilter = filter, reload = TRUE))
        if(!"analyte" %in% colnames(dt)){
          if("analyte_name" %in% colnames(dt)){
            dt <- dt[, analyte := analyte_name]
          } else{
            dt <- dt[, analyte := ""]
          }
        }
      }
      
      # Datasets
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
        }
      dt <- dt[, response := ifelse(value_reported <0, 0, value_reported)]
      if(logT){
        dt <- dt[, response := mean(log2(response+1), na.rm = TRUE),
                 by = "arm_name,subject_accession,analyte,study_time_collected"]
      } else{
        dt <- dt[, response := mean(response, na.rm = TRUE),
                 by = "arm_name,subject_accession,analyte,study_time_collected"]
      }
      setnames(dt, c("gender", "age_reported", "race"), addPar)
      setkeyv(dt, toKeep)
      dt <- unique(dt[, toKeep, with = FALSE])
      
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
      p <- .qpHeatmap(dt, normalize_to_baseline, legend, text_size)
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



#' @importFrom pheatmap pheatmap
#' @importFrom reshape2 acast
.qpHeatmap = function(dt, normalize_to_baseline, legend, text_size){
  contrast <- "study_time_collected"
  annoCols <- c("arm_name", "subject_accession", contrast, "Gender", "Age", "Race")
  palette <- ISpalette(20)
  
  expr <- parse(text = paste0(contrast, ":=as.factor(", contrast, ")"))
  dt <- dt[, eval(expr)]
  #No need to order by legend. This should be done after.
  if(!is.null(legend)){
    dt <- dt[order(arm_name, study_time_collected, get(legend))]
  } else{
    dt <- dt[order(arm_name, study_time_collected)]
  }
  form <- as.formula(paste("analyte ~ arm_name +", contrast, "+ subject_accession"))
  mat <- acast(data = dt, formula = form, value.var = "response") #drop = FALSE yields NAs
  if(ncol(mat) > 2 & nrow(mat) > 1){
    mat <- mat[rowSums(apply(mat, 2, is.na)) < ncol(mat),, drop = FALSE]
  }
  
  # Annotations:
  anno <- data.frame(unique(dt[, annoCols, with = FALSE]))
  rownames(anno) <- paste(anno$arm_name, anno[, contrast], anno$subject_accession, sep = "_")
  expr <- parse(text = c(rev(legend), contrast, "arm_name"))
  anno <- anno[with(anno, order(eval(expr))),]
  anno <- anno[, c(rev(legend), contrast, "arm_name")] #Select and order the annotation rows
  anno[, contrast] <- as.factor(anno[, contrast])
  anno_color <- colorpanel(n = length(levels(anno[,contrast])), low = "white", high = "black")
  names(anno_color) <- levels(anno[, contrast])
  anno_color <- list(anno_color)
  if(contrast == "study_time_collected"){
    setnames(anno, c("arm_name", contrast), c("Arm Name", "Time"))
    contrast <- "Time"
  }
  names(anno_color) <- contrast
  if("Age" %in% legend){
    anno_color$Age <- c("yellow", "red")
  }
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
  if(type == "violin"){
    geom_type <- geom_violin() #+ stat_summary(fun.y="median", geom="point")
  } else{
    geom_type <- geom_boxplot(outlier.size = 0)
  }
  print(head(dt))
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