# library(ggplot2)
# library(pheatmap)
# library(reshape2)
# library(Rlabkey)
# palette <- rev(brewer.pal(n = 11, name = "RdYlBu"))
# addPar <- c("gender", "age_reported", "race")
# annoCols <- c("name", "subject_accession", "study_time_collected", addPar)
# toKeep <- c("response", "analyte", annoCols)
# 
# e <- try({
# dt <- con$getDataset(dataset, reload = TRUE)
# if(length(grep("analyte",colnames(dt)))==0){
#   dt <- dt[, analyte := ""]
# }
# if(dataset == "elispot"){
#   dt <- dt[, value_reported := (spot_number_reported +1) / cell_number_reported]
# } else if(dataset == "pcr"){
#   if(all(is.na(dt[, threshold_cycle]))){
#     stop("PCR results cannot be displayed for studies that do not use threshold cycles")
#   }
#   dt <- dt[, analyte := entrez_gene_id]
# }
# dt <- dt[, response := mean(log2(value_reported), na.rm = TRUE),
#          by = "name,subject_accession,analyte,study_time_collected"]
# dt <- unique(dt[, toKeep, with = FALSE])
# 
# if(normalize_to_baseline){
#   dt <- dt[,response:=response-response[study_time_collected==0],
#            by="name,subject_accession,analyte"][study_time_collected!=0]
#   ylab <- "Response normalized to baseline"
# } else{
#   ylab <- "Response (log2)"
# }
# })
# 
# if(inherits(e, "try-eror")){
#   type <- "error"
#   error_string <- attr(e, "condition")$message
# }
# 
# #dt: analyte, subject, timepoint, response, arm + addPar
# if(type == "heatmap"){
#   mat <- acast(dt, analyte ~ name + study_time_collected + subject_accession, value.var = "response")
#   anno <- data.frame(unique(dt[, annoCols, with = FALSE]))
#   rownames(anno) <- paste(anno$name, anno$study_time_collected, anno$subject_accession, sep = "_")
#   anno <- anno[, c("study_time_collected", "name")]
#   max <- max(abs(mat))
#   show_rnames <- TRUE
#   cluster_rows <- ifelse(nrow(mat) > 2, TRUE, FALSE)
#   p <- pheatmap(mat, annotation = anno, show_colnames = FALSE,
#                 show_rownames = show_rnames, cluster_cols = FALSE,
#                 cluster_rows = cluster_rows, color = palette,
#                 breaks = seq(-max, max, length.out = length(palette) + 1))
# } else if(type == "boxplot"){
#   p <- qplot(as.factor(study_time_collected), response, data = dt,
#              facets = analyte~name, geom = c("boxplot", "jitter"),
#            xlab = "time", ylab = ylab, color = study_time_collected)
#   print(p)
# } else if(type = "error"){
#   #grid.text(error_string)
#   grid.newpage()
#   grid.text("The datset you are trying to visualize does not follow HIPC standards.
#             Please use the built-in LabKey charts for this data.")
#   #ggplot() + geom_text(aes(x = 1, y = 1, label = error_string))
# }