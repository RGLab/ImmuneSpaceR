## ----knitr, echo = FALSE-------------------------------------------------
library(knitr)
opts_chunk$set(message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8)

## ------------------------------------------------------------------------
library(ImmuneSpaceR)
library(data.table)
library(ggplot2)
library(ggthemr)

## ------------------------------------------------------------------------
con <- CreateConnection("SDY144")
flow <- con$getDataset("fcs_analyzed_result")
hai  <- con$getDataset("hai")
vn   <- con$getDataset("neut_ab_titer")

## ----subset--------------------------------------------------------------
pb <- flow[population_name_reported %in% c("Plasma cells,Freq. of,B lym CD27+",
                                           "Plasmablast,Freq. of,Q3: CD19+, CD20-")]
pb <- pb[, population_cell_number := as.numeric(population_cell_number)]
pb <- pb[study_time_collected == 7 & study_time_collected_unit == "Days"] #13 subjects
pb <- pb[, list(subject_accession, population_cell_number, population_name_reported)]

## ----FC------------------------------------------------------------------
# HAI
hai <- hai[,response:=value_reported/value_reported[study_time_collected==0],
                 by="virus_strain,arm_name,subject_accession"][study_time_collected==30]
hai <- hai[, list(subject_accession, virus_strain, response)]
dat_hai <- merge(hai, pb, by = "subject_accession", allow.cartesian = TRUE)
# VN
vn <- vn[, response:=value_reported/value_reported[study_time_collected==0],
                 by="virus_strain,arm_name,subject_accession"][study_time_collected==30]
vn <- vn[, list(subject_accession, virus_strain, response)]
dat_vn <- merge(vn, pb, by = "subject_accession", allow.cartesian = TRUE)

## ------------------------------------------------------------------------
ggthemr('solarized')

## ----HAI, dev='CairoPNG'-------------------------------------------------
ggplot(dat_hai, aes(x = population_cell_number, y = response)) +
  geom_point() + geom_smooth(method = "lm") +
  facet_grid(virus_strain~population_name_reported, scale = "free") +
  xlab("Number of cells") + ylab("HI fold-increase Day 30 vs. baseline")

## ----VN, dev='CairoPNG'--------------------------------------------------
ggplot(dat_vn, aes(x = population_cell_number, y = response)) +
  geom_point() + geom_smooth(method = "lm") +
  facet_grid(virus_strain~population_name_reported, scale = "free") +
  xlab("Number of cells") + ylab("VN fold-increase Day 30 vs. baseline")

