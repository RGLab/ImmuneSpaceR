## ----knitr-opts, echo = FALSE, message = FALSE, cache = FALSE------------
library(knitr)
opts_chunk$set(cache=FALSE, echo=TRUE, message=FALSE, warning=FALSE,
               fig.width=8, fig.height=4, dpi=100, fig.align="center")

## ----libraries, cache=FALSE----------------------------------------------
library(ImmuneSpaceR)
library(ggplot2)
library(data.table)
library(ggthemr)
ggthemr('solarized')

## ----connection----------------------------------------------------------
study <- CreateConnection(c("SDY269"))
dt_hai <- study$getDataset("hai", reload=TRUE)
dt_fcs <- study$getDataset("fcs_analyzed_result", reload=TRUE)
dt_elispot <- study$getDataset("elispot", reload=TRUE)

## ----data-subset---------------------------------------------------------
# Compute max fold change for HAI, and remove time zero
dt_hai <- dt_hai[,hai_response:=value_reported/value_reported[study_time_collected==0],
                 by="virus_strain,arm_name,subject_accession"][study_time_collected==28]
dt_hai <- dt_hai[,list(hai_response=max(hai_response)),by="arm_name,subject_accession"]

# Define variable for ELISPOT, keep only the IgG class
dt_elispot <- dt_elispot[,elispot_response:=spot_number_reported+1][study_time_collected==7 & analyte=="IgG"]
# Compute % plasmablasts
dt_fcs <- dt_fcs[,fcs_response:=(as.double(population_cell_number)+1) /
                   as.double(base_parent_population)][study_time_collected==7]

## ----merging-------------------------------------------------------------
# Let's key the different datasets
setkeyv(dt_hai, c("subject_accession"))
setkeyv(dt_fcs, c("subject_accession"))
setkeyv(dt_elispot, c("subject_accession"))
dt_all <- dt_hai[dt_fcs, nomatch=0][dt_elispot, nomatch=0]

## ----plot1, dev='CairoPNG'-----------------------------------------------
ggplot(dt_all, aes(x=as.double(fcs_response), y=elispot_response, color=arm_name)) +
  geom_point() + scale_y_log10() + scale_x_log10() + geom_smooth(method="lm") +
  xlab("Total plasmablasts (%)") + ylab("Influenza specific cells\n (per 10^6 PBMCs)")

## ----plot2, dev='CairoPNG'-----------------------------------------------
ggplot(dt_all, aes(x=as.double(hai_response), y=elispot_response, color=arm_name)) +
  geom_point() + scale_x_continuous(trans="log2") + scale_y_log10() +
  geom_smooth(method="lm") + xlab("HAI fold") + ylab("Influenza specific cells\n (per 10^6 PBMCs)")

