## ----knitr-opts, echo = FALSE, message = FALSE, cache = FALSE------------
library(knitr)
opts_chunk$set(cache=FALSE, echo=TRUE, message=FALSE, warning=FALSE,
               fig.width=10, fig.height=14, dpi=100, fig.align="center")

## ----libraries, cache=FALSE----------------------------------------------
library(ImmuneSpaceR)
library(ggplot2)
library(data.table)
library(ggthemr)
ggthemr('solarized')

## ----connection----------------------------------------------------------
study <- CreateConnection(c("SDY180"))
dt_fcs <- study$getDataset("fcs_analyzed_result", reload=TRUE)

## ----data-subset---------------------------------------------------------
dt_fcs19 <- dt_fcs[population_name_reported%like%"Plasma"]
dt_fcs19 <- dt_fcs19[,arm_name:=gsub("Study g", "G", arm_name),]


## ----data-summary--------------------------------------------------------
dt_fcs19_median <- dt_fcs19[,.(median_cell_reported=median(as.double(population_cell_number)+1,na.rm=TRUE)),by=.(arm_name,study_time_collected,population_name_reported)]

## ----, dev='CairoPNG'----------------------------------------------------
ggplot(dt_fcs19, aes(x=as.factor(study_time_collected), y=as.double(population_cell_number)+1))+geom_boxplot()+geom_jitter()+scale_y_log10()+facet_grid(arm_name~population_name_reported, scale="free")+xlab("Time")+ylab(expression(paste("Number of cells/",mu,"l")))+geom_line(data=dt_fcs19_median,aes(x=as.factor(study_time_collected), y=as.double(median_cell_reported),group=1),color="black",size=1.2)+labs(title="Plasma cell abundance after vaccination")

