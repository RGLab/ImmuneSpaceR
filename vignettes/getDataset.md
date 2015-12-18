---
title: "getDataset"
author: "Renan Sauteraud"
date: "2015-12-14"
output: html_document
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Downloading tables with getDataset}
---



# getDataset
This vignette shows detailed examples for all functionalities of the `getDataset`
function.

## <a id="contents"></a>Contents
1. [Connections](#cons)
1. [List the datasets](#datasets)
2. [Download a dataset](#download)
3. [Filters](#filters)
4. [Views](#views)
5. [Caching](#caching)
6. [sessionInfo](#sessioninfo)


## <a id="cons"></a>Connections
As explained into the user guide vignette, datasets must be downloaded from
`ImmuneSpaceConnection` objects. We must first instantiate a connection to 
the study or studies of interest. Throughout this vignette, we will use two
connections, one to a single study, and one to to all available data.


```r
library(ImmuneSpaceR)
```

```
## Loading required package: ggthemr
## Loading required package: ggplot2
```

```r
sdy269 <- CreateConnection("SDY269")
all <- CreateConnection("")
```


## <a id="datasets"></a>List the datasets
Now that the connections have been instantiated, we can start downloading from
them. But we need to figure out which datasets are available within our chosen
studies.

Printing the connections will, among other information, list the datasets
availables. The `listDatasets` method will display only the information we are
looking for.


```r
sdy269$listDatasets()
```

```
## datasets
## 	fcs_sample_files
## 	mbaa
## 	fcs_analyzed_result
## 	gene_expression_files
## 	neut_ab_titer
## 	fcs_control_files
## 	elisa
## 	hai
## 	cohort_membership
## 	elispot
## 	pcr
## 	demographics
## 	hla_typing
## 	kir_typing
## Expression Matrices
## 	LAIV_2008
## 	TIV_2008
```

```r
all$listDatasets()
```

```
## datasets
## 	fcs_sample_files
## 	mbaa
## 	fcs_analyzed_result
## 	gene_expression_files
## 	neut_ab_titer
## 	fcs_control_files
## 	elisa
## 	hai
## 	cohort_membership
## 	elispot
## 	pcr
## 	demographics
## 	hla_typing
## 	kir_typing
## Expression Matrices
## 	TIV_2010
## 	VLplus
## 	VLminus
## 	Cohort2_old
## 	Cohort1_young
## 	Saline_group2
## 	Saline_group1
## 	Pneumovax23_group2
## 	Pneumovax23_group1
## 	Fluzone_group2
## 	Fluzone_group1
## 	TIV_2011
## 	TIV_2007
## 	TIV_2008
## 	LAIV_2008
```

Naturally, `all` contains every dataset available on ImmuneSpace as it combines
all available studies.

Additionaly, when creating connection object with `verbose = TRUE`, a call to
the `getDataset` function with an invalid dataset name will return the list of
valid datasets.

[Back to top](#contents)

## <a id="download"></a>Download

Calling `getDataset` returns a selected dataset as it is displayed on ImmuneSpace. 

```r
hai_269 <- sdy269$getDataset("hai")
hai_all <- all$getDataset("hai")

print(head(hai_269))
```

```
##    participant_id age_reported gender  race          cohort
## 1:  SUB112841.269           28 Female White  TIV Group 2008
## 2:  SUB112834.269           27   Male White  TIV Group 2008
## 3:  SUB112868.269           37   Male White LAIV group 2008
## 4:  SUB112836.269           28 Female White LAIV group 2008
## 5:  SUB112859.269           32   Male White LAIV group 2008
## 6:  SUB112865.269           25 Female White LAIV group 2008
##    study_time_collected study_time_collected_unit
## 1:                    0                      Days
## 2:                   28                      Days
## 3:                    0                      Days
## 4:                   28                      Days
## 5:                    0                      Days
## 6:                    0                      Days
##                              virus_strain value_reported
## 1: A/Uruguay/716/2007  NYMC X-175C (H3N2)              5
## 2:              A/Brisbane/59/2007 (H1N1)            160
## 3:           A/South Dakota/6/2007 (H1N1)             20
## 4: A/Uruguay/716/2007  NYMC X-175C (H3N2)              5
## 5:           A/South Dakota/6/2007 (H1N1)              5
## 6:           A/South Dakota/6/2007 (H1N1)              5
```
Because some datasets such as flow cytometry results can contain a large number
of rows, the function returns `data.table` objects to improve performance. This 
is especially important with multi-study connections.

[Back to top](#contents)

## <a id="filters"></a>Filters

The datasets can be filtered before download. Filters should be created using 
`Rlabkey`'s `makeFilter` function.

Each filter is composed of three part:
 - A column name or column label
 - An operator
 - A value or array of values separated by a semi-colon

```r
library(Rlabkey)
```

```
## Loading required package: RCurl
## Loading required package: bitops
## Loading required package: rjson
```

```r
# Get participants under age of 30
young_filter <- makeFilter(c("age_reported", "LESS_THAN", 30))
# Get a specific list of two participants
pid_filter <- makeFilter(c("participantid", "IN", "SUB112841.269;SUB112834.269"))
```
For a list of available operators, see `?makeFilter`.


```r
# HAI data for participants of study SDY269 under age of 30
hai_young <- sdy269$getDataset("hai", colFilter = young_filter)
# List of participants under age 30
demo_young <- all$getDataset("demographics", colFilter = young_filter)
# ELISPOT assay results for two participants
mbaa_pid2 <- all$getDataset("elispot", colFilter = pid_filter)
```

Note that filtering is done before download. When performance is a concern, it
is faster to do the filtering via the `colfFilter` argument than on the returned 
table.


## <a id="views"></a>Views
Any dataset grid on ImmuneSpace offers the possibility to switch views between 
'Default' and 'Full'. The Default view contains information that is directly
relevant to the user. Sample description and results are joined with basic 
demographic.
However, this is not the way data is organized in the database. The 'Full' view 
is a representation of the data as it is stored on 
[ImmPort](http://www.immport.org/immport-open/public/home/home). The accession 
columns are used under the hood for join operations. They will be useful to
developers and user writing reports to be displayed in ImmuneSpace studies.


![](./img/getDataset-views.png)
Screen capture of the button bar of a dataset grid on ImmuneSpace

The `original_view` argument decides which view is downloaded. If set to `TRUE`,
the 'Full' view is returned.

```r
full_hai <- sdy269$getDataset("hai", original_view = TRUE)
print(colnames(full_hai))
```

```
##  [1] "participant_id"            "arm_accession"            
##  [3] "biosample_accession"       "expsample_accession"      
##  [5] "experiment_accession"      "study_accession"          
##  [7] "study_time_collected"      "study_time_collected_unit"
##  [9] "virus_strain"              "value_reported"           
## [11] "value_preferred"           "unit_reported"            
## [13] "unit_preferred"
```


For additional information, refer to the 'Working with tabular data' video 
tutorial available in the menu bar on any page of the portal.

[Back to top](#contents)

## <a id="caching"></a>Caching

As explained in the user guide, the `ImmuneSpaceConnection` class is a Reference
class. It means its objects have fields accessed by reference. As a consequence,
they can be modified without making a copy of the entire object.
ImmuneSpaceR uses this feature to store downloaded datasets and expression 
matrices. Subsequent calls to `getDataset` with the same input will be faster.

See `?setRefClass` for more information about reference classes.

We can see the data currently cached using the `data_cache` field. This is not
intended to be used for data manipulation and only displayed here to explain
what gets cached and.

```r
pcr <- sdy269$getDataset("pcr")
names(sdy269$data_cache)
```

```
## [1] "GE_matrices"           "filter_state_hai"      "hai"                  
## [4] "filter_state_hai_full" "hai_full"              "filter_state_pcr"     
## [7] "pcr"
```

Different [views](#views) are saved separately.

```r
pcr_ori <- sdy269$getDataset("pcr", original_view = TRUE)
names(sdy269$data_cache)
```

```
## [1] "GE_matrices"           "filter_state_hai"      "hai"                  
## [4] "filter_state_hai_full" "hai_full"              "filter_state_pcr"     
## [7] "pcr"                   "filter_state_pcr_full" "pcr_full"
```

Because of the infinite number of filters and combinations of filters, we do not
cache filtered datasets.

If, for any reason, a specific dataset needs to be redownloaded, the `reload` 
argument will clear the cache for that specific `getDataset` call and download
the table again.

```r
hai_269 <- sdy269$getDataset("hai", reload = TRUE)
```

Finally, it is possible to clear every cached dataset (and expression matrix).

```r
sdy269$clear_cache()
names(sdy269$data_cache)
```

```
## [1] "GE_matrices"
```

Again, the `data.cache` field should never be modified manually. When in doubt,
simply reload the dataset.

[Back to top](#contents)

## <a id="sessioninfo"></a> sessionInfo()

```r
sessionInfo()
```

```
## R Under development (unstable) (2015-11-03 r69594)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.3 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] Rlabkey_2.1.128      rjson_0.2.15         RCurl_1.95-4.7      
##  [4] bitops_1.0-6         ImmuneSpaceR_0.2.41  ggthemr_1.0.1       
##  [7] ggplot2_1.0.1        knitr_1.11           data.table_1.9.6    
## [10] devtools_1.9.1       BiocInstaller_1.21.2
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.2           formatR_1.2.1         RColorBrewer_1.1-2   
##  [4] plyr_1.8.3            tools_3.3.0           zlibbioc_1.17.0      
##  [7] digest_0.6.8          evaluate_0.8          memoise_0.2.1        
## [10] preprocessCore_1.33.0 gtable_0.1.2          parallel_3.3.0       
## [13] proto_0.3-10          stringr_1.0.0         gtools_3.5.0         
## [16] caTools_1.17.1        grid_3.3.0            Biobase_2.31.0       
## [19] pheatmap_1.0.7        gdata_2.17.0          reshape2_1.4.1       
## [22] magrittr_1.5          scales_0.3.0          gplots_2.17.0        
## [25] MASS_7.3-44           BiocGenerics_0.17.1   colorspace_1.2-6     
## [28] KernSmooth_2.23-15    stringi_1.0-1         affy_1.49.0          
## [31] munsell_0.4.2         chron_2.3-47          affyio_1.41.0
```
