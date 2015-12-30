---
title: "Expression matrices with ImmuneSpaceR"
author: "Renan Sauteraud"
date: "2015-12-18"
output: html_document
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Handling expression matrices with ImmuneSpaceR}
---





# Downloading expression data using ImmuneSpaceR
This vignette shows detailed examples for downloading expression matrices,

## <a id="contents"></a>Contents
1. [Connections](#cons)
1. [List the expression matrices](#matrices)
2. [Download an expression matrix](#download)
3. [Summarized matrices](#summary)
4. [Combining matrices](#multi)
5. [Caching](#caching)
6. [sessionInfo](#sessioninfo)


## <a id="cons"></a>Connections
As explained into the user guide vignette, datasets must be downloaded from
`ImmuneSpaceConnection` objects. We must first instantiate a connection to 
the study or studies of interest. Throughout this vignette, we will use two
connections, one to a single study, and one to to all available data.


```r
library(ImmuneSpaceR)
sdy269 <- CreateConnection("SDY269")
all <- CreateConnection("")
```


## <a id="matrices"></a>List the expression matrices
Now that the connections have been instantiated, we can start downloading from
them. But we need to figure out which processed matrices are available within 
our chosen studies.

On www.immunespace.org, in the study of interest or at the project level, the
**Gene expression matrices** table will show the available runs.

Printing the connections will, among other information, list the datasets
availables. The `listDatasets` method will only display the downloadable data.
looking for. With `which = "expression"`, the datasets wont be printed.


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
Using `wich = "expression"`, we can remove the datasets from the output.

```r
all$listDatasets(which = "expression")
```

```
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

Naturally, `all` contains every processed matrices available on ImmuneSpace as 
it combines all available studies.

[Back to top](#contents)

## <a id="download"></a>Download

### By run name

The `getGEMatrix` function will accept any of the run names listed in the 
connection.

```r
TIV_2008 <- sdy269$getGEMatrix("TIV_2008")
```

```
## Downloading matrix..
## Downloading Features..
## Constructing ExpressionSet
```

```r
TIV_2011 <- all$getGEMatrix(x = "TIV_2011")
```

```
## Downloading matrix..
## Downloading Features..
## Constructing ExpressionSet
```
The matrices are returned as `ExpressionSet` where the phenoData slot contains
basic demographic information and the featureData slot shows a mapping of probe
to official gene symbols.

```r
TIV_2008
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 54715 features, 80 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: BS586131 BS586187 ... BS586267 (80 total)
##   varLabels: biosample_accession participant_id ...
##     study_time_collected_unit (5 total)
##   varMetadata: labelDescription
## featureData
##   featureNames: 1007_PM_s_at 1053_PM_at ... AFFX-r2-TagQ-5_at
##     (54715 total)
##   fvarLabels: FeatureId gene_symbol
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation:
```

### By cohort

The `cohort` argument can be used in place of the run name (`x`). Likewise, the
list of valid cohorts can be found in the Gene expression matrices table.

```r
LAIV_2008 <- sdy269$getGEMatrix(cohort = "LAIV group 2008")
```

```
## Downloading matrix..
## Constructing ExpressionSet
```
Note that when cohort is used, `x` is ignored.

[Back to top](#contents)


## <a id="summary"></a>Summarized matrices

By default, the returned `ExpressionSet`s have probe names as features (or rows).
However, multiple probes often match the same gene and merging experiments from
different arrays is impossible at feature level.
When they are available, the `summary` argument allows to return the matrices 
with gene symbols instead of probes.

```r
TIV_2008_sum <- sdy269$getGEMatrix("TIV_2008", summary = TRUE)
```

```
## Downloading matrix..
## Constructing ExpressionSet
```
Probes that do not map to a unique gene are removed and expression is averaged 
by gene.

```r
TIV_2008_sum
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 16910 features, 80 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: BS586131 BS586187 ... BS586267 (80 total)
##   varLabels: biosample_accession participant_id ...
##     study_time_collected_unit (5 total)
##   varMetadata: labelDescription
## featureData
##   featureNames: DDR1 RFC2 ... NUS1P3 (16910 total)
##   fvarLabels: FeatureId gene_symbol
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation:
```

[Back to top](#contents)


## <a id="multi"></a>Combining matrices

With the summarized version of the matrices, it is easy to combine multiple 
experiments, even accross different types of array. This is done internally in
`getGEMatrix` when multiple runs or multiple cohorts are specified.


```r
# Within a study
em269 <- sdy269$getGEMatrix(c("TIV_2008", "LAIV_2008"))
```

```
## Downloading matrix..
## Downloading matrix..
## Downloading Features..
## Constructing ExpressionSet
## Constructing ExpressionSet
```

```r
# Combining accross studies
TIV_seasons <- all$getGEMatrix(c("TIV_2008", "TIV_2010"))
```

```
## Downloading matrix..
## Downloading matrix..
## Downloading Features..
## Downloading Features..
## Constructing ExpressionSet
## Constructing ExpressionSet
```



## <a id="caching"></a>Caching

As explained in the user guide, the `ImmuneSpaceConnection` class is a Reference
class. It means its objects have fields accessed by reference. As a consequence,
they can be modified without making a copy of the entire object.
ImmuneSpaceR uses this feature to store downloaded datasets and expression 
matrices. Subsequent calls to `getGEMatrix` with the same input will be faster.

See `?setRefClass` for more information about reference classes.

We can see the data currently cached using the `data_cache` field. This is not
intended to be used for data manipulation and only displayed here to explain
what gets cached and.

```r
pcr <- sdy269$getDataset("pcr")
names(sdy269$data_cache)
```

```
## [1] "GE_matrices"      "TIV_2008"         "featureset_18"   
## [4] "filter_state_pcr" "pcr"
```

If, for any reason, a specific marix needs to be redownloaded, the `reload` 
argument will clear the cache for that specific `getGEMatrix` call and download
the file and metadata again.

```r
TIV_2008 <- sdy269$getGEMatrix("TIV_2008", reload = TRUE)
```

```
## Downloading matrix..
## Constructing ExpressionSet
```

Finally, it is possible to clear every cached expression matrix (and dataset).

```r
sdy269$clear_cache()
names(sdy269$data_cache)
```

```
## [1] "GE_matrices"
```

Again, the `data.cache` field should never be modified manually. When in doubt,
simply reload the expression matrix.

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
## [1] ImmuneSpaceR_0.2.42  ggthemr_1.0.1        ggplot2_1.0.1       
## [4] knitr_1.11           data.table_1.9.6     devtools_1.9.1      
## [7] BiocInstaller_1.21.2
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.2           RColorBrewer_1.1-2    formatR_1.2.1        
##  [4] plyr_1.8.3            highr_0.5.1           bitops_1.0-6         
##  [7] tools_3.3.0           zlibbioc_1.17.0       digest_0.6.8         
## [10] evaluate_0.8          memoise_0.2.1         preprocessCore_1.33.0
## [13] gtable_0.1.2          parallel_3.3.0        proto_0.3-10         
## [16] stringr_1.0.0         gtools_3.5.0          caTools_1.17.1       
## [19] grid_3.3.0            Biobase_2.31.0        Rlabkey_2.1.128      
## [22] pheatmap_1.0.7        gdata_2.17.0          reshape2_1.4.1       
## [25] magrittr_1.5          scales_0.3.0          gplots_2.17.0        
## [28] codetools_0.2-14      MASS_7.3-44           BiocGenerics_0.17.1  
## [31] colorspace_1.2-6      KernSmooth_2.23-15    stringi_1.0-1        
## [34] affy_1.49.0           RCurl_1.95-4.7        munsell_0.4.2        
## [37] chron_2.3-47          rjson_0.2.15          affyio_1.41.0
```
