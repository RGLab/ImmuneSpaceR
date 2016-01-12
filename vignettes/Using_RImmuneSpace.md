---
date: "2016-01-08"
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{An introduction to using the ImmuneSpaceR package}
---



# A simple introduction on using the ImmuneSpaceR package

This package provides a *thin* wrapper around `Rlabkey` and connects to the **ImmuneSpace* database, making it easier to fetch *datasets*, including gene expression data, hai, and so forth, from specific studies. 

## Contents
1. [Configuration](#configuration)
2. [Connections](#connections)
3. [Datasets](#datasets)
4. [Gene expression](#ge)
5. [Plots](#quick_plot)

## <a id="configuration"></a>Configuration

In order to connect to ImmuneSpace, you will need a `.netrc` file in your home 
directory that will contain a `machine` name (hostname of ImmuneSpace), and 
`login` and `password`. See [here](https://www.labkey.org/wiki/home/Documentation/page.view?name=netrc) for more information.

A netrc file may look like this:
```
machine www.immunespace.org
login myuser@domain.com
password supersecretpassword
```

### Set up your netrc file now
Put it in your home directory. 
If you type:
``` 
ls ~/.netrc
```
at the command prompt, you should see it there. If it's not there, create one 
now. Make sure you have a valid login and password. If you don't have one, go to
[ImmuneSpace](http://www.immunespace.org) now and set yourself up with an 
account. 

## <a id="connections"></a>Instantiate a connection

We'll be looking at study `SDY269`. If you want to use a different study, change
that string. The connections have state, so you can instantiate multiple 
connections to different studies simultaneously.


```r
library(ImmuneSpaceR)
sdy269 <- CreateConnection(study = "SDY269")
sdy269
```

```
## Immunespace Connection to study SDY269
## URL: https://test.immunespace.org/Studies/SDY269
## User: unknown_user at not_a_domain.com
## Available datasets
## 	mbaa
## 	fcs_analyzed_result
## 	gene_expression_files
## 	neut_ab_titer
## 	fcs_control_files
## 	hai
## 	cohort_membership
## 	elispot
## 	pcr
## 	demographics
## 	hla_typing
## 	kir_typing
## 	elisa
## 	fcs_sample_files
## 	MATT_Snapshot_2
## Expression Matrices
## 	LAIV_2008
## 	TIV_2008
```

The call to `CreateConnection` instantiates the connection Printing the object 
shows where it's connected, to what study, and the available data sets and gene 
expression matrices.

Note that when a script is running on ImmuneSpace, some variables set in the 
global environments will automatically indicate which study should be used and 
the `study` argument can be skipped.

To combine data across studies, create a connection that contain all available
studies by leaving the `study` as an empty string.

```r
all <- CreateConnection("")
```

## <a id="datasets"></a>Fetching data sets

We can grab any of the datasets listed in the connection.


```r
sdy269$getDataset("hai")
```

```
##      participant_id age_reported gender  race          cohort
##   1:  SUB112841.269           28 Female White  TIV Group 2008
##   2:  SUB112834.269           27   Male White  TIV Group 2008
##   3:  SUB112868.269           37   Male White LAIV group 2008
##   4:  SUB112836.269           28 Female White LAIV group 2008
##   5:  SUB112859.269           32   Male White LAIV group 2008
##  ---                                                         
## 332:  SUB112883.269           23 Female Asian LAIV group 2008
## 333:  SUB112878.269           28 Female White  TIV Group 2008
## 334:  SUB112834.269           27   Male White  TIV Group 2008
## 335:  SUB112863.269           29 Female White  TIV Group 2008
## 336:  SUB112877.269           27   Male White  TIV Group 2008
##      study_time_collected study_time_collected_unit
##   1:                    0                      Days
##   2:                   28                      Days
##   3:                    0                      Days
##   4:                   28                      Days
##   5:                    0                      Days
##  ---                                               
## 332:                   28                      Days
## 333:                    0                      Days
## 334:                   28                      Days
## 335:                   28                      Days
## 336:                    0                      Days
##                                virus_strain value_reported
##   1: A/Uruguay/716/2007  NYMC X-175C (H3N2)              5
##   2:              A/Brisbane/59/2007 (H1N1)            160
##   3:           A/South Dakota/6/2007 (H1N1)             20
##   4: A/Uruguay/716/2007  NYMC X-175C (H3N2)              5
##   5:           A/South Dakota/6/2007 (H1N1)              5
##  ---                                                      
## 332:           A/South Dakota/6/2007 (H1N1)             20
## 333: A/Uruguay/716/2007  NYMC X-175C (H3N2)              5
## 334:                      B/Brisbane/3/2007             40
## 335: A/Uruguay/716/2007  NYMC X-175C (H3N2)             40
## 336: A/Uruguay/716/2007  NYMC X-175C (H3N2)              5
```


The *sdy269* object is an **R5** class, so it behaves like a true object. 
Functions (like `getDataset`) are members of the object, thus the `$` semantics 
to access member functions.

The first time you retrieve a data set, it will contact the database. The data 
is cached locally, so the next time you call `getDataset` on the same dataset, 
it will retrieve the cached local copy. This is much faster. 


To get only a subset of the data and speed up the download, filters can be 
passed to `getDataset`. The filters are created using the `makeFilter` function
of the `Rlabkey` package.

```r
library(Rlabkey)
myFilter <- makeFilter(c("gender", "EQUAL", "Female"))
hai <- sdy269$getDataset("hai", colFilter = myFilter)
```
See `?makeFilter` for more information on the syntax.

For more information about `getDataset`'s options, refer to the dedicated vignette.


## <a id="datasets"></a>Fetching expression matrices
We can also grab a gene expression matrix


```r
sdy269$getGEMatrix("LAIV_2008")
```

```
## Downloading matrix..
## Downloading Features..
## Constructing ExpressionSet
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 54715 features, 83 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: BS586216 BS586160 ... BS586111 (83 total)
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

The object contacts the DB and downloads the matrix file. This is stored and 
cached locally as a `data.table`. The next time you access it, it will be much 
faster since it won't need to contact the database again.

It is also possible to call this function using multiple matrix names. In this
case, all the matrices are downloaded and combined into a single `ExpressionSet`.

```r
sdy269$getGEMatrix(c("TIV_2008", "LAIV_2008"))
```

```
## Downloading matrix..
## Downloading matrix..
## Constructing ExpressionSet
## Constructing ExpressionSet
## Combining ExpressionSets
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 54715 features, 163 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: BS586131 BS586187 ... BS586111 (163 total)
##   varLabels: biosample_accession participant_id ...
##     study_time_collected_unit (5 total)
##   varMetadata: labelDescription
## featureData
##   featureNames: 1 2 ... 54715 (54715 total)
##   fvarLabels: FeatureId gene_symbol
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation:
```

Finally, the summary argument will let you download the matrix with gene symbols
in place of priobe ids.

```r
gs <- sdy269$getGEMatrix("TIV_2008", summary = TRUE)
```

```
## Downloading matrix..
## Constructing ExpressionSet
```


If the connection was created with `verbose = TRUE`, some functions will display
additional informations such as the valid dataset names.

## Quick plots
A quick plot of a data set can be generated using the `quick_plot` function.

`quick_plot` automatically chooses the type of plot depending on the selected 
dataset.


```r
sdy269$quick_plot("hai")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png) 

```r
sdy269$quick_plot("elisa")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-2.png) 

However, the `type` argument can be used to manually select from "boxplot",
"heatmap", "violin" and "line".

## <a id="crossstudy"></a>Cross study connections
To fetch data from multiple studies, simply create a connection at the project level.


```r
con <- CreateConnection("")
```

This will instantiate a connection at the `Studies` level. Most functions work
cross study connections just like they do on single studies.

You can get a list of datasets and gene expression matrices available accross 
all studies.

```r
con
```

```
## Immunespace Connection to study Studies
## URL: https://www.immunespace.org/Studies/
## User: unknown_user at not_a_domain.com
## Available datasets
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

Likewise, `quick_plot` will plot accross studies. Note that in most cases the
datasets will have too many cohorts/subjects, making the filtering of the data
a necessity. The `colFilter` argument can be used here, as described in the 
`getDataset` sectionm.

```r
plotFilter <- makeFilter(c("cohort", "IN", "TIV 2010;TIV Group 2008"))
con$quick_plot("elispot", filter = plotFilter)
```

```
## Warning: `axis.ticks.margin` is deprecated. Please set `margin` property of
## `axis.text` instead
```

```
## Warning: New theme missing the following elements: panel.ontop,
## strip.switch.pad.grid, strip.switch.pad.wrap
```

```
## Error: No geom called GeomJitter.
```
The figure above shows the ELISPOT results for two different years of TIV 
vaccine cohorts.

