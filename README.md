# ImmuneSpaceR <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->
[![R build status](https://github.com/RGLab/ImmuneSpaceR/workflows/R-CMD-check/badge.svg)](https://github.com/RGLab/ImmuneSpaceR/actions)
[![Codecov test coverage](https://codecov.io/gh/RGLab/ImmuneSpaceR/branch/master/graph/badge.svg)](https://codecov.io/gh/RGLab/ImmuneSpaceR?branch=master)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

[![Years in BioC](http://www.bioconductor.org/shields/years-in-bioc/ImmuneSpaceR.svg)](http://bioconductor.org/packages/release/bioc/html/ImmuneSpaceR.html)
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/ImmuneSpaceR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/ImmuneSpaceR)
[![Downloads](http://www.bioconductor.org/shields/downloads/ImmuneSpaceR.svg)](https://bioconductor.org/packages/stats/bioc/ImmuneSpaceR/)
[![Updated](https://bioconductor.org/shields/lastcommit/release/bioc/ImmuneSpaceR.svg)](http://bioconductor.org/packages/release/bioc/news/ImmuneSpaceR/NEWS)
<!-- badges: end -->


A thin wrapper around Rlabkey to access the [ImmuneSpace](https://www.immunespace.org) database from R.

This package simplifies access to the [HIPC](https://www.immuneprofiling.org/) ImmuneSpace database for R programmers. It takes advantage of the standardization of the database to hide all the [`Rlabkey`](https://cran.r-project.org/web/packages/Rlabkey/index.html) specific code away from the user. The study-specific datasets can be accessed via an object-oriented paradigm.



## Installation

Install from [Bioconductor](http://bioconductor.org/packages/release/bioc/html/ImmuneSpaceR.html):

``` r
install.packages("BiocManager")
BiocManager::install("ImmuneSpaceR")
```

Or install the latest development version via [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html):

``` r
install.packages("remotes")
remotes::install_github("RGLab/ImmuneSpaceR")
```



## Configuration

The ImmuneSpace Portal can be accessed via `ImmuneSpaceR` with the user's credentials. A `.netrc` file storing login and password information is required.

1. [Register](https://www.immunespace.org/login/home/register.view?)
1. Create a netrc file with your ImmuneSpace credetntials using `interactive_netrc()` function in R:

``` r
library(ImmuneSpaceR)
interactive_netrc()
```

If you're familiar with the command-line interface, see [the introductory vignette](https://rglab.github.io/ImmuneSpaceR/articles/Intro_to_ImmuneSpaceR.html).



## Usage

### Create a connection

The general idea is that the user creates an instance of an `ImmuneSpaceConnection` class. The instance configures itself to connect to a specific study, and datasets and gene expression matrices can be retrieved by name.

For example:

``` r
library(ImmuneSpaceR)
con <- CreateConnection("SDY269")
```

will create an instance of SDY269.


### List datasets

Datasets can be listed by:

``` r
con$listDatasets()
```

which will print names of available datasets and gene expression matrices.


### Retrieve datasets

Gene expression matrices or datasets can be retreived by:

``` r
LAIV2008 <- con$getGEMatrix("SDY269_PBMC_LAIV_Geo")
elisa <- con$getDataset("elisa")
```

The connection object *caches* data, so once it is retrieved, the next time you access it, it will use the local cached copy. The package uses a [R6](https://cran.r-project.org/web/packages/R6/index.html) class system to represent the connection to a study and get around some of R's copy-on-change behaviour.


### Visualize

The `plot` method uses [`ggplot2`](https://cran.r-project.org/web/packages/ggplot2/index.html) functions to generate visualizations of datasets, leveraging the standardized dataset tables.

``` r
con$plot("hai")
```



## Examples & Documentation

For more advanced examples and detailed documentation, see [the package vignettes](http://rglab.github.io/ImmuneSpaceR/articles/) and the reports available on [ImmuneSpace](https://www.immunespace.org/).



## Contributing

If you'd like to report bugs/issues/feature requests or contribute to the package, please see [the contributing guidelines](./CONTRIBUTING.md) and join [our Slack workspace](https://immunespace.herokuapp.com/).
